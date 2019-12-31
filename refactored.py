import math
from numba import jit, njit
from numba.typed import List, Dict
from time import time
import numpy as np
import os

start = time()

yr2sec = 60*60*24*365.25          #nr of seconds in a year

dx              = 1e3             # width of each node/column (m)
imax            = 100             # number of nodes
tmax            = 50000*yr2sec    # total amount of time to be modelled in seconds [in years = (x*yr2sec)]
dtout           = tmax/100        # nr of years between write output
dtout_progress  = 10*yr2sec       # nr of years between progress bar update

        
nrOfGrainSizes  = 2

#k   = list(range(nrOfGrainSizes))
#q0  = list(range(nrOfGrainSizes))

k = np.zeros(shape=(2))
q0 = np.zeros(shape=(2))
#print(k)
q0[0]= 4*1.e-7                    # Gravel input (m2/s)  Attention: will be overwritten in setBoudnaryCondtitionValues if declared there!
q0[1]= 3*1.e-7                    # Sand input (m2/s)  Attention: will be overwritten in setBoudnaryCondtitionValues if declared there!
    
k[0]= 1*3.2e-4                    # Gravel diffusivity (m2/s)  Attention: will be overwritten in setBoudnaryCondtitionValues if declared there!
k[1]= 2*3.2e-4                    # Sand diffusivity (m2/s)  Attention: will be overwritten in setBoudnaryCondtitionValues if declared there!

subsidenceRate = 0#2.e-6


def makeList(mylist):
    l = List()
    for i in (mylist):
        l.append(i)
    return l
  
def nbZeros(shape):
    if not hasattr(shape, '__len__'):
        result = _nbZeros_1d(shape)
    else:
        result = _nbZeros_nd(shape)
    return result
      
  
@njit
def _nbZeros_nd(shape):
    l = List()
    for i in range(len(shape)):
        k = List()
        for j in range(shape[i]):
            k.append(0.)
        l.append(k)
    return l

@njit
def _nbZeros_1d(length):
    l = List()
    for i in range(length):
        l.append(0.)
    return l
  
@njit
def nbRange(length):
    l = List()
    for i in range(length):
        l.append(i)
    return l


@njit
def nbStackedLists(shape):
    #import pdb; pdb.set_trace()
    l = List()
    for i in range(shape[0]):
        #k = List()
        for j in range(1,len(shape)):
            m = List()
            for j in range(shape[j]):
                m.append(0.)
            #k.append(m)
        l.append(m)
    return l

@njit
def nbmin(values):
    minVal = values[0]
    for val in values:
        if (val < minVal): minVal = val
    return minVal

@njit
def nbmax(values):
    maxVal = values[0]
    for val in values:
        if (val > maxVal): maxVal = val
    return maxVal

@njit
def numbaSumList(mylist):
    total = 0
    for value in mylist:
        total += value
    return total

#import pdb; pdb.set_trace()

totalSedContent = nbStackedLists((imax+1, nrOfGrainSizes))

heights = Dict()
heights['totalHeight'] = nbZeros(imax+1)
heights['bedrockHeight'] = nbZeros(imax+1)

for i in range(nrOfGrainSizes):
    heights['totalSedContent_'+str(i)] = nbZeros(imax+1)

nodes = Dict()
nodes['nodeFill'] = nbStackedLists((imax+1, 1))
nodes['depTime'] = nbStackedLists((imax+1, 1))
nodes['sedContent'] = nbStackedLists((2, imax+1))

#import pdb; pdb.set_trace()

#for i in range(imax):
    #l = List()
    #for f in range(nrOfGrainSizes):
        #l.append(0.)
    #nodes['sedContent'].append(l)

totalInput= 0
totalOutput= 0
InputPerGrainSize = np.zeros(shape=(2))
OutputPerGrainSize = np.zeros(shape=(2))

dy = 1                            # vertical grid spacing (m)  Important: MUST BE 1(as integer) for code to work (in this version)!
tout = 0.                         # threshold to write next output in years (increased every write action by dtout)
tout_progress = 2*dtout_progress    # threshold to update progress bar, will be increased whenever it is reached
outputTimestep = 0

#for i in range(imax+1):
    #heights['totalHeight'][i] = 0
    ##newHeight[i] = 0
    #heights['bedrockHeight'][i] = 0
    ##columns[i] = []
    #for f in range(nrOfGrainSizes):
        #heights['totalSedContent_'+str(f)][i][f] = 0
        #heights['totalSedHeight'][i][f] = 0

@njit
def setTimestep(heights):
    #import pdb; pdb.set_trace()
    """ Calculating the maximum timestep (dt) as limited by the active layer and the FTCS scheme. """
    dtFTCS = 0.9*(dx*dx)/(2.e0*((k[0]+k[1]))) ##FTCS
    maxHeightDiff = heights['totalHeight'][0]+(q0[0]+q0[1])*dtFTCS/dx-heights['totalHeight'][1] ## dtFTCS is used here as an estimation on the high side. It is better to estimate more sediment input resulting in a high slope then it is to estimate on the low side. A too low estimate of the sediment input can lead to an amount of sediment being removed somewhere that is more the the depth of the active layer. Higher estimates of q0 result in lowe dt and thus a slower run speed.
    for i in range(1,imax):
        localSlope = (heights['totalHeight'][i-1] - heights['totalHeight'][i]) 
        if (localSlope > maxHeightDiff): maxHeightDiff = localSlope
        
    #Obtain minimum dt
    if (maxHeightDiff != 0): ## Cannot devide by 0
        dtHeightLimit = 0.9*(dx*dx  )/( (k[0]+k[1]) * maxHeightDiff )
        dt = dtFTCS 
        if (dtFTCS > dtHeightLimit): dt = dtHeightLimit 
    else:
        dt = 0.1*(dx*dx)/(2.e0*((k[0]+k[1]))) ##FTCS
    return dt

@njit
def sedContentInActiveLayer(nodes, totalSedContent, i, f):
    return totalSedContent[i][f]
    

@njit
def FTCS(dt, dx, q0, k, heights, totalSedContent, nodes, totalInput, totalOutput, InputPerGrainSize, OutputPerGrainSize, subsidenceRate):
    """ Calculate the new topography profile and keep track of the total displacement in- and output of sediment in the meantime. """
    sedIn = nbStackedLists((2, imax+1))
    sedOut = nbStackedLists((2, imax+1))
    newHeight = _nbZeros_1d(imax+1)
    
    totalSedIn = 0
    totalSedOut = 0
    heights['totalHeight'][imax] = 0
    
    ## Determine in- and output per column
    for i in range(imax+1): 
        
        ## The case i==0 is treated as a special case in this loop and i==imax+1 remains untouched so it stays 0
        for f in range(nrOfGrainSizes):
            #Dph= ( D(i+1)+D(i) )/2.e0
            #Dmh= ( D(i)+D(i-1) )/2.e0
            Dph = k[f]
            Dmh = k[f]
            #import pdb; pdb.set_trace()
            if (i==0): 
                #import pdb; pdb.set_trace()
                if (heights['totalHeight'][0] >= heights['totalHeight'][1]): ## Slope goes down to the right
                    sedIn[f][0] = q0[f]*dt/dx
                    totalSedIn -= q0[f]*dt/dx
                else: ## Slope goes down to the left
                    sedIn[f][0] += q0[f]*dt/dx
                    totalSedIn -= q0[f]*dt/dx
                    sedOut[f][0] = 0
            else:
                if (heights['totalHeight'][i-1] > heights['totalHeight'][i]): ## Slope goes down to the right
                    transport = ( (Dph*dt)/(dx*dx) )*( (heights['totalHeight'][i-1]) - ((heights['totalHeight'][i])) )
                    #transport = nbmin((transport, heights['totalSedContent_'+str(f)][i-1]))
                    transport = nbmin((transport, sedContentInActiveLayer(nodes, totalSedContent, i, f)))
                    sedOut[f][i-1] += transport
                    sedIn[f][i] += transport

                elif (heights['totalHeight'][i-1] < heights['totalHeight'][i]): ## Slope goes down to the left
                    transport = ( (Dph*dt)/(dx*dx) )*( (heights['totalHeight'][i]) - ((heights['totalHeight'][i-1])) )
                    transport = nbmin((transport, sedContentInActiveLayer(nodes, totalSedContent, i, f)))
                    sedOut[f][i] += transport
                    sedIn[f][i-1] += transport
    
    newSedContent = nbStackedLists((imax+1, nrOfGrainSizes))
    ## Set new height
    for i in range(imax+1): 
        for f in range(nrOfGrainSizes):
            newSedContent[i][f] = totalSedContent[i][f] + sedIn[f][i]-sedOut[f][i]
            heights['totalHeight'][i] += sedIn[f][i]-sedOut[f][i]
        if (heights['totalHeight'][i] < 0):
            heights['totalHeight'][i] = 0
        
    ## Set boundary condition at distal the end
    heights['totalHeight'][imax] = 0
    for f in range(nrOfGrainSizes):
        sedContentInActiveLayer(nodes, totalSedContent, i, f)
        
        ## Keep track of total in- and output 
        totalOutput += sedIn[f][imax]*dx
        OutputPerGrainSize[f] += sedIn[f][imax]*dx
        totalInput += q0[f]*dt
        InputPerGrainSize[f] += q0[f]*dt
        
    heights = subsidence (heights, imax, subsidenceRate, nrOfGrainSizes)
    return heights, nodes, InputPerGrainSize, OutputPerGrainSize, totalInput, totalOutput

@njit
def subsidence (heights, imax, subsidenceRate, nrOfGrainSizes):
    for i in range(imax+1):
        heights['bedrockHeight'][i] -= subsidenceRate*(imax-i)
        #columns[i]["oldHeight"] = nbmax( columns[i]["oldHeight"]-subsidenceRate*(imax-i) , heights['bedrockHeight'][i])
        heights['totalHeight'][i] = nbmax(( heights['totalHeight'][i]-subsidenceRate*(imax-i) , heights['bedrockHeight'][i] )) 
  
    return heights

## Make time in seconds readable:
def printElapsedTime(elapsedTime):
    if (elapsedTime < 60):
        if (elapsedTime < 1):
            print("This run took less then a second. You might want to take a look at the input, mate.")
        else:
            print("This run took "+str(math.ceil(elapsedTime))+" seconds")
    elif (elapsedTime < 60*60):
        timeInMinutes = math.floor(elapsedTime/60)
        if (timeInMinutes == 1):
            minuteText = " minute"
        else:
            minuteText = " minutes"
            
        if (math.ceil(elapsedTime-timeInMinutes*60) < 1):
            secondText = " second"
        else:
            secondText = " seconds"
        print("This run took "+str(timeInMinutes)+ minuteText +" and "+str(math.ceil(elapsedTime-timeInMinutes*60))+ secondText)
    else:
        timeInMinutes = math.floor(elapsedTime/60)
        hours = math.floor(timeInMinutes/60)
        
        if (hours == 1):
            hourText = " hour"
        else:
            hourText = " hours"
        minutes = int(elapsedTime/60 - hours*60)
        
        if (minutes == 1):
            minuteText = " minute"
        else:
            minuteText = " minutes"
        seconds = int(elapsedTime - hours*60*60 - minutes*60)
        
        if (seconds == 1):
            secondText = " second"
        else:
            secondText = " seconds"
        print("This run took "+str(hours)+ hourText+ ", " +str(minutes)+ minuteText +" and "+str(seconds)+ secondText)
        
        
def makeDirectories():
    if (not os.path.isdir("ISMolD_outputdata")):
        os.mkdir("ISMolD_outputdata")
        
    if (not os.path.isdir("ISMolD_outputdata/relief")):
        os.mkdir("ISMolD_outputdata/relief")

    if (not os.path.isdir("ISMolD_outputdata/nodes")):
        os.mkdir("ISMolD_outputdata/nodes")   

def writeOutput(heights):
    
    # makeTimeNodeDirectory(nodeOutputTimestep)
        
    # for i in range(len(x)-1): ## -1 for the last column is always empty (by design). Therefore there is no need to create a file for it.
        # f = open("ISMolD_outputdata/nodes/time"+str(nodeOutputTimestep)+"/column"+str(i)+".txt", "w")
        # jrange = len(columns[i]["nodes"])
        # for j in range(jrange): 
            # writeline = str(j)+" "+ str(columns[i]["heights['totalHeight']"])
            # writeline += " "+str(columns[i]["heights['bedrockHeight']"])
            # for p in range(nrOfGrainSizes):
                # writeline += " "+str(columns[i]["nodes"][j]["nodeSedContent"][p])
            # writeline += " "+ str(columns[i]["nodes"][j]["depositionTimeInYears"]) 
            # writeline += "\n"
            # f.write(writeline)
        # f.close()
    
    # nodeOutputTimestep += 1
    
    n= int(tout/dtout)
    myfile = open("ISMolD_outputdata/relief/topography"+str(n)+".txt", "w")
    for i in range(len(heights['totalHeight'])):
        writeline = str(i)+" "+str(heights['totalHeight'][i])
        writeline += " "+str(heights['bedrockHeight'][i])
        for f in range(nrOfGrainSizes):
            sedContent = heights['totalSedContent_'+str(f)][i]
            writeline += " "+str(sedContent + heights['bedrockHeight'][i]) ## Note that heights['bedrockHeight'] is generally negative
        writeline += "\n"
        myfile.write(writeline)
    myfile.close()
    
def setNodes(heights, nodes, newTotalHeight, newBedrockHeight):
    
    trace = True
    import pdb; pdb.set_trace()
    for i in range(len(newTotalHeight)):
        
        for f in range(nrOfGrainSizes):
            newSedContent = 1
        
        #if ( (newTotalHeight[i] - heights['totalHeight'][i] == 0): ## Nothing happens
            #return nodes
            
        ### Deposition ##
        #elif all( (newSedContent[q]-column["oldSedContent"][q])>=0 for q in range(nrOfGrainSizes) ): 
            #if (trace): print("ONLY DEPOSITION", i)
            
        ### Erosion ##
        #elif ( (newHeight - column["oldHeight"]) < 0 and all( (newSedContent[q]-column["oldSedContent"][q])<=0 for q in range(nrOfGrainSizes)) ):  ## Erosion
            #if (trace): print("ONLY EROSION", i)
        
        ### Both Deposition and Erosion ##
        #else: 
            #if (trace): print("Both Deposition and Erosion", i)
        
        
    #return nodes


if __name__=="__main__":
    
    print("Initializing...")
    makeDirectories()
    ## Remove old data files
    nr_topo_files = len(os.listdir("ISMolD_outputdata/relief")) #-1 since file starts at 'output0', add -1 for any subdirectory
    if (nr_topo_files > 1):
        os.system('rm ISMolD_outputdata/relief/topography*.txt')
        os.system('rm -r -f ISMolD_outputdata/nodes/*')
    ## Clear or create forcing file
    file = open("ISMolD_outputdata/forcing.txt", "w")
    file.write("")
    file.close()
    
    # import pdb; pdb.set_trace()
    
    t= 0
    while (t <= tmax):
        if (sum(k) == 0):
            print("\n\nFatal Error, there is no sediment transport at all this timestep. \nThis is not allowed, for dt approaches infinity if sum(k)=0. \nTime of error:", t/yr2sec, "yr")
            exit()
        
        for i in range(nrOfGrainSizes):
            if ( k[i] < 0):
                print("Error, negative diffusion. t=", str(t/yr2sec)+"yr", "k["+str(i)+"]="+str(k[i]))
                exit()
            if ( q0[i] < 0):
                print("Error, negative input. t=", str(t/yr2sec)+"yr", "q0["+str(i)+"]="+str(q0[i]))
                exit()
        
        dt = setTimestep(heights)
        
        #Call FTCS to obtain the new height profile
        newHeights, nodes, InputPerGrainSize, OutputPerGrainSize, totalInput, totalOutput = FTCS(dt, dx, q0, k, heights, totalSedContent, nodes, totalInput, totalOutput, InputPerGrainSize, OutputPerGrainSize, subsidenceRate)
        #setNodes(heights, newTotalHeight, newBedrockHeight, nodes)
        heights['totalHeight'] = newHeights['totalHeight']
        heights['bedrockHeight'] = newHeights['bedrockHeight']
        
        
        if (t >= tout):
            writeOutput(heights)
            tout += dtout
            # tout = writeOutput(imax, nrOfGrainSizes, heights['totalHeight'], heights['totalSedContent_'+str(f)], totalInput, InputPerGrainSize, OutputPerGrainSize, q0, sedOut, tout, dtout, outputTimestep)
            outputTimestep += 1
        t += dt
        
        if t > tout_progress:
            print("      "+str( (math.ceil(100000*t/tmax))/1000 )+"%", end="\r") ##Track progress
            tout_progress += dtout_progress
    
    # for i in range(5):
        # printNodes(i, nodes, heights['totalHeight'][i], heights['totalSedContent_'+str(f)][i])
    printElapsedTime(time()-start)
    










