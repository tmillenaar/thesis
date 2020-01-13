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

k = np.zeros(shape=(2))
q0 = np.zeros(shape=(2))

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

def nbZeros(*args):
    assert all([type(val).__name__ == 'int' for val in args]), 'Please only unput integers into nbZeros(), got {}'.format([type(val).__name__ for val in args])
    if len(args) > 1: 
        return _nbZeros(args)
    else: 
        return _nbZeros_1d(args[0])
    
    
# @njit
# def _nbZeros(shape):
    # ar1 = List()
    # ar2 = List()
    # for i in range(shape[-1]):
        # ar1.append(0.)
    # for i in range(len(shape)-2,-1,-1):
        # ar2 = 0
        # ar2 = ar1.copy()
        # ar1 = List()
        # for j in range(shape[i]):   
            # dummyarray = List()
            # for k in ar2:
                # dummyarray.append(k)
            # ar1.append(dummyarray)
            # dummyarray = 0
    # return ar1

# def zeros(*shape):
    # from copy import deepcopy
    # ar1 = []
    # ar2 = []
    # for i in range(shape[-1]):
        # ar1.append(0.)
    # for i in range(len(shape)-2,-1,-1):
        # ar2 = 0
        # ar2 = ar1.copy()
        # ar1 = []
        # for j in range(shape[i]):   
            # dummyarray = deepcopy(ar2)
            # ar1.append(dummyarray)
            # dummyarray = 0
    # return ar1

def _nbZeros(shape):
    l = List()
    if (len(shape) > 1):
        for i in range(shape[0]):
            l.append(_nbZeros(shape[1:]))
    else: 
        l = _nbZeros_1d(shape[0])
    return l

@njit
def _nbZeros_1d(length):
    l = List()
    for i in range(length):
        l.append(0.)
    return l

@njit
def _nbZeros_2d(shape):
    l = List()
    for i in range(shape[0]):
        l.append(_nbZeros_1d(shape[1]))
    return l
    
# nbZeros(3,2,2)
# nbZeros(3)
# import pdb; pdb.set_trace()
# def nbZeros(shape):
    # if not hasattr(shape, '__len__'):
        # result = _nbZeros_1d(shape)
    # else:
        # result = _nbZeros_nd(shape)
    # return result
      
  
# @njit
# def _nbZeros_nd(shape):
    # l = List()
    # for i in range(len(shape)):
        # k = List()
        # for j in range(shape[i]):
            # k.append(0.)
        # l.append(k)
    # return l

# @njit
# def _nbZeros_1d(length):
    # l = List()
    # for i in range(length):
        # l.append(0.)
    # return l
  
@njit
def nbRange(length):
    l = List()
    for i in range(length):
        l.append(i)
    return l

    
def nbStackedLists(*args):
    assert all([type(val).__name__ == 'int' for val in args]), 'Please only unput integers into nbStackedLists(), got {}'.format([type(val).__name__ for val in args])
    return _nbStackedLists(args)

@njit
def _nbStackedLists(shape):
    #import pdb; pdb.set_trace()
    l = List()
    for i in shape:
        k = List()
        for j in range(1,len(shape)):
            m = List()
            for j in range(shape[j]):
                m.append(0.)
            k.append(m)
        l.append(k)
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

totalSedContent = nbZeros(imax+1, nrOfGrainSizes)

# totalHeight = nbZeros(imax+1)
# bedrockHeight = nbZeros(imax+1)

totalHeight = np.zeros(shape=(imax+1))
bedrockHeight = np.zeros(shape=(imax+1))

nodeFill = nbZeros(imax+1, 1)
depTime = nbZeros(imax+1, 1)
nodeSedContent = nbZeros(imax+1, 1, 2)

totalInput= 0
totalOutput= 0
InputPerGrainSize = np.zeros(shape=(2))
OutputPerGrainSize = np.zeros(shape=(2))

dy = 1                            # vertical grid spacing (m)  Important: MUST BE 1(as integer) for code to work (in this version)!
tout = 0.                         # threshold to write next output in years (increased every write action by dtout)
tout_progress = 2*dtout_progress    # threshold to update progress bar, will be increased whenever it is reached
outputTimestep = 0

@njit
def setTimestep(totalHeight):
    """ Calculating the maximum timestep (dt) as limited by the active layer and the FTCS scheme. """
    dtFTCS = 0.9*(dx*dx)/(2.e0*((k[0]+k[1]))) ##FTCS
    maxHeightDiff = totalHeight[0]+(q0[0]+q0[1])*dtFTCS/dx-totalHeight[1] ## dtFTCS is used here as an estimation on the high side. It is better to estimate more sediment input resulting in a high slope then it is to estimate on the low side. A too low estimate of the sediment input can lead to an amount of sediment being removed somewhere that is more the the depth of the active layer. Higher estimates of q0 result in lowe dt and thus a slower run speed.
    for i in range(1,imax):
        localSlope = (totalHeight[i-1] - totalHeight[i]) 
        if (localSlope > maxHeightDiff): maxHeightDiff = localSlope
        
    ## Obtain minimum dt
    if (maxHeightDiff != 0): ## Cannot devide by 0
        dtHeightLimit = 0.9*(dx*dx  )/( (k[0]+k[1]) * maxHeightDiff )
        dt = dtFTCS 
        if (dtFTCS > dtHeightLimit): dt = dtHeightLimit 
    else:
        dt = 0.1*(dx*dx)/(2.e0*((k[0]+k[1]))) ##FTCS
    return dt

@njit
def sedContentInActiveLayer(nodeFill, totalSedContent, i, f):
    return totalSedContent[i][f]
    

@njit
def FTCS(totalHeight, bedrockHeight, totalSedContent, nodeFill, nodeSedContent):#dt, dx, q0, k, 
    """ Calculate the new topography profile and keep track of the total displacement in- and output of sediment in the meantime. """
    sedIn = _nbZeros_2d((imax+1, nrOfGrainSizes))
    sedOut = _nbZeros_2d((imax+1, nrOfGrainSizes))
    newHeight = np.zeros(shape=(imax+1))
    
    totalSedIn = 0
    totalSedOut = 0
    totalHeight[imax] = 0
    
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
                if (totalHeight[0] >= totalHeight[1]): ## Slope goes down to the right
                    sedIn[0][f] = q0[f]*dt/dx
                    # totalSedIn -= q0[f]*dt/dx
                else: ## Slope goes down to the left
                    sedIn[0][f] += q0[f]*dt/dx
                    # totalSedIn -= q0[f]*dt/dx
                    sedOut[0][f] = 0
            else:
                if (totalHeight[i-1] > totalHeight[i]): ## Slope goes down to the right
                    transport = ( (Dph*dt)/(dx*dx) )*( totalHeight[i-1] - (totalHeight[i]) )
                    #transport = nbmin((transport, heights['totalSedContent_'+str(f)][i-1]))
                    transport = nbmin((transport, sedContentInActiveLayer(nodeFill, totalSedContent, i-1, f)))
                    transport = min(transport, totalSedContent[i-1][f])
                    sedOut[i-1][f] += transport
                    sedIn[i][f] += transport

                elif (totalHeight[i-1] < totalHeight[i]): ## Slope goes down to the left
                    transport = ( (Dph*dt)/(dx*dx) )*( (totalHeight[i]) - (totalHeight[i-1]) )
                    transport = nbmin((transport, sedContentInActiveLayer(nodeFill, totalSedContent, i-1, f)))
                    sedOut[i][f] += transport
                    sedIn[i-1][f] += transport
    
    newSedContent = _nbZeros_2d((imax+1, nrOfGrainSizes))
    
    newTotalHeight = np.zeros(shape=(imax+1))
    newBedrockHeight = np.zeros(shape=(imax+1))
    
    ## Set new height
    for i in range(imax+1): 
        newTotalHeight[i] = totalHeight[i]
        for f in range(nrOfGrainSizes):
            newSedContent[i][f] = totalSedContent[i][f] + sedIn[i][f]-sedOut[i][f]
            newTotalHeight[i] += sedIn[i][f]-sedOut[i][f]
        if (newTotalHeight[i] < 0):
            newTotalHeight[i] = 0
    
    ## Set boundary condition at distal the end
    newTotalHeight[imax] = 0
    for f in range(nrOfGrainSizes):
        newSedContent[imax][f] = 0
        # sedContentInActiveLayer(nodeFill, totalSedContent, i, f)
        
        ## Keep track of total in- and output 
        # totalOutput += sedIn[f][imax]*dx
        # OutputPerGrainSize[f] += sedIn[f][imax]*dx
        # totalInput += q0[f]*dt
        # InputPerGrainSize[f] += q0[f]*dt
        
    # newHeights= subsidence (newHeights, imax, subsidenceRate, nrOfGrainSizes)
    # import pdb; pdb.set_trace()
    return newTotalHeight, newBedrockHeight, newSedContent

@njit
def subsidence (totalHeight, bedrockHeight):
    for i in range(imax+1):
        bedrockHeight[i] -= subsidenceRate*(imax-i)
        #columns[i]["oldHeight"] = nbmax( columns[i]["oldHeight"]-subsidenceRate*(imax-i) , bedrockHeight[i])
        totalHeight[i] = nbmax(( totalHeight[i]-subsidenceRate*(imax-i) , bedrockHeight[i] )) 
  
    return totalHeight, bedrockHeight

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

def writeOutput(totalHeight, bedrockHeight, totalSedContent):
    
    # makeTimeNodeDirectory(nodeOutputTimestep)
        
    # for i in range(len(x)-1): ## -1 for the last column is always empty (by design). Therefore there is no need to create a file for it.
        # f = open("ISMolD_outputdata/nodes/time"+str(nodeOutputTimestep)+"/column"+str(i)+".txt", "w")
        # jrange = len(columns[i]["nodes"])
        # for j in range(jrange): 
            # writeline = str(j)+" "+ str(columns[i]["totalHeight"])
            # writeline += " "+str(columns[i]["bedrockHeight"])
            # for p in range(nrOfGrainSizes):
                # writeline += " "+str(columns[i]["nodes"][j]["nodeSedContent"][p])
            # writeline += " "+ str(columns[i]["nodes"][j]["depositionTimeInYears"]) 
            # writeline += "\n"
            # f.write(writeline)
        # f.close()
    
    # nodeOutputTimestep += 1
    
    n= int(tout/dtout)
    myfile = open("ISMolD_outputdata/relief/topography"+str(n)+".txt", "w")
    for i in range(len(totalHeight)):
        writeline = str(i)+" "+str(totalHeight[i])
        writeline += " "+str(bedrockHeight[i])
        for f in range(nrOfGrainSizes):
            #sedContent = heights['totalSedContent_'+str(f)][i]
            writeline += " "+str(totalSedContent[i][f] + bedrockHeight[i]) ## Note that bedrockHeight is generally negative
        writeline += "\n"
        myfile.write(writeline)
    myfile.close()
    
# @njit
def setNodes(totalHeight, bedrockHeight, newHeights, totalSedContent, newSedContent, nodes):
    
    trace = True
    #import pdb; pdb.set_trace()
    for i in range(len(newTotalHeight)):
        
        #for f in range(nrOfGrainSizes):
            #newSedContent = 1
        
        if ( (newTotalHeight[i] - totalHeight[i]) == 0): ## Nothing happens
            return nodes
        
        onlyDeposition = True
        onlyErosion = True
        for f in range(nrOfGrainSizes):
            if (newSedContent[i][f] < totalSedContent[i][f]): onlyDeposition = False
            if (newSedContent[i][f] > totalSedContent[i][f]): onlyErosion = False
            
        ## Deposition ##
        if (onlyDeposition):
            if (trace): print("ONLY DEPOSITION", i)
            flowFractions = np.zeros(shape=nrOfGrainSizes)
            # sedContentChange = nbStackedLists((imax+1, nrOfGrainSizes))
            sedContentChange = np.zeros(shape=(imax+1, nrOfGrainSizes))
            
            for f in range (nrOfGrainSizes):
                sedContentChange[i][f] = newSedContent[i][f] - totalSedContent[i][f]
                
            for f in range (nrOfGrainSizes):    
                if (np.sum(sedContentChange) != 0):
                    flowFractions[f] = sedContentChange[i][f]/np.sum(sedContentChange[i])
                else:
                    flowFractions[f] = 0
            import pdb; pdb.set_trace()
            toBeDeposited = np.sum(sedContentChange)
            while (toBeDeposited > 0):
                currentNodeDensity = nodes['nodeFill'][i][-1]
                currentNodeSedContent = np.zeros(shape=nrOfGrainSizes)
                for f in range(nrOfGrainSizes):
                    currentNodeSedContent[f] = nodes['sedContent'][i][-1][f]
                
        ## Erosion ##
        elif (onlyErosion):  ## Erosion
            if (trace): print("ONLY EROSION", i)
        
        ## Both Deposition and Erosion ##
        else: 
            if (trace): print("Both Deposition and Erosion", i)
        
        
    return nodes


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
        
        dt = setTimestep(totalHeight)
        
        #Call FTCS to obtain the new height profile
        newTotalHeight, newBedrockHeight, newTotalSedContent = FTCS(totalHeight, bedrockHeight, totalSedContent, nodeFill, nodeSedContent)
        # import pdb; pdb.set_trace()
        newTotalHeight, newBedrockHeight = subsidence(newTotalHeight, newBedrockHeight)
        # nodes = setNodes(heights, newHeights, totalSedContent, newTotalSedContent,  nodes)
        
        totalHeight = newTotalHeight
        bedrockHeight = newBedrockHeight
        totalSedContent = newTotalSedContent
        
        
        if (t >= tout):
            writeOutput(totalHeight, bedrockHeight, totalSedContent)
            tout += dtout
            # tout = writeOutput(imax, nrOfGrainSizes, totalHeight, heights['totalSedContent_'+str(f)], totalInput, InputPerGrainSize, OutputPerGrainSize, q0, sedOut, tout, dtout, outputTimestep)
            outputTimestep += 1
        t += dt
        
        if t > tout_progress:
            print("      "+str( (math.ceil(100000*t/tmax))/1000 )+"%", end="\r") ##Track progress
            tout_progress += dtout_progress
    
    # for i in range(5):
        # printNodes(i, nodes, totalHeight[i], heights['totalSedContent_'+str(f)][i])
    printElapsedTime(time()-start)
    










