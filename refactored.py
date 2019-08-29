import math
from numba import jit, njit
from time import time
import numpy as np
import os

start = time()

yr2sec = 60*60*24*365.25          #nr of seconds in a year

dx              = 1e3             # width of each node/column (m)
imax            = 100             # number of nodes
tmax            = 500000*yr2sec    # total amount of time to be modelled in seconds [in years = (x*yr2sec)]
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
    
k[0]= 0.1*3.2e-4                    # Gravel diffusivity (m2/s)  Attention: will be overwritten in setBoudnaryCondtitionValues if declared there!
k[1]= 0.2*3.2e-4                    # Sand diffusivity (m2/s)  Attention: will be overwritten in setBoudnaryCondtitionValues if declared there!

x               = np.zeros(shape=(imax+1))
totalHeight     = np.zeros(shape=(imax+1))
bedrockHeight   = np.zeros(shape=(imax+1))
totalSedContent = np.zeros(shape=(2, imax+1))
totalSedHeight  = np.zeros(shape=(2, imax+1))
nodes         = np.zeros(shape=(nrOfGrainSizes+1,5,10)) # For the first index: 0=density, 1=first grain size, 2=second grain size, etc..  The second index indicates the row and the third indicates the column

totalInput= 0
totalOutput= 0
InputPerGrainSize = np.zeros(shape=(2))
OutputPerGrainSize = np.zeros(shape=(2))

dy = 1                            # vertical grid spacing (m)  Important: MUST BE 1(as integer) for code to work (in this version)!
tout = 0.                         # threshold to write next output in years (increased every write action by dtout)
tout_progress = 0.                # threshold to update progress bar

for i in range(imax+1):
    x[i]=i
    totalHeight[i] = 0
    #newHeight[i] = 0
    bedrockHeight[i] = 0
    #columns[i] = []
    for f in range(nrOfGrainSizes):
        totalSedContent[f][i] = 0
        totalSedHeight[f][i] = 0

@njit
def setTimestep():
    """ Calculating the macimum timestep (dt) as limited by the active layer and the FTCS scheme. """
    dtFTCS = 0.9*(dx*dx)/(2.e0*((k[0]+k[1]))) ##FTCS
    maxHeightDiff = totalHeight[0]+(q0[0]+q0[1])*dtFTCS/dx-totalHeight[1] ## dtFTCS is used here as an estimation on the high side. It is better to estimate more sediment input resulting in a high slope then it is to estimate on the low side. A too low estimate of the sediment input can lead to an amount of sediment being removed somewhere that is more the the depth of the active layer. Higher estimates of q0 result in lowe dt and thus a slower run speed.
    for i in range(1,imax):
        maxHeightDiff = max(maxHeightDiff, abs(totalHeight[i-1] - totalHeight[i]) )
        
    if (maxHeightDiff != 0): ## Cannot devide by 0
        dtHeightLimit = 0.9*(dx*dx)/( (k[0]+k[1]) * maxHeightDiff )
        dt = min(dtFTCS, dtHeightLimit)
    else:
        dt = 0.1*(dx*dx)/(2.e0*((k[0]+k[1]))) ##FTCS
    return dt

@njit
def FTCS(dt, dx, q0, k, totalSedContent, totalInput, totalOutput, InputPerGrainSize, OutputPerGrainSize, totalHeight):
    """ Calculate the new topography profile and keep track of the total displacement in- and output of sediment in the meantime. """
    sedIn = np.zeros(shape=(2,imax+1))
    sedOut = np.zeros(shape=(2,imax+1))
    newHeight = np.zeros(shape=(imax+1))
    ## Loop through columns (for FTCS, density calculations, etc):
    totalSedIn = 0
    totalSedOut = 0
    totalHeight[imax] = 0
    
    for i in range(imax+1): 
        ## The case i==0 is treated as a special case in this loop and i==imax+1 remains untouched so it stays 0
        
        for f in range(nrOfGrainSizes):
            #Dph= ( D(i+1)+D(i) )/2.e0
            #Dmh= ( D(i)+D(i-1) )/2.e0
            Dph = k[f]
            Dmh = k[f]
            
            if (i==0): 
                if (totalHeight[0] >= totalHeight[1]): ## Slope goes down to the right
                    sedIn[f,0] = q0[f]*dt/dx
                    totalSedIn -= q0[f]*dt/dx
                else: ## Slope goes down to the left
                    sedIn[f,0] += q0[f]*dt/dx
                    totalSedIn -= q0[f]*dt/dx
                    sedOut[f,0] = 0
            else:
                
                if (totalHeight[i-1] > totalHeight[i]): ## Slope goes down to the right
                    transport = ( (Dph*dt)/(dx*dx) )*( (totalHeight[i-1]) - ((totalHeight[i])) )
                    #transport = min(transport, sedContentInActiveLayer[f,i-1])
                    transport = min(transport, totalSedContent[f][i-1])
                    sedOut[f,i-1] += transport
                    sedIn[f,i] += transport

                elif (totalHeight[i-1] < totalHeight[i]): ## Slope goes down to the left
                    transport = ( (Dph*dt)/(dx*dx) )*( (totalHeight[i]) - ((totalHeight[i-1])) )
                    #transport = min(transport, sedContentInActiveLayer[f,i])
                    transport = min(transport, totalSedContent[f][i])
                    sedOut[f,i] += transport
                    sedIn[f,i-1] += transport
            
    for i in range(imax+1): 
        for f in range(nrOfGrainSizes):
            totalSedContent[f][i] = totalSedContent[f][i] + sedIn[f][i]-sedOut[f][i]
            totalHeight[i] += sedIn[f][i]-sedOut[f][i]
        if (totalHeight[i] < 0):
            totalHeight[i] = 0
            
    totalHeight[imax] = 0
    for f in range(nrOfGrainSizes):
        totalSedContent[f][imax] = 0
        totalOutput += sedIn[f,imax]*dx
        OutputPerGrainSize[f] += sedIn[f,imax]*dx
        totalInput += q0[f]*dt
        InputPerGrainSize[f] += q0[f]*dt
                        
        #for i in range(imax+1): 
            #for f in range(nrOfGrainSizes):            
                #totalSedIn += sedIn[f,i]
                #totalSedOut += sedOut[f,i]
    return totalSedContent, InputPerGrainSize, OutputPerGrainSize, totalInput, totalOutput, totalHeight, sedOut#, totalSedIn, totalSedOut

def writeOutput(nrOfGrainSizes, totalHeight, totalSedContent, totalInput, InputPerGrainSize, OutputPerGrainSize, q0, sedOut, tout, dtout):
    
    #makeTimeNodeDirectory(nodeOutputTimestep)
        
    #for i in range(len(x)-1): ## -1 for the last column is always empty (by design). Therefore there is no need to create a file for it.
        #file = open("ISMolD_outputdata/nodes/time"+str(nodeOutputTimestep)+"/column"+str(i)+".txt", "w")
        #jrange = len(columns[i]["nodes"])
        #for j in range(jrange): 
            #writeline = str(j)+" "+ str(columns[i]["totalHeight"])
            #writeline += " "+str(columns[i]["bedrockHeight"])
            #for p in range(nrOfGrainSizes):
                #writeline += " "+str(columns[i]["nodes"][j]["nodeSedContent"][p])
            #writeline += " "+ str(columns[i]["nodes"][j]["depositionTimeInYears"]) 
            #writeline += "\n"
            #file.write(writeline)
        #file.close()
    
    #nodeOutputTimestep += 1
    
    n= int(tout/dtout)
    file = open("ISMolD_outputdata/relief/topography"+str(n)+".txt", "w")
    for i in range(len(x)):
        writeline = str(x[i])+" "+str(totalHeight[i])
        #writeline += " "+str(columns[i]["bedrockHeight"])
        writeline += " "+str(0)
        for f in range(nrOfGrainSizes):
            #writeline += " "+str(columns[i]["totalSedContent"][p] + columns[i]["bedrockHeight"]) ## Note that bedrockHeight is generally negative
            writeline += " "+str(totalSedContent[f][i])
        writeline += "\n"
        file.write(writeline)
    file.close()
    
    ## Write out forcing
    file = open("ISMolD_outputdata/forcing.txt", "a")
    writeline = str(t)+" "+str(totalInput)
    for p in range(nrOfGrainSizes):
        writeline += " "+str(InputPerGrainSize[p])
    writeline += " "+str(sum(q0))
    for p in range(nrOfGrainSizes):
        writeline += " "+str(q0[p])
    writeline += " "+str(totalOutput)
    for p in range(nrOfGrainSizes):
        writeline += " "+str(OutputPerGrainSize[p])
    writeline += " "+str(sum(sedOut[:,imax-1])/dt)
    for p in range(nrOfGrainSizes):
        writeline += " "+str(sedOut[p,imax-1]/dt)
    writeline += " "+str(sum(k))
    for p in range(nrOfGrainSizes):
        writeline += " "+str(k[p])  
    #writeline += " "+str(subsidenceRate)
    writeline += "\n"
    file.write(writeline)
    file.close()
    
    
    tout = tout + dtout  #make sure the next output is written in another "dtout" years
    return tout










#def timeloop(tmax, yr2sec, nrOfGrainSizes, totalSedContent, k, q0, tout_progress):
if __name__=="__main__":
    
    ## Remove old data files
    nr_topo_files = len(os.listdir("ISMolD_outputdata/relief")) #-1 since file starts at 'output0', add -1 for amy subdirectory
    if (nr_topo_files > 1):
        os.system('rm ISMolD_outputdata/relief/topography*.txt')
        os.system('rm -r -f ISMolD_outputdata/nodes/*')
    ## Clear or create forcing file
    file = open("ISMolD_outputdata/forcing.txt", "w")
    file.write("")
    file.close()
    
    t= 0
    while (t <= tmax):
        
        #if (sum(k) == 0):
            #print("\n\nFatal Error, there is no sediment transport at all this timestep. \nThis is not allowed, for dt approaches infinity if sum(k)=0. \nTime of error:", t/yr2sec, "yr")
            #exit()
        
        #for i in range(nrOfGrainSizes):
            #if ( k[i] < 0):
                #print("Error, negative diffusion. t=", str(t/yr2sec)+"yr", "k["+str(i)+"]="+str(k[i]))
                #exit()
            #if ( q0[i] < 0):
                #print("Error, negative input. t=", str(t/yr2sec)+"yr", "q0["+str(i)+"]="+str(q0[i]))
                #exit()
        
        dt = setTimestep()
        
        #Call FTCS to obtain the new height profile
        totalSedContent, InputPerGrainSize, OutputPerGrainSize, totalInput, totalOutput, totalHeight, sedOut = FTCS(dt, dx, q0, k, totalSedContent, totalInput, totalOutput, InputPerGrainSize, OutputPerGrainSize, totalHeight)
        
        if (t >= tout):
            tout = writeOutput(nrOfGrainSizes, totalHeight, totalSedContent, totalInput, InputPerGrainSize, OutputPerGrainSize, q0, sedOut, tout, dtout)
            
        t += dt
        
        if t > tout_progress:
            print("      "+str( (math.ceil(100000*t/tmax))/1000 )+"%", end="\r") ##Track progress
            tout_progress += dtout_progress
    
    
#jit()(timeloop(tmax, yr2sec, nrOfGrainSizes, k, q0, tout_progress, dtout_progress))
#timeloop(tmax, yr2sec, nrOfGrainSizes, totalSedContent, k, q0, tout_progress)
    
    
    
end = time()
runTime = end-start
print("This run took {} seconds".format(runTime))
    
    
    
    
