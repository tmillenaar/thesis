## Program ISMoLD (Imperfect Sorting Model of Linear Diffusion)
## Designed and created by Timo Millenaar (tmillenaar@gmail.com)

import math
import os
import pdb ## pdb.set_trace()
import matplotlib.pyplot as plt
import numpy as np
from os import listdir

from functions_for_ISMoLD import *

plot = False
spikeTest = False
trace = False

##Uncomment if not desired:
plot = True 
#spikeTest = True
#trace = True

yr2sec = 60*60*24*365.25    #nr of seconds in a year

dx= 1.e3                    # lateral grid spacing (m)
imax= 101                   # number of nodes
tmax= 10000*yr2sec          # total amount of time to be modelled in seconds [in years = (x*yr2sec)]
dtout = tmax/100            # nr of years between write output
dtout_progress = 10*yr2sec  # nr of years between progress bar update
nrOfGrainSizes = 2
k= list(range(nrOfGrainSizes))
k[0]= 3.2e-4                # Gravel diffusivity (m2/s)
k[1]= 3*3.2e-4              # Sand diffusivity (m2/s)
q0= list(range(nrOfGrainSizes))
q0[0]= 1*1.e-6                # Gravel input (m2/s)
q0[1]= 1.e-6                # Sand input (m2/s)

rho0 = 2700

subsidence_rate = 0.#4e-4*dt/3 #m/yr

## Initialize:
inputFractions = list(range(nrOfGrainSizes))
totalHeight= list(range(imax+1))
newHeight= list(range(imax+1))
x= list(range(imax+1))
sedContentInActiveLayer= np.zeros(shape=(2,imax+1))
sedIncrease = np.zeros(shape=(2,imax+1))
sedIn = np.zeros(shape=(2,imax+1))
sedOut = np.zeros(shape=(2,imax+1))
newSedContent = np.zeros(shape=(2,imax+1))

dy = 1         # vertical grid spacing (m)  Important: MUST BE 1(as integer)!

tout = 0.          # threshold to write next output in years (increased every write action by dtout)
tout_progress = 0.          # threshold to update progress bar
tprogress = 0.

totalInput= 0
totalOutput= 0
InputPerGrainSize = [0,0]
OutputPerGrainSize = [0,0]
maxVolumeLoss = 0

columns = {}
for i in range(imax+1):
    x[i]=i
    totalHeight[i] = 0
    columns[i]= {"totalHeight":0,
                 "totalSedContent":list(range(nrOfGrainSizes)),
                 "nodes":{},
                 }
    for p in range(nrOfGrainSizes):
        columns[i]["totalSedContent"][p] = 0
        newSedContent[p,i]= 0
        

## Remove old data files
#     Note: "listdir" returns the number of files in the folder. ...
#           ... One file returns 1 but the names start at "output0", thus the -1:
outputfilesNr = len(listdir("ISMolD_outputdata/"))-1 
if (outputfilesNr >= 0):
    os.system('rm ISMolD_outputdata/topography*.txt')

## Spike test:
#  Note: the code cannot properly deal with intense spikes of 1 (except for i==0), thus the the spike is 2*dx wide.
if (spikeTest):
    spikeHeight = 100
    spikeLocation = 40
    totalInput+= spikeHeight*dx
    for p in range(nrOfGrainSizes):
        InputPerGrainSize[p] += dx*spikeHeight/2
        q0[p] = 0
    columns[spikeLocation]["totalHeight"] = spikeHeight
    columns[spikeLocation]["totalSedContent"][0] = spikeHeight/2
    columns[spikeLocation]["totalSedContent"][1] = spikeHeight/2
    for q in range(spikeHeight):
        columns[spikeLocation]["nodes"].update({q:{"density":rho0,
                                                   "nodeSedContent":[0.5,0.5]}})

    totalInput+= spikeHeight*dx
    for p in range(nrOfGrainSizes):
        InputPerGrainSize[p] += dx*spikeHeight/2
    columns[spikeLocation+1]["totalHeight"] = spikeHeight
    columns[spikeLocation+1]["totalSedContent"][0] = spikeHeight/2
    columns[spikeLocation+1]["totalSedContent"][1] = spikeHeight/2
    for q in range(spikeHeight):
        columns[spikeLocation+1]["nodes"].update({q:{"density":rho0,
                                                     "nodeSedContent":[0.5,0.5]}})
    #print(columns[spikeLocation])

## ## ## ## ## ## ## ##
##  Main time loop:  ##
## ## ## ## ## ## ## ##
t= 0
while (t < tmax):
    
    ## Varying sediment input or diffusivity through time can be set here:
    #q0[0] = max(0, 2.e-6*math.sin(3.1415 * 6* t/tmax))
    
    ## Calculating max dt as limited by the active layer and the FTCS scheme
    dtFTCS = 0.9*(dx*dx)/(2.e0*(sum(k))) ##FTCS
    heightDiff = columns[0]["totalHeight"]+sum(q0)*dtFTCS/dx-columns[1]["totalHeight"] ## dtFTCS is used here as an estimation on the high side. It is better to estimate more sediment input resulting in a high slope then it is to estimate on the low side. A too low estimate of the sediment input can lead to an amount of sediment being removed somewhere that is more the the depth of the active layer. Higher estimates of q0 result in lowe dt and thus a slower run speed.
    for i in range(1,imax):
        #heightDiff = max(heightDiff, abs(columns[i+1]["totalHeight"] - 2*columns[i]["totalHeight"] + columns[i-1]["totalHeight"]) )
        heightDiff = max(heightDiff, abs(columns[i-1]["totalHeight"] - columns[i]["totalHeight"]) )
    dtHeightLimit = 0
    if (heightDiff != 0): ## Cannot devide by 0
        dtHeightLimit += 0.9*(dx*dx)/( sum(k) * heightDiff )
        dt = min(dtFTCS, dtHeightLimit)
    else:
        dt = 0.1*(dx*dx)/(2.e0*(sum(k))) ##FTCS
    
    if t > tout_progress:
        print("      "+str( (math.ceil(100000*t/tmax))/1000 )+"%", end="\r") ##Track progress
        tout_progress += dtout_progress
    
    ## Set boundary condition at proximal end of the basin:
    for i in range(imax+1):
        ## "oldHeight" and "oldSedContent" are needed for the setNodes function at i==0. "totalHeight" would have been sufficient for the other nodes as the "old" height. The "oldHeight" is given to all nodes here to keep the setNodes function clear for a human to read. Similarly for "oldSedContent":
        columns[i]["oldHeight"] = columns[i]["totalHeight"]
        columns[i]["oldSedContent"] = columns[i]["totalSedContent"].copy()
        
    ## Keep track of volume balance:
    for p in range(nrOfGrainSizes):
        totalInput+= q0[p]*dt
        InputPerGrainSize[p] += q0[p]*dt
        totalOutput+= ( (k[p]*dt)/(dx*dx) )*( columns[imax-1]["totalHeight"] - columns[imax]["totalHeight"] )*dx ## Volume leaving = height*dx at i=imax, which is set to 0
        OutputPerGrainSize[p] += ( (k[p]*dt)/(dx*dx) )*( columns[imax-1]["totalHeight"] - columns[imax]["totalHeight"] )*dx
    
    ## Set active layer properties
    for i in range(imax+1):
        maxNode = len(columns[i]["nodes"]) - 1 ## -1 since node count starts at 0 and len() starts at 1
        if (maxNode > 0):
            for p in range(nrOfGrainSizes):
                sedContentInActiveLayer[p,i] = columns[i]["nodes"][maxNode]["nodeSedContent"][p]
            if (sum(sedContentInActiveLayer[:,i]) < 1): ## which is always the case so long as dy==active layer depth
                remainder = 1-sum(sedContentInActiveLayer[:,i])
                for p in range(nrOfGrainSizes):
                    sedContentInActiveLayer[p,i] += columns[i]["nodes"][maxNode-1]["nodeSedContent"][p]*remainder 
        else:
            for p in range(nrOfGrainSizes):
                try:
                    sedContentInActiveLayer[p,i] = columns[i]["nodes"][0]["nodeSedContent"][p]
                except:
                    sedContentInActiveLayer[p,i] = 0
    #print("sedContentInActiveLayer 0",sedContentInActiveLayer[:,0], columns[0]["totalSedContent"])
    
        
    ## Loop through columns (for FTCS, density calculations, etc):
    for i in range(imax): 
        ## i==0 is treated as a special case in this loop and i==imax+1 remains untouched, having it stay 0
        
        newHeight[i]= columns[i]["totalHeight"] ##Start from current height
        for p in range(nrOfGrainSizes):
            
            newSedContent[p,i] = columns[i]["totalSedContent"][p] ##Start from current sediment content
            ## FTCS:
            #Dph= ( D(i+1)+D(i) )/2.e0
            #Dmh= ( D(i)+D(i-1) )/2.e0
            Dph = k[p]
            Dmh = k[p]
            
            if (i==0): 
                if (columns[0]["totalHeight"] > columns[1]["totalHeight"]): ## Slope goes down to the right
                    sedIn[p,0] = q0[p]*dt/dx
                    sedOut[p,0] = ( (k[p]*dt)/(dx*dx) )*( columns[0]["totalHeight"] - columns[1]["totalHeight"] )
                    sedOut[p,0] = min(sedOut[p,0], sedContentInActiveLayer[p,0])
                else: ## Slope goes down to the left
                    sedIn[p,0] = ( (k[p]*dt)/(dx*dx) )*( columns[1]["totalHeight"] - columns[0]["totalHeight"] )
                    sedIn[p,0] = min(sedIn[p,0], sedContentInActiveLayer[p,1])
                    sedIn[p,0] += q0[p]*dt/dx
                    sedOut[p,0] = 0
            else: 
                if (columns[i-1]["totalHeight"] > columns[i+1]["totalHeight"]): ## Slope goes down to the right
                    sedIn[p,i] = ( (Dph*dt)/(dx*dx) )*( columns[i-1]["totalHeight"] - columns[i]["totalHeight"] )
                    sedIn[p,i] = min(sedIn[p,i], sedContentInActiveLayer[p,i-1])
                    sedOut[p,i] = ( (Dmh*dt)/(dx*dx) )*( columns[i]["totalHeight"] - columns[i+1]["totalHeight"] )
                    sedOut[p,i] = min(sedOut[p,i], sedContentInActiveLayer[p,i])
                else: ## Slope goes down to the left
                    sedIn[p,i] = ( (Dph*dt)/(dx*dx) )*( columns[i+1]["totalHeight"] - columns[i]["totalHeight"] )
                    sedIn[p,i] = min(sedIn[p,i], sedContentInActiveLayer[p,i+1])
                    sedOut[p,i] = ( (Dmh*dt)/(dx*dx) )*( columns[i]["totalHeight"] - columns[i-1]["totalHeight"] )
                    sedOut[p,i] = min(sedOut[p,i], sedContentInActiveLayer[p,i])
            
            newSedContent[p,i] += sedIn[p,i]-sedOut[p,i]
            newHeight[i] += sedIn[p,i]-sedOut[p,i]
            if (abs(sedIn[p,i]-sedOut[p,i]) > 1):
                print("active Layer > 1",i , dt, sedIn[p,i]-sedOut[p,i])
        
        ## Update the nodes to suite newHeight and newSedContent:
        columns[i] = setNodes(i, k, newHeight[i], columns[i], newSedContent[:,i], dt, dx, dy, rho0)  
        
    ## End column loop (i)
    
    ## Overwrite old the topography with the new profile:
    for i in range(0,imax):
        columns[i]["totalHeight"] = newHeight[i]
        for p in range(nrOfGrainSizes):
            columns[i]["totalSedContent"][p] = newSedContent[p,i]
    
    if (t >= tout):
        n= int(tout/dtout)
        f = open("ISMolD_outputdata/topography"+str(n)+".txt", "w")
        for i in range(len(x)):
            f.write(str(x[i])+" "+str(columns[i]["totalHeight"])+" "+str(columns[i]["totalSedContent"][0])+" "+str(columns[i]["totalSedContent"][1])+"\n")
        f.close()
        tout = tout + dtout
        
    t += dt
    
## End time loop

## ## ## ## ## ## ## ##
##    Wrapping up:   ##
## ## ## ## ## ## ## ##

totalDepositVolume= 0
totalNodeVolume= 0
DepositVolPerGrainSize = [0,0]
totalGrainSizeFraction = [0,0]
for i in range(imax+1):
    totalHeight[i]= columns[i]["totalHeight"]
    totalDepositVolume+= columns[i]["totalHeight"]*dx
    for p in range(nrOfGrainSizes):
        DepositVolPerGrainSize[p] += columns[i]["totalSedContent"][p]*dx
    for j in range(0,len(columns[i]["nodes"])+1):
        try:
            totalNodeVolume+= dx*dy * columns[i]["nodes"][j]["density"]/rho0
        except:
            pass
        for p in range(nrOfGrainSizes):
            try:
                totalGrainSizeFraction[p] = totalGrainSizeFraction[p] + columns[i]["nodes"][j]["nodeSedContent"][p]
            except:
                pass

totalGrainSizeFraction_sum = sum(totalGrainSizeFraction)
for p in range(nrOfGrainSizes):
    totalGrainSizeFraction[p] = totalGrainSizeFraction[p]/totalGrainSizeFraction_sum 

if (trace):
    pdb.set_trace()
    
#print(columns[0]["nodes"])
#print("")

print("totalInput:", str(totalInput)+"m^3")
print("InputPerGrainSize:", str(InputPerGrainSize)+"m^3")
print("totalOutput:", str(totalOutput)+"m^3")
print("OutputPerGrainSize:", str(OutputPerGrainSize)+"m^3")
print("")
print("totalDepositVolume:", str(totalDepositVolume)+"m^3")
print("totalNodeVolume:", str(totalNodeVolume)+"m^3")
print("DepositVolPerGrainSize:", str(DepositVolPerGrainSize)+"m^3", "  sum:", sum(DepositVolPerGrainSize))
print("totalGrainSizeFraction", totalGrainSizeFraction)
print("")
print("Volume error:", str(totalDepositVolume+totalOutput-totalInput)+"m^3")
print("Error in %:", str(100*(totalDepositVolume+totalOutput-totalInput)/totalInput)+"%")
print("")
for p in range(nrOfGrainSizes):
    print("Volume error per grain size ("+str(p)+"):", str(DepositVolPerGrainSize[p]+OutputPerGrainSize[p]-InputPerGrainSize[p])+"m^3")
    print("Error in %:", str(100*(DepositVolPerGrainSize[p]+OutputPerGrainSize[p]-InputPerGrainSize[p])/InputPerGrainSize[p])+"%")
print("")
print("Node volume error:", str(totalNodeVolume+totalOutput-totalInput)+"m^3")
print("Node volume error in %:", str(100*(totalNodeVolume+totalOutput-totalInput)/totalInput)+"%")
print("")


if (plot):
    print("Plotting...")
    os.system("python3 animate_ISMoLD.py")





