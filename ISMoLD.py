## Program ISMoLD (Imperfect Sorting Model of Linear Diffusion)
## Designed and created by Timo Millenaar (tmillenaar@gmail.com)

import math
import os
import matplotlib.pyplot as plt
import numpy as np
from os import listdir

from functions_for_ISMoLD import *

plot = False
spikeTest = False

##Uncomment if not desired
plot = True 
#spikeTest = True

yr2sec = 60*60*24*365.25
dx= 1.e3       # lateral grid spacing (m)
dy = 1         # vertical grid spacing (m)  Important: MUST BE 1(as integer)!
imax= 101      # number of nodes
totalHeight= list(range(imax+1))
newHeight= list(range(imax+1))
x=list(range(imax+1))
sedIncrease = np.zeros(shape=(2,imax+1))
sedIn = np.zeros(shape=(2,imax+1))
sedOut = np.zeros(shape=(2,imax+1))
newSedContent = np.zeros(shape=(2,imax+1))
t= 0
dt= 1          # time step in (yr)
tmax= 10000*yr2sec     # max number of years in seconds (x*yr2sec)
dtout = tmax/100     # nr of years between write output
dtout_progress = 10*yr2sec     # nr of years between progress bar update
tout = 0.          # threshold to write next output in years (increased every write action by dtout)
tout_progress = 0.          # threshold to update progress bar
tprogress = 0.
dtprogress = tmax/1000. 
k= list(range(2))
k[0]= 3.2e-4       # Gravel diffusivity (m2/s)
k[1]= 4*3.2e-4       # Sand diffusivity (m2/s)
q0= list(range(2))
q0[0]= 2.e-6      # Gravel input (m2/s)
q0[1]= 1.e-6      # Sand input (m2/s)

rho0 = dy * 2700

totalInput= 0
totalOutput= 0
InputPerGrainSize = [0,0]
OutputPerGrainSize = [0,0]
subsidence_rate = 0.#4e-4*dt/3 #m/yr
maxVolumeLoss = 0


## Remove old data files
#     Note: "listdir" returns the number of files in the folder. ...
#           ... One file returns 1 but the names start at "output0", thus the -1:
outputfilesNr = len(listdir("ISMolD_outputdata/"))-1 
if (outputfilesNr >= 0):
    os.system('rm ISMolD_outputdata/topography*.txt')

## Initialize:
columns = {}
for i in range(imax+1):
    x[i]=i
    totalHeight[i] = 0
    columns[i]= {"totalHeight":0,
                 "TotalSedContent":list(range(2)),
                 "nodes":{},
                 }
    for p in range(2):
        columns[i]["TotalSedContent"][p] = 0
        newSedContent[p,i]= 0

## Spike test:
#  Note: the code cannot properly deal with intense spikes of 1 (except for i==0), thus the the spike is 2*dx wide.
if (spikeTest):
        
    totalInput+= 100*dx
    for p in range(2):
        InputPerGrainSize[p] += 50*dx
        q0[p] = 0
    columns[50]["totalHeight"] = 100
    columns[50]["TotalSedContent"][0] = 50
    columns[50]["TotalSedContent"][1] = 50
    for q in range(100):
        columns[50]["nodes"].update({q:{"density":rho0}})

    totalInput+= 100*dx
    for p in range(2):
        InputPerGrainSize[p] += 50*dx
    columns[51]["totalHeight"] = 100
    columns[51]["TotalSedContent"][0] = 50
    columns[51]["TotalSedContent"][1] = 50
    for q in range(100):
        columns[51]["nodes"].update({q:{"density":rho0}})


## ## ## ## ## ## ## ##
##  Main time loop:  ##
## ## ## ## ## ## ## ##
while (t < tmax):
    
    ## Varying sediment input or diffusivity through time can be set here:
    #q0[0] = max(0, 2.e-6*math.sin(3.1415 * 6* t/tmax))
    
    ## Calculating max dt as limited by the active layer and the FTCS scheme
    dtFTCS = 0.9*(dx*dx)/(2.e0*(sum(k))) ##FTCS
    heightDiff = columns[0]["totalHeight"]+sum(q0)*dt/dx-columns[1]["totalHeight"]
    for i in range(1,imax):
        heightDiff = max(heightDiff, abs(columns[i+1]["totalHeight"] - 2*columns[i]["totalHeight"] + columns[i-1]["totalHeight"]) )
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
        columns[i].update({"oldHeight":columns[i]["totalHeight"]}) ## Needed for the setNodes function at i==0. "totalHeight" would have been sufficient for the other nodes at the "old" height. The "oldHeight" is given to all nodes here to keep the setNodes function clear for a human to read.
    for p in range(2):
        columns[0]["TotalSedContent"][p] += q0[p]*dt/dx
        columns[0]["totalHeight"] += q0[p]*dt/dx
    
    ## The right end of basin always remains 0 for it is not altered by the FTCS scheme. Therefore there is no need to specify that boundary condition here.
    
    ## Keep track of volume balance:
    for p in range(2):
        totalInput+= q0[p]*dt
        InputPerGrainSize[p] += q0[p]*dt
        totalOutput+= ( (k[p]*dt)/(dx*dx) )*( columns[imax-1]["totalHeight"] - columns[imax]["totalHeight"] )*dx ## Volume leaving = height*dx at i=imax, which is set to 0
        OutputPerGrainSize[p] += ( (k[p]*dt)/(dx*dx) )*( columns[imax-1]["totalHeight"] - columns[imax]["totalHeight"] )*dx
    
    ## Loop through columns (for FTCS, density calculations, etc):
    for i in range(1,imax):
        newHeight[i]= columns[i]["totalHeight"] ##Start from current height
        for p in range(2):
            newSedContent[p,i] = columns[i]["TotalSedContent"][p] ##Start from current sediment content
            ## FTCS:
            #Dph= ( D(i+1)+D(i) )/2.e0
            #Dmh= ( D(i)+D(i-1) )/2.e0
            Dph = k[p]
            Dmh = k[p]
            if (columns[i-1]["totalHeight"] > columns[i+1]["totalHeight"]): ## Slope goes down to the right
                sedIn[p,i] = ( (Dph*dt)/(dx*dx) )*( columns[i-1]["totalHeight"] - columns[i]["totalHeight"] )
                sedIn[p,i] = min(sedIn[p,i], columns[i-1]["TotalSedContent"][p])
                sedOut[p,i] = ( (Dmh*dt)/(dx*dx) )*( columns[i]["totalHeight"] - columns[i+1]["totalHeight"] )
                sedOut[p,i] = min(sedOut[p,i], columns[i]["TotalSedContent"][p])
            else: ## Slope goes down to the left
                sedIn[p,i] = ( (Dph*dt)/(dx*dx) )*( columns[i+1]["totalHeight"] - columns[i]["totalHeight"] )
                sedIn[p,i] = min(sedIn[p,i], columns[i+1]["TotalSedContent"][p])
                sedOut[p,i] = ( (Dmh*dt)/(dx*dx) )*( columns[i]["totalHeight"] - columns[i-1]["totalHeight"] )
                sedOut[p,i] = min(sedOut[p,i], columns[i]["TotalSedContent"][p])
            
            newSedContent[p,i] += sedIn[p,i]-sedOut[p,i]
            newHeight[i] += sedIn[p,i]-sedOut[p,i]
            if (abs(sedIn[p,i]-sedOut[p,i]) > 1):
                print("active Layer > 1",i , dt, sedIn[p,i]-sedOut[p,i])
        
        columns[i] = setNodes(i, k, newHeight[i], columns[i], dt, dx, dy, rho0)  ## Update density, (todo grain size fractions and porosity)
        
    ## End column loop (i)
    
    ## Update in proximal boundary column:
    newHeight[0] = columns[0]["totalHeight"]
    for p in range(2):
        sedOut[p,0] = ( (k[p]*dt)/(dx*dx) )*( columns[0]["totalHeight"] - columns[1]["totalHeight"] )
        sedOut[p,0] = min(sedOut[p,0], columns[0]["TotalSedContent"][p]) 
        newSedContent[p,0] = columns[0]["TotalSedContent"][p] - sedOut[p,0]
        newHeight[0] -= sedOut[p,0]
        
    columns[0] = setNodes(0, k, newHeight[0], columns[0], dt, dx, dy, rho0)
    ## Overwrite old the topography with the new profile:
    for i in range(0,imax):
        columns[i]["totalHeight"] = newHeight[i]
        for p in range(2):
            columns[i]["TotalSedContent"][p] = newSedContent[p,i]
    
    
    if (t >= tout):
        n= int(tout/dtout)
        f = open("ISMolD_outputdata/topography"+str(n)+".txt", "w")
        for i in range(len(x)):
            f.write(str(x[i])+" "+str(columns[i]["totalHeight"])+" "+str(columns[i]["TotalSedContent"][0])+" "+str(columns[i]["TotalSedContent"][1])+"\n")
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
for i in range(imax+1):
    totalHeight[i]= columns[i]["totalHeight"]
    totalDepositVolume+= columns[i]["totalHeight"]*dx
    for p in range(2):
        DepositVolPerGrainSize[p] += columns[i]["TotalSedContent"][p]*dx
    for j in range(0,len(columns[i]["nodes"])):
        try:
            totalNodeVolume+= dx*dy * columns[i]["nodes"][j]["density"]/rho0
        except:
            pass

print("totalInput:", str(totalInput)+"m^3")
print("InputPerGrainSize:", str(InputPerGrainSize)+"m^3")
print("totalOutput:", str(totalOutput)+"m^3")
print("OutputPerGrainSize:", str(OutputPerGrainSize)+"m^3")
print("")
print("totalDepositVolume:", str(totalDepositVolume)+"m^3")
print("totalNodeVolume:", str(totalNodeVolume)+"m^3")
print("DepositVolPerGrainSize:", str(DepositVolPerGrainSize)+"m^3", "  sum:", sum(DepositVolPerGrainSize))
print("")
print("Volume error:", str(totalDepositVolume+totalOutput-totalInput)+"m^3")
print("Error in %:", str(100*(totalDepositVolume+totalOutput-totalInput)/totalInput)+"%")
print("")
for p in range(2):
    print("Volume error per grain size ("+str(p)+"):", str(DepositVolPerGrainSize[p]+OutputPerGrainSize[p]-InputPerGrainSize[p])+"m^3")
    print("Error in %:", str(100*(DepositVolPerGrainSize[p]+OutputPerGrainSize[p]-InputPerGrainSize[p])/InputPerGrainSize[p])+"%")
print("")
print("Node volume error:", str(totalNodeVolume+totalOutput-totalInput)+"m^3")
print("Node volume error in %:", str(100*(totalNodeVolume+totalOutput-totalInput)/totalInput)+"%")
print("")


if (plot):
    print("Plotting...")
    os.system("python3 animate_ISMoLD.py")


    





#columns[1]["nodes"].update({i:{"fractions":[0.3, 0.3, 0.4]}})

