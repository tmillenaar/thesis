## Program ISMoLD (Imperfect Sorting Model of Linear Diffusion)
## Designed and created by Timo Millenaar (tmillenaar@gmail.com)

import math
import os
import matplotlib.pyplot as plt
import numpy as np
from os import listdir

from functions_for_ISMoLD import *

outputfilesNr = len(listdir("ISMolD_outputdata/"))-1 #-1 since file starts at 'output0'
if (outputfilesNr >= 0):
    os.system('rm ISMolD_outputdata/topography*.txt') ##remove old data files

yr2sec = 60*60*24*365.25
dx= 1.e3       # lateral grid spacing (m)
dy = 1         # vertical grid spacing (m)  Important: MUST BE 1(as integer)!
imax= 101      # number of nodes
totalHeight= list(range(imax+1))
newHeight= list(range(imax+1))
x=list(range(imax+1))
sedIncrease = np.zeros(shape=(2,imax+1))
newSedContent = np.zeros(shape=(2,imax+1))
t= 0
dt= 1          # time step in (yr)
tmax= 10000*yr2sec     # max number of years
dtout = tmax/100     # nr of years between write output
dtout_progress = 10*yr2sec     # nr of years between progress bar update
tout = 0.          # threshold to write next output in years (increased every write action by dtout)
tout_progress = 0.          # threshold to update progress bar
tprogress = 0.
dtprogress = tmax/1000. 
k= list(range(2))
k[0]= 5*3.2e-5       # Gravel diffusivity (m2/s)
k[1]= 3.2e-4       # Sand diffusivity (m2/s)
q0= list(range(2))
q0[0]= 4.e-6      # Gravel input (m2/s)
q0[1]= 2.e-6      # Sand input (m2/s)

rho0 = dy * 2700

totalInput= 0
totalOutput= 0
subsidence_rate = 0.#4e-4*dt/3 #m/yr
maxVolumeLoss = 0

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


## ## ## ## ## ## ## ##
##  Main time loop:  ##
## ## ## ## ## ## ## ##
while (t < tmax):
    #dt = 0.9*(dx*dx)/(2.e0*(max(k)))
    dt = 0.9*(dx*dx)/(2.e0*(sum(k)))
    #q0[0] = max(0, 2.e-6*math.sin(3.1415 * 6* t/tmax))
    
    #dt = min(tmax/1000, 0.9*(dx*dx)/(2.e0*(k[0])) )
    if t > tout_progress:
        print("      "+str( (math.ceil(100000*t/tmax))/1000 )+"%", end="\r") ##Track progress
        tout_progress += dtout_progress
    
    ## Set boundary condition at proximal end of the basin:
    
    ## Important note! If columns[0]["totalHeight"] is changed here to accommodate the boundary condition, no exception has to be made in ...
    ## ... the FTCS method for i==1 and i==0. BUT the old columns[0]["totalHeight"] is required in setNodes. So it would need some sort ...
    ## ... of columns[0]["OLD_height"] attribute.
    #newHeight[0] = columns[0]["totalHeight"]
    for p in range(2):
        #newHeight[0] = newHeight[0]+q0[p]*dt/dx
        sedIncrease[p,0] = q0[p]*dt/dx
        columns[0]["TotalSedContent"][p] += sedIncrease[p,0]
    
    ## The right end of basin always remains 0 for it is not altered by the FTCS scheme. Therefore there is no need to specify that boundary condition here.
    
    ## Keep track of volume balance:
    for p in range(2):
        totalInput+= q0[p]*dt
        totalOutput+= ( (k[p]*dt)/(dx*dx) )*( columns[imax-1]["totalHeight"] - columns[imax]["totalHeight"] )*dx ## Volume leaving = height*dx at i=imax, which is set to 0
    ## Loop through columns (for FTCS, density calculations, etc):
    for i in range(1,imax):
        newHeight[i]= columns[i]["totalHeight"] ##Start from current height
        newSedContent[p,i] = columns[i]["TotalSedContent"][p] ##Start from current sediment content
        for p in range(2):
            ## FTCS:
            #Dph= ( D(i+1)+D(i) )/2.e0
            #Dmh= ( D(i)+D(i-1) )/2.e0
            Dph = k[p]
            Dmh = k[p]
            sedIncrease[p,i] = ( (Dph*dt)/(dx*dx) )*( columns[i+1]["totalHeight"] - columns[i]["totalHeight"] )
            #newHeight[i]= newHeight[i] + ( (Dph*dt)/(dx*dx) )*( columns[i+1]["totalHeight"] - columns[i]["totalHeight"] )
            if (i==1):
                ##Consider effect of proximal boundary condition (i.e. i==0):
                sedIncrease[p,1] = sedIncrease[p,1] - ( (Dmh*dt)/(dx*dx) )*( columns[1]["totalHeight"] - columns[0]["totalHeight"]-q0[p]*dt/dx ) 
                #newHeight[1]= newHeight[1] - ( (Dmh*dt)/(dx*dx) )*( columns[1]["totalHeight"] - columns[0]["totalHeight"]-q0[p]*dt/dx ) 
            else:
                sedIncrease[p,i] = sedIncrease[p,i] - ( (Dmh*dt)/(dx*dx) )*( columns[i]["totalHeight"] - columns[i-1]["totalHeight"] )
            sedIncrease[p,i] = min(sedIncrease[p,i], columns[i-1]["TotalSedContent"][p])
            newSedContent[p,i] += sedIncrease[p,i]
            newHeight[i] += sedIncrease[p,i]
        
        ##Initial spike test:
        #if (i==70 and t==0): 
            #newHeight[70] = 100
            #for q in range(100):
                #columns[70]["nodes"].update({q:{"density":rho0}})
            #print(columns[70])
        columns[i] = setNodes(i, k, newHeight[i], columns[i], dt, dx, dy, rho0)  ## Update density, (todo grain size fractions and porosity)
        
    ## End column loop (i)
    
    ## Set nodes in proximal boundary column:
    for p in range(2):
        sedIncrease[p,0] = sedIncrease[p,0] - ( (k[p]*dt)/(dx*dx) )*( columns[0]["totalHeight"]+q0[p]*dt/dx - columns[1]["totalHeight"] )
        sedIncrease[p,0] = max(sedIncrease[p,0], -columns[0]["TotalSedContent"][p])
        newSedContent[p,0] += sedIncrease[p,0] 
        #columns[0]["TotalSedContent"][p] += sedIncrease[p,0]- q0[p]*dt/dx ## q0 is added before FTCS to both sedIncrease and "TotalSedContent". q0 should not be added twice so it is subtracted once here.
        newHeight[0] += sedIncrease[p,0]
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
        #call writeOutput(n, imax, x, hn, bedrock, totalHeight)
        tout = tout + dtout
    
    
    t += dt
    
## End time loop

## ## ## ## ## ## ## ##
##    Wrapping up:   ##
## ## ## ## ## ## ## ##

totalDepositVolume= 0
totalNodeVolume= 0
for i in range(imax+1):
    totalHeight[i]= columns[i]["totalHeight"]
    totalDepositVolume+= columns[i]["totalHeight"]*dx
    for j in range(0,len(columns[i]["nodes"])):
        try:
            totalNodeVolume+= dx*dy * columns[i]["nodes"][j]["density"]/rho0
        except:
            pass
    
    
print(columns[0]["TotalSedContent"], sum(columns[0]["TotalSedContent"]), totalHeight[0])
print("totalInput:", str(totalInput)+"m^3")
print("totalOutput:", str(totalOutput)+"m^3")
print("totalDepositVolume:", str(totalDepositVolume)+"m^3")
print("totalNodeVolume:", str(totalNodeVolume)+"m^3")
print("")
print("Volume error:", str(totalDepositVolume+totalOutput-totalInput)+"m^3")
print("Error in %:", str(100*(totalDepositVolume+totalOutput-totalInput)/totalInput)+"%")
print("")
print("Node volume error:", str(totalNodeVolume+totalOutput-totalInput)+"m^3")
print("Node volume error in %:", str(100*(totalNodeVolume+totalOutput-totalInput)/totalInput)+"%")
print("")

#plt.plot(x, totalHeight)
#plt.show()

os.system("python3 animate_ISMoLD.py")


    





#columns[1]["nodes"].update({i:{"fractions":[0.3, 0.3, 0.4]}})

