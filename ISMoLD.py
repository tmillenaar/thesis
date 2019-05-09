## Program ISMoLD (Imperfect Sorting Model of Linear Diffusion)
## Designed and created by Timo Millenaar (tmillenaar@gmail.com)

import math
import os
import matplotlib.pyplot as plt

from functions_for_ISMoLD import *

try:
    os.system('rm ISMolD_outputdata/topography*.txt') ##remove old data files
except:
    pass

yr2sec = 60*60*24*365.25
dx= 1.e3       # lateral grid spacing (m)
dy = 1         # vertical grid spacing (m)
imax= 101      # number of nodes
totalHeight= list(range(imax+1))
newHeight= list(range(imax+1))
x=list(range(imax+1))
t= 0
dt= 1          # time step in (yr)
tmax= 100000*yr2sec     # max number of years
dtout = tmax/100     # nr of years between write output
dtout_progress = 10*yr2sec     # nr of years between progress bar update
tout = 0.          # threshold to write next output in years (increased every write action by dtout)
tout_progress = 0.          # threshold to update progress bar
tprogress = 0.
dtprogress = tmax/1000. 
k= list(range(3))
k[0]= 3.2e-4       # Gravel diffusivity (m2/s)
k[1]= 3.2e-3       # Sand diffusivity (m2/s)
k[2]= 3.2e-3       # Silt diffusivity (m2/s)
q0= list(range(3))
q0[0]= 2.e-6      # Gravel input (m2/s)
q0[1]= 1.e-10      # Sand input (m2/s)
q0[2]= 0           # silt input (m2/s)

rho0 = 2700

totalInput= 0
totalOutput= 0
subsidence_rate = 0.#4e-4*dt/3 #m/yr
maxVolumeLoss = 0



## Initialize:
columns = {}
for i in range(imax+1):
    x[i]=i
    totalHeight[i] = 0
    columns[i]= {"height":0,
                 "nodes":{}
                 }

## ## ## ## ## ## ## ##
##  Main time loop:  ##
## ## ## ## ## ## ## ##
while (t < tmax):
    
    dt = 0.9*(dx*dx)/(2.e0*(max(k)))
    q0[0] = max(0, 2.e-6*math.sin(3.1415 * 6* t/tmax))
    
    #dt = min(tmax/1000, 0.9*(dx*dx)/(2.e0*(k[0])) )
    if t > tout_progress:
        print("      "+str( (math.ceil(100000*t/tmax))/1000 )+"%", end="\r") ##Track progress
        tout_progress += dtout_progress
    
    ## Set boundary condition at proximal end of the basin:
    
    ## Important note! If columns[0]["height"] is changed here to accommodate the boundary condition, no exception has to be made in ...
    ## ... the FTCS method for i==1 and i==0. BUT the old columns[0]["height"] is required in setNodes. So it would need some sort ...
    ## ... of columns[0]["OLD_height"] attribute.
    newHeight[0] = columns[0]["height"]+q0[0]*dt/dx
    
    ## The right end of basin always remains 0 for it is not altered by the FTCS scheme. Therefore there is no need to specify that boundary condition here.
    
    ## Keep track of volume balance:
    totalInput+= q0[0]*dt
    totalOutput+= ( (k[0]*dt)/(dx*dx) )*( columns[imax-1]["height"] - columns[imax]["height"] )*dx ## Volume leaving = height*dx at i=imax
    
    ## Loop through columns (for FTCS, density calculations, etc):
    for i in range(1,imax):
        #Dph= ( D(i+1)+D(i) )/2.e0
        #Dmh= ( D(i)+D(i-1) )/2.e0
        Dph = k[0]
        Dmh = k[0]
        newHeight[i]= columns[i]["height"] + ( (Dph*dt)/(dx*dx) )*( columns[i+1]["height"] - columns[i]["height"] )
        if (i==1):
            ##Consider effect of proximal boundary condition (i.e. i==0):
            newHeight[1]= newHeight[1] - ( (Dmh*dt)/(dx*dx) )*( columns[1]["height"] - columns[0]["height"]-q0[0]*dt/dx ) 
        else:
            newHeight[i]= newHeight[i] - ( (Dmh*dt)/(dx*dx) )*( columns[i]["height"] - columns[i-1]["height"] )
        
        ##Initial spike test:
        #if (i==70 and t==0): 
            #newHeight[70] = 100
            #for p in range(100):
                #columns[70]["nodes"].update({p:{"density":rho0}})
            #print(columns[70])
        columns[i] = setNodes(i, k, newHeight[i], columns[i], dt, dx, dy, rho0)
        
    ## Set nodes in proximal boundary column:
    newHeight[0]= newHeight[0] - ( (Dph*dt)/(dx*dx) )*( columns[0]["height"]+q0[0]*dt/dx - columns[1]["height"] )
    columns[0] = setNodes(0, k, newHeight[0], columns[0], dt, dx, dy, rho0)
    
    ## Overwrite old the topography with the new profile:
    for i in range(0,imax):
        columns[i]["height"] = newHeight[i]
    
    
    if (t >= tout):
        n= int(tout/dtout)
        f = open("ISMolD_outputdata/topography"+str(n)+".txt", "w")
        for i in range(len(x)):
            f.write(str(x[i])+" "+str(columns[i]["height"])+"\n")
        f.close()
        #call writeOutput(n, imax, x, hn, bedrock, totalHeight)
        tout = tout + dtout
    
    
    t += dt
    

## ## ## ## ## ## ## ##
##    Wrapping up:   ##
## ## ## ## ## ## ## ##

totalDepositVolume= 0
totalNodeVolume= 0
for i in range(imax+1):
    totalHeight[i]= columns[i]["height"]
    totalDepositVolume+= columns[i]["height"]*dx
    for j in range(0,len(columns[i]["nodes"])):
        try:
            totalNodeVolume+= dx*dy * columns[i]["nodes"][j]["density"]/rho0
        except:
            pass
    
    
#print(columns[70])
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
print("Plotting...")
#plt.plot(x, totalHeight)
#plt.show()

os.system("python3 animate_ISMoLD.py")


    





#columns[1]["nodes"].update({i:{"fractions":[0.3, 0.3, 0.4]}})

