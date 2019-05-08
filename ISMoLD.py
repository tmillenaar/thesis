## Program ISMoLD (Imperfect Sorting Model of Linear Diffusion)
## Designed and created by Timo Millenaar (tmillenaar@gmail.com)

def setNodes(i, k, newHeight, column, dt, dx, dy, rho0):
    if (newHeight[i] >= column["height"]):  ## Deposition
        dispose_density = ( newHeight[i] - column["height"] ) * dx
        j = math.floor( column["height"] )
        
        for j in range (math.floor( column["height"] ), math.ceil( newHeight[i] )):
            #print(j, math.ceil( newHeight[i]))
            if (j == math.floor( column["height"]) ):
                try: 
                    starting_box_density = dy * column["nodes"][j]["density"]
                except:
                    starting_box_density= 0
                
                if (j==math.ceil( newHeight[i] )-1): ## In the usual case that not more then one dy is added in height this timestep
                    if (starting_box_density+ rho0 *(newHeight[i]-column["height"]) > rho0):
                        #if i==1:
                            #print(starting_box_density)
                        column["nodes"].update({j-1:{"density":rho0}})
                        remainder = starting_box_density+ rho0 *(newHeight[i]-column["height"]) -rho0
                        column["nodes"].update({j:{"density":( remainder )}})
                    else:
                        remainder = starting_box_density+ rho0 *(newHeight[i]-column["height"])
                        column["nodes"].update({j:{"density":remainder}})
                    
            elif (j==math.ceil( newHeight[i] )-1):
                if (starting_box_density+ rho0 *(newHeight[i]-column["height"]) > rho0):
                    #print("This should not have happened either!!", j)
                    column["nodes"].update({j-1:{"density":rho0}})
                    #print("starting_box_density "+str(starting_box_density))
                    #print("nr"+str(starting_box_density/rho0 + j-newHeight[i])+" "+str(starting_box_density/rho0)+" "+str(j-newHeight[i]))
                    remainder = starting_box_density+ rho0 *(newHeight[i]-column["height"]) -rho0
                    column["nodes"].update({j:{"density":( remainder )}})
                else:
                    print("oh no!")
                    remainder = starting_box_density+ rho0 *(newHeight[i]-column["height"])
                    column["nodes"].update({j:{"density":( remainder )}})
            else:
                column["nodes"].update({j:{"density":rho0}})
                
    else:  ## Erosion
        for j in range(math.floor(newHeight[i]), math.ceil(column["height"])+1):
            del column["nodes"][j]
    return column

import math
import matplotlib.pyplot as plt

yr2sec = 60*60*24*365.25
dx= 1.e3       # lateral grid spacing (m)
dy = 1         # vertical grid spacing (m)
imax= 101      # number of nodes
totalHeight= list(range(imax+1))
newHeight= list(range(imax+1))
t= 0
dt= 1          # time step in (yr)
tmax= 100000*yr2sec     # max number of years
dtout = 10*yr2sec     # nr of years between write output
tout = 0.          # threshold to write next output in years (increased every write action by dtout)
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
    totalHeight[i] = 0
    columns[i]= {"height":0,
                 "nodes":{}
                 }

## ## ## ## ## ## ## ##
##  Main time loop:  ##
## ## ## ## ## ## ## ##
while (t < tmax):
    
    dt = 0.9*(dx*dx)/(2.e0*(max(k)))
    #dt = min(tmax/1000, 0.9*(dx*dx)/(2.e0*(k[0])) )
    if t > tout:
        print("      "+str( (math.ceil(100000*t/tmax))/1000 )+"%", end="\r") ##Track progress
        tout += dtout
    
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
            newHeight[i]= newHeight[i] - ( (Dmh*dt)/(dx*dx) )*( columns[i]["height"] - columns[i-1]["height"]-q0[0]*dt/dx ) ##Proximal boundary condition
        else:
            newHeight[i]= newHeight[i] - ( (Dmh*dt)/(dx*dx) )*( columns[i]["height"] - columns[i-1]["height"] )
        
        columns[i] = setNodes(i, k, newHeight, columns[i], dt, dx, dy, rho0)
        
    ## Set nodes in proximal boundary column:
    newHeight[0]= newHeight[0] - ( (Dph*dt)/(dx*dx) )*( columns[0]["height"]+q0[0]*dt/dx - columns[1]["height"] )
    columns[0] = setNodes(0, k, newHeight, columns[0], dt, dx, dy, rho0)
    
    ## Overwrite old the topography with the new profile:
    for i in range(0,imax):
        columns[i]["height"] = newHeight[i]
    
    t += dt
    

## ## ## ## ## ## ## ##
##    Wrapping up:   ##
## ## ## ## ## ## ## ##

totalDepositVolume= 0
totalNodeVolume= 0
x=list(range(imax+1))
for i in range(imax+1):
    #print(i, columns[i]["height"])
    #columns[i]["height"]= columns[i]["height"]
    x[i]=i
    totalHeight[i]= columns[i]["height"]
    totalDepositVolume+= columns[i]["height"]*dx
    for j in range(0,len(columns[i]["nodes"])):
        try:
            totalNodeVolume+= dx*dy * columns[i]["nodes"][j]["density"]/rho0
        except:
            pass
    
    
print(columns[0])
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
plt.plot(x, totalHeight)
plt.show()


    





#columns[1]["nodes"].update({i:{"fractions":[0.3, 0.3, 0.4]}})

#print(columns[1]["nodes"])



#print(columns)
