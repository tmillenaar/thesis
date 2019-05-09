## Program ISMoLD (Imperfect Sorting Model of Linear Diffusion)
## Designed and created by Timo Millenaar (tmillenaar@gmail.com)

def setNodes(i, k, newHeight, column, dt, dx, dy, rho0):
    
    if (newHeight[i] >= column["height"]):  ## Deposition
        dispose_density = ( newHeight[i] - column["height"] ) * dx
        j = math.floor( column["height"] )
        
        for j in range (math.floor( column["height"] ), math.ceil( newHeight[i] )):
            if (j == math.floor( column["height"]) ):
                if (j==math.ceil( newHeight[i] )-1): ## In the usual case that the first j is also the last
                    dh = newHeight[i] - column["height"]
                else:
                    dh = math.ceil( column["height"]) - column["height"]
            else:
                if (j==math.ceil( newHeight[i] )-1): ## If current j is final j
                    dh = newHeight[i]-j
                else:
                    dh = dy
            
            try:
                current_box_density = dy * column["nodes"][j]["density"]
            except:
                current_box_density = 0
                
            if (j==math.ceil( newHeight[i] )-1): ## Final j
                if (current_box_density+ rho0 *(dh) >= rho0):
                    column["nodes"].update({j-1:{"density":rho0}})
                    remainder = current_box_density+ rho0 *(dh) -rho0
                    print(1,remainder)
                    column["nodes"].update({j:{"density":( remainder )}})
                else:
                    remainder = current_box_density+ rho0 *(dh)
                    print(2,remainder)
                    column["nodes"].update({j:{"density":remainder}})
            else:
                print('yep')
                column["nodes"].update({j:{"density":rho0}})
                
    else:  ## Erosion
        #starting_box_density = dy * column["nodes"][ math.ceil(column["height"]) ]["density"]
        #print(math.ceil(column["height"]), math.floor(newHeight[i]))
        print("start")
        for j in range(math.floor(column["height"]), math.floor(newHeight[i]), -1):
            print(j)
            if (j == math.ceil( column["height"] )):
                if (j==math.ceil( newHeight[i] )-1): ## In the usual case that the first j is also the last
                    dh = column["height"] - newHeight[i]
                else:
                    dh = column["height"] - math.floor( column["height"] )
            else:
                if (j==math.ceil( newHeight[i] )-1): ## If current j is final j
                    dh = j - newHeight[i] 
                else:
                    dh = dy
                    
            try:
                current_box_density = dy * column["nodes"][j]["density"]
            except:
                current_box_density = 0
                    
            if (j==math.floor( newHeight[i] )+1): ## In the usual case that not more then one dy is removed in height this timestep
                if (rho0 *(dh) >= current_box_density):
                    try:
                        del column["nodes"][j]
                    except:
                        pass
                    remainder = rho0*( dh ) - current_box_density
                    #print(1, remainder, rho0*( j-newHeight[i] ), -current_box_density)
                    column["nodes"].update({j-1:{"density":remainder}})
                else:
                    remainder = current_box_density - rho0*( dh )
                    #print(2, remainder, rho0*( j-newHeight[i] ), -current_box_density)
                    column["nodes"].update({j:{"density":remainder}})
            else:
                try:
                    del column["nodes"][j]
                except:
                    pass
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
tmax= 1000*yr2sec     # max number of years
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

totalInput= 100*dx
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
        
        if (i==70 and t==0): 
            newHeight[70] = 100
            for p in range(100):
                columns[70]["nodes"].update({p:{"density":rho0}})
            #print(columns[70])
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
    
    
print(columns[70])
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
#plt.show()


    





#columns[1]["nodes"].update({i:{"fractions":[0.3, 0.3, 0.4]}})

#print(columns[1]["nodes"])



#print(columns)
