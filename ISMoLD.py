## Program ISMoLD (Imperfect Sorting Model of Linear Diffusion)
## Designed and created by Timo Millenaar (tmillenaar@gmail.com)

import math
import os
import pdb ## pdb.set_trace()
import matplotlib.pyplot as plt
import numpy as np
from os import listdir
from decimal import Decimal

from functions_for_ISMoLD import *

makeDirectories()

enableSubsidence = False
animateHeight = False
plotNodes = False
spikeTest = False
pbdtrace = False
testSedTransport = False

##Uncomment if not desired:
enableSubsidence = True

#plotHeigt = True 
#plotNodes = True

##spikeTest = True
#pbdtrace = True
#testSedTransport = True

yr2sec = 60*60*24*365.25      #nr of seconds in a year

dx= 1e3                      # lateral grid spacing (m)
imax= 100                     # number of nodes
tmax= 100000*yr2sec             # total amount of time to be modelled in seconds [in years = (x*yr2sec)]
dtout = tmax/100              # nr of years between write output
dtout_progress = 10*yr2sec    # nr of years between progress bar update
nrOfGrainSizes = 2
k= list(range(nrOfGrainSizes))
k[0]= 1*3.2e-4                # Gravel diffusivity (m2/s)
k[1]= 2*3.2e-4                # Sand diffusivity (m2/s)
q0= list(range(nrOfGrainSizes))
q0[0]= 3*1.e-6                # Gravel input (m2/s)
q0[1]= 2*1.e-6                # Sand input (m2/s)
subsidenceRate = 1e-6         # Subsidence rate at the proximal end of the basin (m/yr), negative values result in uplift
spikeLocation = 50             # Only relevant if spikeTest == True

transportDensity = 2700
#transportPorosity = 0.3

## Initialize:
totalSedInputThroughTime0 = list(range(nrOfGrainSizes))
totalSedInputThroughTime1 = list(range(nrOfGrainSizes))
inputFractions = list(range(nrOfGrainSizes))
totalHeight= list(range(imax+1))
newHeight= list(range(imax+1))
x= list(range(imax+1))
sedContentInActiveLayer = np.zeros(shape=(2,imax+1))
totalSedcontenti0 = []
totalSedcontenti1 = []
activeLayeri0 = []
activeLayeri1 = []
time = []
sedIncrease = np.zeros(shape=(2,imax+1))
sedIn = np.zeros(shape=(2,imax+1))
sedOut = np.zeros(shape=(2,imax+1))

dy = 1         # vertical grid spacing (m)  Important: MUST BE 1(as integer) for code to work (in this version)!

tout = 0.          # threshold to write next output in years (increased every write action by dtout)
tout_progress = 0.          # threshold to update progress bar
tprogress = 0.

totalInput= 0
totalOutput= 0
InputPerGrainSize = [0,0]
OutputPerGrainSize = [0,0]
maxVolumeLoss = 0
nodeOutputTimestep = 0

totalSedInputThroughTime0 = []
totalSedInputThroughTime1 = []
totalSedOutputThroughTime0 = []
totalSedOutputThroughTime1 = []
columns = {}
for i in range(imax+1):
    x[i]=i
    totalHeight[i] = 0
    columns[i]= {"totalHeight":0,
                 "bedrockHeight": 0,
                 "totalSedContent":list(range(nrOfGrainSizes)),
                 "nodes":{},
                 "newSedContent": list(range(nrOfGrainSizes)),
                 "totalSedHeight": list(range(nrOfGrainSizes)),
                 }
    for p in range(nrOfGrainSizes):
        columns[i]["totalSedContent"][p] = 0
        columns[i]["newSedContent"][p] = 0
        columns[i]["totalSedHeight"][p] = 0
        

## Remove old data files
nr_topo_files = len(listdir("ISMolD_outputdata/relief")) #-1 since file starts at 'output0', add -1 for amy subdirectory
if (nr_topo_files > 1):
    os.system('rm ISMolD_outputdata/relief/topography*.txt')
    os.system('rm -r -f ISMolD_outputdata/nodes/*')

## Spike test:
#  Note: the code cannot properly deal with intense spikes with a width of 1 dx (except for i==0), thus the the spike is 2*dx wide.
if (spikeTest):
    spikeHeight = 100
    #spikeLocation = 3
    totalInput+= spikeHeight*dx
    for p in range(nrOfGrainSizes):
        InputPerGrainSize[p] += dx*spikeHeight/2
        q0[p] = 0
    columns[spikeLocation]["totalHeight"] = spikeHeight
    columns[spikeLocation]["totalSedContent"][0] = spikeHeight/2
    columns[spikeLocation]["totalSedContent"][1] = spikeHeight/2
    for q in range(spikeHeight):
        columns[spikeLocation]["nodes"].update({q:{"density":transportDensity,
                                                   "nodeSedContent":[0.5,0.5]}})

    totalInput+= spikeHeight*dx
    for p in range(nrOfGrainSizes):
        InputPerGrainSize[p] += dx*spikeHeight/2
    columns[spikeLocation+1]["totalHeight"] = spikeHeight
    columns[spikeLocation+1]["totalSedContent"][0] = spikeHeight/2
    columns[spikeLocation+1]["totalSedContent"][1] = spikeHeight/2
    for q in range(spikeHeight):
        columns[spikeLocation+1]["nodes"].update({q:{"density":transportDensity,
                                                     "nodeSedContent":[0.5,0.5]}})

## ## ## ## ## ## ## ##
##  Main time loop:  ##
## ## ## ## ## ## ## ##
t= 0
while (t < tmax):
    
    # Varying sediment input or diffusivity through time can be set here:
    q0[0] = max(0, 0.5e-6*math.sin(3.1415 * 4* t/tmax))
    q0[1] = max(0, 0.5e-6*math.sin(3.1415 * 4* t/tmax))
    
    #if(t>0.4*tmax): 
        #for p in range(nrOfGrainSizes):
            #q0[p] = 0.e-7
    
    #if(t<0.15*tmax): 
        #q0[0] = 1.e-6
        #q0[1] = 0
    #elif(t<0.3*tmax ): 
        #q0[0] = 0
        #q0[1] = 1.e-6
    #elif (t>0.3*tmax):
        #for p in range(nrOfGrainSizes):
            #q0[p] = 0.e-7
    
    if (spikeTest == True):
        for p in range(nrOfGrainSizes):
            q0[p] = 0
    ## Calculating max dt as limited by the active layer and the FTCS scheme
    dtFTCS = 0.9*(dx*dx)/(2.e0*(sum(k))) ##FTCS
    maxHeightDiff = columns[0]["totalHeight"]+sum(q0)*dtFTCS/dx-columns[1]["totalHeight"] ## dtFTCS is used here as an estimation on the high side. It is better to estimate more sediment input resulting in a high slope then it is to estimate on the low side. A too low estimate of the sediment input can lead to an amount of sediment being removed somewhere that is more the the depth of the active layer. Higher estimates of q0 result in lowe dt and thus a slower run speed.
    for i in range(1,imax):
        #maxHeightDiff = max(maxHeightDiff, abs(columns[i+1]["totalHeight"] - 2*columns[i]["totalHeight"] + columns[i-1]["totalHeight"]) )
        maxHeightDiff = max(maxHeightDiff, abs(columns[i-1]["totalHeight"] - columns[i]["totalHeight"]) )
    dtHeightLimit = 0
    if (maxHeightDiff != 0): ## Cannot devide by 0
        dtHeightLimit += 0.9*(dx*dx)/( sum(k) * maxHeightDiff )
        dt = min(dtFTCS, dtHeightLimit)
    else:
        dt = 0.1*(dx*dx)/(2.e0*(sum(k))) ##FTCS
    #dt = 0.1*dt
    
    if t > tout_progress:
        print("      "+str( (math.ceil(100000*t/tmax))/1000 )+"%", end="\r") ##Track progress
        tout_progress += dtout_progress
    
    ## Set boundary condition at proximal end of the basin:
    for i in range(imax+1):
        ## "oldHeight" and "oldSedContent" are needed for the setNodes function at i==0. "totalHeight" would have been sufficient for the other nodes as the "old" height. The "oldHeight" is given to all nodes here to keep the setNodes function clear for a human to read. Similarly for "oldSedContent":
        columns[i]["oldHeight"] = columns[i]["totalHeight"]
        columns[i]["oldSedContent"] = columns[i]["totalSedContent"].copy()
        
    if (enableSubsidence): subsidence (columns, imax, dt*subsidenceRate/yr2sec, nrOfGrainSizes)
    
    #if (columns[imax-8]["totalHeight"] > columns[imax-6]["totalHeight"]): ## Slope goes down to the right
        #print("SLOPE TO RIGHT")
        #print(columns[imax-1]["totalHeight"]-columns[imax-1]["bedrockHeight"])
        
    ## Keep track of volume balance:
    time.append(t)
    for p in range(nrOfGrainSizes):
        totalInput+= q0[p]*dt
        InputPerGrainSize[p] += q0[p]*dt
        if (p == 0):
            totalSedInputThroughTime0.append(InputPerGrainSize[p])
        if (p == 1):
            totalSedInputThroughTime1.append(InputPerGrainSize[p])
            
            
        if (columns[imax-1]["totalHeight"] > columns[imax]["totalHeight"]): ## Slope goes down to the right
            output = ( (k[p]*dt)/(dx*dx) )*( (columns[imax-1]["totalHeight"]) - (columns[imax]["totalHeight"]) )*dx ## Volume leaving = height*dx at i=imax, which is set to 0
            output = min(output, sedContentInActiveLayer[p,i-1])
            #totalOutput += output
            #OutputPerGrainSize[p] += ( (k[p]*dt)/(dx*dx) )*( columns[imax-1]["totalHeight"] - (columns[imax]["totalHeight"]) )*dx
            if (p==0):
                totalSedOutputThroughTime0.append(OutputPerGrainSize[p])
            if (p==1):
                totalSedOutputThroughTime1.append(OutputPerGrainSize[p])
        else:
            totalOutput+= 0
            OutputPerGrainSize[p] += 0
            if (p==0):
                totalSedOutputThroughTime0.append(OutputPerGrainSize[p])
            if (p==1):
                totalSedOutputThroughTime1.append(OutputPerGrainSize[p])
            
    
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
        for p in range(nrOfGrainSizes):
            if (sedContentInActiveLayer[p,i] < 0):
                print("sedContentInActiveLayer < 0:", sedContentInActiveLayer[:,i], columns[i])
                exit()
    
        
    ## Loop through columns (for FTCS, density calculations, etc):
    totalSedIn = 0
    totalSedOut = 0
    for i in range(imax+1): 
        ## The case i==0 is treated as a special case in this loop and i==imax+1 remains untouched, having it stay 0
        
        newHeight[i]= columns[i]["totalHeight"] ##Start from current height
        for p in range(nrOfGrainSizes):
            sedIn[p,i] = 0
            sedOut[p,i] = 0
            columns[i]["newSedContent"][p] = columns[i]["totalSedContent"][p] ##Start from current sediment content
            ## FTCS:
            #Dph= ( D(i+1)+D(i) )/2.e0
            #Dmh= ( D(i)+D(i-1) )/2.e0
            Dph = k[p]
            Dmh = k[p]
            
            if (i==0): 
                if (columns[0]["totalHeight"] >= columns[1]["totalHeight"]): ## Slope goes down to the right
                    sedIn[p,0] = q0[p]*dt/dx
                    totalSedIn -= q0[p]*dt/dx
                else: ## Slope goes down to the left
                    sedIn[p,0] += q0[p]*dt/dx
                    totalSedIn -= q0[p]*dt/dx
                    sedOut[p,0] = 0
            else:
                
                if (columns[i-1]["totalHeight"] > columns[i]["totalHeight"]): ## Slope goes down to the right
                    transport = ( (Dph*dt)/(dx*dx) )*( (columns[i-1]["totalHeight"]) - ((columns[i]["totalHeight"])) )
                    transport = min(transport, sedContentInActiveLayer[p,i-1])
                    sedOut[p,i-1] += transport
                    sedIn[p,i] += transport

                elif (columns[i-1]["totalHeight"] < columns[i]["totalHeight"]): ## Slope goes down to the left
                    transport = ( (Dph*dt)/(dx*dx) )*( (columns[i]["totalHeight"]) - ((columns[i-1]["totalHeight"])) )
                    transport = min(transport, sedContentInActiveLayer[p,i])
                    sedOut[p,i] += transport
                    sedIn[p,i-1] += transport
                        
            
    ## Loop through columns for a second time to fix the minima and maxima:
    for i in range(imax+1): 
        originalTotalSedOut = 0
        fractionToLeft = list(range(nrOfGrainSizes))
        fractionToRight = list(range(nrOfGrainSizes))
        for p in range(nrOfGrainSizes):            
            fractionToLeft[p] = 0
            fractionToRight[p] = 0
            originalTotalSedOut += sedOut[p,i]
        for p in range(nrOfGrainSizes):
            totalSedIn += sedIn[p,i]
            totalSedOut += sedOut[p,i]
    
    if (testSedTransport):
        print("") 
        print(totalSedIn, totalSedOut, totalSedIn - totalSedOut)
        for i in range(imax+1):
            print(i, sedIn[:,i], sedOut[:,i], sedContentInActiveLayer[p,i])
        print("") 
    if (totalSedIn - totalSedOut > 1e-14):
        print("Error, sedIn != sedOut", totalSedIn, totalSedOut, totalSedIn - totalSedOut)
        for i in range(imax+1):
            print(i, sedIn[:,i], sedOut[:,i], sedContentInActiveLayer[p,i])
        exit()
           
    for p in range(nrOfGrainSizes):
        totalOutput += sedIn[p,imax]*dx
        OutputPerGrainSize[p] += sedIn[p,imax]*dx
        
    ## Loop through columns for a second time to fix the minima and maxima:
    for i in range(imax): 
        ## The case i==0 is treated as a special case in this loop and i==imax+1 remains untouched, having it stay 0
        
        newHeight[i]= columns[i]["totalHeight"] ##Start from current height
        for p in range(nrOfGrainSizes):
            
            columns[i]["newSedContent"][p] = columns[i]["newSedContent"][p] + sedIn[p,i]-sedOut[p,i]
            #print("HERE", p, sedIn[p,0], columns[0]["totalSedContent"][p], columns[0]["newSedContent"][p], sedIn[p,i], sedOut[p,i], sedContentInActiveLayer[p,i])
            if (columns[i]["newSedContent"][p] < columns[i]["bedrockHeight"]-1e-14 and subsidenceRate > 0): 
                print("Attention: newSedContent["+str(p)+"] (column "+str(i)+") was set to "+str(columns[i]["bedrockHeight"])+", originally:", columns[i]["newSedContent"][p])
                columns[i]["newSedContent"][p] = columns[i]["bedrockHeight"]
                print(p, columns[i]["newSedContent"][p], sedIn[p,i], sedOut[p,i], "   Active layer:", sedContentInActiveLayer[p,i])
                print("")
                print("")
                print("")
                printColumn(columns[i], i, newHeight[i], columns[i]["newSedContent"])
                exit()
                
            newHeight[i] += sedIn[p,i]-sedOut[p,i]
            if (abs(sedIn[p,i]-sedOut[p,i]) < -1): ## Sediment input at the proximal end (i=0) is determined by the boundary condition and may be more then 1m in height
                print("Error, active Layer < -1",i , dt, sedIn[p,i], -sedOut[p,i], sedIn[p,i]-sedOut[p,i], sedContentInActiveLayer[:,i])
                exit()
                
        if (i==15): 
            totalSedcontenti0.append(columns[15]["newSedContent"][0])
            totalSedcontenti1.append(columns[15]["newSedContent"][1])
            activeLayeri0.append(sedContentInActiveLayer[0,15])
            activeLayeri1.append(sedContentInActiveLayer[1,15])
            #time.append(t)
        
        ## Update the nodes to suite newHeight and newSedContent:
        columns[i] = setNodes(i, k, newHeight[i], columns[i], columns[i]["newSedContent"], dt, dx, dy, transportDensity)  
        
    ## End column loop (i)
    
    ## Overwrite old the topography with the new profile:
    for i in range(0,imax):
        columns[i]["totalHeight"] = newHeight[i]
        for p in range(nrOfGrainSizes):
            columns[i]["totalSedContent"][p] = columns[i]["newSedContent"][p]
    
    #print("type:",newHeight.dtype, newSedContent.dtype, columns[0]["totalHeight"].dtype)
    #print("type:",newHeight[0], columns[0]["totalHeight"], newSedContent[0,0], columns[0]["nodes"][0]["nodeSedContent"])
    #print(np.dtype(newSedContent[0,0]))
    #print("")
    if (t >= tout):
        
        makeTimeNodeDirectory(nodeOutputTimestep)
        
        for i in range(len(x)-1): ## -1 for the last column is always empty (by design). Therefore there is no need to create a file for it.
            f = open("ISMolD_outputdata/nodes/time"+str(nodeOutputTimestep)+"/column"+str(i)+".txt", "w")
            for j in range(len(columns[i]["nodes"])):
                writeline = str(j)+" "+ str(columns[i]["totalHeight"])
                writeline += " "+str(columns[i]["bedrockHeight"])
                for p in range(nrOfGrainSizes):
                    writeline += " "+str(columns[i]["nodes"][j]["nodeSedContent"][p])
                writeline += "\n"
                f.write(writeline)
            f.close()
        
        nodeOutputTimestep += 1
        
        n= int(tout/dtout)
        f = open("ISMolD_outputdata/relief/topography"+str(n)+".txt", "w")
        for i in range(len(x)):
            writeline = str(x[i])+" "+str(columns[i]["totalHeight"])
            writeline += " "+str(columns[i]["bedrockHeight"])
            for p in range(nrOfGrainSizes):
                writeline += " "+str(columns[i]["totalSedContent"][p] + columns[i]["bedrockHeight"]) ## Note that bedrockHeight is generally negative
            writeline += "\n"
            f.write(writeline)
        f.close()
        tout = tout + dtout
        
    t += dt
    
## End time loop

## ## ## ## ## ## ## ##
##    Wrapping up:   ##
## ## ## ## ## ## ## ##

totalDepositVolume= 0
totalNodeVolume= 0
totalVolumeBelow0= 0
depositVolumeBelow0= 0
DepositVolPerGrainSize = [0,0]
totalGrainSizeFraction = [0,0]
for i in range(imax+1):
    #totalHeight[i]= columns[i]["totalHeight"]
    totalVolumeBelow0 += -columns[i]["bedrockHeight"]*dx
    if(columns[i]["totalHeight"]-columns[i]["bedrockHeight"] == 0):
        depositVolumeBelow0 += 0
    elif(columns[i]["totalHeight"] > 0):
        print("Too much", i, columns[i]["totalHeight"], columns[i]["bedrockHeight"])
        depositVolumeBelow0 += -columns[i]["bedrockHeight"]*dx ## "-" for bedrockHeight value is negative during subsidence
    else:
        print("Too little", i, columns[i]["bedrockHeight"] - columns[i]["totalHeight"], columns[i]["bedrockHeight"])
        depositVolumeBelow0 += (columns[i]["totalHeight"] - columns[i]["bedrockHeight"])*dx
    totalDepositVolume+= (columns[i]["totalHeight"]-columns[i]["bedrockHeight"])*dx
    for p in range(nrOfGrainSizes):
        DepositVolPerGrainSize[p] += columns[i]["totalSedContent"][p]*dx
    for j in range(0,len(columns[i]["nodes"])+1):
        try:
            totalNodeVolume+= dx*dy * columns[i]["nodes"][j]["density"]/transportDensity
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

if (pbdtrace):
    pdb.set_trace()

    
#printColumn(columns[2], 2, newHeight[2], columns[2]["newSedContent"])
#for i in range(imax):
    #printColumn(columns[i], i, newHeight[i], columns[i]["newSedContent"])
    #print("")

print("totalInput:", str(totalInput)+"m^3")
print("InputPerGrainSize:", str(InputPerGrainSize)+"m^3")
print("totalOutput:", str(totalOutput)+"m^3")
print("OutputPerGrainSize:", str(OutputPerGrainSize)+"m^3")
print("")
if(columns[0]["bedrockHeight"] < 0):
    print("totalVolumeBelow 0:", str(totalVolumeBelow0)+"m^3")
    print("depositVolumeBelow 0:", str(depositVolumeBelow0)+"m^3")
print("totalDepositVolume:", str(totalDepositVolume)+"m^3")
print("totalNodeVolume:", str(totalNodeVolume)+"m^3")
print("DepositVolPerGrainSize:", str(DepositVolPerGrainSize)+"m^3", "  sum:", sum(DepositVolPerGrainSize))
print("totalGrainSizeFraction", totalGrainSizeFraction)
print("")
if (totalDepositVolume+totalOutput-totalInput > 0):
    print("Volume error:", str(totalDepositVolume+totalOutput-totalInput)+"m^3", " (positive => excess sediment in body)")
elif(totalDepositVolume+totalOutput-totalInput == 0):
    print("Volume error:", str(totalDepositVolume+totalOutput-totalInput)+"m^3", " (Really, 0 error? Impressive and unexpected)")
else:
    print("Volume error:", str(totalDepositVolume+totalOutput-totalInput)+"m^3", " (negative => missing sediment in body)")
print("Error in %:", str(100*(totalDepositVolume+totalOutput-totalInput)/totalInput)+"%")
print("")
for p in range(nrOfGrainSizes):
    print("Volume error per grain size ("+str(p)+"):", str(DepositVolPerGrainSize[p]+OutputPerGrainSize[p]-InputPerGrainSize[p])+"m^3")
    print("Error in %:", str(100*(DepositVolPerGrainSize[p]+OutputPerGrainSize[p]-InputPerGrainSize[p])/InputPerGrainSize[p])+"%")
print("")
print("Node volume error:", str(totalNodeVolume+totalOutput-totalInput)+"m^3")
print("Node volume error in %:", str(100*(totalNodeVolume+totalOutput-totalInput)/totalInput)+"%")
print("totalDepositVolume-totalNodeVolume:", Decimal(totalDepositVolume), Decimal(totalNodeVolume), totalDepositVolume-totalNodeVolume)
print(Decimal(totalDepositVolume))
print(Decimal(totalNodeVolume))
print("")


#fig, ax1 = plt.subplots()

#ax1.plot(time, totalSedcontenti0, label="Gravel")
#ax1.plot(time, totalSedcontenti1, label="Sand")

#ax2 = ax1.twinx()
#ax2.plot(time, activeLayeri0, linestyle="--", label="Gravel in active layer")
#ax2.plot(time, activeLayeri1, linestyle="--", label="Sand in active layer")

#plt.title("Total sed content in column 15 through time")
#fig.legend()
##plt.show()

fig, ax1 = plt.subplots()

ax1.plot(time, totalSedInputThroughTime0, label="Gravel Input")
ax1.plot(time, totalSedInputThroughTime1, label="Sand Input")

ax2 = ax1.twinx()
ax2.plot(time, totalSedOutputThroughTime0, linestyle="--", label="Gravel output")
ax2.plot(time, totalSedOutputThroughTime1, linestyle="--", label="Sand output")

plt.title("Total sediment input and output through time")
fig.legend()
plt.show()

if (animateHeight):
    print("Plotting...")
    os.system("python3 animate_ISMoLD.py")
    
if (plotNodes):
    print("Plotting...")
    os.system("python3 nodeAnimation.py")
    




