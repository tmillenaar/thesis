import math
from numba import jit, njit
from time import time
import numpy as np
import os

start = time()

yr2sec = 60*60*24*365.25          #nr of seconds in a year

dx              = 1e3             # width of each node/column (m)
imax            = 100             # number of nodes
tmax            = 5000*yr2sec    # total amount of time to be modelled in seconds [in years = (x*yr2sec)]
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

x               = np.zeros(shape=(imax+1))
totalHeight     = np.zeros(shape=(imax+1))
bedrockHeight   = np.zeros(shape=(imax+1))
totalSedContent = np.zeros(shape=(imax+1, 2))
totalSedHeight  = np.zeros(shape=(imax+1, 2))
nodes           = np.zeros(shape=(nrOfGrainSizes+2,imax+1,1)) # For the first index: 0=density, 1=first grain size, 2=second grain size(, etc..) and the last one is depositionTimeInYears. The second index indicates the row and the third indicates the column.

totalInput= 0
totalOutput= 0
InputPerGrainSize = np.zeros(shape=(2))
OutputPerGrainSize = np.zeros(shape=(2))

dy = 1                            # vertical grid spacing (m)  Important: MUST BE 1(as integer) for code to work (in this version)!
tout = 0.                         # threshold to write next output in years (increased every write action by dtout)
tout_progress = 2*dtout_progress    # threshold to update progress bar, will be increased whenever it is reached
outputTimestep = 0

for i in range(imax+1):
    x[i]=i
    totalHeight[i] = 0
    #newHeight[i] = 0
    bedrockHeight[i] = 0
    #columns[i] = []
    for f in range(nrOfGrainSizes):
        totalSedContent[i][f] = 0
        totalSedHeight[i][f] = 0

@njit
def setTimestep():
    """ Calculating the maximum timestep (dt) as limited by the active layer and the FTCS scheme. """
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
def FTCS(dt, dx, q0, k, totalSedContent, totalInput, totalOutput, InputPerGrainSize, OutputPerGrainSize, totalHeight, bedrockHeight, subsidenceRate):
    """ Calculate the new topography profile and keep track of the total displacement in- and output of sediment in the meantime. """
    sedIn = np.zeros(shape=(2,imax+1))
    sedOut = np.zeros(shape=(2,imax+1))
    newHeight = np.zeros(shape=(imax+1))
    ## Loop through columns (for FTCS, density calculations, etc):
    totalSedIn = 0
    totalSedOut = 0
    totalHeight[imax] = 0
    for i in range(imax+1): 
        # if i==1: print("before:", totalHeight[i])
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
                    # print(i)
                    transport = ( (Dph*dt)/(dx*dx) )*( (totalHeight[i-1]) - ((totalHeight[i])) )
                    #transport = min(transport, sedContentInActiveLayer[f,i-1])
                    transport = min(transport, totalSedContent[i-1][f])
                    # print("transport",transport, totalSedContent[i-1])
                    sedOut[f,i-1] += transport
                    sedIn[f,i] += transport

                elif (totalHeight[i-1] < totalHeight[i]): ## Slope goes down to the left
                    transport = ( (Dph*dt)/(dx*dx) )*( (totalHeight[i]) - ((totalHeight[i-1])) )
                    #transport = min(transport, sedContentInActiveLayer[f,i])
                    transport = min(transport, totalSedContent[i][f])
                    sedOut[f,i] += transport
                    sedIn[f,i-1] += transport
            
        # if i==1: print(sedIn[:,1], sedOut[:,1])    
    for i in range(imax+1): 
        for f in range(nrOfGrainSizes):
            totalSedContent[i][f] = totalSedContent[i][f] + sedIn[f][i]-sedOut[f][i]
            totalHeight[i] += sedIn[f][i]-sedOut[f][i]
        if (totalHeight[i] < 0):
            totalHeight[i] = 0
        # if i==1: print("After:", totalHeight[i],"\n")
    totalHeight[imax] = 0
    for f in range(nrOfGrainSizes):
        totalSedContent[imax][f] = 0
        totalOutput += sedIn[f,imax]*dx
        OutputPerGrainSize[f] += sedIn[f,imax]*dx
        totalInput += q0[f]*dt
        InputPerGrainSize[f] += q0[f]*dt
                        
        #for i in range(imax+1): 
            #for f in range(nrOfGrainSizes):            
                #totalSedIn += sedIn[f,i]
                #totalSedOut += sedOut[f,i]
    
    subsidence (totalHeight, bedrockHeight, imax, subsidenceRate, nrOfGrainSizes)
    return totalSedContent, InputPerGrainSize, OutputPerGrainSize, totalInput, totalOutput, totalHeight, bedrockHeight, sedOut#, totalSedIn, totalSedOut

@njit
def subsidence (totalHeight, bedrockHeight, imax, subsidenceRate, nrOfGrainSizes):
    for i in range(imax+1):
        bedrockHeight[i] -= subsidenceRate*(imax-i)
        #columns[i]["oldHeight"] = max( columns[i]["oldHeight"]-subsidenceRate*(imax-i) , bedrockHeight[i])
        totalHeight[i] = max( totalHeight[i]-subsidenceRate*(imax-i) , bedrockHeight[i]) 
  
    return totalHeight, bedrockHeight


def writeOutput(imax, nrOfGrainSizes, totalHeight, totalSedContent, totalInput, InputPerGrainSize, OutputPerGrainSize, q0, sedOut, tout, dtout, outputTimestep):
    
    #Make data directory for the current timestep
    if (not os.path.isdir("ISMolD_outputdata/nodes/time"+str(outputTimestep))):
        os.mkdir("ISMolD_outputdata/nodes/time"+str(outputTimestep))
        
    for i in range(imax): 
        with open("ISMolD_outputdata/nodes/time"+str(outputTimestep)+"/column"+str(i)+".txt", "w") as file:
        #file = open("ISMolD_outputdata/nodes/time"+str(outputTimestep)+"/column"+str(i)+".txt", "w")
            jrange = math.floor(totalHeight[i]) #Set jrange to be the node at the surface
            for j in range(jrange): 
                writeline = str(j)+" "+ str(totalHeight[i])
                writeline += " "+str(bedrockHeight[i])
                for f in range(nrOfGrainSizes):
                    writeline += " "+str(nodes[f+1,i,j])
                writeline += " "+ str(nodes[len(nodes)-1,i,j])
                writeline += "\n"
                file.write(writeline)
        #file.close()
    
    outputTimestep += 1
    
    n= int(tout/dtout)
    file = open("ISMolD_outputdata/relief/topography"+str(n)+".txt", "w")
    for i in range(imax):
        writeline = str(x[i])+" "+str(totalHeight[i])
        writeline += " "+str(bedrockHeight[i])
        for f in range(nrOfGrainSizes):
            writeline += " "+str(totalSedContent[i][f] + bedrockHeight[i]) ## Note that bedrockHeight is generally negative
            #writeline += " "+str(totalSedContent[i][f])
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



def printNodes(i, nodes, totalHeight, totalSedContent):
    
    if (totalHeight != 0):
        summedGrainSizes = 0
        if (hasattr(totalSedContent, "__len__")):
            nrOfGrainSizes = len(totalSedContent)
            for f in range(nrOfGrainSizes):
                summedGrainSizes += totalSedContent[f]
        else:
            summedGrainSizes = totalSedContent
        print("Column {}, total fill:".format(i), nodes[0][i], 'summed:', summedGrainSizes)
        if (hasattr(totalSedContent, "__len__")):
            for f in range(nrOfGrainSizes):
                print("Column {}, f={}:".format(i, f), nodes[f+1][i])
            print('\n')
        else:
            for f in range(len(nodes)-1): ## -2 for fill and deposittime columns
                print("Column {}, f={}:".format(i, 'unknown'), nodes[f+1][i])
    #exit()
    
def checkNodes(nodes, dy):
    for i in range(len(nodes[0])):
        for j in range(len(nodes[0,i])):
            if ( j != len(nodes[0,i])-1 and nodes[0][i][j] < dy ): #-1 for if len(nodes[0,i])==1, j==0
                
                if (j+1 <= len(nodes[0][i])-1):
                    nextNodeContent = nodes[0][i][j+1]
                else:
                    nextNodeContent = 0
                    
                if (j-1 < 0):
                    prevNodeContent = 1.
                else:
                    prevNodeContent = nodes[0][i][j-1]    
                    
                if ( nextNodeContent != 0. ): #prevNodeContent != 1. or
                    # print('Column {}, total fill:'.format(i),nodes[0][i])
                    nextNodeFractions = list(range(nrOfGrainSizes))
                    for f in range(nrOfGrainSizes):
                        # print('Column {}, f:{}:'.format(i,f),nodes[f][i])
                        nextNodeFractions[f] = nodes[f+1][i][j+1] / nodes[0][i][j+1] 
                    # print('nextNodeFractions', nextNodeFractions)
                    deficit = dy - nodes[0][i][j]
                    # print('deficit', deficit)
                    for f in range(nrOfGrainSizes):
                        # print("Hit", deficit * nextNodeFractions[f])
                        # print('Now:', nodes[0][i][j], 'expected:', nodes[0][i][j] + deficit * nextNodeFractions[f])
                        nodes[0][i][j] += deficit * nextNodeFractions[f]
                        # print('Got:', nodes[0][i][j])
                        nodes[f+1][i][j] += deficit * nextNodeFractions[f]
                        nodes[f+1][i][j+1] -= deficit * nextNodeFractions[f]
                        # print('During: Column {}, total fill:'.format(i),nodes[0][i])
                        # print('During: Column {}, f:{}:'.format(i,f),nodes[f][i])
                    # print(nextNodeContent, nodes[0][i][j],  prevNodeContent)
                    # raise ValueError('Column {} is not properly filled'.format(i))
                    # print('After: Column {}, total fill:'.format(i),nodes[0][i])
                    # for f in range(nrOfGrainSizes):
                        # print('After: Column {}, f:{}:'.format(i,f),nodes[f][i])
    return nodes
def erodeNodes(oldSedContent, totalSedContent, bedrockHeight, oldHeight, nodes):
    yetToBeEroded = oldSedContent - (totalSedContent-bedrockHeight)
    j = math.floor(oldHeight)
    print("////", yetToBeEroded, oldSedContent, totalSedContent)
    while (yetToBeEroded > 0):
        
        currentNodeSedContent = nodes[f+1,i,j]
        print("yetToBeEroded", yetToBeEroded, currentNodeSedContent)
        if (currentNodeSedContent < yetToBeEroded): ## Remaining change exceeds one node
            nodes[f+1,i,j] = 0
            nodes[0,i,j] -= currentNodeSedContent
            yetToBeEroded -= currentNodeSedContent
            j -= 1
            if (j<0):
                raise valueError('j cannot be <0', yetToBeEroded)
        else:  ## Remaining change stays within one node
            nodes[f+1,i,j] -= yetToBeEroded
            nodes[0,i,j] -= yetToBeEroded
            yetToBeEroded = 0
    print("////", yetToBeEroded, oldSedContent, totalSedContent)    
    return nodes
		
    
    
def depositNodes(totalSedContent, bedrockHeight, currentNodeFill, oldSedContent, oldHeight, nodes, flowFraction, f, i, j, dy):
    yr2sec = 60*60*24*365.25          #nr of seconds in a year
    yetToBeDeposited = (totalSedContent - oldSedContent) - bedrockHeight
    # print(yetToBeDeposited)
    j = math.floor(oldHeight)
    while (yetToBeDeposited > 0):
        # print('j:',j, yetToBeDeposited, currentNodeFill, dy, flowFraction, (dy-currentNodeFill)*flowFraction)
        # print(nodes[0][i])
        # print(f,i,j, len(nodes[0][i]), oldHeight, totalSedContent)
        # currentNodeFill = nodes[0,i,j]currentNodeFill = nodes[0,i,j]
        print("---", j, nodes[0][i][j], currentNodeFill, dy-currentNodeFill, flowFraction, (dy-currentNodeFill)*flowFraction, yetToBeDeposited, nodes[0][i])
        nodes[len(nodes)-1,i,j] = t/yr2sec   #set time of deposition
        if ((dy-currentNodeFill)*flowFraction < yetToBeDeposited): #Remaining change exceeds one node
            nodes[0,i,j] += (dy-currentNodeFill)*flowFraction
            nodes[f+1,i,j] += (dy-currentNodeFill)*flowFraction
            yetToBeDeposited -= (dy-currentNodeFill)*flowFraction
            j+=1
            currentNodeFill = 0 #fill of next node
        else: #Remaining change stays in one node
            nodes[0,i,j] += yetToBeDeposited
            nodes[f+1,i,j] += yetToBeDeposited
            yetToBeDeposited = 0
    if nodes[0,i,j] > dy:
        print("Error, node bigger then dy")
        print(nodes[0,i])
        printNodes(i, nodes, totalHeight[i], totalSedContent)
    return nodes
        

#@njit
def setNodes(i, k, nodes, totalHeight, bedrockHeight, nrOfGrainSizes, totalSedContent, dt, dx, dy, t):
    yr2sec = 60*60*24*365.25      #nr of seconds in a year
    
    trace = True
    
    # Infer the height of the column before alteration from how filled the nodes are
    oldHeight = 0
    oldSedContent = np.zeros(shape=(nrOfGrainSizes))
    for j in range(len(nodes[0,i,:])):
        oldHeight += nodes[0,i,j]
        for f in range(nrOfGrainSizes):
            oldSedContent[f] += nodes[f+1,i,j]
            
            
    onlyErosion = True            
    onlyDeposition = True
    for f in range(nrOfGrainSizes):
        if ( totalSedContent[f]-oldSedContent[f] < 0 ): #erosion occurs
            onlyDeposition = False
        if ( totalSedContent[f]-oldSedContent[f] > 0 ): #deposition occurs
            onlyErosion = False
        
    # # Increase array size of nodes if necessary 
        # # if (int(totalHeight-bedrockHeight) is not 0): print((totalHeight-bedrockHeight), len(nodes[0,i]))
        # if ( (totalHeight-bedrockHeight) >= len(nodes[0,i])):
            # deficite = math.ceil(totalHeight-bedrockHeight) - len(nodes[0,i])
            # # print('deficite', deficite, len(nodes[0,i]), math.ceil(totalHeight-bedrockHeight))
            # newArrayShape = np.zeros(shape=(len(nodes), len(nodes[0]), deficite))# len(nodes[0,i]) + math.ceil(yetToBeDeposited)))
            # # print(np.shape(nodes), np.shape(newArrayShape))
            # nodes = np.concatenate((nodes,newArrayShape), axis=2)
            # # print('After concatenation:', np.shape(nodes))
        
    # if ((totalHeight - oldHeight) != 0): print("check tussendoor:", totalSedContent)
    if ( (totalHeight - oldHeight) == 0): ## Nothing happens
        pass
    ##------------##
    ## Deposition ##
    ##------------##
    elif ( (totalHeight-bedrockHeight - oldHeight) > 0 and onlyDeposition == True ): 
        if (trace): print("Only Deposition, column:", i)
        flowFractions = list(range(nrOfGrainSizes))
        sedContentChange  = list(range(nrOfGrainSizes))
        for f in range (nrOfGrainSizes):
            sedContentChange[f] = totalSedContent[f] - oldSedContent[f]
            
        for f in range (nrOfGrainSizes):    
            if (sum(sedContentChange) != 0):
                flowFractions[f] = sedContentChange[f]/sum(sedContentChange)
            else:
                flowFractions[f] = 0
                
        # print("Before deposit nodes:")
        # print(oldHeight, totalHeight, len(nodes[0][i]), nodes[1,i], nodes[2,i])
        currentNodeFill = nodes[0,i,math.floor(oldHeight)]
        for f in range(nrOfGrainSizes):
            nodes = depositNodes(totalSedContent[f], bedrockHeight, currentNodeFill, oldSedContent[f], oldHeight, nodes, flowFractions[f], f, i, j, dy)
        # print("After deposit nodes:")
        # print(oldHeight, totalHeight, len(nodes[0][i]), nodes[1,i], nodes[2,i])
        
        # yetToBeDeposited = ((totalHeight-bedrockHeight) - oldHeight)  ## Height in m to be added to existing column
        # flowFractions = np.zeros(shape=(nrOfGrainSizes))
        # for f in range (nrOfGrainSizes):    
            # if (yetToBeDeposited != 0):
                # flowFractions[f] = (totalSedContent[f]-oldSedContent[f])/yetToBeDeposited
            # else:
                # flowFractions[f] = 0
        
        # currentNodeSedContent = np.zeros(shape=(nrOfGrainSizes))
        # j = math.floor(oldHeight) #Set j to be the node at the surface
        # while (yetToBeDeposited > 0):
            # currentNodeFill = nodes[0,i,j]
            # nodes[len(nodes)-1,i,j] = t/yr2sec
            # for f in range(nrOfGrainSizes):
                # currentNodeSedContent[f] = nodes[f+1,i,j]
            # if ((dy-currentNodeFill) <= yetToBeDeposited): #Remaining change exceeds at least one node
                # nodes[0,i,j] = dy
                # for f in range(nrOfGrainSizes):
                    # nodes[f+1,i,j] = (dy - currentNodeSedContent[f]) * flowFractions[f]
                # yetToBeDeposited -= (dy-currentNodeFill)
                # j += 1
            # else: #Remaining change stays within the node
                # nodes[0,i,j] = currentNodeFill + yetToBeDeposited
                # nodes[len(nodes)-1,i,j] = t/yr2sec
                # for f in range(nrOfGrainSizes):
                    # nodes[f+1,i,j] = currentNodeSedContent[f] + yetToBeDeposited * flowFractions[f]
                # yetToBeDeposited = 0
            
        
    
    ###---------##
    ### Erosion ##
    ###---------##
    elif ( (totalHeight - bedrockHeight - oldHeight) < 0 and onlyErosion == True):
        if (trace): print("Only Erosion, column:", i)
        for f in range(nrOfGrainSizes):
            nodes = erodeNodes(oldSedContent[f], totalSedContent[f], bedrockHeight, oldHeight, nodes)
        
    
    ##-----------------------------##
    ## Both Deposition and Erosion ##
    ##-----------------------------##
    else: 
        if (trace): print("\n\nBoth Deposition and Erosion, column", i)
        flowFractions = list(range(nrOfGrainSizes))
        sedContentChange  = list(range(nrOfGrainSizes))
        for f in range (nrOfGrainSizes):
            sedContentChange[f] = totalSedContent[f] - oldSedContent[f]
        
        print("Before node erosion", nodes[0][i], nodes[1][i], nodes[2][i])
        for f in range (nrOfGrainSizes):    
            if (sedContentChange[f] < 0):
                nodes = erodeNodes(oldSedContent[f], totalSedContent[f], bedrockHeight, oldHeight, nodes)
        print("After node erosion", nodes[0][i], nodes[1][i], nodes[2][i])
        
        ## Make sure nodes are correctly filled (not tested)
        # for j in range(nodes[0][i].__length__):
            # if (j < dy and j != 0  and j != nodes[0][i].__length__-1):
                # missing = dy - nodes[0][i][j]
                # nextNodeFraction = list(range(nrOfGrainSizes))
                # for f in range(nrOfGrainSizes):
                    # nextNodeFraction[f] = nodes[f+1][i][j+1]/nodes[0][i][j+1]
                    # nodes[f+1][i][j+j] -= missing*nextNodeFraction[f]
                    # nodes[f+1][i][j] += missing*Height[f]
                    # nodes[0][i][j] += missing*nextNodeFraction[f]
        
        
        j = math.floor(oldHeight)
        currentNodeFill = nodes[0][i][j]
        print("currentNodeFill", currentNodeFill)
        for f in range (nrOfGrainSizes):            
            if (sedContentChange[f] > 0):
                flowFraction = 1.0 ## If nrOfGrainSizes > 2, this is no longer valid
                print("f:", f, flowFraction)
                nodes = depositNodes(totalSedContent[f], bedrockHeight, currentNodeFill, oldSedContent[f], oldHeight, nodes, flowFraction, f, i, j, dy)
        print("___  ___", nodes[0][i], totalHeight, totalSedContent)
    nodesum = 0
    for j in range(len(nodes[0][i])):
        nodesum += nodes[0][i][j]
    if (not np.allclose(nodesum, totalHeight, atol=1e-14)):
        print("Height mismatch", totalHeight, nodesum, totalHeight-nodesum, "totalSedContent", sum(totalSedContent))
        printNodes(i, nodes, totalHeight, totalSedContent)
        exit()
    return nodes

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


#def timeloop(tmax, yr2sec, nrOfGrainSizes, totalSedContent, k, q0, tout_progress):
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
        
        dt = setTimestep()
        
        #Call FTCS to obtain the new height profile
        totalSedContent, InputPerGrainSize, OutputPerGrainSize, totalInput, totalOutput, totalHeight, bedrockHeight, sedOut = FTCS(dt, dx, q0, k, totalSedContent, totalInput, totalOutput, InputPerGrainSize, OutputPerGrainSize, totalHeight, bedrockHeight, subsidenceRate)
        
        #Increase array size of nodes if necessary 
        #if (int(totalHeight-bedrockHeight) is not 0): print((totalHeight-bedrockHeight), len(nodes[0,i]))
        if ( (max(totalHeight)-min(bedrockHeight)+1) >= len(nodes[0,i])): ## +1 is to allow some room for inefficient filling, which can occur during deposition but is corrected for in checkNodes()
            deficite = math.ceil(max(totalHeight)-min(bedrockHeight)+1) - len(nodes[0,i])
            #print('deficite', deficite, len(nodes[0,i]), math.ceil(totalHeight-bedrockHeight))
            newArrayShape = np.zeros(shape=(len(nodes), len(nodes[0]), deficite))# len(nodes[0,i]) + math.ceil(yetToBeDeposited)))
            #print(np.shape(nodes), np.shape(newArrayShape))
            nodes = np.concatenate((nodes,newArrayShape), axis=2)
            #print('After concatenation:', np.shape(nodes))
    
        
        #Update the node array to match the new topography profile
        for i in range(imax+1):
            # if (totalHeight[i] != 0):
                # print('Before nodeset:')
                # printNodes(i, nodes, totalHeight[i], totalSedContent[i])
            nodes = setNodes(i, k, nodes, totalHeight[i], bedrockHeight[i], nrOfGrainSizes, totalSedContent[i], dt, dx, dy, t)
            # if (totalHeight[i] != 0):
                # print('After nodeset:')
                # printNodes(i, nodes, totalHeight[i], totalSedContent[i])
        
        nodes = checkNodes(nodes, dy)
        if (t >= tout):
            tout = writeOutput(imax, nrOfGrainSizes, totalHeight, totalSedContent, totalInput, InputPerGrainSize, OutputPerGrainSize, q0, sedOut, tout, dtout, outputTimestep)
            
        t += dt
        
        if t > tout_progress:
            print("      "+str( (math.ceil(100000*t/tmax))/1000 )+"%", end="\r") ##Track progress
            tout_progress += dtout_progress
    
    for i in range(5):
        printNodes(i, nodes, totalHeight[i], totalSedContent[i])
    printElapsedTime(time()-start)
    
    
    
    
    
## todo: check erosion, some nodes have negative sedcontent
## implement mass loss check FTCS
## implement mass loss check setNodes and children










