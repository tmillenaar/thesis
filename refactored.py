import math
from numba import jit, njit
from numba.typed import List, Dict
from time import time
import numpy as np
import math
import os

start = time()

yr2sec = 60*60*24*365.25          #nr of seconds in a year

dx              = 1e3             # width of each node/column (m)
imax            = 100             # number of nodes
tmax            = 100000*yr2sec    # total amount of time to be modelled in seconds [in years = (x*yr2sec)]
dtout           = tmax/25        # nr of years between write output
dtout_progress  = 100*yr2sec       # nr of years between progress bar update

        
nrOfGrainSizes  = 2

k = np.zeros(shape=(2))
q0 = np.zeros(shape=(2))

q0[0]= 4*1.e-7                    # Gravel input (m2/s)  Attention: will be overwritten in setBoudnaryCondtitionValues if declared there!
q0[1]= 3*1.e-7                    # Sand input (m2/s)  Attention: will be overwritten in setBoudnaryCondtitionValues if declared there!
    
k[0]= 1*3.2e-4                    # Gravel diffusivity (m2/s)  Attention: will be overwritten in setBoudnaryCondtitionValues if declared there!
k[1]= 2*3.2e-4                    # Sand diffusivity (m2/s)  Attention: will be overwritten in setBoudnaryCondtitionValues if declared there!

subsidenceRate = 0#2.e-6

nodeOutputTimestep = 0            # Used in writing node output files

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
    
    
# @njit(cache=True)
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

@njit(cache=True)
def _nbZeros_1d(length):
    l = List()
    for i in range(length):
        l.append(0.)
    return l

@njit(cache=True)
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
      
  
# @njit(cache=True)
# def _nbZeros_nd(shape):
    # l = List()
    # for i in range(len(shape)):
        # k = List()
        # for j in range(shape[i]):
            # k.append(0.)
        # l.append(k)
    # return l

# @njit(cache=True)
# def _nbZeros_1d(length):
    # l = List()
    # for i in range(length):
        # l.append(0.)
    # return l
  
@njit(cache=True)
def nbRange(length):
    l = List()
    for i in range(length):
        l.append(i)
    return l

    
def nbStackedLists(*args):
    assert all([type(val).__name__ == 'int' for val in args]), 'Please only unput integers into nbStackedLists(), got {}'.format([type(val).__name__ for val in args])
    return _nbStackedLists(args)

@njit(cache=True)
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

@njit(cache=True)
def nbmin(values):
    minVal = values[0]
    for val in values:
        if (val < minVal): minVal = val
    return minVal

@njit(cache=True)
def nbmax(values):
    maxVal = values[0]
    for val in values:
        if (val > maxVal): maxVal = val
    return maxVal

@njit(cache=True)
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

@njit(cache=True)
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

@njit(cache=True)
def sedContentInActiveLayer(nodeFill, nodeSedContent, i):
    
    maxNode = len(nodeFill[i]) - 1 ## -1 since node count starts at 0 and len() starts at 1
    fillInActiveLayer = 0
    sedFractions = np.zeros(shape=(nrOfGrainSizes))
    sedContentInActiveLayer = np.zeros(shape=(nrOfGrainSizes))
    
    for j in range(maxNode, -1, -1):
        remainder = 1-fillInActiveLayer
        if (remainder  <= 0): break
        
        if (nodeFill[i][j] > 0):
            for f in range(nrOfGrainSizes):
                sedFractions[f] = nodeSedContent[i][j][f]/nodeFill[i][j]
            
            for f in range(nrOfGrainSizes):
                toBeAdded = nbmin((remainder*sedFractions[f], nodeSedContent[i][j][f]))
                sedContentInActiveLayer[f] += toBeAdded
                fillInActiveLayer += toBeAdded
        # import pdb; pdb.set_trace()
    return sedContentInActiveLayer

@njit(cache=True)
def FTCS(totalHeight, bedrockHeight, totalSedContent, nodeFill, nodeSedContent, q0, subsidenceRate):#dt, dx, q0, k, 
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
                    transport = nbmin((transport, sedContentInActiveLayer(nodeFill, nodeSedContent, i-1)[f]))
                    transport = min(transport, totalSedContent[i-1][f])
                    sedOut[i-1][f] += transport
                    sedIn[i][f] += transport

                elif (totalHeight[i-1] < totalHeight[i]): ## Slope goes down to the left
                    transport = ( (Dph*dt)/(dx*dx) )*( (totalHeight[i]) - (totalHeight[i-1]) )
                    transport = nbmin((transport, sedContentInActiveLayer(nodeFill, nodeSedContent, i-1)[f]))
                    # transport = nbmin((transport, totalSedContent[i-1][f]))
                    sedOut[i][f] += transport
                    sedIn[i-1][f] += transport
    
    newSedContent = _nbZeros_2d((imax+1, nrOfGrainSizes))
    
    newTotalHeight = np.zeros(shape=(imax+1))
    # newBedrockHeight = np.zeros(shape=(imax+1))
    
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
        
    # newTotalHeight, newBedrockHeight= subsidence (newTotalHeight, newBedrockHeight, subsidenceRate)
    # import pdb; pdb.set_trace()
    return newTotalHeight, newSedContent

@njit(cache=True)
def subsidence (totalHeight, bedrockHeight, subsidenceRate):
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


def writeOutput(totalHeight, bedrockHeight, totalSedContent, nodeSedContent, nodeOutputTimestep, outputTimestep):
    
    print("      "+str( (math.ceil(100000*t/tmax))/1000 )+"%       Writing output {}".format(outputTimestep), end="\r") ##Update progress tracker with writer
    
    if (not os.path.isdir("ISMolD_outputdata/nodes/time"+str(nodeOutputTimestep))):
        os.mkdir("ISMolD_outputdata/nodes/time"+str(nodeOutputTimestep))
        
    for i in range(imax): ## -1 for the last column is always empty (by design). Therefore there is no need to create a file for it.
        outputFile = open("ISMolD_outputdata/nodes/time"+str(nodeOutputTimestep)+"/column"+str(i)+".txt", "w")
        jrange = len(nodeFill[i])
        for j in range(jrange): 
            writeline = str(j)+" "+ str(totalHeight[i])
            writeline += " "+str(bedrockHeight[i])
            for f in range(nrOfGrainSizes):
                writeline += " "+str(nodeSedContent[i][j][f])
            # writeline += " "+ str(columns[i]["nodes"][j]["depositionTimeInYears"]) 
            writeline += "\n"
            outputFile.write(writeline)
        outputFile.close()
    
    nodeOutputTimestep += 1
    
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
    
    return nodeOutputTimestep
    
@njit(cache=True)
def depositNodes(sedContentChange, i, f, flowFraction, nodeFill, nodeSedContent, oldNodeFills):
    trace = False
    j = math.floor(totalHeight[i]-bedrockHeight[i])
    toBeDeposited = sedContentChange
    while (toBeDeposited > 0):
        currentNodeFill = oldNodeFills[j]
        currentNodeSedContent = nodeSedContent[i][j][f]
        remainingNodeFill = 1.-currentNodeFill
        reservedFill = remainingNodeFill - remainingNodeFill*flowFraction # Reserved space in node j for other grain sizes 
        # import pdb; pdb.set_trace()
            
            ## Todo depositTime.append(_nbZeros_1d(1))
        
        ## Note: newNodeSedContent must be set before newNodeDensity, for the variable depositDensity is altered in setting newNodeDensity and must be used unaltered when setting newNodeSedContent
        ## Set new nodeSedContent:
        # if i ==0: import pdb; pdb.set_trace()
        if (currentNodeFill + toBeDeposited > 1-reservedFill): ## Cross a node boundary
            if (trace): print("ONLY DEPOSITION, cross boundary. Column:", i)
            # import pdb; pdb.set_trace()
            # nodeSedContent[i].append(_nbZeros_1d(nrOfGrainSizes))  ## TODO move one layer up. Compare length to totalHeight
            # nodeFill[i].append(0.)  ## TODO move one layer up. Compare length to totalHeight
            nodeSedContent[i][j][f] += remainingNodeFill * flowFraction
            nodeFill[i][j] += remainingNodeFill * flowFraction
            toBeDeposited -= remainingNodeFill * flowFraction
            
                # column["nodes"][j].update({"depositionTimeInYears":t/yr2sec})
        else: ## Change stays within one node:
            if (trace): print("ONLY DEPOSITION, inter node", i)
            nodeSedContent[i][j][f] += toBeDeposited# * flowFraction
            nodeFill[i][j] += toBeDeposited# * flowFraction
            toBeDeposited = 0
                
            # if(sum(current_nodeSedContent) == 0):
                # column["nodes"][j].update({"depositionTimeInYears":t/yr2sec})
        j+=1
    return nodeFill, nodeSedContent
    
@njit(cache=True)
def erodeNodes(toBeEroded, i, f, nodeFill, nodeSedContent, maxNode):
    # j = math.floor(totalHeight[i]-bedrockHeight[i])
    j = maxNode
    # if totalHeight[i] > 1: 
    # print('New maxNode:', j, 'i:', i, totalHeight[i], totalHeight[i], totalHeight[i]-bedrockHeight[i])
    while (toBeEroded > 0):
        if totalHeight[i] > 1: print('i:', i, 'j:', j, 'f:', f, 't: ', t)
        currentNodeSedContent = nodeSedContent[i][j][f]
        # import pdb; pdb.set_trace()
        if (toBeEroded > currentNodeSedContent): ## Cross a node boundary
            # if (trace): print("ONLY DEPOSITION, cross boundary. Column:", i)
            nodeSedContent[i][j][f] -= currentNodeSedContent
            nodeFill[i][j] -= currentNodeSedContent
            toBeEroded -= currentNodeSedContent
            nodeFill, nodeSedContent = checkNodes(i, nodeFill, nodeSedContent)
            # if (currentNodeSedContent < 1e-17):
                # toBeEroded = 0
            # if j==0: 
            if totalHeight[i] > 1: print('HIT', 'i:', i, 'j:', j, toBeEroded, nodeSedContent[i][j][f], currentNodeSedContent)
            
        else: ## Change stays within one node:
            # if (trace): print("ONLY DEPOSITION, inter node", i)
            nodeSedContent[i][j][f] -= toBeEroded
            nodeFill[i][j] -= toBeEroded
            break
            
        j-=1
        # if (j < 0 and toBeEroded!=0):
        if (j < 0):
            print(toBeEroded > 0)
            print("j<0, mass not eroded:", toBeEroded, 'j, i, f:', j, i, f)
            print(nodeSedContent[i])
            break
    return nodeFill, nodeSedContent

@njit(cache=True)
def checkNodes(i, nodeFill, nodeSedContent):

    ## Make sure only the top node is not fully filled
    maxNode = len(nodeFill[i]) - 1 ## -1 since node count starts at 0 and len() starts at 1
    if (maxNode == -1): maxNode = 0
    # j=maxNode
    nextNodeFraction = np.zeros(shape=(nrOfGrainSizes,))
    for j in range(maxNode):
        nextNodeContentSum = 0
        for f in range(nrOfGrainSizes):
            nextNodeContentSum += nodeSedContent[i][j+1][f] #= np.sum(nodeSedContent[i][j+1])
            
        if (nextNodeContentSum == 0): 
            return nodeFill, nodeSedContent
            # print("Error! nextNodeContentSum = 0", column) ## Should not be the case and cannot divide by 0 in next mention
            # exit()
        for f in range (nrOfGrainSizes):
            nextNodeFraction[f] = nodeSedContent[i][j+1][f]/nextNodeContentSum
        if (nodeFill[i][j] < 1 and j < maxNode):
            # nodeFill, nodeSedContent = checkNodes(i, nodeFill, nodeSedContent)
            # currentDensity = column["nodes"][j]["density"]
            # missingDensity = rho0 - column["nodes"][j]["density"]
            # print('Before: ', j, nodeSedContent[i])
            missingFill = 1 - nodeFill[i][j]
            redistributedFill = nbmin((missingFill, nodeFill[i][j+1]))
            nodeFill[i][j] += redistributedFill
            nodeFill[i][j+1] -= redistributedFill
            for f in range(nrOfGrainSizes):
                nodeSedContent[i][j][f] += redistributedFill * nextNodeFraction[f]
                nodeSedContent[i][j+1][f] -= redistributedFill * nextNodeFraction[f]
            # print('After: ', j, nodeSedContent[i])
    return nodeFill, nodeSedContent
    
@njit(cache=True)
def setNodes(totalHeight, bedrockHeight, newTotalHeight, newBedrockHeight, totalSedContent, newSedContent, nodeFill, nodeSedContent):
    trace = False
    #import pdb; pdb.set_trace()
    for i in range(len(newTotalHeight)):
        
        #for f in range(nrOfGrainSizes):
            #newSedContent = 1
        
        if ( (newTotalHeight[i] - totalHeight[i]) == 0): ## Nothing happens
            return nodeFill, nodeSedContent ## Note, this is probably false in some edge case where no change occurs halfway the slope, though this seems improbable
        
        onlyDeposition = True
        onlyErosion = True
        for f in range(nrOfGrainSizes):
            if (newSedContent[i][f] < totalSedContent[i][f]): onlyDeposition = False
            if (newSedContent[i][f] > totalSedContent[i][f]): onlyErosion = False
            
        ## Deposition ##
        if (onlyDeposition):
            
            flowFractions = np.zeros(shape=nrOfGrainSizes)
            # sedContentChange = nbStackedLists((imax+1, nrOfGrainSizes))
            sedContentChange = np.zeros(shape=(nrOfGrainSizes))
            
            for f in range (nrOfGrainSizes):
                sedContentChange[f] = (newSedContent[i][f]-newBedrockHeight[i]) - (totalSedContent[i][f]-bedrockHeight[i])
                
            for f in range (nrOfGrainSizes):    
                if (np.sum(sedContentChange) != 0):
                    flowFractions[f] = sedContentChange[f]/np.sum(sedContentChange)
                else:
                    flowFractions[f] = 0
            
            # import pdb; pdb.set_trace()
            nrNodes = len(nodeFill[i])-1
            if (nrNodes < newTotalHeight[i]-newBedrockHeight[i]):
                lackingNodeHeight = newTotalHeight[i] - newBedrockHeight[i] - nrNodes
                for j in range(math.ceil(lackingNodeHeight)):
                    nodeSedContent[i].append(_nbZeros_1d(nrOfGrainSizes))
                    nodeFill[i].append(0.)
                    
            oldNodeFills = nodeFill[i].copy()
            for f in range(nrOfGrainSizes):
                nodeFill, nodeSedContent = depositNodes(sedContentChange[f], i, f, flowFractions[f], nodeFill, nodeSedContent, oldNodeFills)
            
            # toBeDeposited = np.sum(sedContentChange)
            # j = len(nodeFill[i])-1
            # while (toBeDeposited > 0):
                # currentNodeFill = nodeFill[i][j]
                # currentNodeSedContent = nodeSedContent[i][j]
                # # import pdb; pdb.set_trace()
                    
                    # ## Todo depositTime.append(_nbZeros_1d(1))
                
                # ## Note: newNodeSedContent must be set before newNodeDensity, for the variable depositDensity is altered in setting newNodeDensity and must be used unaltered when setting newNodeSedContent
                # ## Set new nodeSedContent:
                # # if i ==0: import pdb; pdb.set_trace()
                # if (currentNodeFill + toBeDeposited > 1): ## Cross a node boundary
                    # nodeSedContent[i].append(_nbZeros_1d(nrOfGrainSizes))
                    # nodeFill[i].append(0.)
                    # remainingNodeFill = 1.-currentNodeFill
                    # for f in range(nrOfGrainSizes):
                        # newNodeSedContent = currentNodeSedContent[f] + remainingNodeFill * flowFractions[f]
                        # nodeSedContent[i][j][f] = newNodeSedContent
                    # nodeFill[i][j] = nodeFill[i][j] + remainingNodeFill
                    # toBeDeposited -= remainingNodeFill
                        
                        # # column["nodes"][j].update({"depositionTimeInYears":t/yr2sec})
                # else: ## Change stays within one node:
                    # for f in range(nrOfGrainSizes):
                        # newNodeSedContent = currentNodeSedContent[f] + toBeDeposited * flowFractions[f]
                        # nodeSedContent[i][j][f] = newNodeSedContent
                    # nodeFill[i][j] += toBeDeposited
                    # toBeDeposited = 0
                        
                    # # if(sum(current_nodeSedContent) == 0):
                        # # column["nodes"][j].update({"depositionTimeInYears":t/yr2sec})
                # j+=1
        ## Erosion ##
        elif (onlyErosion):  ## Erosion
            if (trace): print("ONLY EROSION", i)
            sedContentChange = np.zeros(shape=(nrOfGrainSizes))
        
            for f in range(nrOfGrainSizes):
                sedContentChange[f] = newSedContent[i][f] - totalSedContent[i][f] ##Negative for erosion
            
            # nodeFill, nodeSedContent = checkNodes(i, nodeFill, nodeSedContent)            
            for f in range(nrOfGrainSizes):
                if sedContentChange[f] < 0:
                    maxNode = math.floor(totalHeight[i]-bedrockHeight[i])
                    nodeFill, nodeSedContent = erodeNodes(-sedContentChange[f], i, f, nodeFill, nodeSedContent, maxNode)
    
        ## Both Deposition and Erosion ##
        else: 
            if (trace): print("Both Deposition and Erosion", i)
            flowFractions = np.zeros(shape=nrOfGrainSizes)
            # sedContentChange = nbStackedLists((imax+1, nrOfGrainSizes))
            sedContentChange = np.zeros(shape=(nrOfGrainSizes))
            
            for f in range(nrOfGrainSizes):
                sedContentChange[f] = (newSedContent[i][f]-newBedrockHeight[i]) - (totalSedContent[i][f]-bedrockHeight[i])
                
            for f in range (nrOfGrainSizes):  ## Note: This only works with two grain sizes        
                if (sedContentChange[f] > 0):
                    flowFractions[f] = 1
                else:
                    flowFractions[f] = 0
            
            ## Erode
            for f in range(nrOfGrainSizes):
                if sedContentChange[f] < 0:
                    maxNode = math.floor(totalHeight[i]-bedrockHeight[i])
                    nodeFill, nodeSedContent = erodeNodes(-sedContentChange[f], i, f, nodeFill, nodeSedContent, maxNode)
            
            # Since material has been eroded, the column has to be checked for non-full nodes
            # nodeFill, nodeSedContent = checkNodes(i, nodeFill, nodeSedContent)
            
            # import pdb; pdb.set_trace()
            nrNodes = len(nodeFill[i])-1
            if (nrNodes < newTotalHeight[i]-newBedrockHeight[i]):
                lackingNodeHeight = newTotalHeight[i] - newBedrockHeight[i] - nrNodes
                for j in range(math.ceil(lackingNodeHeight)):
                    nodeSedContent[i].append(_nbZeros_1d(nrOfGrainSizes))
                    nodeFill[i].append(0.)
            
            oldNodeFills = nodeFill[i].copy()
            for f in range(nrOfGrainSizes):
                if sedContentChange[f] > 0:
                    nodeFill, nodeSedContent = depositNodes(sedContentChange[f], i, f, flowFractions[f], nodeFill, nodeSedContent, oldNodeFills)
            
        # nodeFill, nodeSedContent = checkNodes(i, nodeFill, nodeSedContent)
    return nodeFill, nodeSedContent

def setPeriodicForcingValues(t, nrOfGrainSizes, periods, amplitudes, averages, minval="NULL"):
    if (type(periods) == list):
        if (len(periods) == 0 or len(amplitudes) == 0 or len(averages) == 0):
            print("Error, input given for setPeriodicForcingsubsidenceRate is not sufficient.")
            exit()
        if (len(periods) != len(amplitudes) or len(periods) != len(averages)):
            print("Error, the number of periods, amplitudes and averages must be the same. You supplied: nr of periods="+str(len(periods))+", nr of amplitudes="+str(len(amplitudes))+"and nr of averages="+str(len(averages))+".")
            exit()
    
    if (type(periods) != list):
        value = amplitudes*math.sin(2*3.1415 * t * 1/periods)+averages
        if (minval!="NULL"): value = max(minval, value)
    else: ## If multiple values are supplied
        value = 0
        for i in range(len(periods)):
            value += amplitudes[i]*math.sin(2*3.1415 * t * 1/periods[i])+averages[i]
        if (minval!="NULL"): value = max(minval, value)
    return value


if __name__=="__main__":
    
    print("Compiling...", end='\r')
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
        
        ###     setPeriodicForcingValues(t, nrOfGrainSizes, periods, amplitudes, averages, minval="NULL")   note: mival is optional!
        q0[0] = setPeriodicForcingValues(t, nrOfGrainSizes, [0.5*tmax], [2.0e-7] , [4.0e-7], 0)
        q0[1] = setPeriodicForcingValues(t, nrOfGrainSizes, [0.5*tmax], [1.0e-7] , [3.0e-7], 0)
    
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
        subsidenceRate = 2e-5
        #Call FTCS to obtain the new height profile
        newTotalHeight, newTotalSedContent = FTCS(totalHeight, bedrockHeight, totalSedContent, nodeFill, nodeSedContent, q0, subsidenceRate)
        
        # import pdb; pdb.set_trace()
        newTotalHeight, newBedrockHeight = subsidence(newTotalHeight, bedrockHeight, subsidenceRate)
        nodeFill, nodeSedContent = setNodes(totalHeight, bedrockHeight, newTotalHeight, newBedrockHeight, totalSedContent, newTotalSedContent, nodeFill, nodeSedContent)
        # import pdb; pdb.set_trace()
        # for i in range(imax):
            # nodeFill, nodeSedContent = checkNodes(i, nodeFill, nodeSedContent)
        totalHeight = newTotalHeight
        bedrockHeight = newBedrockHeight
        totalSedContent = newTotalSedContent
        
        if (t >= tout):
            nodeOutputTimestep = writeOutput(totalHeight, bedrockHeight, totalSedContent, nodeSedContent, nodeOutputTimestep, outputTimestep)
            tout += dtout
            # tout = writeOutput(imax, nrOfGrainSizes, totalHeight, heights['totalSedContent_'+str(f)], totalInput, InputPerGrainSize, OutputPerGrainSize, q0, sedOut, tout, dtout, outputTimestep)
            outputTimestep += 1
        t += dt
        
        if t > tout_progress:
            print("      "+str( (math.ceil(100000*t/tmax))/1000 )+"%".ljust(30), end="\r") ##Track progress
            tout_progress += dtout_progress
    
    # for i in range(5):
        # printNodes(i, nodes, totalHeight[i], heights['totalSedContent_'+str(f)][i])
    printElapsedTime(time()-start)
    import pdb; pdb.set_trace()
    










