import math
import os 
from decimal import Decimal

trace = False
#trace = True ## Keep track in the terminal wether there was deposition, erosion or both
    
def makeDirectories():
    if (not os.path.isdir("ISMolD_outputdata")):
        os.mkdir("ISMolD_outputdata")
        
    if (not os.path.isdir("ISMolD_outputdata/relief")):
        os.mkdir("ISMolD_outputdata/relief")

    if (not os.path.isdir("ISMolD_outputdata/nodes")):
        os.mkdir("ISMolD_outputdata/nodes")
        
def makeTimeNodeDirectory(timestep) :
    if (not os.path.isdir("ISMolD_outputdata/nodes/time"+str(timestep))):
        os.mkdir("ISMolD_outputdata/nodes/time"+str(timestep))
        
def printColumn(column, i, newHeight, newSedContent):
    
    maxNode = len(column["nodes"]) - 1 ## -1 since node count starts at 0 and len() starts at 1
    if (maxNode == -1): maxNode = 0
        
    columnDensity = 0
    newTotalSedContent = [0,0] ## Change if nrOfGrainSizes is changed, now hardcoded
    for j in range(maxNode+1):
        columnDensity += column["nodes"][j]["density"]
        for p in range(2): ## hardcoded
            newTotalSedContent[p] += column["nodes"][j]["nodeSedContent"][p]
    newTotalHeight = columnDensity/2700 ## hardcoded
    
    print("Column "+str(i)+":")
    print("totalHeight:", column["totalHeight"])
    print("totalSedContent:", column["totalSedContent"])
    for j in range(len(column["nodes"])):
        print(j, column["nodes"][j])
    print("Summed:  ", Decimal(newTotalHeight), Decimal(newTotalSedContent[0]), Decimal(newTotalSedContent[1]))
    print("Expected:", Decimal(newHeight), Decimal(newSedContent[0]), Decimal(newSedContent[1]))
    print("")






def setNodes(i, k, newHeight, column, newSedContent, dt, dx, dy, rho0):
    
    #if (i==5): 
        #trace = False
    #else:
        #trace = False
    
    nrOfGrainSizes = len(newSedContent)
    
    current_nodeSedContent = list(range(nrOfGrainSizes))
    sedContentChange = list(range(nrOfGrainSizes))
    
    maxNode = len(column["nodes"]) - 1 ## -1 since node count starts at 0 and len() starts at 1
    totalSedContentBefore = list(range(nrOfGrainSizes))
    for p in range(nrOfGrainSizes):
        totalSedContentBefore[p] = 0
        
    if (maxNode == -1): 
        maxNode = 0
    else:
        for j in range(maxNode+1):
            for p in range(nrOfGrainSizes):
                totalSedContentBefore[p] += column["nodes"][j]["nodeSedContent"][p]
            
    if ( (newHeight - column["oldHeight"]) == 0): ## Nothing happens
        return column
    
    
    
    
    
    ##------------##
    ## Deposition ##
    ##------------##
    elif ( (newHeight - column["oldHeight"]) > 0 and all( (newSedContent[q]-column["oldSedContent"][q])>=0 for q in range(nrOfGrainSizes)) ): 
        if (trace): print("ONLY DEPOSITION", i)
        
        flowFractions = list(range(nrOfGrainSizes))
        
        for p in range (nrOfGrainSizes):
            sedContentChange[p] = newSedContent[p] - column["oldSedContent"][p]
            
        for p in range (nrOfGrainSizes):    
            if (sum(sedContentChange) != 0):
                flowFractions[p] = sedContentChange[p]/sum(sedContentChange)
            else:
                flowFractions[p] = 0
        
        maxNode = len(column["nodes"]) - 1 ## -1 since node count starts at 0 and len() starts at 1
        if (maxNode < 0): 
            j=0
        else:
            j=maxNode
        
        depositDensity = rho0 * (newHeight - column["oldHeight"])
        
        while (depositDensity > 0):
            
            ## Obtain current values:
            try:
                current_node_density = column["nodes"][j]["density"]
            except:
                current_node_density = 0
                
            for p in range(nrOfGrainSizes):
                try:
                    current_nodeSedContent[p] = column["nodes"][j]["nodeSedContent"][p]
                except:
                    current_nodeSedContent[p] = 0
            
            if (current_node_density > rho0):
                print("Error, node too dense!", "j:", j, "i:", i)
                printColumn(column, i, newHeight, newSedContent)
                exit()
                    
            ## Note: newNodeSedContent must be set before newNodeDensity, for the variable depositDensity is altered after being used in newNodeDensity and must be used unaltered in newNodeSedContent
            ## Set new nodeSedContent:
            if (current_node_density + depositDensity > rho0): ## Cross a node boundary
                remainingDensity = rho0-current_node_density
                remainingHeight = remainingDensity/rho0
                for p in range(nrOfGrainSizes):
                    newNodeSedContent = current_nodeSedContent[p] + remainingHeight * flowFractions[p]
                    try:
                        column["nodes"][j]["nodeSedContent"][p] = newNodeSedContent
                    except:
                        try:
                            column["nodes"][j]["nodeSedContent"].update({p: newNodeSedContent })
                        except:
                            try:
                                column["nodes"][j].update({"nodeSedContent":{p: newNodeSedContent }})
                            except:
                                column["nodes"].update({j:{"nodeSedContent":{p: newNodeSedContent }}})
                    if (newNodeSedContent < 0):
                        print("newNodeSedContent < 0:", i, p, newNodeSedContent, current_nodeSedContent[p], remainingHeight, flowFractions[p], current_node_density)
                        exit()
            else: ## Change stays within one node:
                newNodeDensity = current_node_density + depositDensity
                for p in range(nrOfGrainSizes):
                    newNodeSedContent = current_nodeSedContent[p] + depositDensity/rho0 * flowFractions[p]
                    try:
                        column["nodes"][j]["nodeSedContent"][p] = newNodeSedContent
                    except:
                        try:
                            column["nodes"][j]["nodeSedContent"].update({p: newNodeSedContent })
                        except:
                            try:
                                column["nodes"][j].update({"nodeSedContent":{p: newNodeSedContent }})
                            except:
                                column["nodes"].update({j:{"nodeSedContent":{p: newNodeSedContent }}})
                    if (newNodeSedContent < 0):
                        print(newNodeSedContent)
                        exit()
            ## Set new node density:
            if (current_node_density + depositDensity > rho0):
                newNodeDensity = rho0
                depositDensity -= remainingDensity
            else:
                newNodeDensity = current_node_density + depositDensity
                depositDensity -= depositDensity
            
            try:
                column["nodes"][j]["density"] = newNodeDensity
            except:
                try:
                    column["nodes"][j].update({"density":newNodeDensity})
                except:
                    column["nodes"].update({j:{"density":newNodeDensity}})
            if (column["nodes"][j]["density"] > rho0):
                print("Error 1, node density too large:")
                printColumn(column, i, newHeight, newSedContent)
                exit()
            j+=1
            
            if (newNodeDensity > rho0):
                print("Error, newNodeDensity too high:", newNodeDensity)
        
        maxNode = len(column["nodes"]) - 1 ## -1 since node count starts at 0 and len() starts at 1
        if (maxNode < 0): maxNode = 0
        for j in range(maxNode):
            if (column["nodes"][j]["density"] > rho0):
                print("Error 2, node density too large:")    
                printColumn(column, i, newHeight, newSedContent)
                exit()
                
                
                
                
                
    ##---------##
    ## Erosion ##
    ##---------##
    elif ( (newHeight - column["oldHeight"]) < 0 and all( (newSedContent[q]-column["oldSedContent"][q])<=0 for q in range(nrOfGrainSizes)) ):  ## Erosion
        if (trace): print("ONLY EROSION", i)
        #erosionFractions = list(range(nrOfGrainSizes))
        erosionContent = list(range(nrOfGrainSizes))
        
        for p in range (nrOfGrainSizes):
            sedContentChange[p] = column["oldSedContent"][p] - newSedContent[p]
        
        erosionContent = sedContentChange ## list of len nrOfGrainSizes
            
        #for p in range (nrOfGrainSizes):
            #if (sum(sedContentChange) != 0):
                #erosionFractions[p] = sedContentChange[p]/sum(sedContentChange)
            #else:
                #erosionFractions[p] = 0
        
        maxNode = len(column["nodes"]) - 1 ## -1 since node count starts at 0 and len() starts at 1
    
        if (maxNode < 0): 
            print("Error, cannot erode nonexistent material.")
            exit()
        else:
            j=maxNode
        
        erosionDensity = rho0 * ( column["oldHeight"] - newHeight ) ## Make this dependent on the nodes!!!!
    
        if (j<0): print("Error, cannot remove node", j)
        ## Obtain current values:
        try:
            current_node_density = column["nodes"][j]["density"]
        except:
            current_node_density = 0
            
        #remainingDensity = erosionDensity-current_node_density
        #remainingHeight = remainingDensity/rho0
        
        ## Set new nodeSedContent:
        #print(j, column)
        #print("")
        for p in range(nrOfGrainSizes):
            erodeNodes(p, column, newSedContent, rho0)
            
        nodeContentSum = 0
        for p in range(nrOfGrainSizes):
            nodeContentSum += column["nodes"][maxNode]["nodeSedContent"][p]
        if ( column["nodes"][maxNode]["density"] <= 0 or nodeContentSum <= 0):
            #print("Deleting", maxNode, "from", column["nodes"])
            del column["nodes"][maxNode]
            
        ## Make sure only the top node is not fully filled
        maxNode = len(column["nodes"]) - 1 ## -1 since node count starts at 0 and len() starts at 1
        #maxNode = math.floor(column["oldHeight"])
        if (maxNode == -1):
            maxNode = 0
        else:
            j=maxNode
        if (maxNode > 0):
            nextNodeFraction = list(range(nrOfGrainSizes))
            for j in range(maxNode):
                nextNodeContentSum = 0
                for p in range(nrOfGrainSizes):
                    nextNodeContentSum += column["nodes"][j+1]["nodeSedContent"][p]
                    
                if (nextNodeContentSum == 0): print("Error! nextNodeContentSum = 0", column)
                for p in range (nrOfGrainSizes):
                    nextNodeFraction[p] = column["nodes"][j+1]["nodeSedContent"][p]/nextNodeContentSum
                if (column["nodes"][j]["density"] < rho0 and j < maxNode):
                    currentDensity = column["nodes"][j]["density"]
                    missingDensity = rho0 - column["nodes"][j]["density"]
                    movedDensity = min(missingDensity, column["nodes"][j+1]["density"])
                    column["nodes"][j]["density"] += movedDensity 
                    column["nodes"][j+1]["density"] -= movedDensity 
                    for p in range(nrOfGrainSizes):
                        column["nodes"][j]["nodeSedContent"][p] += (movedDensity/rho0) * nextNodeFraction[p]
                        column["nodes"][j+1]["nodeSedContent"][p] -= (movedDensity/rho0) * nextNodeFraction[p]
            if (len(column["nodes"]) > 0):
                if (column["nodes"][maxNode]["density"] <= 0):
                    del column["nodes"][maxNode]
            else:
                print("Error, node "+str(maxNode)+" in column "+str(i)+" cannot be deleted for it does not exist", column)
            
            
            
            
            
    ##-----------------------------##
    ## Both Deposition and Erosion ##
    ##-----------------------------##
    else: 
        if (trace): print("Both Deposition and Erosion", i)

        ## First erode all material
        erosionContent = list(range(nrOfGrainSizes))
        
        erosionDensity = rho0 * ( column["oldHeight"] - newHeight ) ## Make this dependent on the nodes!!!!
        
        for p in range (nrOfGrainSizes):
            sedContentChange[p] = newSedContent[p] - column["oldSedContent"][p]
                
        for p in range(nrOfGrainSizes):
            if( newSedContent[p] - column["oldSedContent"][p] < 0 ): ## If erosion for this grain size (p)
                erodeNodes(p, column, newSedContent, rho0)
            
        ## Remove empty nodes:
        maxNode = len(column["nodes"]) - 1 ## -1 since node count starts at 0 and len() starts at 1
        if (maxNode == -1): maxNode = 0
            
        nodeContentSum = 0
        for p in range(nrOfGrainSizes):
            nodeContentSum += column["nodes"][maxNode]["nodeSedContent"][p]
        if ( column["nodes"][maxNode]["density"] <= 0 or nodeContentSum <= 0):
            del column["nodes"][maxNode]
            
            
        
        ## Make sure only the top node is not fully filled
        maxNode = len(column["nodes"]) - 1 ## -1 since node count starts at 0 and len() starts at 1
        if (maxNode == -1): maxNode = 0
        j=maxNode
        
        nextNodeFraction = list(range(nrOfGrainSizes))
        for j in range(maxNode):
            nextNodeContentSum = 0
            for p in range(nrOfGrainSizes):
                nextNodeContentSum += column["nodes"][j+1]["nodeSedContent"][p]
                
            if (nextNodeContentSum == 0): 
                print("Error! nextNodeContentSum = 0", column) ## Should not be the case and cannot divide by 0 in next mention
                exit()
            for p in range (nrOfGrainSizes):
                nextNodeFraction[p] = column["nodes"][j+1]["nodeSedContent"][p]/nextNodeContentSum
            if (column["nodes"][j]["density"] < rho0 and j < maxNode):
                currentDensity = column["nodes"][j]["density"]
                missingDensity = rho0 - column["nodes"][j]["density"]
                movedDensity = min(missingDensity, column["nodes"][j+1]["density"])
                column["nodes"][j]["density"] += movedDensity 
                column["nodes"][j+1]["density"] -= movedDensity 
                for p in range(nrOfGrainSizes):
                    column["nodes"][j]["nodeSedContent"][p] += (movedDensity/rho0) * nextNodeFraction[p]
                    column["nodes"][j+1]["nodeSedContent"][p] -= (movedDensity/rho0) * nextNodeFraction[p]
        
        ## Remove empty nodes:
        nodeContentSum = 0
        for p in range(nrOfGrainSizes):
            nodeContentSum += column["nodes"][maxNode]["nodeSedContent"][p]
        if ( column["nodes"][maxNode]["density"] <= 0 or nodeContentSum <= 0):
            del column["nodes"][maxNode]
                
                
                
        ## Deposit remaining material
        depositionContent = list(range(nrOfGrainSizes))
        current_nodeSedContent = list(range(nrOfGrainSizes))
        flowFractions = list(range(nrOfGrainSizes))
        
        maxNode = len(column["nodes"]) - 1 ## -1 since node count starts at 0 and len() starts at 1
        if (maxNode == -1): maxNode = 0
        try:
            current_node_density = column["nodes"][maxNode]["density"]
        except:
            current_node_density = 0
        
        for p in range(nrOfGrainSizes):
            try:
                current_nodeSedContent[p] = column["nodes"][maxNode]["nodeSedContent"][p]
            except:
                current_nodeSedContent[p] = 0
            if (sedContentChange[p] > 0):
                depositionContent[p] = sedContentChange[p] 
            else:
                depositionContent[p] = 0
            
        initialDepositionContent = depositionContent.copy() ## Is list with length nrOfGrainSizes
        
        for p in range (nrOfGrainSizes):
            if (sedContentChange[p] > 0 and sum(depositionContent) != 0):
                flowFractions[p] = sedContentChange[p]/sum(depositionContent)
            else:
                flowFractions[p] = 0
                
        maxNode = len(column["nodes"]) - 1 ## Start maxNode before grain size loop (index p) for if the first iteration makes a new node, it should not directly become the new maxNode for it would leave the "old" maxNode underfilled.
        if (maxNode == -1): maxNode = 0
        
        for p in range(nrOfGrainSizes):
            if (current_node_density + sum(initialDepositionContent)*rho0 > rho0 ): ## Cross a node boundary
                #print("Cross a node boundary")
                
                ## Fill existing node:
                remainder = 1 - current_node_density/rho0
                newNodeSedContent = current_nodeSedContent[p] + remainder * flowFractions[p]
                depositionContent[p] -= remainder * flowFractions[p]
                try:
                    column["nodes"][maxNode]["nodeSedContent"][p] = newNodeSedContent
                except:
                    try:
                        column["nodes"][maxNode]["nodeSedContent"].update({p: newNodeSedContent })
                    except:
                        try:
                            column["nodes"][maxNode].update({"nodeSedContent":{p: newNodeSedContent }})
                        except:
                            column["nodes"].update({maxNode:{"nodeSedContent":{p: newNodeSedContent }}})
                
                try:
                    column["nodes"][maxNode]["density"] = rho0
                except:
                    column["nodes"][maxNode].update({"density":rho0})
                
                ## Set new node:
                newNodeSedContent = depositionContent[p]
                try:
                    column["nodes"][maxNode+1]["nodeSedContent"][p] = newNodeSedContent
                except:
                    try:
                        column["nodes"][maxNode+1]["nodeSedContent"].update({p: newNodeSedContent })
                    except:
                        try:
                            column["nodes"][maxNode+1].update({"nodeSedContent":{p: newNodeSedContent }})
                        except:
                            column["nodes"].update({maxNode+1:{"nodeSedContent":{p: newNodeSedContent }}})
                
                try:
                    column["nodes"][maxNode+1]["density"] += depositionContent[p]*rho0
                except:
                    try:
                        column["nodes"][maxNode+1].update({"density":depositionContent[p]*rho0})
                    except:
                        column["nodes"].update({maxNode+1:{"density":depositionContent[p]*rho0}})
                    
            else: ## Change stays within one node:
                #print("Change within 1 layer!")
                newNodeSedContent = current_nodeSedContent[p] + depositionContent[p]
                try:
                    column["nodes"][maxNode]["nodeSedContent"][p] = newNodeSedContent
                except:
                    print("Error! This node should exist.", "p:", p, "maxNode", maxNode, column)
                    exit()

                try:
                    column["nodes"][maxNode]["density"] += depositionContent[p]*rho0
                except:
                    column["nodes"][maxNode].update({"density":depositionContent[p]*rho0})
                if (column["nodes"][maxNode]["density"] > rho0):
                    print("Error 3, node density too large:")
                    printColumn(column, i, newHeight, newSedContent)
                    exit()
    
        maxNode = len(column["nodes"]) - 1 ## Start maxNode before p loop for if the first iteration made a new node, it may not directly become the new maxNode for it would leave the "old" maxNode underfilled.
        if (maxNode == -1): maxNode = 0
        
        if (column["nodes"][maxNode]["density"] > rho0):
            print("Error 4, node density too large:", column)
            printColumn(column, i, newHeight, newSedContent)
            exit()
    
    
    ## Wrapping up
    maxNode = len(column["nodes"]) - 1 ## -1 since node count starts at 0 and len() starts at 1
    totalSedContentAfter = list(range(nrOfGrainSizes))
    for p in range(nrOfGrainSizes):
        totalSedContentAfter[p] = 0
        
    if (maxNode == -1): 
        maxNode = 0
    else:
        for j in range(maxNode+1):
            for p in range(nrOfGrainSizes):
                totalSedContentAfter[p] += column["nodes"][j]["nodeSedContent"][p]
            
    actualSedContentChange = list(range(nrOfGrainSizes))
    difference = list(range(nrOfGrainSizes))
    for p in range(nrOfGrainSizes):
        actualSedContentChange[p] = totalSedContentAfter[p]-totalSedContentBefore[p]
        difference[p] = actualSedContentChange[p] - sedContentChange[p]
        
        
    if (sum(difference) > 1e-12): 
        print(totalSedContentBefore, totalSedContentAfter)
        print("i: ",i, "Total sed content change: Actual:", actualSedContentChange, "intended:", sedContentChange, "difference", difference)
        exit()
    totalDeposited = 0
    totalEroded = 0
    for p in range(nrOfGrainSizes):
        if (sedContentChange[p] >= 0):
            totalDeposited += sedContentChange[p]
        else:
            totalEroded += sedContentChange[p]
    
    columnDensity = 0
    newTotalSedContent = list(range(nrOfGrainSizes))
    for p in range(nrOfGrainSizes):
        newTotalSedContent[p] = 0
        
    maxNode = len(column["nodes"])-1
    if (maxNode < 0): maxNode = 0
    
    if (len(column["nodes"]) > 0):
        for j in range(maxNode+1):
            columnDensity += column["nodes"][j]["density"]
            if(column["nodes"][j]["density"] < rho0 and j < maxNode): 
                print("Error, column not properly filled:")
                printColumn(column, i, newHeight, newSedContent)
                exit()
            for p in range(nrOfGrainSizes):
                newTotalSedContent[p] += column["nodes"][j]["nodeSedContent"][p]
        newTotalHeight = columnDensity/rho0
        
        #if (i==5): print("Height Mismatch column "+str(i)+":", newHeight - newTotalHeight)
        if (newHeight-column["bedrockHeight"] - newTotalHeight < -1e-10 or newHeight-column["bedrockHeight"] - newTotalHeight > 1e-10):
            print("")
            print("totalHeight Mismatch:", newHeight - newTotalHeight, i, sedContentChange)
            printColumn(column, i, newHeight, newSedContent)
            exit()
        elif(newSedContent[0] - newTotalSedContent[0] < -1e-10 or newSedContent[0] - newTotalSedContent[0] > 1e-10):
            print("")
            print("sedContent 0 Mismatch:", newSedContent[0] - newTotalSedContent[0], sedContentChange, i)
            exit()
        elif(newSedContent[1] - newTotalSedContent[1] < -1e-10 or newSedContent[1] - newTotalSedContent[1] > 1e-10):
            print("")
            print("sedContent 1 Mismatch:", newSedContent[1] - newTotalSedContent[1], sedContentChange, i)
            exit()
        
    maxNode = len(column["nodes"]) - 1 
    if (maxNode == -1): maxNode = 0
    
    sedContentSum = 0
    for j in range(maxNode):
        sedContentSum += column["nodes"][j]["density"]
        if(column["nodes"][j]["density"] < rho0 and j < maxNode): 
            print("Error, column not properly filled:")
            printColumn(column, i, newHeight, newSedContent)
            #exit()
        if (column["nodes"][j]["density"] > rho0):
            print("Error 5, node density too large:")
            printColumn(column, i, newHeight, newSedContent)
            exit()
            
        for p in range(nrOfGrainSizes):
            if(column["nodes"][j]["nodeSedContent"][p] < 0):
                print("Error, something < 0")
                printColumn(column, i, newHeight, newSedContent)
                exit()
                
    if(len(column["nodes"]) == 0):
        return column
    elif (len(column["nodes"][maxNode]["nodeSedContent"]) < 2):
        print("maxNode", maxNode)
        print("Major error, not enough grain sizes in nodeSedContent!:", column)
        exit()
    else:
        return column
    




def erodeNodes(p, column, newSedContent, rho0):
    sedContentChange = column["oldSedContent"][p] - newSedContent[p]
    erosionContent = sedContentChange 
    
    before = 0
    for j in range(len(column["nodes"])):
        before += column["nodes"][j]["density"]/2700
    
    maxNode = len(column["nodes"]) - 1 ## -1 since node count starts at 0 and len() starts at 1

    if (maxNode < 0): 
        print("Error, cannot erode nonexistent material.")
        exit()
    else:
        j=maxNode
        
    #print("")
    while (erosionContent  > 0):
        if (j < 0):
            print("Error, j cannot be negative. j: "+str(j), "erosionContent:", erosionContent, column)
            exit()
        #print("erosionContent", erosionContent[p], p, j)
        try:
            current_nodeSedContent = column["nodes"][j]["nodeSedContent"][p]
        except:
            current_nodeSedContent = 0
            print("Error, while eroding there should not be nodes without nodeSedContent.", "j: "+str(j), column)
            exit()
            
        if (erosionContent > current_nodeSedContent ): ## Cross a node boundary
            if (trace): print("Erosion crossing boundary")
            erosionContent -= current_nodeSedContent 
            column["nodes"][j]["density"] -= current_nodeSedContent * rho0
            if (erosionContent < 1e-20): erosionContent = 0
            column["nodes"][j]["nodeSedContent"][p] = 0
        else: ## Change stays within one node:
            if (trace): print("Erosion within 1 node")
            newNodeSedContent = current_nodeSedContent - erosionContent
            #if (trace): print("Erosion within 1 node", "newNodeSedContent:", Decimal(newNodeSedContent))
            #column["nodes"][j]["nodeSedContent"][p] = newNodeSedContent
            column["nodes"][j]["nodeSedContent"][p] -= erosionContent
            column["nodes"][j]["density"] -= erosionContent*rho0
            if (newNodeSedContent < 0): 
                print("Error, newNodeSedContent is negative:", newNodeSedContent)
                exit()
            erosionContent = 0
        j-=1
        
    after = 0
    for j in range(len(column["nodes"])):
        after += column["nodes"][j]["density"]/2700
    #print("Actual change:", Decimal(before-after), "Inteded:", Decimal(sedContentChange), "Difference:", (before-after)-sedContentChange)
    return column





def subsidence (columns, imax, subsidenceRate, nrOfGrainSizes):
    for i in range(imax+1):
        columns[i]["bedrockHeight"] = columns[i]["bedrockHeight"]-subsidenceRate*(imax-i)
        columns[i]["oldHeight"] = max( columns[i]["oldHeight"]-subsidenceRate*(imax-i) , columns[i]["bedrockHeight"])
        columns[i]["totalHeight"] = max( columns[i]["totalHeight"]-subsidenceRate*(imax-i) , columns[i]["bedrockHeight"]) 
  
    return columns

