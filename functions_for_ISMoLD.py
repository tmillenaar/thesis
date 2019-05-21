import math
import os 

def makeDirectories():
    if (not os.path.isdir("ISMolD_outputdata/relief")):
        os.mkdir("ISMolD_outputdata/relief")

    if (not os.path.isdir("ISMolD_outputdata/nodes")):
        os.mkdir("ISMolD_outputdata/nodes")


def setNodes(i, k, newHeight, column, newSedContent, dt, dx, dy, rho0):
    nrOfGrainSizes = len(newSedContent)
    
    current_nodeSedContent = list(range(nrOfGrainSizes))
    sedContentChange = list(range(nrOfGrainSizes))
    
    if ( (newHeight - column["oldHeight"]) == 0): ## Nothing happens
        return column
    
    ##------------##
    ## Deposition ##
    ##------------##
    elif ( (newHeight - column["oldHeight"]) > 0 and all( (newSedContent[q]-column["oldSedContent"][q])>0 for q in range(nrOfGrainSizes)) ): 
        
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
                
            remainingDensity = rho0-current_node_density
            remainingHeight = remainingDensity/rho0
            
            for p in range(nrOfGrainSizes):
                try:
                    current_nodeSedContent[p] = column["nodes"][j]["nodeSedContent"][p]
                except:
                    current_nodeSedContent[p] = 0
                    
            ## Note: newNodeSedContent must be set before newNodeDensity, for the variable depositDensity is altered after being used in newNodeDensity and must be used unaltered in newNodeSedContent
            ## Set new nodeSedContent:
            if (current_node_density + depositDensity > rho0): ## Cross a node boundary
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
                
            j+=1
            
            
    ##---------##
    ## Erosion ##
    ##---------##
    elif ( (newHeight - column["oldHeight"]) < 0 and all( (newSedContent[q]-column["oldSedContent"][q])<0 for q in range(nrOfGrainSizes)) ):  ## Erosion
            
        erosionFractions = list(range(nrOfGrainSizes))
        erosionContent = list(range(nrOfGrainSizes))
        
        for p in range (nrOfGrainSizes):
            sedContentChange[p] = column["oldSedContent"][p] - newSedContent[p]
        
        erosionContent = sedContentChange
            
        for p in range (nrOfGrainSizes):
            if (sum(sedContentChange) != 0):
                erosionFractions[p] = sedContentChange[p]/sum(sedContentChange)
            else:
                erosionFractions[p] = 0
        
        maxNode = len(column["nodes"]) - 1 ## -1 since node count starts at 0 and len() starts at 1
    
        if (maxNode < 0): 
            print("Error, cannot erode nonexistent material.")
            exit()
        else:
            j=maxNode
        
        erosionDensity = rho0 * ( column["oldHeight"] - newHeight ) ## Make this dependent on the nodes!!!!
        
        while (erosionDensity  > 0):
            if (j<0): print("Error, cannot remove node", j)
            ## Obtain current values:
            try:
                current_node_density = column["nodes"][j]["density"]
            except:
                current_node_density = 0
                
            remainingDensity = erosionDensity-current_node_density
            remainingHeight = remainingDensity/rho0
            
            for p in range(nrOfGrainSizes):
                try:
                    current_nodeSedContent[p] = column["nodes"][j]["nodeSedContent"][p]
                except:
                    current_nodeSedContent[p] = 0
                    print("Error, while eroding there should not be nodes without nodeSedContent.")
                    
            ## Set new nodeSedContent:
            if (current_node_density < erosionDensity ): ## Cross a node boundary
                for p in range (nrOfGrainSizes):
                    erosionContent[p] -= column["nodes"][j]["nodeSedContent"][p]
                del column["nodes"][j] 
                erosionDensity -= current_node_density
                
            else: ## Change stays within one node:
                newNodeDensity = current_node_density - erosionDensity
                column["nodes"][j]["density"] = newNodeDensity
                for p in range(nrOfGrainSizes):
                    newNodeSedContent = current_nodeSedContent[p] - erosionContent[p]
                    column["nodes"][j]["nodeSedContent"][p] = newNodeSedContent
                erosionDensity -= erosionDensity
                
            j-=1
        
    else: ## Both Deposition and Erosion
        
        ## First erode all material
        for p in range(nrOfGrainSizes):
            
            if( newSedContent[p] - column["oldSedContent"][p] < 0 ): ## If erosion 
                
                erosionDensity = abs(newSedContent[p] - column["oldSedContent"][p]) *rho0
                
                maxNode = len(column["nodes"]) - 1 ## -1 since node count starts at 0 and len() starts at 1
                if (maxNode < 0): 
                    j = 0
                else:
                    j=maxNode
                    
                while (erosionDensity > 0):
                    if (j < 0): print("Error, j cannot be negative. j:", j)
                    currentSedContent = column["nodes"][j]["nodeSedContent"][p]
                    if (currentSedContent*rho0 < erosionDensity): ## Change exeeds one node boundary
                        column["nodes"][j]["nodeSedContent"][p] = 0
                        column["nodes"][j]["density"] -= min(erosionDensity, currentSedContent*rho0)
                        erosionDensity -= min(erosionDensity, currentSedContent*rho0)
                    else: ## Change stays within same node
                        column["nodes"][j]["nodeSedContent"][p] = currentSedContent - erosionDensity/rho0
                        column["nodes"][j]["density"] -= erosionDensity
                        erosionDensity = 0
                    j-=1
            
        ## Make sure only the top node is not fully filled
        maxNode = len(column["nodes"]) - 1 ## -1 since node count starts at 0 and len() starts at 1
        if (maxNode == -1):
            maxNode = 0
        else:
            j=maxNode
        
        nextNodeFraction = list(range(nrOfGrainSizes))
        for j in range(maxNode):
            nextNodeContentSum = 0
            for p in range(nrOfGrainSizes):
                nextNodeContentSum += column["nodes"][j+1]["nodeSedContent"][p]
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
            
        ## Remaining deposition
        flowFractions = list(range(nrOfGrainSizes))
        
        for p in range (nrOfGrainSizes):
            sedContentChange[p] = newSedContent[p] - column["oldSedContent"][p]
            
        totalDepositionDensity = 0
        for p in range (nrOfGrainSizes):
            if (sedContentChange[p] > 0):
                totalDepositionDensity += sedContentChange[p]*rho0
                
        for p in range (nrOfGrainSizes):    
            if (sedContentChange[p] > 0):
                flowFractions[p] = sedContentChange[p]*rho0/totalDepositionDensity
            else:
                flowFractions[p] = 0
        
        maxNode = len(column["nodes"]) - 1 ## -1 since node count starts at 0 and len() starts at 1
        if (maxNode == -1):
            maxNode = 0
        else:
            j=maxNode
            
        depositDensity = totalDepositionDensity
        while (depositDensity > 0):
            
            ## Obtain current values:
            try:
                current_node_density = column["nodes"][j]["density"]
            except:
                current_node_density = 0
                
            remainingDensity = rho0-current_node_density
            remainingHeight = remainingDensity/rho0
            
            for p in range(nrOfGrainSizes):
                try:
                    current_nodeSedContent[p] = column["nodes"][j]["nodeSedContent"][p]
                except:
                    current_nodeSedContent[p] = 0
            if (current_node_density + depositDensity > rho0): ## Cross a node boundary
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
                    
            j+=1
        
        
    return column


