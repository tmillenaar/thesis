import math

def setNodes(i, k, newHeight, column, newSedContent, dt, dx, dy, rho0):
    nrOfGrainSizes = len(newSedContent)
    flowFractions = list(range(nrOfGrainSizes))
    
    #dSedContent = list(range(nrOfGrainSizes))
    
    if (newHeight >= column["oldHeight"]):  ## Deposition
        
        for p in range (nrOfGrainSizes):
            if (sum(newSedContent) != 0):
                flowFractions[p] = newSedContent[p]/sum(newSedContent)
            else:
                flowFractions[p] = 0
        #if (sum(flowFractions) != 0): print(flowFractions)
                
        depositDensity = rho0 * (newHeight - column["oldHeight"])
        
        for j in range (math.floor( column["oldHeight"] ), math.ceil( newHeight )):
            if (j == math.floor( column["oldHeight"]) ): ## Start j
                if (j==math.ceil( newHeight )-1): ## In the usual case that the first j is also the last
                    dh = newHeight - column["oldHeight"]
                    #for p in range(nrOfGrainSizes):
                        #dSedContent = column["oldSedContent"][p] - newSedContent[p]
                else:
                    dh = math.ceil( column["oldHeight"]) - column["oldHeight"]
                    #dSedContent = math.ceil( column["oldSedContent"][p] - newSedContent[p] )
            else:
                if (j==math.ceil( newHeight )-1): ## If current j is final j
                    dh = newHeight-j
                    #dSedContent = newSedContent[p] - j
                else:
                    dh = dy 
            
            try:
                current_box_density = column["nodes"][j]["density"]
                #if (current_box_density < 0): print("error, current_box_density = ", current_box_density, column["oldHeight"], newHeight )
            except:
                current_box_density = 0
                
            if (j==math.ceil( newHeight )-1): ## Final j
                #if (current_box_density + rho0 *(dh) >= rho0): ##Final height change exeeds one node boundary:
                #if (j==math.floor( column["oldHeight"] )): ##Height change is within the same node:
                remainder = current_box_density + depositDensity
                try:
                    column["nodes"][j]["density"] = remainder
                except:
                    column["nodes"].update({j:{"density":remainder}})
                for p in range(nrOfGrainSizes):
                    try:
                        remainder = column["nodes"][j]["nodeSedContent"][p] 
                    except:
                        remainder = 0
                    try:
                        column["nodes"][j]["nodeSedContent"][p] = remainder+ dh * flowFractions[p]
                    except:
                        column["nodes"][j].update({"nodeSedContent":{p:remainder+ dh * flowFractions[p]}})
                        
            else: ##Fill whole node from 0 to full:
                
                column["nodes"][j]["density"] = rho0
                depositDensity -= (rho0 - current_box_density)
                for p in range(nrOfGrainSizes):
                    try:
                        remainder = column["nodes"][j]["nodeSedContent"][p] 
                    except:
                        remainder = 0
                    try:
                        column["nodes"][j]["nodeSedContent"][p] = remainder+ dh * flowFractions[p]
                    except:
                        column["nodes"][j].update({"nodeSedContent":{p:remainder+ dh * flowFractions[p]}})
                
    else:  ## Erosion
        for j in range(math.floor(column["oldHeight"]), math.floor(newHeight)-1, -1): ## loop from top down
            
            if (j == math.floor( column["oldHeight"] )): ## Start j
                if (j==math.floor( newHeight )): ## In the usual case that the first j is also the last
                    dh = column["oldHeight"] - newHeight
                else:
                    dh = column["oldHeight"] - math.floor( column["totalHeight"] )
            else:
                if (j==math.floor( newHeight )): ## If current j is final j
                    dh = math.ceil( newHeight ) - newHeight
                else:
                    dh = dy
                    
            try:
                current_box_density = column["nodes"][j]["density"]
            except:
                current_box_density = 0
                    
            if (j==math.floor( newHeight )): ## Final j
                depositFractions = list(range(nrOfGrainSizes))
                for p in range(nrOfGrainSizes):
                    depositFractions[p] = column["nodes"][j]["nodeSedContent"][p]/sum(column["nodes"][j]["nodeSedContent"])
                remainder = current_box_density - rho0*( dh )
                column["nodes"][j]["density"] = remainder
                for p in range(nrOfGrainSizes):
                    try:
                        remainder = column["nodes"][j]["nodeSedContent"][p] 
                        column["nodes"][j]["nodeSedContent"][p] = remainder- dh * depositFractions[p]
                    except:
                        pass
            else: ## Node is above the newHeight thus is to be eroded
                try:
                    del column["nodes"][j] 
                except:
                    pass
    return column

