import math

def setNodes(i, k, newHeight, column, dt, dx, dy, rho0):
    
    if (newHeight >= column["totalHeight"]):  ## Deposition
        
        for j in range (math.floor( column["totalHeight"] ), math.ceil( newHeight )):
            if (j == math.floor( column["totalHeight"]) ): ## Start j
                if (j==math.ceil( newHeight )-1): ## In the usual case that the first j is also the last
                    dh = newHeight - column["totalHeight"]
                else:
                    dh = math.ceil( column["totalHeight"]) - column["totalHeight"]
            else:
                if (j==math.ceil( newHeight )-1): ## If current j is final j
                    dh = newHeight-j
                else:
                    dh = dy ## not nesecary actually, but this is effectively the case when a whole node is set from 0 to rho0
            
            try:
                current_box_density = column["nodes"][j]["density"]
            except:
                current_box_density = 0
                
            if (j==math.ceil( newHeight )-1): ## Final j
                if (current_box_density+ rho0 *(dh) >= rho0): ##Final height change exeeds one node boundary:
                    column["nodes"].update({j-1:{"density":rho0}})
                    remainder = current_box_density+ rho0 *(dh) -rho0
                    column["nodes"].update({j:{"density":( remainder )}})
                else: ##Height change is within the same node:
                    remainder = current_box_density+ rho0 *(dh)
                    column["nodes"].update({j:{"density":remainder}})
            else: ##Fill whole node from 0 to full:
                column["nodes"].update({j:{"density":rho0}})
                
    else:  ## Erosion

        for j in range(math.floor(column["totalHeight"]), math.floor(newHeight)-1, -1): ## loop from top down
            
            if (j == math.floor( column["totalHeight"] )): ## Start j
                if (j==math.floor( newHeight )): ## In the usual case that the first j is also the last
                    dh = column["totalHeight"] - newHeight
                else:
                    dh = column["totalHeight"] - math.floor( column["totalHeight"] )
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
                remainder = current_box_density - rho0*( dh )
                column["nodes"].update({j:{"density":remainder}})
            else: ## Node is above the newHeight thus is to be eroded
                try:
                    del column["nodes"][j] 
                except:
                    pass
    return column
