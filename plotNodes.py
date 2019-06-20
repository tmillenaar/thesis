import os
import math         ## Used for floor() and ceil() to obtain integers from floats
import warnings     ## Used to filter numpy warnings when opening empty text/data files
import numpy as np  
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation 
from os import listdir

saveImage = False
showAnimation = False

#saveImage = True
showAnimation = True


detailFactor = 20 ## Mesh size for the plot. Higher detailFactor creates a smoother bedrock slope. Example: if detailFactor = 10, each mesh cell is divided into 10x10 subcells


if (detailFactor <= 0): detailFactor = 1 ## Though I cannot imagine someone wants negative detail ^^' 

if (not os.path.isdir("ISMolD_outputdata/nodes/time0")):
    print("Error, ISMolD_outputdata/nodes/time0 does not exist.")
    exit()

print("Determining the extent of the data...", end="\r")
upperBound = 0
lowerBound = 0
max_timestep = len(os.listdir("ISMolD_outputdata/relief"))-1 #-1 since file starts at 'time0'
nrColumns = len(listdir("ISMolD_outputdata/nodes/time"+str(max_timestep)))-1 #-1 since file starts at 'output0', add -1 for amy subdirectory
for i in range (nrColumns):
    with warnings.catch_warnings(): 
        warnings.simplefilter("ignore") ## Ignores the warnings when numpy opens an empty file (most columns are empty/unfilled at the start of the run)
        data = np.loadtxt("ISMolD_outputdata/nodes/time"+str(max_timestep)+"/column"+str(i)+".txt")
    upperBound = max(upperBound, int(data.size/5))
    try:
        lowerBound = min(lowerBound, min(data[:,2]))
    except:
        pass

lowerBound = math.floor(lowerBound)        

heatmapData = np.zeros(shape=(detailFactor*upperBound-detailFactor*lowerBound, detailFactor*nrColumns))

if (upperBound == 0):
    print("Error, no data")
    exit()

print("Preparing the image...                          ", end="\r") ## These spaces are to overwrite the whole of the previous line

for i in range(nrColumns):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        columnData = np.loadtxt("ISMolD_outputdata/nodes/time"+str(max_timestep)+"/column"+str(i)+".txt") ## nodeID, totalHeight, bedrockHeight, gravel, sand
    for j in range(len(columnData)):
        for l in range(detailFactor):
            xmodifier = detailFactor*i+l
            for k in range(detailFactor):
                ymodifier = detailFactor*j+k-detailFactor*lowerBound
                if (int(columnData.size/4) > 1):
                    heatmapData[ymodifier+math.ceil(detailFactor*columnData[j,2]),xmodifier] = columnData[j,3]
                elif(int(columnData.size/4) == 1): 
                    heatmapData[k-detailFactor*lowerBound+math.ceil(detailFactor*columnData[2]),xmodifier] = columnData[3]
            
fig, ax = plt.subplots()
im = ax.imshow(heatmapData, vmin=0, vmax=1, aspect="equal", extent=[0,nrColumns,upperBound,lowerBound]) ## , extent=[80,120,32,upperBound]
fig.colorbar(im).set_label("Gravel fraction")
ax.set_title("Sediment content per node")
ax.invert_yaxis()

#if (saveImage):
    ##anim.save('sedimentContentInNodes.gif', dpi=200, writer='imagemagick')
    #os.system("xdg-open sedimentContentInNodes.gif")
    
print("Done", end="\r")
if (showAnimation):
    print("Showing figure...    ")
    plt.show()
