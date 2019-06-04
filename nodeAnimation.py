import os
import math         ## Used for floor() and ceil() to obtain integers from floats
import warnings     ## Used to filter numpy warnings when opening empty text/data files
import numpy as np  
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation 
from os import listdir

saveImage = False
showAnimation = False

saveImage = True
#showAnimation = True



detailFactor = 20
animationInterval = 40 ## Time interval between frames in ms


if (detailFactor <= 0): detailFactor = 1 ## Though I cannot imagine someone wants negative detail ^^' 

if (not os.path.isdir("ISMolD_outputdata/nodes/time0")):
    print("Error, directory ISMolD_outputdata/nodes/time0 does not exist.")
    exit()

print("Determining the extent of the data...", end="\r")
upperBound = 0
lowerBound = 0
max_timestep = len(os.listdir("ISMolD_outputdata/relief"))-1 #-1 since file starts at 'time0'
for t in range(max_timestep):
    nrColumns = len(listdir("ISMolD_outputdata/nodes/time"+str(t)))-1 #-1 since file starts at 'output0', add -1 for amy subdirectory

    for i in range (nrColumns):
        with warnings.catch_warnings(): 
            warnings.simplefilter("ignore") ## Ignores the warnings when numpy opens an empty file (most columns are empty/unfilled at the start of the run)
            data = np.loadtxt("ISMolD_outputdata/nodes/time"+str(t)+"/column"+str(i)+".txt")
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
        columnData = np.loadtxt("ISMolD_outputdata/nodes/time0/column"+str(i)+".txt") ## nodeID, totalHeight, bedrockHeight, gravel, sand
    
    for j in range(len(columnData)):
        for l in range(detailFactor):
            if (int(columnData.size/5) > 1):
                heatmapData[j-detailFactor*lowerBound+math.ceil(columnData[j,2]),i] = columnData[j,3]
            elif(int(columnData.size/5) == 1): 
                heatmapData[detailFactor*lowerBound+math.ceil(columnData[2]),i] = columnData[3]
            
fig, ax = plt.subplots()
im = ax.imshow(heatmapData, vmin=0, vmax=1, aspect="equal", extent=[0,nrColumns,upperBound,lowerBound]) ## , extent=[80,120,32,upperBound]
fig.colorbar(im).set_label("Gravel fraction")
ax.set_title("Sediment content per node")
ax.invert_yaxis()

def update(t):
    label = 't = {0}0 kyr'.format(t)
    #print(label)
    
    # Update the line and the axes (with a new xlabel). Return a tuple of
    # "artists" that have to be redrawn for this frame.
    nrColumns = len(listdir("ISMolD_outputdata/nodes/time"+str(t)))-1 #-1 since file starts at 'output0', add -1 for amy subdirectory

    for i in range (nrColumns):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            data = np.loadtxt("ISMolD_outputdata/nodes/time"+str(t)+"/column"+str(i)+".txt")
    heatmapData = np.zeros(shape=(detailFactor*upperBound-detailFactor*lowerBound, detailFactor*nrColumns))
    
    for i in range(nrColumns):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            columnData = np.loadtxt("ISMolD_outputdata/nodes/time"+str(t)+"/column"+str(i)+".txt") ## nodeID, totalHeight, gravel, sand
        for j in range(len(columnData)):
            for l in range(detailFactor):
                xmodifier = detailFactor*i+l
                for k in range(detailFactor):
                    ymodifier = detailFactor*j+k-detailFactor*lowerBound
                    if (int(columnData.size/4) > 1):
                        heatmapData[ymodifier+math.ceil(detailFactor*columnData[j,2]),xmodifier] = columnData[j,3]
                    elif(int(columnData.size/4) == 1): 
                        heatmapData[k-detailFactor*lowerBound+math.ceil(detailFactor*columnData[2]),xmodifier] = columnData[3]
    im.set_data(heatmapData)
    im.set_clim(vmin=0, vmax=1)
    return im

# FuncAnimation will call the 'update' function for each frame; here
# animating over 10 frames, with an interval of 200ms between frames.
max_timestep = len(os.listdir("ISMolD_outputdata/relief"))-1 #-1 since file starts at 'time0'

anim = FuncAnimation(fig, update, frames=np.arange(0, max_timestep), interval=animationInterval)
if (saveImage):
    print("Making and Saving animation (This may take a while, maybe even several minutes) ...    ", end="\r")
    anim.save('sedimentContentInNodes.gif', dpi=80, writer='imagemagick')

if (showAnimation):
    print("Showing figure...    ")
    plt.show()
