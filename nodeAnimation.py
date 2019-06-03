import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation 
from os import listdir
import warnings

if (not os.path.isdir("ISMolD_outputdata/nodes/time0")):
    print("Error, directory ISMolD_outputdata/nodes/time0 does not exist.")
    exit()

print("Determining the extend of the data...", end="\r")
upperBound = 0
max_timestep = len(os.listdir("ISMolD_outputdata/relief"))-1 #-1 since file starts at 'time0'
for t in range(max_timestep):
    nrColumns = len(listdir("ISMolD_outputdata/nodes/time"+str(t)))-1 #-1 since file starts at 'output0', add -1 for amy subdirectory

    for i in range (nrColumns):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            data = np.loadtxt("ISMolD_outputdata/nodes/time"+str(t)+"/column"+str(i)+".txt")
        upperBound = max(upperBound, int(data.size/4))
        
heatmapData = np.zeros(shape=(upperBound, nrColumns))

if (upperBound == 0):
    print("Error, no data")
    exit()

print("Handeling data...                          ", end="\r") ## These spaces are to overwrite the whole of the previous line

for i in range(nrColumns):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        columnData = np.loadtxt("ISMolD_outputdata/nodes/time0/column"+str(i)+".txt") ## nodeID, totalHeight, gravel, sand
    
    for j in range(len(columnData)):
        if (int(columnData.size/4) > 1):
            heatmapData[j,i] = columnData[j,2]
        elif(int(columnData.size/4) == 1): 
            heatmapData[0,i] = columnData[2]
            
fig, ax = plt.subplots()
im = ax.imshow(heatmapData, vmin=0, vmax=1, aspect="equal")
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
    heatmapData = np.zeros(shape=(upperBound, nrColumns))
    
    for i in range(nrColumns):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            columnData = np.loadtxt("ISMolD_outputdata/nodes/time"+str(t)+"/column"+str(i)+".txt") ## nodeID, totalHeight, gravel, sand
        for j in range(len(columnData)):
            if (int(columnData.size/4) > 1):
                heatmapData[j,i] = columnData[j,2]
            elif(int(columnData.size/4) == 1): 
                heatmapData[0,i] = columnData[2]
    im.set_data(heatmapData)
    im.set_clim(vmin=0, vmax=1)
    return im

# FuncAnimation will call the 'update' function for each frame; here
# animating over 10 frames, with an interval of 200ms between frames.
max_timestep = len(os.listdir("ISMolD_outputdata/relief"))-1 #-1 since file starts at 'time0'

anim = FuncAnimation(fig, update, frames=np.arange(0, max_timestep), interval=10)
print("Saving figure...    ", end="\r")
anim.save('sedimentContentInNodes.gif', dpi=80, writer='imagemagick')
print("Showing figure...    ")
plt.show()
