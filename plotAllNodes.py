import os
import math         ## Used for floor() and ceil() to obtain integers from floats
import warnings     ## Used to filter numpy warnings when opening empty text/data files
import numpy as np  
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation 
from matplotlib.colors import LinearSegmentedColormap
from os import listdir


detailFactor = 20 ## Mesh size for the plot. Higher detailFactor creates a smoother bedrock slope. Example: if detailFactor = 10, each mesh cell is divided into 10x10 subcells


## Set custom color scale:
        
## blue-green-yellow
mycdict = {'red':  ((0.0, 1.0, 1.0),
                   (0.01, 0.0, 0.0),
                   (0.25, 0.0, 0.0),
                   (0.5, 0.0, 0.0),
                   (0.75, 0.2, 0.2),
                   (0.95, 0.8, 0.8),
                   (1.0, 1.0, 1.0)),

         'green': ((0.0, 1.0, 1.0),
                   (0.01, 0.0, 0.0),
                   (0.25, 0.2, 0.2),
                   (0.5, 0.7, 0.7),
                   (0.75, 0.7, 0.7),
                   (0.95, 0.8, 0.8),
                   (1.0, 0.8, 0.8)),

         'blue':  ((0.0, 1.0, 1.0),
                   (0.01, 0.4, 0.4),
                   (0.25, 0.7, 0.7),
                   (0.5, 0.7, 0.7),
                   (0.75, 0.2, 0.2),
                   (0.95, 0.0, 0.0),
                   (1.0, 0.0, 0.0))
        }
         


if (detailFactor <= 0): detailFactor = 1 ## Though I cannot imagine someone wants negative detail ^^' 

if (not os.path.isdir("ISMolD_outputdata/nodes/time0")):
    print("Error, ISMolD_outputdata/nodes/time0 does not exist.")
    exit()

print("Determining the extent of the data...", end="\r")
upperBound = 0
lowerBound = 0
max_timestep = len(os.listdir("ISMolD_outputdata/relief"))-1 #-1 since file starts at 'time0'
nrColumns = len(listdir("ISMolD_outputdata/nodes/time"+str(max_timestep)))-1 #-1 since file starts at 'output0', add -1 for amy subdirectory

if (not os.path.isdir("allNodes")):
    os.mkdir("allNodes")
else:
    if (len(listdir("ISMolD_outputdata/nodes/time"+str(max_timestep))) != 0 ):
        os.system("rm allNodes/nodesPlot*.png" )

totalElapsedTime = 0
for i in range(nrColumns):
    data = np.loadtxt("ISMolD_outputdata/nodes/time"+str(max_timestep)+"/column"+str(i)+".txt") 
    try:
        totalElapsedTime = max( max(data[:,5]), totalElapsedTime)
    except: ## Triggers if the data file has only one line
        totalElapsedTime = max( data[5], totalElapsedTime)

for t in range(max_timestep):

    nrColumns = len(listdir("ISMolD_outputdata/nodes/time"+str(t)))-1 #-1 since file starts at 'output0', add -1 for amy subdirectory
    for i in range (nrColumns):
        with warnings.catch_warnings(): 
            warnings.simplefilter("ignore") ## Ignores the warnings when numpy opens an empty file (most columns are empty/unfilled at the start of the run)
            data = np.loadtxt("ISMolD_outputdata/nodes/time"+str(t)+"/column"+str(i)+".txt")
        try:
            upperBound = max(upperBound, max(data[:,1]))
        except:
            pass
        try:
            lowerBound = min(lowerBound, min(data[:,2]))
        except:
            pass

upperBound = math.ceil(upperBound+0.1*(upperBound-lowerBound))
lowerBound = math.floor(lowerBound-0.1*(upperBound-lowerBound))     

heatmapData = np.zeros(shape=(detailFactor*upperBound-detailFactor*lowerBound, detailFactor*nrColumns))

if (upperBound == 0):
    print("Error, no data")
    exit()
        
for t in range(max_timestep):
    print("Making image "+str(t)+"/"+str(max_timestep)+"                     ", end="\r")
    heatmapData = np.zeros(shape=(detailFactor*upperBound-detailFactor*lowerBound, detailFactor*nrColumns))
    for i in range(nrColumns):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            columnData = np.loadtxt("ISMolD_outputdata/nodes/time"+str(t)+"/column"+str(i)+".txt") ## nodeID, totalHeight, bedrockHeight, gravel, sand, depositTime
            
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
    im = ax.imshow(heatmapData, vmin=0, vmax=1, aspect="auto", extent=[0,nrColumns,upperBound,lowerBound], cmap=LinearSegmentedColormap('CustomColorSet', mycdict)) ## , extent=[80,120,32,upperBound]
    fig.colorbar(im).set_label("Gravel fraction")
    ax.set_title("Sediment content per node. Time: "+str(int(t*totalElapsedTime/(1000*max_timestep)))+"kyr")
    ax.invert_yaxis()
    ax.set_xlabel('Distance form source [km]')
    ax.set_ylabel('Height [m]')
    
    plt.tight_layout(rect=[0,0,0.85,1])
    plt.savefig("allNodes/nodesPlot"+str(t)+".png")
    plt.close(fig)

print("Done. Images saved in /allNodes")

