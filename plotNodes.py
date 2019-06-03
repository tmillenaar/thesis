import os
import numpy as np
import matplotlib.pyplot as plt

if (not os.path.isdir("ISMolD_outputdata/nodes")):
    print("Error, directory ISMolD_outputdata/nodes does not exist.")
    exit()

nr_frames = len(os.listdir("ISMolD_outputdata/relief"))-1 #-1 since file starts at 'output0', add -1 for any subdirectory

max_timestep = len(os.listdir("ISMolD_outputdata/relief"))-1 #-1 since file starts at 'time0'
print("nr_timesteps ", nr_timesteps )

upperBound = 0
for i in range (nr_frames):
    data = np.loadtxt('ISMolD_outputdata/nodes/column'+str(i)+'.txt')
    upperBound = max(upperBound, int(data.size/4))

heatmapData = np.zeros(shape=(upperBound, nr_frames))

if (upperBound == 0):
    print("Error, no data")
    exit()

for i in range(nr_frames):
    columnData = np.loadtxt("ISMolD_outputdata/nodes/column"+str(i)+".txt") ## nodeID, totalHeight, gravel, sand
    print("")
    print("Column", i,columnData)
    for j in range(len(columnData)):
        if (int(columnData.size/4) > 1):
            heatmapData[j,i] = columnData[j,2]
        elif(int(columnData.size/4) == 1): 
            heatmapData[0,i] = columnData[2]
    #data = np.append(data, [[columnData]])
    #data.append(np.loadtxt("ISMolD_outputdata/nodes/column"+str(i)+".txt"))
    
fig, ax = plt.subplots()
im = ax.imshow(heatmapData)
ax.invert_yaxis()

#print(data)
print(heatmapData)
#plt.plot(data[:,0],data[:,2], linewidth = 3, label='Relief')  
#plt.plot(data[:,0],data[:,1], '--', linewidth = 4, label='silt')
#plt.plot(data[:,0],data[:,2], linewidth = 1, label='sand')
#plt.plot(data[:,0],data[:,3], linewidth = 1, label='gravel')
#plt.plot(data[:,0],data[:,4], linewidth = 3, label='bedrock')
#plt.legend(loc=1, borderaxespad=0.)
plt.show()
