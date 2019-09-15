import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation 
from os import listdir


if (not os.path.isdir("ISMolD_outputdata/relief")):
    print("Error, directory ISMolD_outputdata/relief does not exist.")
    exit()

nr_frames = len(listdir("ISMolD_outputdata/relief"))-1 #-1 since file starts at 'output0', add -1 for amy subdirectory

lowerBound = 0
upperBound = 0
for i in range (nr_frames):
    data = np.loadtxt('ISMolD_outputdata/relief/topography'+str(i)+'.txt')
    lowerBound = min(lowerBound, min(data[:,2]-10))
    upperBound = max(upperBound, max(data[:,1]+20))
    
#data = np.delete(data, (0), axis=0)

fig, ax = plt.subplots()
fig.set_tight_layout(True)
plt.ylim(lowerBound, upperBound)
plt.grid(color='black', linestyle=(0, (1, 10)), linewidth=1)
plt.xlabel('Distance from source in km')
plt.ylabel('Height in m')

# Query the figure's on-screen size and DPI. Note that when saving the figure to
# a file, we need to provide a DPI for that separately.

# Plot a scatter that persists (isn't redrawn) and the initial line.
totalHeight, = ax.plot(data[:,0], data[:,1], 'b-', linewidth=3, label='Relief')
bedrock, = ax.plot(data[:,0], data[:,2], 'r-', linewidth=3, label='Basement')
gravel, = ax.plot(data[:,0], data[:,3], 'green', linewidth=2, label='Gravel')
sand, = ax.plot(data[:,0], data[:,4], 'orange', linewidth=2, label='Sand')
#silt, = ax.plot(data[:,0], data[:,1], 'purple', linewidth=2)

def update(i):
    label = 't = {0}0 kyr'.format(i)
    #print(label)
    
    # Update the line and the axes (with a new xlabel). Return a tuple of
    # "artists" that have to be redrawn for this frame.
    data = np.loadtxt('ISMolD_outputdata/relief/topography'+str(i)+'.txt')
    totalHeight.set_ydata(data[:,1])
    bedrock.set_ydata(data[:,2])
    gravel.set_ydata(data[:,3])
    sand.set_ydata(data[:,4])
    #bedrock.set_ydata(data[:,4])
    #silt.set_ydata(data[:,1])
    ax.set_title(label)
    return totalHeight, bedrock, gravel, sand, ax

if __name__ == '__main__':
    # FuncAnimation will call the 'update' function for each frame; here
    # animating over 10 frames, with an interval of 200ms between frames.
    anim = FuncAnimation(fig, update, frames=np.arange(0, nr_frames), interval=50)
    plt.legend(loc=1, borderaxespad=0.)
    # anim.save('test.gif', dpi=80, writer='imagemagick')
    plt.show()
