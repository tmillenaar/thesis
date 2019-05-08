#import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
#import matplotlib.animation as anim
#from matplotlib.animation import Animation 
from os import listdir


nr_frames = len(listdir("outputData/"))-2 #-1 since file starts at 'output0' and -1 for the 'old' directory

data = np.loadtxt('outputData/output'+str(nr_frames)+'.xy')

#data = np.delete(data, (0), axis=0)

fig, ax = plt.subplots()
fig.set_tight_layout(True)
plt.ylim(min(data[:,4]-50),max(data[:,5]+100))
plt.xlabel('Distance from source in km')
plt.ylabel('Height in m')

# Query the figure's on-screen size and DPI. Note that when saving the figure to
# a file, we need to provide a DPI for that separately.

# Plot a scatter that persists (isn't redrawn) and the initial line.

totalHeight, = ax.plot(data[:,0], data[:,5], 'b-', linewidth=3, label='Relief')
bedrock, = ax.plot(data[:,0], data[:,4], 'r-', linewidth=3, label='Basement')
gravel, = ax.plot(data[:,0], data[:,3], 'green', linewidth=2, label='Gravel')
sand, = ax.plot(data[:,0], data[:,2], 'orange', linewidth=2, label='Sand')
#silt, = ax.plot(data[:,0], data[:,1], 'purple', linewidth=2)

def update(i):
    label = 't = {0}0 kyr'.format(i)
    #print(label)
    
    # Update the line and the axes (with a new xlabel). Return a tuple of
    # "artists" that have to be redrawn for this frame.
    data = np.loadtxt('outputData/output'+str(i)+'.xy')
    totalHeight.set_ydata(data[:,5])
    gravel.set_ydata(data[:,3])
    sand.set_ydata(data[:,2])
    bedrock.set_ydata(data[:,4])
    #silt.set_ydata(data[:,1])
    ax.set_title(label)
    return totalHeight, bedrock, gravel, sand, ax

if __name__ == '__main__':
    # FuncAnimation will call the 'update' function for each frame; here
    # animating over 10 frames, with an interval of 200ms between frames.
    print('Plotting')
    anim = FuncAnimation(fig, update, frames=np.arange(0, nr_frames), interval=20)
    #if len(sys.argv) > 1 and sys.argv[1] == 'save':
        #anim.save('line.gif', dpi=80, writer='imagemagick')
    #else:
        ## plt.show() will just loop the animation forever.
    plt.legend(loc=1, borderaxespad=0.)
    anim.save('test.gif', dpi=80, writer='imagemagick')
    plt.show()
