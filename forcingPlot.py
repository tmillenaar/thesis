import os
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("ISMolD_outputdata/forcing.txt")
#data = np.delete(data, (0), axis=0)

fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)

ax1.set_title("Input")
ax1.plot(data[:,0],data[:,1], linewidth = 1, label='Total Input')  
ax1.plot(data[:,0],data[:,2], linewidth = 1, label='Gravel Input')
ax1.plot(data[:,0],data[:,3], linewidth = 1, label='Sand Input')

ax2.set_title("Output")
ax2.plot(data[:,0],data[:,4], '--', linewidth = 1, label='Total Output')
ax2.plot(data[:,0],data[:,5], '--', linewidth = 1, label='Gravel Output')
ax2.plot(data[:,0],data[:,6], '--', linewidth = 1, label='Sand Output')

ax3.set_title("Diffusivity")
ax3.plot(data[:,0],data[:,7], linewidth = 1, label='Total diffusivity')  
ax3.plot(data[:,0],data[:,8], linewidth = 1, label='Gravel diffusivity')
ax3.plot(data[:,0],data[:,9], linewidth = 1, label='Sand diffusivity')

ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout(rect=[0,0,0.75,1])
plt.show()
