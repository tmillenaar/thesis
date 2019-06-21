import os
import numpy as np
import matplotlib.pyplot as plt

yr2sec = 60*60*24*365.25      #nr of seconds in a year

data = np.loadtxt("ISMolD_outputdata/forcing.txt")
#data = np.delete(data, (0), axis=0)

fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)

ax1.set_title("Input")
ax1.plot(data[:,0]/yr2sec,data[:,1], linewidth = 1, label='Cumulative Combined ')  
ax1.plot(data[:,0]/yr2sec,data[:,2], linewidth = 1, label='Cumulative Gravel ')
ax1.plot(data[:,0]/yr2sec,data[:,3], linewidth = 1, label='Cumulative Sand ')
ax11 = ax1.twinx()
ax11.plot(data[:,0]/yr2sec,data[:,4], '--', linewidth = 1, label='Current Combined ')  
ax11.plot(data[:,0]/yr2sec,data[:,5], '--', linewidth = 1, label='Current Gravel ')
ax11.plot(data[:,0]/yr2sec,data[:,6], '--', linewidth = 1, label='Current Sand ')


ax2.set_title("Output")
ax2.plot(data[:,0]/yr2sec,data[:,7], linewidth = 1, label='Cumulative Combined ')
ax2.plot(data[:,0]/yr2sec,data[:,8], linewidth = 1, label='Cumulative Gravel ')
ax2.plot(data[:,0]/yr2sec,data[:,9], linewidth = 1, label='Cumulative Sand ')
ax22 = ax2.twinx()
ax22.plot(data[:,0]/yr2sec,data[:,10], '--', linewidth = 1, label='Combined Input')  
ax22.plot(data[:,0]/yr2sec,data[:,11], '--', linewidth = 1, label='Gravel Input')
ax22.plot(data[:,0]/yr2sec,data[:,12], '--', linewidth = 1, label='Sand Input')

ax3.set_title("Diffusivity and Subsidence")
ax3.plot(data[:,0]/yr2sec,data[:,13], linewidth = 1, label='Total diffusivity')  
ax3.plot(data[:,0]/yr2sec,data[:,14], linewidth = 1, label='Gravel diffusivity')
ax3.plot(data[:,0]/yr2sec,data[:,15], linewidth = 1, label='Sand diffusivity')
ax33 = ax3.twinx()
ax33.plot(data[:,0]/yr2sec,data[:,16], '--', linewidth = 1, color="purple", label='Subsidence Rate') 

ax1.set_ylabel('m^2')
ax11.set_ylabel('m^2/s')
ax2.set_ylabel('m^2')
ax22.set_ylabel('m^2/s')
ax3.set_ylabel('m^2')
ax33.set_ylabel('m^2/s')

## Use scientific notation for axes:
ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax11.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax22.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax3.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax33.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ax3.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

ax3.tick_params(labelrotation=20)

ax1.legend(loc='center left', bbox_to_anchor=(1.15, 0.8))
#ax2.legend(loc='center left', bbox_to_anchor=(1.3, 0.8))
#ax3.legend(loc='center left', bbox_to_anchor=(1.3, 0.8))
ax11.legend(loc='center left', bbox_to_anchor=(1.15, 0.1))
#ax22.legend(loc='center left', bbox_to_anchor=(1.3, 0.2))
ax33.legend(loc='center left', bbox_to_anchor=(1.15, 2.3))
plt.tight_layout(rect=[0,0,0.80,1])
plt.show()
