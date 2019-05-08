#import os
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('outputData/output9.xy')
data = np.delete(data, (0), axis=0)

    
plt.plot(data[:,0],data[:,5], linewidth = 3, label='Relief')  
plt.plot(data[:,0],data[:,1], '--', linewidth = 4, label='silt')
plt.plot(data[:,0],data[:,2], linewidth = 1, label='sand')
plt.plot(data[:,0],data[:,3], linewidth = 1, label='gravel')
plt.plot(data[:,0],data[:,4], linewidth = 3, label='bedrock')
plt.legend(loc=1, borderaxespad=0.)
plt.show()
