import os
import warnings     ## Used to filter numpy warnings when opening empty text/data files
import numpy as np
import matplotlib.pyplot as plt
import math
import sys


max_timestep = len(os.listdir("ISMolD_outputdata/relief"))-1 #-1 since file starts at 'time0'

## Assign True or False
timeInKyr = True     ## if False, plots in yr
    
xticks = [0, .2, .4, .6, .8, 1]

if (len(sys.argv) ==1 ):
    columns = ["", 10, 30, 50, 70, 90]
else:
    columns = sys.argv

## Find latest deposition time:
## Note: in the model, time 0 at the start of the model. When making logs, time=0 is considered to be the present. Therefore, the final timestep of the model reperesents present time.
nrColumns = len(os.listdir("ISMolD_outputdata/nodes/time"+str(max_timestep)))-1 #-1 since file starts at 'output0', add -1 for any subdirectory

totalElapsedTime = 0
for i in range(nrColumns):
    with warnings.catch_warnings(): 
        warnings.simplefilter("ignore") ## Ignores the warnings when numpy opens an empty file (most columns are empty/unfilled at the start of the run)
        data = np.loadtxt("ISMolD_outputdata/nodes/time"+str(max_timestep)+"/column"+str(i)+".txt") 
    try:
        totalElapsedTime = max( max(data[:,5]), totalElapsedTime)
    except: 
        try: ## Triggers if the data file has only one line
            totalElapsedTime = max( data[5], totalElapsedTime)
        except: ## Triggers if the data file is empty
            pass
print ("totalElapsedTime",int(totalElapsedTime),"yrs")

## Find longest column of those desired:
yticks = []
y_tick_labels = []
for j in range(len(columns)-1):
    logColumn = columns[j+1] ## The first argument "[0]" is the name of the program, which is the argument we want

    ## nodeID, totalHeight, bedrockHeight, gravel, sand, depositTime
    with warnings.catch_warnings(): 
        warnings.simplefilter("ignore") ## Ignores the warnings when numpy opens an empty file (most columns are empty/unfilled at the start of the run)
        data = np.loadtxt("ISMolD_outputdata/nodes/time"+str(max_timestep)+"/column"+str(logColumn)+".txt") 
        
    if (data.size>6):
        if(j==0 or max(data[:,0])>max(yticks)):
            new_tick_locations = list(range(len(data[:,0])))
            stepsize = math.ceil(len(new_tick_locations)/30)
            nrofticks = math.ceil(len(new_tick_locations) / stepsize)
            for i in range(nrofticks+1):
                yticks.append( i*stepsize )
                y_tick_labels.append( i*stepsize )
    else:
        print("")
        print("Error, column "+str(columns[j])+" (and likely more distal ones too) does not have enough data to draw a log. Two or more nodes are needed.")
        exit()

## Plot figures
if (len(columns) > 2):
    fig, ax = plt.subplots(ncols= len(columns)-1 )
    for j in range(len(columns)-1):
        logColumn = columns[j+1] ## The first argument "[0]" is the name of the program, which is the argument we want

        ## nodeID, totalHeight, bedrockHeight, gravel, sand, depositTime
        with warnings.catch_warnings(): 
            warnings.simplefilter("ignore") ## Ignores the warnings when numpy opens an empty file (most columns are empty/unfilled at the start of the run)
            data = np.loadtxt("ISMolD_outputdata/nodes/time"+str(max_timestep)+"/column"+str(logColumn)+".txt") 
        timeticks = []
        time_tick_labels = []
        new_tick_locations = list(range(len(data[:,0])))
        nrofticks = math.ceil(len(new_tick_locations) / stepsize)
        yScaleRatio = nrofticks/max(yticks)
        
        for i in range(nrofticks+2):
            timeticks.append(i*stepsize)
            if (timeInKyr): 
                if ((len(data[:,0])-1)-i*stepsize >= 0):
                    time_tick_labels.append(int( (totalElapsedTime - data[(len(data[:,0])-1)-i*stepsize,5])/100)/10 )
                else:
                    time_tick_labels.append("")
            else:
                time_tick_labels.append(int( data[(nrofticks-i-1)*stepsize,5]) )
                
        ax[j].set_yticks(yticks)
        ax[j].set_ylim([-5,max(yticks)+5])
        ax[j].plot(data[:,3],max(data[:,0])-data[:,0], linewidth = 1, label='Gravel fraction')
        #ax[j].plot(data[:,3],data[:,0], linewidth = 1, label='Gravel fraction')
        ax1 = ax[j].twinx()
        #ax1.set_ylim([-5*yScaleRatio, yScaleRatio*max(timeticks)+5*yScaleRatio])
        #ax1.set_ylim([0,max(yticks)])
        ax1.plot(data[:,3],max(data[:,0])-data[:,0], linewidth = 1, label='Gravel fraction')
        #ax1.plot(data[:,3],data[:,0], linewidth = 1, label='Gravel fraction')
        ax1.set_yticks(timeticks)
        ax1.set_ylim([-5,max(yticks)+5])
        #ax[j].set_yticklabels(y_tick_labels)
        ax[j].set_xticks(xticks)
        ax1.set_xticks(xticks)
        ax[j].set_xticklabels(xticks, rotation = 90)
        ax1.set_yticklabels(time_tick_labels)
        ax[j].invert_yaxis()
        ax1.invert_yaxis()
        
        ax[j].set_title("Column "+str(columns[j+1]))
        ax[j].grid()
    ax[0].set_ylabel("Depth [m]")
        
else: ## If only one column is requested
    fig, ax = plt.subplots(ncols= len(columns)-1 )
    logColumn = columns[j+1] ## The first argument "[0]" is the name of the program, which is the argument we want

    ## nodeID, totalHeight, bedrockHeight, gravel, sand, depositTime
    with warnings.catch_warnings(): 
        warnings.simplefilter("ignore") ## Ignores the warnings when numpy opens an empty file (most columns are empty/unfilled at the start of the run)
        data = np.loadtxt("ISMolD_outputdata/nodes/time"+str(max_timestep)+"/column"+str(logColumn)+".txt") 
    timeticks = []
    time_tick_labels = []
    new_tick_locations = list(range(len(data[:,0])))
    nrofticks = math.ceil(len(new_tick_locations) / stepsize)
    yScaleRatio = nrofticks/max(yticks)
    
    for i in range(nrofticks+2):
        timeticks.append(i*stepsize)
        if (timeInKyr): 
            if ((len(data[:,0])-1)-i*stepsize >= 0):
                time_tick_labels.append(int( (totalElapsedTime - data[(len(data[:,0])-1)-i*stepsize,5])/100)/10 )
            else:
                time_tick_labels.append("")
        else:
            time_tick_labels.append(int( data[(nrofticks-i-1)*stepsize,5]) )
            
    ax.set_yticks(yticks)
    ax.set_ylim([-5,max(yticks)+5])
    ax.plot(data[:,3],max(data[:,0])-data[:,0], linewidth = 1, label='Gravel fraction')
    #ax[j].plot(data[:,3],data[:,0], linewidth = 1, label='Gravel fraction')
    ax1 = ax.twinx()
    #ax1.set_ylim([-5*yScaleRatio, yScaleRatio*max(timeticks)+5*yScaleRatio])
    #ax1.set_ylim([0,max(yticks)])
    ax1.plot(data[:,3],max(data[:,0])-data[:,0], linewidth = 1, label='Gravel fraction')
    #ax1.plot(data[:,3],data[:,0], linewidth = 1, label='Gravel fraction')
    ax1.set_yticks(timeticks)
    ax1.set_ylim([-5,max(yticks)+5])
    #ax[j].set_yticklabels(y_tick_labels)
    ax.set_xticks(xticks)
    ax1.set_xticks(xticks)
    ax.set_xticklabels(xticks, rotation = 90)
    ax1.set_yticklabels(time_tick_labels)
    ax.invert_yaxis()
    ax1.invert_yaxis()
    
    ax.set_title("Column "+str(columns[j+1]))
    ax.grid()  
    ax.set_ylabel("Depth [m]")

if (timeInKyr):
    ax1.set_ylabel("Time of deposition [ka]")
else:
    ax1.set_ylabel("Time of deposition [yr ago]")
    
    

plt.tight_layout()
plt.show()
