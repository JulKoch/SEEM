###
# Author: Julian Koch (juko@geus.dk)
# Date: June/2016
###
# This is an example for how to use the Connectivity analysis to quantify spatial 
# similarity between two datasets. In this example a reference dataset is
# referred to as "obs" and a "model" refers to the dataset to be evaluated


import numpy as np
import scipy
import matplotlib.pyplot as plt
from SEEM import connectivity

# Where is the example data?     
Path='../Data/' 

# Do you want to plot the Connectivity results?
plot=True

# read data
mydata=scipy.io.loadmat((Path+'ET_data_2008.mat'))
myET=mydata.get('ET_data_out')
mydata=scipy.io.loadmat((Path+'LST_data_2008'))
myLST=mydata.get('LST_data_out')

print ('read data done!')

# Which variable do you want to investigate?
# Evapotranspiration=myET or land-surface-temperature=myLST

myVar=myET

# Which scenario do you want to investigate?
# Reduced spatial detail in vegetation parameters =1 or
# Reduced spatial detail in climatic forcing = 2

myScen=1

size=myVar.shape

Con_high=np.empty([size[2]])
Con_low=np.empty([size[2]])

# call the connectivity function for each day; 365. This will take a couple of mins
for i in range(0,size[2]):

    temp=connectivity(myVar[:,:,i,0],(myVar[:,:,i,myScen]))
    Con_high[i]=temp[8] # final connectivity metric for high phase 
    Con_low[i]=temp[9] # final connectivity metric for low phase
    print (('Day '+str(i)+' done'))
    
    
#### The rest of the script is for plotting the Connectivity output
if plot==True:
    
    x = [0,0+121,0+242,365]
    xlabels = ['jan08','may08','sep08','jan09']
    # plot timeseries of the pattern similairty score based on the EOF analysis
    plt.figure(0)
    plt.figure(figsize=(6,3))    
    plt.plot(Con_high,'r',label='Pattern Score')
    plt.xticks(x,xlabels,rotation=45)
    lim=np.nanmax(Con_high)    
    plt.axis((0,365,0,lim))
    plt.ylabel('Pattern Similarity Score')
    plt.title('Connectivity High Phase')

    plt.figure(1)
    plt.figure(figsize=(6,3))    
    plt.plot(Con_low,'b',label='Pattern Score')
    plt.xticks(x,xlabels,rotation=45)
    lim=np.nanmax(Con_low)    
    plt.axis((0,365,0,lim))
    plt.ylabel('Pattern Similarity Score')
    plt.title('Connectivity Low Phase')

