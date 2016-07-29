###
# Author: Julian Koch (juko@geus.dk)
# Date: June/2016
###
# This is an example for how to use the FSS function to quantify spatial 
# similarity between two datasets. In this example a reference dataset is
# referred to as "obs" and a "model" refers to the dataset to be evaluated

import scipy
import numpy as np
from SEEM import count_neigh
from SEEM import search_neigh
import matplotlib.pyplot as plt 


# Where is the example data?     
Path='../Data/' 

# Do you want to plot the EOF results?
plot=True

# read data
mydata=scipy.io.loadmat((Path+'ET_data_2008.mat'))
myET=mydata.get('ET_data_out')
mydata=scipy.io.loadmat((Path+'LST_data_2008'))
myLST=mydata.get('LST_data_out')

print ('read data done!')

# Which variable do you want to investigate?
# Evapotranspiration=myET or land-surface-temperature=myLST
myVar=myLST

# Which scenario do you want to investigate?
# Reduced spatial detail in vegetation parameters =1 or
# Reduced spatial detail in climatic forcing = 2
myScen=1


size=myVar.shape

# Here you need to define a combination of threshold percentiles and scales
# Define the arrays in equal size, because only direct combinations are considered, e.g. perc_FSS[0] and scales_FSS[0], perc_FSS[1] and scales_FSS[1], etc
# This concept of critical scale helps to reduce the computational time, but requires thought on how to define the critical combinations of scale and threshold.
# In our example, we allow higher uncertainty to the extreme percentiles and moderate percentiles are addressed  at finer scale, because placement errors are less tolerated here than for extreme percentiles. 
perc_FSS=(2,10,20,80,90,98)
scales_FSS=(25,15,5,5,15,25) # prefer uneven number

# to save the results
MSE=np.empty([len(scales_FSS),size[2]]) # mean squared error is an intermedaite result
FSS=np.empty([len(scales_FSS),size[2]]) # final FSS score
 

# Compute FSS for each combination of scale and threshold at each of the 365 days. This will take a couple mins
for i in range(0,size[2]):
    for a in range(0,len(perc_FSS)):     
        # Compute the fractions for model and observed. Fractions relate to occurrences of values exceeding a certain threshold at a given window size (scale).          
        # Fractions are computed from dividing the search neighbor function with the count neighbor function.          
        frac_model=search_neigh(myVar[:,:,i,myScen],scales_FSS[a],perc_FSS[a])/count_neigh(myVar[:,:,i,myScen],scales_FSS[a])
        frac_obs=search_neigh(myVar[:,:,i,0],scales_FSS[a],perc_FSS[a])/count_neigh(myVar[:,:,i,0],scales_FSS[a])
        # The mean-squared-error is computed for each combination of scale and threshold.          
        MSE[a,i]=np.nanmean(np.square(frac_model-frac_obs))
        # FSS is the MSE normalized by a worst case MSE that reflects the condition of zero match between the observed and modeled fractions; simply adding the squared fractions. 
        FSS[a,i]=1-MSE[a,i]/(np.nanmean(np.square(frac_model))+np.nanmean(np.square(frac_obs)))
    print (('day: '+str(i)+' done!'))

FSS_out_low=np.mean(FSS[0:3,:],axis=0)
FSS_out_high=np.mean(FSS[3:6,:],axis=0)
#### The rest of the script is for plotting the FSS output
if plot==True:
    x = [0,0+121,0+242,365]
    xlabels = ['jan08','may08','sep08','jan09']
    plt.figure(0) 
    plt.figure(figsize=(10,4))
    plt.plot(FSS[0,:],label='2nd @ 25km scale',color = '0.3')
    plt.plot(FSS[1,:],label='10th @ 15km scale',color = '0.5')
    plt.plot(FSS[2,:],label='20th @ 5km scale',color = '0.7')
    plt.plot(FSS_out_low,label='mean',color='r')
    plt.ylabel('FSS')
    plt.xticks(x,xlabels,rotation=45)
    plt.title('FSS for low percentiles')
    plt.legend(loc=4)
    plt.axis((0,365,0,1))

    plt.figure(1) 
    plt.figure(figsize=(10,4))
    plt.plot(FSS[5,:],label='98th @ 25km scale',color = '0.3')
    plt.plot(FSS[4,:],label='90th @ 15km scale',color = '0.5')
    plt.plot(FSS[3,:],label='80th @ 5km scale',color = '0.7')
    plt.plot(FSS_out_high,label='mean',color='r')
    plt.ylabel('FSS')
    plt.xticks(x,xlabels,rotation=45)
    plt.title('FSS for high percentiles')
    plt.legend(loc=4)
    plt.axis((0,365,0,1))