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

# Which day of the year do you want to investigate
day=10



size=myVar.shape

# Here you need to define a combination of threshold percentiles and scales
perc_FSS=(5,10,20,80,90,95)
scales_FSS=(1,3,5,9,15,25,50,99) # prefer uneven number

# to save the results
MSE=np.empty([len(scales_FSS),len(perc_FSS)]) # mean squared error is an intermedaite result
FSS=np.empty([len(scales_FSS),len(perc_FSS)]) # final FSS score

# Compute FSS for each combination of scale and threshold
for a in range(0,len(scales_FSS)):
    for b in range(0,len(perc_FSS)):     
        # Compute the fractions for model and observed. Fractions relate to occurrences of values exceeding a certain threshold at a given window size (scale).          
        # Fractions are computed from dividing the search neighbor function with the count neighbor function.          
        frac_model=search_neigh(myVar[:,:,day,myScen],scales_FSS[a],perc_FSS[b])/count_neigh(myVar[:,:,day,myScen],scales_FSS[a])
        frac_obs=search_neigh(myVar[:,:,day,0],scales_FSS[a],perc_FSS[b])/count_neigh(myVar[:,:,day,0],scales_FSS[a])
        # The mean-squared-error is computed for each combination of scale and threshold.          
        MSE[a,b]=np.nanmean(np.square(frac_model-frac_obs))
        # FSS is the MSE normalized by a worst case MSE that reflects the condition of zero match between the observed and modeled fractions; simply adding the squared fractions. 
        FSS[a,b]=1-MSE[a,b]/(np.nanmean(np.square(frac_model))+np.nanmean(np.square(frac_obs)))
    print (('scale: '+str(a)+' done!'))
#### The rest of the script is for plotting the FSS output
if plot==True:

    for i in range(0,len(perc_FSS)):  
        plt.plot(scales_FSS,FSS[:,i],label=('threshold: '+str(perc_FSS[i])+'%'))
        plt.ylabel('FSS')
        plt.xlabel('scale [km]')
        plt.title('FSS for various threshold percentiles')
        plt.legend(loc=4)

