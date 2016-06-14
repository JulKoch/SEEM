###
# Author: Julian Koch (juko@geus.dk)
# Date: June/2016
###
# This is an example of how to use the EOF function to quantify spatial 
# similarity between two datasets. In this example a reference dataset is
# referred to as "obs" and a "model" refers to the dataset to be evaluated

import scipy
import numpy as np
from SEEM import EOF_similarity
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
EOF=np.empty([size[2]])

# Prepare data for the EOF analysis. The new dimensions reflect the positions 
# in space as rows and the timesteps as columns. All empty rows are deleted.
# 2410 locations and 365 timesteps  
 
obs=np.reshape(myVar[:,:,:,0],(size[0]*size[1],size[2]))
mask=obs[:,0] # mask is used later to reshape the EOF output to the orginal size
obs=obs[~np.isnan(obs).any(axis=1)]
model=np.reshape(myVar[:,:,:,myScen],(size[0]*size[1],size[2]))
model=model[~np.isnan(model).any(axis=1)]    

# EOF - Analysis
a=EOF_similarity(obs,model)

## output
EOF[:]=a[0] # Skill Score
eof1=a[1] # EOF1
eof2=a[2] # EOF2
eof3=a[3] # EOF3
var=a[4] # explained variance
load_obs=a[5] # loadings obs
load_model=a[6] # loadings model
print(('Scenario '+str(myScen)+' EOF done'))
##        

#### The rest of the script is for plotting the EOF output
if plot==True:

# The EOF maps have to be reshaped from 1D to 2D using the 1D mask saved before.  
    eof1_map=np.empty([size[0],size[1]])
    eof2_map=0*eof1_map
    eof3_map=0*eof1_map
    k=0
    l=0
    for i in range(0,size[0]):
        for j in range(0,size[1]):
            if np.isfinite(mask[k]):
                eof1_map[i,j]=eof1[l]
                eof2_map[i,j]=eof2[l]
                eof3_map[i,j]=eof3[l]
                l=l+1
                k=k+1
            else:
                eof1_map[i,j]=np.nan
                eof2_map[i,j]=np.nan
                eof3_map[i,j]=np.nan
                k=k+1
                    
    eof1_map=np.flipud(eof1_map)            
    eof2_map=np.flipud(eof2_map)     
    eof3_map=np.flipud(eof3_map) 
## Finished reshaping the maps.    

# Plotting details for X-axis    
    x = [0,0+121,0+242,365]
    xlabels = ['jan08','may08','sep08','jan09']
    
    plt.figure(figsize=(12,6))
    plt.subplots_adjust(hspace=0.2,wspace=0.2)    
    # plot EOF1-3 maps    
     
    #EOF1
    plt.subplot(2,3,1) 
    lim=np.max([np.nanmax(eof1_map),np.abs(np.nanmin(eof1_map))])    
    plt.imshow(eof1_map,cmap='RdGy',interpolation='none')
    plt.title('EOF1 - '+str(var[0]*100)[0:4]+'% explained variance')
    plt.axis('off')
    plt.colorbar(orientation ='vertical',ticks=[-lim,0,lim],fraction=0.03)       
    plt.clim(-lim,lim) 
    
    # EOF2    
    plt.subplot(2,3,2)
    lim=np.max([np.nanmax(eof2_map),np.abs(np.nanmin(eof2_map))])
    plt.imshow(eof2_map,cmap='RdGy',interpolation='none')
    plt.title('EOF2 - '+str(var[1]*100)[0:4]+'% explained variance')
    plt.axis('off')
    plt.colorbar(orientation ='vertical',ticks=[-lim,0,lim],fraction=0.03)       
    plt.clim(-lim,lim)
    
    # EOF3    
    plt.subplot(2,3,3)
    lim=np.max([np.nanmax(eof3_map),np.abs(np.nanmin(eof3_map))])
    plt.imshow(eof3_map,cmap='RdGy',interpolation='none')
    plt.title('EOF3 - '+str(var[2]*100)[0:4]+'% explained variance')
    plt.axis('off')
    plt.colorbar(orientation ='vertical',ticks=[-lim,0,lim],fraction=0.03)       
    plt.clim(-lim,lim)
    
    
    # plot loadings
    # loadings associated to EOF1    
    plt.subplot(2,3,4)
    plt.plot(load_obs[:,0],'k')
    plt.plot(load_model[:,0],'r')
    plt.xticks(x,xlabels,rotation=45)
    lim=np.max([np.nanmax(load_obs[:,0]),np.abs(np.nanmin(load_obs[:,0]))])    
    plt.axis((0,365,-lim,lim))
    plt.yticks([-lim,0,lim])
    plt.ylabel('EOF Loadings')
    # loadings associated to EOF2  
    plt.subplot(2,3,5)
    plt.plot(load_obs[:,1],'k')
    plt.plot(load_model[:,1],'r')
    plt.xticks(x,xlabels,rotation=45)
    lim=np.max([np.nanmax(load_obs[:,1]),np.abs(np.nanmin(load_obs[:,1]))])    
    plt.axis((0,365,-lim,lim))
    plt.yticks([-lim,0,lim])
    # loadings associated to EOF3      
    plt.subplot(2,3,6)
    plt.plot(load_obs[:,2],'k',label='baseline')
    plt.plot(load_model[:,2],'r',label='scenario')
    plt.xticks(x,xlabels,rotation=45)
    lim=np.max([np.nanmax(load_obs[:,2]),np.abs(np.nanmin(load_obs[:,2]))])    
    plt.axis((0,365,-lim,lim))
    plt.yticks([-lim,0,lim])
    plt.legend(loc=4,fontsize=10,frameon=False,ncol=2,columnspacing=0.2,labelspacing=0.1,handletextpad=0.1)         
    
   