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
from matplotlib import colors 
from SEEM import connectivity

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

myScen=2

size=myVar.shape

# Which DOY to use as an example?
j=97

temp=connectivity(myVar[:,:,j,0],(myVar[:,:,j,myScen]))
con1=temp[0]
con2=temp[1]
con3=temp[2]
con4=temp[3]


cmap = colors.ListedColormap(['0.9', '0.9'])

a=temp[4]
a[a==0]=np.nan
b=temp[5]
b[b==0]=np.nan
c=temp[6]
c[c==0]=np.nan
d=temp[7]
d[d==0]=np.nan
from parula import parula_cm

cmi=np.nanmin(myVar[:,:,j,0])
cma=np.nanmax(myVar[:,:,j,0])

mask=np.isfinite(np.flipud(myVar[:,:,j,1]))*1
mask=mask.astype('float')
mask[mask==0]=np.nan


plt.figure(figsize=(8,12))    
plt.subplots_adjust(hspace=0.1,wspace=0) 
mi=np.nanmin(myVar[:,:,j,0])
ma=np.nanmax(myVar[:,:,j,0])
plt.subplot(4,2,1)
plt.imshow(np.flipud(myVar[:,:,j,0]),cmap=parula_cm,interpolation='none',vmin=mi, vmax=ma)
plt.clim(mi,ma)

plt.title('baseline')
plt.axis('off')


plt.subplot(4,2,2)
plt.imshow(np.flipud(myVar[:,:,j,myScen]),cmap=parula_cm,interpolation='none',vmin=mi, vmax=ma)
plt.clim(mi,ma)
plt.colorbar(fraction=0.04)
plt.title('scenario')
plt.axis('off')


plt.subplot(4,2,3)
plt.imshow(mask,interpolation='none',cmap=cmap)
plt.imshow(np.flipud(c[:,:,20]),cmap=plt.get_cmap('Dark2'),interpolation='none')
plt.axis('off')
plt.title('Cluster-20th percentile')


plt.subplot(4,2,4)
plt.imshow(mask,interpolation='none',cmap=cmap)
plt.imshow(np.flipud(d[:,:,20]),cmap=plt.get_cmap('Dark2'),interpolation='none')
plt.axis('off')
plt.title('Cluster-20th percentile')


plt.subplot(4,2,5)
plt.imshow(mask,interpolation='none',cmap=cmap)
plt.imshow(np.flipud(a[:,:,80]),cmap=plt.get_cmap('Dark2'),interpolation='none')
plt.axis('off')
plt.title('Cluster-80th percentile')

plt.subplot(4,2,6)
plt.imshow(mask,interpolation='none',cmap=cmap)
plt.imshow(np.flipud(b[:,:,80]),cmap=plt.get_cmap('Dark2'),interpolation='none')
plt.axis('off')
plt.title('Cluster-80th percentile')

plt.subplot(4,1,4)
plt.plot(con1,'k--',label='baseline - high')
plt.plot(con2,'r--',label='scenario - high')
plt.plot(con3,'k',label='baseline - low')
plt.plot(con4,'r',label='scenario - low')
plt.ylabel('Connectivity')
plt.legend(loc=2,fontsize=10,frameon=True,ncol=1,columnspacing=0.2,labelspacing=0.3,handletextpad=0.1) 
plt.xlabel('Low                Percentile                High')



