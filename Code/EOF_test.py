import scipy
import numpy as np
from SEEM import EOF_similarity
import matplotlib.pyplot as plt 
     
Path='../Data/' 

# read data
mydata=scipy.io.loadmat((Path+'ET_data_2008.mat'))
myET=mydata.get('ET_data_out')
mydata=scipy.io.loadmat((Path+'LST_data_2008'))
myLST=mydata.get('LST_data_out')

print ('read data done!')

myVar=myET
myScen=1

size=myVar.shape
EOF=np.empty([size[2]])

#EOF-prepare data   
obs=np.reshape(myVar[:,:,:,0],(size[0]*size[1],size[2]))
mask=obs[:,0]
obs=obs[~np.isnan(obs).any(axis=1)]
model=np.reshape(myVar[:,:,:,myScen],(size[0]*size[1],size[2]))
model=model[~np.isnan(model).any(axis=1)]    

# EOF
a=EOF_similarity(obs,model)
EOF[:]=a[0]
eof1=a[1]
eof2=a[2]
eof3=a[3]
var=a[4]
load_obs=a[5]
load_model=a[6]
print(('Scenario '+str(myScen)+' EOF done'))
        

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


x = [0,0+121,0+242,365]
xlabels = ['jan08','may08','sep08','jan09']

plt.figure(figsize=(12,6))    
    
plt.subplots_adjust(hspace=0.2,wspace=0.2) 
plt.subplot(2,3,1)
plt.imshow(eof1_map,cmap='RdGy',interpolation='none')
plt.title('EOF1')
plt.axis('off')
plt.colorbar(orientation ='vertical',ticks=[-1.5,0,1.5],fraction=0.03)
plt.clim(-1.5,1.5)
plt.subplot(2,3,2)
plt.imshow(eof2_map,cmap='RdGy',interpolation='none')
plt.title('EOF2')
plt.axis('off')
plt.colorbar(orientation ='vertical',ticks=[-0.8,0,0.8],fraction=0.03)
plt.clim(-0.8,0.8)
plt.subplot(2,3,3)
plt.imshow(eof3_map,cmap='RdGy',interpolation='none')
plt.title('EOF3')
plt.axis('off')
plt.colorbar(orientation ='vertical',ticks=[-1,0,1],fraction=0.03)
plt.clim(-1,1)

plt.subplot(2,3,4)
plt.plot(load_obs[:,0],'k')
plt.plot(load_model[:,0],'r')
plt.yticks([-0.15,-0.075,0,0.075,0.15])
plt.xticks(x,xlabels,rotation=45)
plt.axis((0,365,-0.15,0.15))
plt.ylabel('EOF Loadings')
plt.subplot(2,3,5)
plt.plot(load_obs[:,1],'k')
plt.plot(load_model[:,1],'r')
plt.yticks([-0.1,-0.05,0,0.05,0.1])
plt.xticks(x,xlabels,rotation=45)
plt.axis((0,365,-0.1,0.1))
plt.clim(-1.2,1.2)
plt.subplot(2,3,6)
plt.plot(load_obs[:,2],'k',label='baseline')
plt.plot(load_model[:,2],'r',label='scenario')
plt.yticks([-0.2,-0.1,0,0.1,0.2])
plt.xticks(x,xlabels,rotation=45)
plt.axis((0,365,-0.2,0.2))
plt.legend(loc=4,fontsize=10,frameon=False,ncol=2,columnspacing=0.2,labelspacing=0.1,handletextpad=0.1)         
    
   