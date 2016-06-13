import os
import numpy as np
import scipy
from scipy import stats
from scipy.ndimage import measurements
from SEEM import EOF_similarity
     
Path='../Data/' 

# read data
mydata=scipy.io.loadmat((Path+'ET_data_2008.mat'))
myET=mydata.get('ET_data_out')
mydata=scipy.io.loadmat((Path+'LST_data_2008'))
myLST=mydata.get('LST_data_out')

print ('read data done!')

myVar=myET

size=myVar.shape
EOF=np.empty([size[2],size[3]-1])


for i in range(0,2):   
        #EOF-prepare data   
        obs=np.reshape(myVar[:,:,:,0],(size[0]*size[1],size[2]))
        obs=obs[~np.isnan(obs).any(axis=1)]
        model=np.reshape(myVar[:,:,:,i+1],(size[0]*size[1],size[2]))
        model=model[~np.isnan(model).any(axis=1)]    

        # EOF
        EOF[:,i]=EOF_similarity(obs,model)

        print(('Scenario '+str(i)+' EOF done'))
        
    
   