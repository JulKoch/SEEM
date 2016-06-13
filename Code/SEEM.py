import os
import numpy as np
import scipy
from scipy import stats
from scipy.ndimage import measurements
 
def EOF_similarity(obs,model):
    array_size=obs.shape
    for i in range(0,array_size[1]):
        obs[:,i]=obs[:,i]-np.mean(obs[:,i])
        model[:,i]=model[:,i]-np.mean(model[:,i])        
    data=np.concatenate((obs,model),axis=1)
    U, s, V = scipy.linalg.svd(data, full_matrices=True)
    #Y = np.dot((data),U);
    var_exp=np.cumsum(s**2/np.sum(s**2));
    #U=-1*U
    load_obs,load_model=np.array_split(np.transpose(V), 2)    
    dif=np.abs(load_obs-load_model)
    var=var_exp[0]
    for l in range(1,obs.shape[1]*2):
        var=np.append(var,(var_exp[l]-var_exp[l-1]))
    skill=np.empty([1,1])   
    for l in range(0,obs.shape[1]):
        skill=np.append(skill,np.sum(dif[l,:]*var))
    skill=np.delete(skill,0,axis=0)        
    return skill 
     
    
#Where am I? 
Path='D:/SpatialSensitivity/analysis/data_final/' 
os.chdir(Path)
# dir files
files=os.listdir(Path)
# number of scenarios, each scenarios comes with 4 files
n_scen=6

# read baseline data
i=0
mydata=scipy.io.loadmat((Path+files[i]))
myET_avg=mydata.get('ETday_avg_1km')
mydata=scipy.io.loadmat((Path+files[i+1]))
myLST_avg=mydata.get('LSTday_avg_1km')
mydata=scipy.io.loadmat((Path+files[i+2]))
myLST_rise=mydata.get('LST_rise_1km')

# extend the arrays to concatenate in 4th dimension
myET_avg=myET_avg[...,np.newaxis]
myLST_avg=myLST_avg[...,np.newaxis]
myLST_rise=myLST_rise[...,np.newaxis]

# start with the first scenario and read data from each of them

for i in [18,27,75,78,81,84]:
    mydata=scipy.io.loadmat((Path+files[i]))
    myET_avg_temp=mydata.get('ETday_avg_1km')
    mydata=scipy.io.loadmat((Path+files[i+1]))
    myLST_avg_temp=mydata.get('LSTday_avg_1km')
    mydata=scipy.io.loadmat((Path+files[i+2]))
    myLST_rise_temp=mydata.get('LST_rise_1km')
  
    
    myET_avg=np.concatenate((myET_avg,myET_avg_temp[...,np.newaxis]),axis=3)
    myLST_avg=np.concatenate((myLST_avg,myLST_avg_temp[...,np.newaxis]),axis=3)
    myLST_rise=np.concatenate((myLST_rise,myLST_rise_temp[...,np.newaxis]),axis=3)

print ('read data done!')

size=myET_avg.shape
EOF=np.empty([size[2],size[3],3])




for i in range(0,n_scen+1):   
        #EOF-prepare data
        #ET_avg    
        obsET_avg=np.reshape(myET_avg[:,:,:,0],(size[0]*size[1],size[2]))
        obsET_avg=obsET_avg[~np.isnan(obsET_avg).any(axis=1)]
        modelET_avg=np.reshape(myET_avg[:,:,:,i],(size[0]*size[1],size[2]))
        modelET_avg=modelET_avg[~np.isnan(modelET_avg).any(axis=1)]    
        #LST_avg    
        obsLST_avg=np.reshape(myLST_avg[:,:,:,0],(size[0]*size[1],size[2]))
        obsLST_avg=obsLST_avg[~np.isnan(obsLST_avg).any(axis=1)]
        modelLST_avg=np.reshape(myLST_avg[:,:,:,i],(size[0]*size[1],size[2]))
        modelLST_avg=modelLST_avg[~np.isnan(modelLST_avg).any(axis=1)]
        #LST_rise  
        obsLST_rise=np.reshape(myLST_rise[:,:,:,0],(size[0]*size[1],size[2]))
        obsLST_rise=obsLST_rise[~np.isnan(obsLST_rise).any(axis=1)]
        modelLST_rise=np.reshape(myLST_rise[:,:,:,i],(size[0]*size[1],size[2]))
        modelLST_rise=modelLST_rise[~np.isnan(modelLST_rise).any(axis=1)]    
        # EOF
        EOF[:,i,0]=EOF_similarity(obsET_avg,modelET_avg)
        EOF[:,i,1]=EOF_similarity(obsLST_avg,modelLST_avg)
        EOF[:,i,2]=EOF_similarity(obsLST_rise,modelLST_rise)
        print(('Scenario '+str(i)+' EOF done'))
        

# smoothing

metric_raw=EOF
metric_smooth=np.empty([size[2],size[3],3])


def movingaverage(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    out=np.convolve(interval, window, 'same')
    out[0]=np.nanmean(interval[0:7])
    out[1]=np.nanmean(interval[0:8]) 
    out[2]=np.nanmean(interval[0:8])
    out[1092]=np.nanmean(interval[1084:1092])
    out[1093]=np.nanmean(interval[1085:1093]) 
    out[1094]=np.nanmean(interval[1087:1094])      
    #out[0:3]=interval[0:3]
    #out[1092:1095]=interval[1092:1095]
    return out

for i in range(0,3):# for each variable
    for j in range(0,size[3]):# for each scenario
        metric_smooth[:,j,i]=movingaverage(metric_raw[:,j,i],7)

np.save('D:/SpatialPatterns/SurveyII/analysis/metrics/out/EOF_smooth.npy', metric_smooth)

np.save('D:/SpatialPatterns/SurveyII/analysis/metrics/out/EOF.npy', metric_raw)




    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
