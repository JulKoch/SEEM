###
# Author: Julian Koch (juko@geus.dk)
# Date: June/2016
###
#########################################################################################
# Empirical-Orthogonal-Functions-Analysis
# Input: two arrays that represent the space-time matrix of a variable. 
# For the spatial EOF-analysis (most common) represent the rows locations in
# space and the columns represent the temporal dimension. The two inputs need
# to have identical dimensions. This functions returns the first three EOFs
# which is usually sufficient, as most of the variance is explained by the 
# first three EOFs. Changes can easily be implemented.
  
import numpy as np
import scipy
from scipy.ndimage import measurements

def EOF_similarity(obs,model):
    array_size=obs.shape
# transform the datasets to its spatial anomalies which is an essential 
# step for the EOF analysis     
    for i in range(0,array_size[1]):
        obs[:,i]=obs[:,i]-np.mean(obs[:,i])
        model[:,i]=model[:,i]-np.mean(model[:,i])        
# Build the integral datamatrix that contains both datasets along the time axis    
    data=np.concatenate((obs,model),axis=1)
# Conduct the EOF via the svd function
    U, s, V = scipy.linalg.svd(data, full_matrices=True)
    Y = np.dot((data),V);
# Compute the amount of explained variance    
    var_exp=np.cumsum(s**2/np.sum(s**2));
# Split the resulting loadings into two    
    load_obs,load_model=np.array_split(np.transpose(V), 2)    
    dif=np.abs(load_obs-load_model)
# reverse cumsum of the explained variance    
    var=var_exp[0]
    for l in range(1,obs.shape[1]*2):
        var=np.append(var,(var_exp[l]-var_exp[l-1]))
# Compute the skill score based on the weighted sum of the absolute loading 
# differences and return skill (0 for perfect agreement). 
    skill=np.empty([1,1])   
    for l in range(0,obs.shape[1]):
        skill=np.append(skill,np.sum(dif[l,:]*var))
    skill=np.delete(skill,0,axis=0)        
# output: skill score, EOF maps 1-3, explained variance, loadings obs, loadings model
    return (skill,Y[:,1],Y[:,2],Y[:,3],var,load_obs,load_model) 
#########################################################################################    
# Connectivity Analysis
# Input: two arrays with the same dimensions at a single time step. 
    
def connectivity(obs,model):   
    # The s array defines the directions of possible connections.     
    s = [[1,1,1],[1,1,1],[1,1,1]] # 8 way connection
    #s = [[0,1,0],[1,1,1],[[0,1,0]] # 4 way connection      
    array_size=obs.shape    
    # The connectivity is computed at each percentile. 
    # First, compute percentiles of each map. Reshape to 1D     
    obs_re=np.reshape(obs,[array_size[0]*array_size[1]])
    model_re=np.reshape(model,[array_size[0]*array_size[1]]) 
    # Delete NANs    
    obs_re=obs_re[~np.isnan(obs_re)]  
    model_re=model_re[~np.isnan(model_re)]
    # Compute percentiles
    per_obs=np.percentile(obs_re,np.linspace(1,100,num=100))
    per_model=np.percentile(model_re,np.linspace(1,100,num=100))
    # The connectivity is calculated for the low and the high phase.
    # The high phase considers values exceeding the threshold percentile.
    # The high phase considers values exceeding the thresholds percentile.     
    #  Initialize: low phase
    # connectivity for each threshold percentile       
    con_obs_low=np.empty([100]) 
    con_model_low=np.empty([100])
    # cluster maps for each threshold percentile
    cl_obs_low=np.empty([array_size[0],array_size[1],100])
    cl_model_low=np.empty([array_size[0],array_size[1],100])
    # Initialize: high phase
    # connectivity for each threshold percentile        
    con_obs_high=np.empty([100])
    con_model_high=np.empty([100])
    # cluster maps for each threshold percentile
    cl_obs_high=np.empty([array_size[0],array_size[1],100])
    cl_model_high=np.empty([array_size[0],array_size[1],100])
    
    # iterate through all percentiles    
    for j in range(0,100):    
        # truncate low phase at percentile j and transfer to binary        
        temp_obs_low=(obs<=per_obs[j]).astype(int) # smaller than threshold
        temp_model_low=(model<=per_model[j]).astype(int)
        # perform the cluster analysis on the binary map.
        # it returns a map of connected clusters where each cluster gets a unique ID
        cl_obs_low[:,:,j], num_obs = measurements.label(temp_obs_low,structure=s)
        cl_model_low[:,:,j], num_obs = measurements.label(temp_model_low,structure=s)
        # this returns the size (number of pixels) of each cluster
        area_obs_low = measurements.sum(temp_obs_low, cl_obs_low[:,:,j], index=np.arange(cl_obs_low[:,:,j].max() + 1))        
        area_model_low = measurements.sum(temp_model_low, cl_model_low[:,:,j], index=np.arange(cl_model_low[:,:,j].max() + 1))
        # the connectivity metric describes the proportion of pairs of cells that are connected among all possible pairs of connected cells        
        con_obs_low[j]=np.sum(np.square(area_obs_low))/np.square(np.sum(area_obs_low))      
        con_model_low[j]=np.sum(np.square(area_model_low))/np.square(np.sum(area_model_low))
        # as a perfromance metric the RMSE is computed for the connectivity of obs and model at each percentile
        out_low=np.sqrt(np.nanmean((con_model_low-con_obs_low)**2))        
        
        # truncate high phase at percentile j and transfer to binary         
        temp_obs_high=(obs>=per_obs[j]).astype(int) # greater than threshold
        temp_model_high=(model>=per_model[j]).astype(int)
        # perform the cluster analysis on the binary map.
        # it returns a map of connected clusters where each cluster gets a unique ID
        cl_obs_high[:,:,j], num_obs = measurements.label(temp_obs_high,structure=s)
        cl_model_high[:,:,j], num_obs = measurements.label(temp_model_high,structure=s)
        # this returns the size (number of pixels) of each cluster
        area_obs_high = measurements.sum(temp_obs_high, cl_obs_high[:,:,j], index=np.arange(cl_obs_high[:,:,j].max() + 1))        
        area_model_high = measurements.sum(temp_model_high, cl_model_high[:,:,j], index=np.arange(cl_model_high[:,:,j].max() + 1))
         # the connectivity metric describes the proportion of pairs of cells that are connected among all possible pairs of connected cells          
        con_obs_high[j]=np.sum(np.square(area_obs_high))/np.square(np.sum(area_obs_high))      
        con_model_high[j]=np.sum(np.square(area_model_high))/np.square(np.sum(area_model_high))
        # as a performance metric the RMSE is computed for the connectivity of obs and model at each percentile
        out_high=np.sqrt(np.nanmean((con_model_high-con_obs_high)**2))
    
    #output: 
    # connectivity at each percentile for obs-high phase
    # connectivity at each percentile for model-high phase
    # connectivity at each percentile for obs-low phase
    # connectivity at each percentile for model-low phase
    # cluster maps at each percentile for obs-high phase
    # cluster maps at each percentile for model-high phase
    # cluster maps at each percentile for obs-low phase
    # cluster maps at each percentile for model-low phase 
    # final connectivity metric for high phase
    # final connectivity metric for low phase
    return (con_obs_high,con_model_high,con_obs_low,con_model_low,cl_obs_high,cl_model_high,cl_obs_low,cl_model_low,out_high,out_low)    
#########################################################################################     