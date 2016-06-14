###
# Author: Julian Koch (juko@geus.dk)
# Date: June/2016
###

# Empirical-Orthogonal-Functions-Analysis
# Input: two arrays that represent the space-time matrix of a variable. 
# For the spatial EOF-analysis (most common) represent the rows locations in
# space and the columns represent the temporal dimension. The two inputs need
# to have identical dimensions. This functions returns the first three EOFs
# which is usually sufficient, as most of the variance is explained by the 
# first three EOFs. Changes can easily be implemented.
  
import numpy as np
import scipy

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
     