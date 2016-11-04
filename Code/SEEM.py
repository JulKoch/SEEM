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
    return (skill,Y[:,0],Y[:,1],Y[:,2],var,load_obs,load_model) 
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
# Fractions Skill Score (FSS) functions

# in order to calculate FSS at large scale it makes it easier if we just extend the original array in all directions, 
# insert the original data and fill the remaining grids with NaN values. 
def extend_array(myarray):
   maxlen=np.max(myarray.shape) 
   newarray=np.empty([myarray.shape[0]+2*maxlen,myarray.shape[1]+2*maxlen])
   newarray[:]=np.nan   
   newarray[maxlen:maxlen+myarray.shape[0],maxlen:maxlen+myarray.shape[1]]=myarray    
   return newarray

# This function counts the number of possible neighbors at a given scale n (window size). 
# It takes into consideration the catchment morphology, grids close to the boundary will have less possible neighbors.     
# This function is most meaningful with uneven scales.
def count_neigh(myarray, scale): 
   newarray=np.empty([myarray.shape[0],myarray.shape[1]])
   newarray[:]=np.nan
   mask=np.isnan(myarray)
   myarray=np.abs(np.isnan(myarray).astype(int)-1) # convert into binary
   if scale==1: # if scale ==1 then there are no neighbors and we can use the original array.
       newarray=myarray
   else:
       myarray=extend_array(myarray) # extend the array       
       maxlen=np.max(newarray.shape)        
       d=(scale-1)/2 # the scale is divided by two and the original grid is placed the center. 
       for o in range(maxlen,maxlen+newarray.shape[0]): # loop through all grids
           for p in range(maxlen,maxlen+newarray.shape[1]):
               sub=myarray[o-d-1:o+d+1,p-d-1:p+d+1]
               newarray[o-maxlen,p-maxlen]=np.nansum(sub) # count the neighbors within +/-d
   newarray = newarray.astype(float)            
   newarray[mask]=np.nan # grids outside the boundary are set to NaN 
   return newarray # the output is an array of the same size as the input array and contains the neighbor count at a scale for each grid   
    
# This function searches for neighbors that belong to a certain threshold percentile (perc) class and fall into a window of size scale.
# If a threshold percentile is below 50 it will automatically use the grids that fall below that value: E.g. perc=20 will focus on the lowest 20%.
# Opposed, perc > 50 will focus on grids exceeding that threshold: E.g. perc=90 focuses on the top 10% of grids.    
def search_neigh(myarray,scale,perc):
    newarray=np.empty([myarray.shape[0],myarray.shape[1]])
    newarray[:]=np.nan
    mask=np.isnan(myarray)    
    t=np.nanpercentile(myarray,perc)
    if perc>50: # check if perc is greater or small than 50   
        myarray=(myarray>=t).astype(int) # # if smaller the function will focus on grids falling below the threshold.  
    else:
        myarray=(myarray<=t).astype(int) # if smaller the function will focus on grids falling below the threshold. 
    if scale==1: # if scale ==1 then there are no neighbors and we can use the original array.
       newarray=myarray
    else:
       myarray=extend_array(myarray) # extend the array
       maxlen=np.max(newarray.shape) 
       d=(scale-1)/2 # the scale is divided by two and the original grid is placed the center.
       for o in range(maxlen,maxlen+newarray.shape[0]): # loop through all grids
           for p in range(maxlen,maxlen+newarray.shape[1]):
               sub=myarray[o-d-1:o+d+1,p-d-1:p+d+1]
               newarray[o-maxlen,p-maxlen]=np.nansum(sub)
    newarray = newarray.astype(float)           
    newarray[mask]=np.nan  # grids outside the boundary are set to NaN 
    return newarray # the output is an array of the same size as the input array and contains the neighbor count at a scale for each grid constrain to a threshold percentile.     
#########################################################################################  
# parula color map
# This is my preferred color map and it is used by e.g. matlab as standard.
# It is not yet part of matplotlib and therefore needs to be imported separately.

from matplotlib.colors import LinearSegmentedColormap

parameters = {'xp': [8.6840042204210732, -21.378780072752193, -13.715204596193757, -94.858944936224646, 53.679179741776352, -9.6],
              'yp': [-36.497267048261136, -44.616177228232843, 10.381246780010315, 84.763008758371967, 6.324059763008762, 41],
              'min_JK': 25.4069767442,
              'max_JK': 95}

cm_data = [[ 0.26710521,  0.03311059,  0.6188155 ],
       [ 0.26493929,  0.04780926,  0.62261795],
       [ 0.26260545,  0.06084214,  0.62619176],
       [ 0.26009691,  0.07264411,  0.62951561],
       [ 0.25740785,  0.08360391,  0.63256745],
       [ 0.25453369,  0.09395358,  0.63532497],
       [ 0.25147146,  0.10384228,  0.6377661 ],
       [ 0.24822014,  0.11337029,  0.6398697 ],
       [ 0.24478105,  0.12260661,  0.64161629],
       [ 0.24115816,  0.131599  ,  0.6429888 ],
       [ 0.23735836,  0.14038009,  0.64397346],
       [ 0.23339166,  0.14897137,  0.64456048],
       [ 0.22927127,  0.15738602,  0.64474476],
       [ 0.22501278,  0.16563165,  0.64452595],
       [ 0.22063349,  0.17371215,  0.64390834],
       [ 0.21616055,  0.18162302,  0.64290515],
       [ 0.21161851,  0.18936156,  0.64153295],
       [ 0.20703353,  0.19692415,  0.63981287],
       [ 0.20243273,  0.20430706,  0.63776986],
       [ 0.19784363,  0.211507  ,  0.63543183],
       [ 0.19329361,  0.21852157,  0.63282872],
       [ 0.18880937,  0.2253495 ,  0.62999156],
       [ 0.18442119,  0.23198815,  0.62695569],
       [ 0.18014936,  0.23844124,  0.62374886],
       [ 0.17601569,  0.24471172,  0.62040016],
       [ 0.17204028,  0.25080356,  0.61693715],
       [ 0.16824123,  0.25672163,  0.6133854 ],
       [ 0.16463462,  0.26247158,  0.60976836],
       [ 0.16123449,  0.26805963,  0.60610723],
       [ 0.15805279,  0.27349243,  0.60242099],
       [ 0.15509948,  0.27877688,  0.59872645],
       [ 0.15238249,  0.28392004,  0.59503836],
       [ 0.14990781,  0.28892902,  0.59136956],
       [ 0.14767951,  0.29381086,  0.58773113],
       [ 0.14569979,  0.29857245,  0.58413255],
       [ 0.1439691 ,  0.30322055,  0.58058191],
       [ 0.14248613,  0.30776167,  0.57708599],
       [ 0.14124797,  0.31220208,  0.57365049],
       [ 0.14025018,  0.31654779,  0.57028011],
       [ 0.13948691,  0.32080454,  0.5669787 ],
       [ 0.13895174,  0.32497744,  0.56375063],
       [ 0.13863958,  0.32907012,  0.56060453],
       [ 0.138537  ,  0.3330895 ,  0.55753513],
       [ 0.13863384,  0.33704026,  0.55454374],
       [ 0.13891931,  0.34092684,  0.55163126],
       [ 0.13938212,  0.34475344,  0.54879827],
       [ 0.14001061,  0.34852402,  0.54604503],
       [ 0.14079292,  0.35224233,  0.54337156],
       [ 0.14172091,  0.35590982,  0.54078769],
       [ 0.14277848,  0.35953205,  0.53828312],
       [ 0.14395358,  0.36311234,  0.53585661],
       [ 0.1452346 ,  0.36665374,  0.5335074 ],
       [ 0.14661019,  0.3701591 ,  0.5312346 ],
       [ 0.14807104,  0.37363011,  0.52904278],
       [ 0.1496059 ,  0.3770697 ,  0.52692951],
       [ 0.15120289,  0.3804813 ,  0.52488853],
       [ 0.15285214,  0.38386729,  0.52291854],
       [ 0.15454421,  0.38722991,  0.52101815],
       [ 0.15627225,  0.39056998,  0.5191937 ],
       [ 0.15802555,  0.39389087,  0.5174364 ],
       [ 0.15979549,  0.39719482,  0.51574311],
       [ 0.16157425,  0.40048375,  0.51411214],
       [ 0.16335571,  0.40375871,  0.51254622],
       [ 0.16513234,  0.40702178,  0.51104174],
       [ 0.1668964 ,  0.41027528,  0.50959299],
       [ 0.16864151,  0.41352084,  0.50819797],
       [ 0.17036277,  0.41675941,  0.50685814],
       [ 0.1720542 ,  0.41999269,  0.50557008],
       [ 0.17370932,  0.42322271,  0.50432818],
       [ 0.17532301,  0.42645082,  0.50313007],
       [ 0.17689176,  0.42967776,  0.50197686],
       [ 0.17841013,  0.43290523,  0.5008633 ],
       [ 0.17987314,  0.43613477,  0.49978492],
       [ 0.18127676,  0.43936752,  0.49873901],
       [ 0.18261885,  0.44260392,  0.49772638],
       [ 0.18389409,  0.44584578,  0.49673978],
       [ 0.18509911,  0.44909409,  0.49577605],
       [ 0.18623135,  0.4523496 ,  0.494833  ],
       [ 0.18728844,  0.45561305,  0.49390803],
       [ 0.18826671,  0.45888565,  0.49299567],
       [ 0.18916393,  0.46216809,  0.49209268],
       [ 0.18997879,  0.46546084,  0.49119678],
       [ 0.19070881,  0.46876472,  0.49030328],
       [ 0.19135221,  0.47208035,  0.48940827],
       [ 0.19190791,  0.47540815,  0.48850845],
       [ 0.19237491,  0.47874852,  0.4876002 ],
       [ 0.19275204,  0.48210192,  0.48667935],
       [ 0.19303899,  0.48546858,  0.48574251],
       [ 0.19323526,  0.48884877,  0.48478573],
       [ 0.19334062,  0.49224271,  0.48380506],
       [ 0.19335574,  0.49565037,  0.4827974 ],
       [ 0.19328143,  0.49907173,  0.48175948],
       [ 0.19311664,  0.50250719,  0.48068559],
       [ 0.192864  ,  0.50595628,  0.47957408],
       [ 0.19252521,  0.50941877,  0.47842186],
       [ 0.19210087,  0.51289469,  0.47722441],
       [ 0.19159194,  0.516384  ,  0.47597744],
       [ 0.19100267,  0.51988593,  0.47467988],
       [ 0.19033595,  0.52340005,  0.47332894],
       [ 0.18959113,  0.5269267 ,  0.47191795],
       [ 0.18877336,  0.530465  ,  0.47044603],
       [ 0.18788765,  0.53401416,  0.46891178],
       [ 0.18693822,  0.53757359,  0.46731272],
       [ 0.18592276,  0.54114404,  0.46563962],
       [ 0.18485204,  0.54472367,  0.46389595],
       [ 0.18373148,  0.5483118 ,  0.46207951],
       [ 0.18256585,  0.55190791,  0.4601871 ],
       [ 0.18135481,  0.55551253,  0.45821002],
       [ 0.18011172,  0.55912361,  0.45615277],
       [ 0.17884392,  0.56274038,  0.45401341],
       [ 0.17755858,  0.56636217,  0.45178933],
       [ 0.17625543,  0.56998972,  0.44946971],
       [ 0.174952  ,  0.57362064,  0.44706119],
       [ 0.17365805,  0.57725408,  0.44456198],
       [ 0.17238403,  0.58088916,  0.4419703 ],
       [ 0.17113321,  0.58452637,  0.43927576],
       [ 0.1699221 ,  0.58816399,  0.43648119],
       [ 0.1687662 ,  0.5918006 ,  0.43358772],
       [ 0.16767908,  0.59543526,  0.43059358],
       [ 0.16667511,  0.59906699,  0.42749697],
       [ 0.16575939,  0.60269653,  0.42428344],
       [ 0.16495764,  0.6063212 ,  0.42096245],
       [ 0.16428695,  0.60993988,  0.41753246],
       [ 0.16376481,  0.61355147,  0.41399151],
       [ 0.16340924,  0.61715487,  0.41033757],
       [ 0.16323549,  0.62074951,  0.40656329],
       [ 0.16326148,  0.62433443,  0.40266378],
       [ 0.16351136,  0.62790748,  0.39864431],
       [ 0.16400433,  0.63146734,  0.39450263],
       [ 0.16475937,  0.63501264,  0.39023638],
       [ 0.16579502,  0.63854196,  0.38584309],
       [ 0.16712921,  0.64205381,  0.38132023],
       [ 0.168779  ,  0.64554661,  0.37666513],
       [ 0.17075915,  0.64901912,  0.37186962],
       [ 0.17308572,  0.65246934,  0.36693299],
       [ 0.1757732 ,  0.65589512,  0.36185643],
       [ 0.17883344,  0.65929449,  0.3566372 ],
       [ 0.18227669,  0.66266536,  0.35127251],
       [ 0.18611159,  0.66600553,  0.34575959],
       [ 0.19034516,  0.66931265,  0.34009571],
       [ 0.19498285,  0.67258423,  0.3342782 ],
       [ 0.20002863,  0.67581761,  0.32830456],
       [ 0.20548509,  0.67900997,  0.3221725 ],
       [ 0.21135348,  0.68215834,  0.31587999],
       [ 0.2176339 ,  0.68525954,  0.30942543],
       [ 0.22432532,  0.68831023,  0.30280771],
       [ 0.23142568,  0.69130688,  0.29602636],
       [ 0.23893914,  0.69424565,  0.28906643],
       [ 0.2468574 ,  0.69712255,  0.28194103],
       [ 0.25517514,  0.69993351,  0.27465372],
       [ 0.26388625,  0.70267437,  0.26720869],
       [ 0.27298333,  0.70534087,  0.25961196],
       [ 0.28246016,  0.70792854,  0.25186761],
       [ 0.29232159,  0.71043184,  0.2439642 ],
       [ 0.30253943,  0.71284765,  0.23594089],
       [ 0.31309875,  0.71517209,  0.22781515],
       [ 0.32399522,  0.71740028,  0.21959115],
       [ 0.33520729,  0.71952906,  0.21129816],
       [ 0.3467003 ,  0.72155723,  0.20298257],
       [ 0.35846225,  0.72348143,  0.19466318],
       [ 0.3704552 ,  0.72530195,  0.18639333],
       [ 0.38264126,  0.72702007,  0.17822762],
       [ 0.39499483,  0.72863609,  0.17020921],
       [ 0.40746591,  0.73015499,  0.1624122 ],
       [ 0.42001969,  0.73158058,  0.15489659],
       [ 0.43261504,  0.73291878,  0.14773267],
       [ 0.44521378,  0.73417623,  0.14099043],
       [ 0.45777768,  0.73536072,  0.13474173],
       [ 0.47028295,  0.73647823,  0.1290455 ],
       [ 0.48268544,  0.73753985,  0.12397794],
       [ 0.49497773,  0.73854983,  0.11957878],
       [ 0.5071369 ,  0.73951621,  0.11589589],
       [ 0.51913764,  0.74044827,  0.11296861],
       [ 0.53098624,  0.74134823,  0.11080237],
       [ 0.5426701 ,  0.74222288,  0.10940411],
       [ 0.55417235,  0.74308049,  0.10876749],
       [ 0.56550904,  0.74392086,  0.10885609],
       [ 0.57667994,  0.74474781,  0.10963233],
       [ 0.58767906,  0.74556676,  0.11105089],
       [ 0.59850723,  0.74638125,  0.1130567 ],
       [ 0.609179  ,  0.74719067,  0.11558918],
       [ 0.61969877,  0.74799703,  0.11859042],
       [ 0.63007148,  0.74880206,  0.12200388],
       [ 0.64030249,  0.74960714,  0.12577596],
       [ 0.65038997,  0.75041586,  0.12985641],
       [ 0.66034774,  0.75122659,  0.1342004 ],
       [ 0.67018264,  0.75203968,  0.13876817],
       [ 0.67990043,  0.75285567,  0.14352456],
       [ 0.68950682,  0.75367492,  0.14843886],
       [ 0.69900745,  0.75449768,  0.15348445],
       [ 0.70840781,  0.75532408,  0.15863839],
       [ 0.71771325,  0.75615416,  0.16388098],
       [ 0.72692898,  0.75698787,  0.1691954 ],
       [ 0.73606001,  0.75782508,  0.17456729],
       [ 0.74511119,  0.75866562,  0.17998443],
       [ 0.75408719,  0.75950924,  0.18543644],
       [ 0.76299247,  0.76035568,  0.19091446],
       [ 0.77183123,  0.76120466,  0.19641095],
       [ 0.78060815,  0.76205561,  0.20191973],
       [ 0.78932717,  0.76290815,  0.20743538],
       [ 0.79799213,  0.76376186,  0.21295324],
       [ 0.8066067 ,  0.76461631,  0.21846931],
       [ 0.81517444,  0.76547101,  0.22398014],
       [ 0.82369877,  0.76632547,  0.2294827 ],
       [ 0.832183  ,  0.7671792 ,  0.2349743 ],
       [ 0.8406303 ,  0.76803167,  0.24045248],
       [ 0.84904371,  0.76888236,  0.24591492],
       [ 0.85742615,  0.76973076,  0.25135935],
       [ 0.86578037,  0.77057636,  0.25678342],
       [ 0.87410891,  0.77141875,  0.2621846 ],
       [ 0.88241406,  0.77225757,  0.26755999],
       [ 0.89070781,  0.77308772,  0.27291122],
       [ 0.89898836,  0.77391069,  0.27823228],
       [ 0.90725475,  0.77472764,  0.28351668],
       [ 0.91550775,  0.77553893,  0.28875751],
       [ 0.92375722,  0.7763404 ,  0.29395046],
       [ 0.9320227 ,  0.77712286,  0.29909267],
       [ 0.94027715,  0.7779011 ,  0.30415428],
       [ 0.94856742,  0.77865213,  0.3091325 ],
       [ 0.95686038,  0.7793949 ,  0.31397459],
       [ 0.965222  ,  0.7800975 ,  0.31864342],
       [ 0.97365189,  0.78076521,  0.32301107],
       [ 0.98227405,  0.78134549,  0.32678728],
       [ 0.99136564,  0.78176999,  0.3281624 ],
       [ 0.99505988,  0.78542889,  0.32106514],
       [ 0.99594185,  0.79046888,  0.31648808],
       [ 0.99646635,  0.79566972,  0.31244662],
       [ 0.99681528,  0.80094905,  0.30858532],
       [ 0.9970578 ,  0.80627441,  0.30479247],
       [ 0.99724883,  0.81161757,  0.30105328],
       [ 0.99736711,  0.81699344,  0.29725528],
       [ 0.99742254,  0.82239736,  0.29337235],
       [ 0.99744736,  0.82781159,  0.28943391],
       [ 0.99744951,  0.83323244,  0.28543062],
       [ 0.9973953 ,  0.83867931,  0.2812767 ],
       [ 0.99727248,  0.84415897,  0.27692897],
       [ 0.99713953,  0.84963903,  0.27248698],
       [ 0.99698641,  0.85512544,  0.26791703],
       [ 0.99673736,  0.86065927,  0.26304767],
       [ 0.99652358,  0.86616957,  0.25813608],
       [ 0.99622774,  0.87171946,  0.25292044],
       [ 0.99590494,  0.87727931,  0.24750009],
       [ 0.99555225,  0.88285068,  0.2418514 ],
       [ 0.99513763,  0.8884501 ,  0.23588062],
       [ 0.99471252,  0.89405076,  0.2296837 ],
       [ 0.99421873,  0.89968246,  0.2230963 ],
       [ 0.99370185,  0.90532165,  0.21619768],
       [ 0.99313786,  0.91098038,  0.2088926 ],
       [ 0.99250707,  0.91666811,  0.20108214],
       [ 0.99187888,  0.92235023,  0.19290417],
       [ 0.99110991,  0.92809686,  0.18387963],
       [ 0.99042108,  0.93379995,  0.17458127],
       [ 0.98958484,  0.93956962,  0.16420166],
       [ 0.98873988,  0.94533859,  0.15303117],
       [ 0.98784836,  0.95112482,  0.14074826],
       [ 0.98680727,  0.95697596,  0.12661626]]

parula_cm = LinearSegmentedColormap.from_list(__file__, cm_data)   
