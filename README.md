# SEEM
## Spatial Evaluation of Environmental Models
***
Keywords: 
- Spatial Patterns
- Model Evaluation
- Remote Sensing
- Spatial Statistics
- Human Perception
- Visual Comparison

***
This repository contains methodologies that can be incorporated as spatial performance metrics to evaluate spatial model output. Simulated patterns are compared with reference patterns (e.g. remote sensing observations) and the spatial performance is quantified by a skill score.

*14/06/2016*

The EOF example is up and running. Detailed documentation will be added. In the meantime please check Koch et al. (2015) in the documentation folder for reference.

*23/06/2016*

The Connectivity analysis example is up and running. Detailed documentation will be added. The first example focus on the illustration of the methodology where the connectivity analysis is applied on a single map comparison. The second example applies the connectivity analysis on a series of maps and gives a similarity score for each timestep. Here the connectivity is assessed individually for the high and the low phase. In the meantime please check Koch et al. (2016) in the documentation folder for reference.
 
*22/07/2016*

Based on Koch et al. (2015), located in the literature folder, we have uploaded the 12 synthetic land-surface-temperature maps and the corresponding reference map. The 12 maps are perturbation from the reference map and a detailed description of the perturbations strategies can be found in the paper. A web-based survey of the human perception was set up (https://kwiksurveys.com/s.asp?sid=05a50l3t1jb4eyp326466#/) to employ the well trained human perception to rate the pattern similarity between the 12 maps and the reference. Data is located under the "HumanPerceptionSurvey_I" folder and the file ranking.txt contains the final results for each map; a value of 1 indicates complete similarity. The ranking can be utilized to evaluate spatial performance in their ability to mimic the human perception. 

*29/07/2016*

The Fractions Skill Score (FSS) analysis has been added to the repository. The code follows the description by Roberts and Lean (2008). The first example FSS_example_single_day computes and plots FSS for various scales and percentile thresholds at a single day. The second example calculates FSS at critical pairs of threshold and scale for every year in a one year period. The calculation is limited to specified pairs in order to reduce computational time.  

*22/02/2017*
We have uploaded the data obtained from our zooniverse project Pattern Perception (https://www.zooniverse.org/projects/jukoch/pattern-perception/home). The data is placed in the folder "HumanPerceptionSurvey_II". For each variable, the maps are made available as .mat files, where the four dimensions refer to X, Y, nday, nscenario. nscenario=1 refers to the baseline map which is used as reference. The resulting similarity scores for the six scenarios (nscenario=2-7) states the spatial similarity of that respective scenario to the baseline given for each of the 365 days.  

## Work in progress!!
