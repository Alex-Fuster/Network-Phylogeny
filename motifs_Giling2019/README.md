## Food web Motif Analysis
Code for Giling et al. "Plant diversity alters the representation of motifs in food webs"

### 01_Giling_et_al_build_and_analyse_networks.R
This script first constructs the matrices of interaction probability among species (plants, consumers and static resources) that co-occurred on an experimental plot at a particular time. These networks are constructed for multiple scenarios of interaction probability by altering the parameters alpha and beta (see Methods section "Food web construction"). 

Subsequently, the code samples these probabilities a specified number of times and calculates a selection of network properties (with the custom 'web.properties' function defined in section 2). This includes the motif_counter function of Borrelli 2015. Each network is also rewired a specified number of times with the functions 'curving' (modified from Borrelli 2015) and curve_ball (Strona et al. 2014) and z-scores are calculated. For motifs, the z-scores are normalised.

There are optional lines of code that can be included to run the habitat complexity scenario (see Supplementary Note 2) or to subset the network to consumers only.

The looping proceedure is implemented for running in parallel with packages foreach and doParallel.


### 02_Giling_et_al_motif_analysis.R
This code uses the network properties generated by code "01_Giling_et_al_build_and_analyse_networks.R" to test the effect of plant species richness on the counts (section 3) and z-score of motifs (section 4). This is performed by the functions "plot.freq" and "plot.zscores", respectively. These functions plot the data and predicted relationships and save the results of the linear mixed models.

The data to accompany this script and reproduce the analyses presented in the manuscript are available on figshare (DOI: 10.6084/m9.figshare.7605242).


### 03_Giling_et_al_motif_analysis_Supplementary.R
This code uses the network properties generated by code "01_Giling_et_al_build_and_analyse_networks.R" to examine the effects of parameters a and b on network properties and motif z score (section 3), analyse trends in consumer richness (section 4), and analyse the results of the scenario designed to test sensitivity to potential changes in habitat complexity (section 5).

The data to accompany this script and reproduce the analyses presented in the manuscript are available on figshare (DOI: 10.6084/m9.figshare.7605242).