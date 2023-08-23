# Code
########################################

(1) Simulations are run for each type of interaction in separated scripts, and parameters to yield the desired results are tested in an additional script for each interaction type:

**simulation_fac_comp.Rmd**
**test_parameters_fac_comp.Rmd**

**simulation_foodweb.Rmd**
**test_parameters_foodweb.Rmd**

**simulation_neutral.Rmd**
**test_parameters_neutral.Rmd**

(2) Computing the phylogenetic signal for all simulations' data generated in (1) is done in the script:

The script **compute_phylosignal_1simulation.Rmd** computes the phylogenetic signal only one of the simulations and contains explanations for each step.

The script **compute_phylosignal_nsimulations.Rmd** computes the phylogenetic signal for all simulations generated from (1).
