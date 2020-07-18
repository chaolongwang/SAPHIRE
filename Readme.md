## Modelling transmission dynamics of COVID-19

This repository includes data and codes for reproducing the results in the manuscript:

Hao X, Cheng S, Wu D, Wu T, Lin X, and Wang C (2020). Reconstruction of the full transmission of COVID-19 in Wuhan. **Nature**, doi: 10.1038/s41586-020-2554-8. 

[Link]: https://www.nature.com/articles/s41586-020-2554-8

In the following sections, we will describe the purposes of each major directory and the scripts in that folder. **One should keep in mind that the scripts and data here are meant for reproduction of the results in the paper only.** Scripts were tested on R-3.6.x.

**We encourage you to raise purely technical questions through Github issues so that we can answer your questions ASAP.**

### Prerequisite third-party R packages

You may need to install the following R packages if you have not done so yet:

- BayesianTools
- coda
- cairoDevice
- vioplot
- readr
- corrplot
- IDPmisc

### Descriptions of folders, scripts and other files

Scripts in folders `scripts_main`, `scripts_resurgence`  and `scripts_sensitivity`  are meant to be run directly, while scripts in other folders are supportive and not meant to be run directly by users.

#### Folder `scripts_main `

R scripts for our main analyses. Since we constantly need to use functions defined in `R` directory and previous outputs located at `output` directory, it is necessary to set the `code_root` variable at the beginning of the script properly. `code_root` should be set to the directory under which this README file is located. `code_root` should ends with `/`. The same applies to other scripts that are meant to be directly run by users. For example, suppose you have git cloned SAPHIRE at `/home/Sarah/SAPHIRE` and this `Readme.md` file is at `/home/Sarah/SAPHIRE`, then `code_root`  variable should be set to `/home/Sarah/SAPHIRE/`.

- `Run_SEIR_main_analysis.R`: R script to reproduce the main analyses (**Fig. 2**). This R script will call `SEIRfitting` function to perform the analysis. **After the run, please first inspect `output/par_traj_run_main_analysis.png`  visually to make sure the MCMC run has converged.** If convergence has not been achieved, rerun the script with different random seeds, or specify a good initial parameter values.

- `confirm_convergence.R`:  Confirm and test convergence of MCMC by comparing three independent runs. This script will reproduce **Supplementary Fig. 10**. **This requires output from `Run_SEIR_main_analysis.R` . Run  `Run_SEIR_main_analysis.R`  first before running this script.**



#### Folder ` scripts_resurgence`

R scripts for risk of resurgence estimations.

- `Run_resurge_simulation.R`: R script to run the risk of resurgence estimations after control measures is lifted using the parameters from the main and sensitivity s8 analyses. **This requires output from `scripts_main/Run_SEIR_main_analysis.R` and `scripts_sensitivity/Run_SEIR_s8.R`. Run  `Run_SEIR_main_analysis.R` and `scripts_sensitivity/Run_SEIR_s8.R` first before running this script.**

- `Run_resurge_plot_fig3_A.R`: R script to reproduce **Fig. 3A**. **This requires output from `scripts_main/Run_SEIR_main_analysis.R` and `Run_resurge_simulation.R`.**

- `Run_resurge_plot_fig3_BC.R`: R script to reproduce **Fig. 3B/C**. **This requires output from `scripts_sensitivity/Run_SEIR_s8.R` and `Run_resurge_simulation.R`.**



#### Folder `R`
This folder contains major R functions used for our model fitting, parameter estimation and producing results figures. 

- `fun_SEIRpred.R`: Function `SEIRpred` evolves the system according to deterministic model specified by Eqs. 1-7.

- `fun_SEIRsimu.R`: Function `SEIRsimu` evolves the system according to the stochastic model specified by Eqs. 10-16.

- `fun_SEIRfitting.R` : SEIR model fitting and parameter estimation by MCMC with the Delayed Rejection Adaptive Metropolis (DRAM) algorithm.

- `fun_SEIRresurge.R`: Resurgence simulation using the SIER model

- `fun_SEIRplot.R`: The plot function for the figures of main and sensitivity analyses

- `fun_R0estimate.R` : The function to calculate the $R_0$ (basic/effective reproduction number)

- `fun_Findzero.R`: Function `Findzero`, given previous parameter estimation results and corresponding initial conditions, will find dates of zero case, i.e., date when I=0 and when E+P+I+A=0.

- `init_cond.R`: Function `generate_init_condi` creates a list containing all parameters, useful constants, initial conditions of the population set according to the parameters, and several functions that need to be accessed in various R functions.

- `correlationPlot_modified.R`: A modified version of `correlationPlot` in package `BayesianTools`.



#### Folder `scripts_sensitivity`

R scripts for out sensitivity analyses (s1-s9).

- `Run_SEIR_s1.R`: R script to run the sensitivity analyses s1: Adjust the reported incidences from January 29 to February 1 to their average. 
- `Run_SEIR_s2.R`: R script to run the sensitivity analyses s2: Assume an incubation period of 4.1 days (lower 95% CI from reference 2) and presymptomatic infectious period of 1.1 days (lower 95% CI from reference 10 is 0.8 days, but our discrete stochastic model requires $D_p>1$), equivalent to setting $D_e=3$ and $D_p=1.1$, and adjust P(0) and E(0) accordingly.
- `Run_SEIR_s3.R`: R script to run the sensitivity analyses s3: Assume an incubation period of 7 days (upper 95% CI from reference 2) and presymptomatic infectious period of 3 days (upper 95% CI from reference 10), equivalent to set $D_e=4$ and $D_p=3$, and adjust P(0) and E(0) accordingly.
- `Run_SEIR_s4.R`: R script to run the sensitivity analyses s4: Assume the transmissibility of the presymptomatic and unascertained cases is $α=0.46$ (lower 95% CI from reference 15) of the ascertained cases. 
- `Run_SEIR_s5.R`: R script to run the sensitivity analyses s5: Assume the transmissibility of the presymptomatic and unascertained cases is $α=0.62$ (upper 95% CI from reference 15) of the ascertained cases. 
- `Run_SEIR_s6.R`: R script to run the sensitivity analyses s6: Assume the initial ascertainment rate is $r_0=0.14$ (lower 95% CI of the estimate using Singapore data) and adjust A(0), P(0), and E(0) accordingly.
- `Run_SEIR_s7.R`: R script to run the sensitivity analyses s7: Assume the initial ascertainment rate is $r_0=0.42$ (upper 95% CI of the estimate using Singapore data) and adjust A(0), P(0), and E(0) accordingly.
- `Run_SEIR_s8.R`: R script to run the sensitivity analyses s8: Assume the initial ascertainment rate is $r_0=1$ (theoretical upper limit) and adjust A(0), P(0), and E(0) accordingly.
- `Run_SEIR_s9.R`: R script to run the sensitivity analyses s9: Assume no unascertained cases by fixing $r_0=r_{12}=r_3=r_4=r_5=1$. 



#### Folder `data`

This folder contains the main data used in the study. Detailed description of the data can be found in the following paper:

Pan A, Liu L, Wang C, Guo H, Hao X, Wang Q, Huang J, He N, Yu H, Lin X, Wei S, Wu T (2020). Association of Public Health Interventions With the Epidemiology of the COVID-19 Outbreak in Wuhan, China. **JAMA**, 323(19):1915-1923.

- `Covid19CasesWH.csv`: This file contains daily counts of laboratory-confirmed cases with onset between Dec 8, 2019 and Mar 8, 2020.  


#### Folder `output `

This folder stores the results of parameters estimations and the fitting plots. 



