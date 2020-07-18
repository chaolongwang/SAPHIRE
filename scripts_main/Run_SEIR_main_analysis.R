## IMPORTANT: Please set code_root variable properly. 
## code_root should be set to the directory where the repository README file is located. 
## For more information, please read the repository README file
code_root="~/SAPHIRE/"

setwd(paste0(code_root, "scripts_main"))
library(BayesianTools)
library(vioplot)
library("corrplot")
library(readr)
library(cairoDevice)

##
source(paste0(code_root, "R/fun_SEIRpred.R"))
source(paste0(code_root, "R/fun_SEIRsimu.R"))
source(paste0(code_root, "R/fun_SEIRfitting.R"))
source(paste0(code_root, "R/init_cond.R"))
source(paste0(code_root, "R/fun_R0estimate.R"))
source(paste0(code_root, "R/correlationPlot_modified.R"))
source(paste0(code_root, "R/fun_SEIRplot.R"))
source(paste0(code_root, "R/fun_Findzero.R"))
##

init_sets_list=get_init_sets_list(r0 = 0.23)

# good initial conditions
# c(1.284, 0.384, 0.174, 0.096, 0.161, -0.046, -0.379, 0.569)

SEIRfitting(init_sets_list, randomize_startValue = T,
            run_id = "main_analysis", output_ret = T, skip_MCMC=F)

## to evaluate convergence, we run another two rounds of this program
SEIRfitting(init_sets_list, randomize_startValue = T, run_id = "main_analysis_rep1", output_ret = T, skip_MCMC=F)
SEIRfitting(init_sets_list, randomize_startValue = T, run_id = "main_analysis_rep2", output_ret = T, skip_MCMC=F)
