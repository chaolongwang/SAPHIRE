rm(list = ls())
## IMPORTANT: Please set code_root variable properly. 
## code_root should be set to the directory where the repository README file is located. 
## For more information, please read the repository README file
code_root="~/SAPHIRE/"

setwd(paste0(code_root, "scripts_main"))
library(BayesianTools)
library(vioplot)
library("corrplot")
library(readr)

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

# flatten daily new cases on Jan 29, 30, 31, Feb 1
init_sets_list$daily_new_case_all[29:32] = round(mean(init_sets_list$daily_new_case_all[29:32]),0)
init_sets_list$daily_new_case[29:32] = round(mean(init_sets_list$daily_new_case[29:32]),0)

SEIRfitting(init_sets_list, randomize_startValue=T, run_id = "s1", output_ret = T)
