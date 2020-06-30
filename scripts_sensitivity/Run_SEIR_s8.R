rm(list = ls())
## IMPORTANT: set code_root properly!
code_root="~/jianguoyun/WangLabAdmin/COVID-19/NatSEIR_Rcode/"

setwd(paste0(code_root, "scripts_sensitivity"))
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

init_sets_list <- generate_init_condi(r0 = 0.42)

SEIRfitting(init_sets_list, randomize_startValue=T, run_id = "s8", output_ret = T, skip_MCMC = F)
