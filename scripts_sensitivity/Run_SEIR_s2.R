rm(list = ls())
## IMPORTANT: set code_root properly before running this script!
code_root="~/jianguoyun/WangLabAdmin/COVID-19/NatSEIR_Rcode/"

setwd(paste0(code_root, "scripts_main"))
library(BayesianTools)
library(vioplot)
library("corrplot")

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

init_sets_list <- generate_init_condi(r0 = 0.23, De=3, Dp=1.1)

SEIRfitting(init_sets_list, randomize_startValue=T, run_id = "s2", output_ret = T)
