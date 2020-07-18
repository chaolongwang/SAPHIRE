rm(list = ls())
## IMPORTANT: Please set code_root variable properly. 
## code_root should be set to the directory where the repository README file is located. 
## For more information, please read the repository README file
code_root="~/SAPHIRE/"

setwd(paste0(code_root, "scripts_sensitivity"))
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

init_sets_list <- generate_init_condi(r0 = 1)

pars_density_allAscertained <- function(pars) {
    d_vec <- rep(NA, 4)
    ##b12, b3, b4, b5
    for(i in c(1:4)) {
      d_vec[i] <- dunif(pars[i], 0, 2, log = T)
    }
    return(sum(d_vec))
}

pars_sampler_allAscertained <- function(n = 1) {
  s_vec <- matrix(NA, n, 4)
  ## b12, b3, b4, b5
  for(i in c(1:4)) {
    s_vec[, i] <- runif(n, 0, 2) 
  }
  return(s_vec)
}
  
transform_var_allAscertained=function(pars) {
  b_vec <- pars[1:4]
  b_vec <- c(b_vec[1], b_vec[1], b_vec[2:4])
  ##
  r12 <- 1
  r3 <- 1
  r4 <- 1
  r5 <- 1
  r_vec <- c(r12,r12,r3,r4,r5)
  
  return(list(b_vec, r_vec))
}

init_sets_list$var_trans_fun=transform_var_allAscertained
init_sets_list$par_lower = c(b12 = 0, b3 = 0, b4 = 0, b5 = 0)
init_sets_list$par_upper = c(b12 = 2, b3 = 2, b4 = 2, b5 = 2)


SEIRfitting(init_sets_list, randomize_startValue=T, run_id = "s9", output_ret = T, skip_MCMC=F,
            pars_density=pars_density_allAscertained, pars_sampler=pars_sampler_allAscertained,
            pars_name=c("b12", "b3", "b4", "b5"))

