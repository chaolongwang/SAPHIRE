rm(list = ls())
##
code_root="C:/Users/xingj/Documents/WangLabAdmin/COVID-19/NatSEIR_Rcode/"
setwd(paste0(code_root, "scripts_resurgence"))
##
source(paste0(code_root, "R/fun_SEIRresurge.R"))
source(paste0(code_root, "R/init_cond.R"))


##   main models
init_sets_list <- get_init_sets_list(r0 = 0.23)
pars_est_dat <- read.table("../output/pars_est_run_main_analysis.txt", header = T)
pars_est_dat <- as.matrix(tail(pars_est_dat, 10000))
##
SEIRresurge(pars = pars_est_dat[1, ], init_settings = init_sets_list, lift_type = 2, zero_days = 7)
SEIRresurge(pars = pars_est_dat[1, ], init_settings = init_sets_list, lift_type = 2, zero_days = 7, outbreak_sign = 150, return.Mat = T)
#
type_1_date_mat <- matrix(NA, 10000, 30)
type_2_date_mat <- matrix(NA, 10000, 30)
colnames(type_1_date_mat) <- colnames(type_2_date_mat) <- 1:30

for(days_num in 1:30) {
  # days_num <- 15
  try({
    simu_outbreak_date1 <- apply(pars_est_dat, 1, function(x) SEIRresurge(pars = x, init_settings = init_sets_list, lift_type = 1, zero_days = days_num))
    simu_outbreak_date2 <- apply(pars_est_dat, 1, function(x) SEIRresurge(pars = x, init_settings = init_sets_list, lift_type = 2, zero_days = days_num))
    type_1_date_mat[, days_num] <- simu_outbreak_date1
    type_2_date_mat[, days_num] <- simu_outbreak_date2
    # sum(!is.na(simu_outbreak_date1)) /  length(simu_outbreak_date1)
    # sum(!is.na(simu_outbreak_date2)) /  length(simu_outbreak_date2)
    cat(days_num, fill = T)
  })
}

write.table(type_1_date_mat, "../output/simu_main_resurgence_date_type_1_lift.txt", quote = F, row.names = F, col.names = T, sep = "\t")
write.table(type_2_date_mat, "../output/simu_main_resurgence_date_type_2_lift.txt", quote = F, row.names = F, col.names = T, sep = "\t")


##  S8 models
init_sets_list <- get_init_sets_list(r0 = 1)
pars_est_dat <- read.table("../output/pars_est_run_s8.txt", header = T)
pars_est_dat <- as.matrix(tail(pars_est_dat, 10000))
##

type_1_date_mat <- matrix(NA, 10000, 30)
type_2_date_mat <- matrix(NA, 10000, 30)
colnames(type_1_date_mat) <- colnames(type_2_date_mat) <- 1:30

for(days_num in 1:30) {
  # days_num <- 14
  try({
    simu_outbreak_date1 <- apply(pars_est_dat, 1, function(x) SEIRresurge(pars = x, init_settings = init_sets_list, lift_type = 1, zero_days = days_num))
    simu_outbreak_date2 <- apply(pars_est_dat, 1, function(x) SEIRresurge(pars = x, init_settings = init_sets_list, lift_type = 2, zero_days = days_num))
    type_1_date_mat[, days_num] <- simu_outbreak_date1
    type_2_date_mat[, days_num] <- simu_outbreak_date2
    # sum(!is.na(simu_outbreak_date1)) /  length(simu_outbreak_date1)
    # sum(!is.na(simu_outbreak_date2)) /  length(simu_outbreak_date2)
    cat(days_num, fill = T)
  })
}

write.table(type_1_date_mat, "../output/simu_s8_resurgence_date_type_1_lift.txt", quote = F, row.names = F, col.names = T, sep = "\t")
write.table(type_2_date_mat, "../output/simu_s8_resurgence_date_type_2_lift.txt", quote = F, row.names = F, col.names = T, sep = "\t")








