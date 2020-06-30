## Find the date when 0 case occurs
## five periods: Jan 1-9 (index 1-9), Jan 10-22 (index 10-22), Jan 23-Feb 1 (index 23-32), Feb 2-16 (index 33-47), Feb 17- (index 48-60)
#' @param pars_estimate           a vetor of parameters: c(b12, b3, b3, b5, r12, delta3, delta4, delta5)
#' @param init_settings           a list of initial values and known parameters
#
Findzero <- function(pars_estimate, init_settings) {
  idx_to_date <- function(d) {
    d <- as.numeric(d)
    mydate <- paste(rep(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep"),c(31, 29, 31, 30, 31, 30, 31, 31, 30)),
                    c(1:31, 1:29, 1:31, 1:30, 1:31, 1:30, 1:31, 1:31, 1:30))
    return(mydate[d])
  }
  ##
  init_settings$days_to_fit <- 1:200
  ##
  est_mat <- apply(pars_estimate, 1, function(x) SEIRsimu(pars = x, init_settings = init_settings, num_periods = 5)[, c("I", "P", "E", "A")])
  est_mat_I <- est_mat[1:200, ]
  est_mat_EPIA <- est_mat[1:200, ] + est_mat[201:400, ] + est_mat[401:600, ] + est_mat[601:800, ]
  ## I = 0
  est <- apply(est_mat_I, 2, function(x) which(x == 0)[1])
  lowdate <- idx_to_date(round(quantile(est, 0.025), 0))
  meandate <- idx_to_date(round(mean(est), 0))
  upperdate <- idx_to_date(round(quantile(est, 0.975), 0))
  intdate1 <- paste(meandate, " (", lowdate, " to ", upperdate,")", sep = "")
  rm(est)
  # E + P + I+ A = 0
  est <- apply(est_mat_EPIA, 2, function(x) tail(which(x != 0), 1) + 1)
  lowdate <- idx_to_date(round(quantile(est, 0.025), 0))
  meandate <- idx_to_date(round(mean(est), 0))
  upperdate <- idx_to_date(round(quantile(est, 0.975), 0))
  intdate2 <- paste(meandate, " (", lowdate, " to ", upperdate,")", sep = "")
  rm(est)
  #
  intdate <- c(intdate1, intdate2)
  names(intdate) <- c("I=0", "E+P+I+A=0")
  return(intdate)
}

## Usage example
#rm(list = ls())
#code_root="C:/Users/xingj/Documents/WangLabAdmin/COVID-19/SEIRcode/"
#setwd(paste0(code_root, "output"))
#source(paste0(code_root, "R/init_cond.R"))
#source(paste0(code_root, "R/fun_SEIRsimu.R"))
##source(paste0(code_root, "R/fun_Findzero.R"))

#init_sets_list <- get_init_sets_list(r0 = 0.23)
#pars_est_dat <- read.table("pars_est_run_main_analysis.txt", header = T)
#pars_est_dat <- as.matrix(tail(pars_est_dat, 10000))

# est_zero_date <- Findzero(pars_estimate = pars_est_dat[1:1000, ], init_settings = init_sets_list)
# 
# est_zero_date





