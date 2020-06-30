## Ro for the five period
## five periods: Jan 1-9 (index 1-9), Jan 10-22 (index 10-22), Jan 23-Feb 1 (index 23-32), Feb 2-16 (index 33-47), Feb 17- (index 48-60)
#' @param pars                    a vetor of parameters: c(b12, b3, b3, b5, r12, delta3, delta4, delta5)
#' @param init_settings           a list of initial values and known parameters
#' @param Dp                      presymptomatic infectious period
#' @param Di                      symptomatic infectious period
#' @param De                      latent period
#' @param Dq                      duration from illness onset to hospitalization
#' @param Dh                      hospitalization period                   
#' @param alpha                   ratio of the transmission rate of unascertained over ascertained case
#' @param N                       population size
#' @param flowN_vec               daily inbound and outbound size during five periods (n)
#' @param init_states             initial c(S, E, P, Is, A, H, R)
#' @param days_to_fit             the days to fit the model
#' @param b                       transmission rate of ascertained cases
#' @param r                       ascertainment rate  
#' @param num_periods             number of periods to simulate
estimate_R <- function(pars, init_settings) {
  tmp_ret=init_settings$var_trans_fun(pars)
  b_vec=tmp_ret[[1]]
  r_vec=tmp_ret[[2]]
  stage_intervals=init_settings$stage_intervals
  n_stage=length(stage_intervals)
  ##
  ##
  Di <- init_settings$Di
  Dp <- init_settings$Dp
  Dq_vec <- init_settings$Dq
  alpha <- init_settings$alpha
  N <- init_settings$N
  flowN_vec <- init_settings$flowN
  #
  R0_est <- rep(NA, n_stage)
  for(i in 1:n_stage) {
    b <- b_vec[i]
    r <- r_vec[i]
    Dq <- Dq_vec[i]
    n <- flowN_vec[i]
    R0_est[i] <- alpha * b / (1 / Dp + n / N) + (1 - r) * alpha * b / (1 / Di + n / N) + r * b / (1 / Di + 1 / Dq)
    rm(i, b, r, Dq, n)
  }
  return(R0_est)
}

