## Deterministic SEIR model
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
#' @param stage_intervals         default corresponding staget intervals (Jan1, Jan9), (Jan10, Jan22), (Jan23, Feb1), (Feb2, Feb16), (Feb17, last day (varying))
#' 
#################################################################################################
SEIRpred <- function(pars, 
                     init_settings) {
  tmp_ret=init_settings$var_trans_fun(pars)
  b_vec=tmp_ret[[1]]
  r_vec=tmp_ret[[2]]
  stage_intervals=init_settings$stage_intervals
  n_stage=length(stage_intervals)
  
  ##
  Di <- init_settings$Di
  Dp <- init_settings$Dp
  De <- init_settings$De
  Dq_vec <- init_settings$Dq
  alpha <- init_settings$alpha
  Dh <- init_settings$Dh
  N <- init_settings$N
  flowN_vec <- init_settings$flowN
  init_states <- init_settings$init_states
  days_to_fit <- init_settings$days_to_fit
  ## ODE function based on deterministic SEIR model
  update_func <- function(stage_pars, states_old) {
    ## stage pars
    b <- stage_pars[1]
    r <- stage_pars[2]
    Dq <- stage_pars[3]
    n <- stage_pars[4]
    ## old states number: c(S, E, P, I, A, H, R)
    S <- states_old[1]
    E <- states_old[2]
    P <- states_old[3]
    I <- states_old[4]
    A <- states_old[5]
    H <- states_old[6]
    R <- states_old[7]
    ## new values
    S_new <- S - b * S * (alpha * P + I + alpha * A) / N + n - n * S / N
    E_new <- E + b * S * (alpha * P + I + alpha * A) / N - E / De - n * E / N
    P_new <- P +  E / De  - P / Dp - n * P / N
    I_new <- I + r * P / Dp - I / Di - I / Dq
    A_new <- A + (1 - r) * P / Dp - A / Di - n * A / N
    H_new <- H + I / Dq - H / Dh
    R_new <- R + H / Dh + (A + I) / Di - n * R / N
    Onset_expect <- r * P / Dp
    ##
    return(c(S_new, E_new, P_new, I_new, A_new, H_new, R_new, Onset_expect))
  }
  ## matrix for results
  states_mat <- matrix(0, length(days_to_fit), length(init_states) + 2)
  states_mat[, 1] <- days_to_fit
  colnames(states_mat) <- c("time", "S", "E", "P", "I", "A", "H", "R", "Onset_expect")
  
  myold_states <- init_states
  
  for (i_stage in 1:n_stage) {
    stage_pars_setings <- c(b = b_vec[i_stage], r = r_vec[i_stage], Dq = Dq_vec[i_stage], n = flowN_vec[i_stage])
    for (d in stage_intervals[[i_stage]][["start"]]:stage_intervals[[i_stage]][["end"]]) {
      states_mat[d, -1] <- update_func(stage_pars = stage_pars_setings, states_old = myold_states)
      myold_states <- states_mat[d, -1]
    }
  }
  return(states_mat)
}

