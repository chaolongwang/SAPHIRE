generate_init_condi <- function(r0,
                                Di = 2.9,
                                Dp = 2.3,
                                De = 2.9,
                                Dq = c(21, 15, 10, 6, 3),
                                alpha = 0.55,
                                Dh = 30,
                                N = 10000000,
                                flowN = c(500000, 800000, 0, 0, 0)
                                ) {
  
  stopifnot(r0>=0 & r0<=1 & Di>=0 & Dp>=0 & De>=0 & all(Dp>=0) & alpha>=0 & alpha<=1 & Dh>=0 & N>=0 & all(flowN>=0))
  
  ## N            : population size
  ## H0           : initial number of hospitalized cases based on the reports
  ## R0           : initial number of removed individuals
  ## De           : latent period
  ## r0           : initial ascertainment rate
  ## realData     : real data from the CDC
  R0 <- 0
  H0 <- 27
  
  realData_all <- read.csv("../data/Covid19CasesWH.csv", row.names = 1)  
  realData <- realData_all[-c(1:24), ] # the 25th row correspond to 1 Jan
  jan1_idx = 25
  
  #cases from Jan 1 to Feb 29
  daily_new_case <- realData[1:60, 1]
  daily_new_case_all <- realData[, 1]
  ##
  E0 <- sum(realData_all[(jan1_idx+round(Dp)):(jan1_idx+round(Dp)+round(De)-1),1]) / r0 ## Jan 3-5 for De=2.9 and Dp=2.3
  # E0 <- (40 + 23 + 47) / r0                                              
  P0 <- sum(realData_all[jan1_idx:(jan1_idx+round(Dp)-1),1]) / r0                     ## Jan 1-2 for Dp=2.3
  # P0 <- (41 + 34) / r0
  I0 <- sum(realData_all[(jan1_idx-round(Di)):(jan1_idx-1),1])                             ## Dec 29-31 for Di=2.9
  # I0 <- 11 + 13 + 10                                     
  A0 <- I0 * (1 - r0) / r0
  S0 <- N - E0 - P0 - I0 - A0 - H0 - R0
  init_states <- round(c(S = S0, E = E0, P = P0, I = I0, A = A0, H = H0, R = R0), 0)
  
  ## helper function
  # transform variables to a form that SEIRpred can use
  # so that SEIRpred can be re-used as much as possible
  transform_var_main_5stage=function(pars) {
    b_vec <- pars[1:4]
    b_vec <- c(b_vec[1], b_vec[1], b_vec[2:4])
    ##
    r12 <- pars[5]
    r3 <- 1 / (1 + (1 - r12) / (r12 * exp(pars[6])))
    r4 <- 1 / (1 + (1 - r3) / (r3 * exp(pars[7])))
    r5 <- 1 / (1 + (1 - r4) / (r4 * exp(pars[8])))
    r_vec <- c(r12,r12,r3,r4,r5)
    
    return(list(b_vec, r_vec))
  }
  
  return(list(Di=Di,
              Dp=Dp,
              De=De,
              Dq=Dq,
              alpha=alpha,
              Dh=Dh,
              N=N,
              flowN=flowN,
              daily_new_case = daily_new_case, 
              daily_new_case_all = daily_new_case_all, 
              init_states = init_states,
              days_to_fit=1:60,
              stage_intervals=list(
                c(start=1, end=9),
                c(start=10, end=22),
                c(start=23, end=32),
                c(start=33, end=47),
                c(start=48, end=60)
              ),
              var_trans_fun=transform_var_main_5stage,
         par_lower = c(b12 = 0, b3 = 0, b4 = 0, b5 = 0, r12 = 0, delta3 = -10, delta4 = -10, delta5 = -10),
         par_upper = c(b12 = 2, b3 = 2, b4 = 2, b5 = 2, r12 = 1, delta3 = 10, delta4 = 10, delta5 = 10)))
  # TODO: please confirm the following:
  # boundaries for delta3-5 will not be used, they are here merely to meet the formality imposed by runMCMC
}

# get_init_sets_list is an alias of generate_init_condi in order not to break exsiting code
get_init_sets_list = generate_init_condi

delta_mean <- 0
delta_sd <- 1
beta_shape1 <- 7.3
beta_shape2 <- 24.6
