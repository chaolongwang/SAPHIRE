library(cairoDevice)

default_pars_density <- function(pars) {
  d_vec <- rep(NA, 8)
  ##b12, b3, b4, b5
  for(i in c(1:4)) {
    d_vec[i] <- dunif(pars[i], 0, 2, log = T)
  }
  ## r12
  r12 <- pars[5]
  d_vec[5] <- dbeta(r12, beta_shape1, beta_shape2, log = T)
  ## r3
  delta_3 <- pars[6]
  d_vec[6] <- dnorm(delta_3, delta_mean, delta_sd, log = T)
  ## r4
  delta_4 <- pars[7]
  d_vec[7] <- dnorm(delta_4, delta_mean, delta_sd, log = T)
  ## r5
  delta_5 <- pars[8]
  d_vec[8] <- dnorm(delta_5, delta_mean, delta_sd, log = T)
  ##
  return(sum(d_vec))
}

default_pars_sampler <- function(n = 1) {
  s_vec <- matrix(NA, n, 8)
  ## b12, b3, b4, b5
  for(i in c(1:4)) {
    s_vec[, i] <- runif(n, 0, 2) 
  }
  ## r12 
  r12 <- rbeta(n, beta_shape1, beta_shape2)
  s_vec[, 5] <- r12
  ## r3
  s_vec[, 6] <- rnorm(n, delta_mean, delta_sd)
  ## r4
  s_vec[, 7] <- rnorm(n, delta_mean, delta_sd)
  ## r5
  s_vec[, 8] <- rnorm(n, delta_mean, delta_sd)
  return(s_vec)
}

## wrapper for the analysis run
## Create BayesianSetup and settings, lower/upper for parameters: b12, b3, b3, b5, r12, delta3, delta4, delta5
#' @param init_set_list   initial settings produced by init_cond.R
#' @param randomize_startValue    this function will randomly generate an initial condtion if this argument is set to T; If you want to specify your own initial condition, set to F
#' @param startValue  If randomize_startValue is set to T, you can your own initial condition using this argument; If randomize_startValue is set to T, this argument will be ignored
#' @param output_ret  Whether to output parameter estimates output by MCMC
#' @param run_id this ID is meant to distinguish different runs (different parameters, random seeds, etc.). Run_ID will be included in the file names of all outputs
#' @param skip_MCMC This is meant for redrawing all results without rerunning MCMC
#' @param panel_B_R_ylim the y limit in panel B of the main result plot
SEIRfitting=function(init_sets_list, 
                     randomize_startValue=F, 
                     startValue=NA, 
                     output_ret=T, 
                     run_id=0, 
                     skip_MCMC=F, 
                     panel_B_R_ylim=4,
                     plot_combined_fig=T,
                     pars_density=default_pars_density,
                     pars_sampler=default_pars_sampler,
                     pars_name=c("b12", "b3", "b4", "b5", "r12", "delta3", "delta4", "delta5"),
                     calc_clearance=T,
                     n_burn_in=4000,
                     n_iterations=180000) {
  if (randomize_startValue & !is.na(startValue)) {
    print("startValue will be ignored since randomize_startValue is set to TRUE!")
  } else if (!randomize_startValue & is.na(startValue)) {
    print("Please specify a startValue since you have set randomize_startValue to FALSE! Exiting!")
    q(save="no")
  }
  
  onset_obs <- init_sets_list$daily_new_case
  init_states <- init_sets_list$init_states
  n_pars = length(pars_name)
  n_stage = length(init_sets_list$stage_intervals)
  
  ## take a try: pars = c(b12, b3, b3, b5, r12, delta3, delta4, delta5)
  # SEIRpred(pars = c(1.4, 0.4, 0.1, 0.1, 0.5, -1, 0, 0), init_settings = init_sets_list)[, "Onset_expect"]
  ################################################################################################################
  ## likelihood function
  loglh_func <- function(pars){
    ypred <- SEIRpred(pars, init_settings = init_sets_list)
    ypred <- ypred[, "Onset_expect"]
    
    # meant to suppress warnings when ypred is negative
    suppressWarnings(p <- dpois(onset_obs, ypred))
    
    if(any(p == 0) || any(is.nan(p))){
      logL <- -Inf
    }else{
      logL <- sum(log10(p))
    }
    return(logL)
  }
  ## take a try
  # loglh_func(pars = c(1.4, 0.4, 0.1, 0.1, 0.5, -1, 0, 0))
  
  ## Create BayesianSetup and settings, lower/upper for parameters: b12, b3, b3, b5, r12, delta3, delta4, delta5
  
  
  ## take a try 
  
  ## take a try 
  # pars_sampler(n = 1)
  pars_prior <- createPrior(density = pars_density, sampler = pars_sampler, 
                            lower = init_sets_list$par_lower, upper = init_sets_list$par_upper)
  
  if (!skip_MCMC) {
    bayesSEIR <- createBayesianSetup(loglh_func, prior = pars_prior)
  
    if (randomize_startValue) {  
      startValue=pars_sampler()
      while (is.infinite(loglh_func(startValue))) {
        startValue=pars_sampler()
      }
    }
    
    ## DRAM: Adaptive MCMC, prior optimization, delayed rejection
    # startValue = c(b12=1.2, b3=0.4, b4=0.2, b5=0.1, r12=0.5, delta3=-1, delta4=0, delta5=0)
    # startValue = c(b12 = 1.359, b3 = 0.537, b4 = 0.203, b5 = 0.196, r12 = 0.305, delta3 = -0.964, delta4 = -0.593, delta5 = -0.309)
    mh_settings = list(startValue = startValue,
                       adapt = T, DRlevels = 2, iterations = n_iterations, thin = 10)
    mh_out <- runMCMC(bayesianSetup = bayesSEIR, sampler = "Metropolis", settings = mh_settings)
    #plot(mh_out)
    mcmc_pars_estimate <- getSample(mh_out, start = n_burn_in+2, thin = 1)  ## set start = 2002 as the burn in period
    mcmc_pars_estimate <- round(mcmc_pars_estimate, 3)
    
    colnames(mcmc_pars_estimate) <- pars_name 
    
    if (output_ret) {
      write.table(mcmc_pars_estimate, paste0("../output/pars_est_run_",run_id,".txt"), quote = F, row.names = F, sep = "\t")
    }
  } else {
    mcmc_pars_estimate = read.table(paste0("../output/pars_est_run_",run_id,".txt"), header = T)
    pars_name = names(mcmc_pars_estimate)
  }
  
  summary_string=paste0(paste(pars_name, collapse = ","), "\n")
  
  par_str=list()
  for (i_par in 1:n_pars) {
    par_str[[i_par]]=paste0(round(mean(mcmc_pars_estimate[,i_par]),2), " (",
                            round(quantile(mcmc_pars_estimate[,i_par],0.025),2)," - " , 
                            round(quantile(mcmc_pars_estimate[,i_par],0.975),2), ")")
  }
  
  summary_string = paste0(summary_string, paste(par_str,collapse = ", "),"\n\n")
  
  estRt_mat <- apply(mcmc_pars_estimate, 1, function(x) estimate_R(pars = x, init_settings = init_sets_list))

  summary_string = paste0(summary_string, paste0("stage",1:n_stage,collapse=","), "\n")
  
  r_str=list()
  
  if (n_stage>1) {
    for (i_stage in 1:n_stage) {
      r_str[[i_stage]]=paste0(round(mean(estRt_mat[i_stage,]),2), " (",
                              round(quantile(estRt_mat[i_stage,],0.025),2)," - " , 
                              round(quantile(estRt_mat[i_stage,],0.975),2), ")")
    }
  } else {
    r_str[[1]]=paste0(round(mean(estRt_mat),2), " (",
                            round(quantile(estRt_mat,0.025),2)," - " , 
                            round(quantile(estRt_mat,0.975),2), ")")
  }
  
  summary_string = paste0(summary_string, paste(r_str,collapse = ", "),"\n\n")
  
  if (calc_clearance) {
    clearance_date = Findzero(mcmc_pars_estimate, init_sets_list)
    
    summary_string = paste0(summary_string, paste(names(clearance_date), collapse = ", "))
    
    summary_string = paste0(summary_string, "\n", paste(clearance_date, collapse = ", "), "\n")
  }
  
  write_file(summary_string, paste0("../output/summary_run_",run_id,".txt"))
  
  cairo_pdf(paste0("../output/par_cor_run_",run_id,".pdf"),width=10,height=10)
  correlationPlot_modified(mcmc_pars_estimate, scaleCorText = F)
  dev.off()
  
  png(paste0("../output/par_hist_run_",run_id,".png"))
  par(mfrow = c(2, 4))
  for(i in 1:n_pars) {
    hist(mcmc_pars_estimate[, i], xlab = pars_name[i], main = "", col = "red")
    rm(i)
  }
  dev.off()
  
  png(paste0("../output/par_traj_run_",run_id,".png"), width=1000, height=500)
  par(mfrow = c(2, 4))
  for(i in 1:n_pars) {
    plot(1:nrow(mcmc_pars_estimate), mcmc_pars_estimate[, i], ylab = pars_name[i], xlab = "iter", main = "", type = "l")
    rm(i)
  }
  dev.off()
  
  if (plot_combined_fig) {
    SEIRplot(pars_estimate = mcmc_pars_estimate, file_name = run_id, init_settings = init_sets_list, panel_B_R_ylim = panel_B_R_ylim)
  }
  
  par(mfrow = c(1, 1))
  # corrplot(cor(mcmc_pars_estimate))
  # pairs(mcmc_pars_estimate)
  
}
