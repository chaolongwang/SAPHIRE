## SEIR model plots for the six panels
## five periods: Jan 1-9 (index 1-9), Jan 10-22 (index 10-22), Jan 23-Feb 1 (index 23-32), Feb 2-16 (index 33-47), Feb 17- (index 48-60)
#' @param pars_estimate           a vetor of parameters: c(b12, b3, b3, b5, r12, delta3, delta4, delta5)
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
#################################################################################################
SEIRplot <- function(pars_estimate, file_name, init_settings, panel_B_R_ylim=4) {
  init_settings$days_to_fit <- 1:68
  library(vioplot)
  ##
  onset_obs_all <- init_settings$daily_new_case_all
  ptime <- 1:length(onset_obs_all)
  mydate <- c(paste("Jan", 1:31), paste("Feb", 1:29), paste("Mar", 1:8))
  #
  pdf(paste0("../output/Figure_", file_name, ".pdf"), width = 9, height = 10)
  par(mar = c(4, 5, 2.5, 1))
  layout(matrix(c(1:6), byrow = T, nrow = 3))
  
  ##########################################   Panel A  ##########################################################
  estN_mat <- apply(pars_estimate, 1, function(x) SEIRsimu(pars = x, init_settings = init_settings, num_periods = 5)[, "Onset_expect"])
  estN_mean <- round(apply(estN_mat, 1, mean), 0)
  estN_up <- round(apply(estN_mat, 1, function(x) quantile(x, 0.975)), 0)
  estN_low <- round(apply(estN_mat, 1, function(x) quantile(x, 0.025)), 0)
  # start A
  plot(ptime, estN_mean, ylim = c(0, max(estN_up, onset_obs_all) * 1.05), xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
  mtext("Onset date (2020)", side = 1, line  = 3, cex = 1.01)
  mtext("No. of ascertained cases", side = 2, line = 3, cex = 1.01)
  axis(1, at = seq(1, 68, 11), labels = mydate[seq(1, 68, 11)])
  #
  abline(v = c(10, 23, 33, 48, 61), lty = 3, lwd = 2, col = "darkgrey")
  text(c(10, 23, 33, 48, 61), par()$usr[4], labels = mydate[c(10, 23, 33, 48, 61)], col = "darkgrey", pos = 3, xpd = T)
  #
  polygon(c(ptime[1:61], rev(ptime[1:61])), c(estN_up[1:61], rev(estN_low[1:61])), col = "#F39B7FB2", border = NA)
  polygon(c(ptime[-c(1:60)], rev(ptime[-c(1:60)])), c(estN_up[-c(1:60)], rev(estN_low[-c(1:60)])), col = "#4DBBD5B2", border = NA)
  #
  points(ptime[1:60], estN_mean[1:60], col = "#BC3C29FF", pch = 16, cex = 0.8)
  points(ptime[-c(1:60)], estN_mean[-c(1:60)], col = "#0072B5FF", pch = 17, cex = 0.8)
  points(ptime, onset_obs_all, col = "black", pch = 4, cex = 0.8)
  #
  legend("topleft", legend = c("Observed", "Fitted",  "Predicted"), col = c("black",  "#BC3C29FF", "#0072B5FF"), pch = c(4, 16, 17), bty = "n")
  text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.12, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "A", xpd = T, cex = 2)
  
  ##########################################   Panel B  ##########################################################
  estRt_mat <- apply(pars_estimate, 1, function(x) estimate_R(pars = x, init_settings = init_settings))
  estRt_mat <- t(estRt_mat)
  ##
  rt_mean <- sprintf("%.2f", round(apply(estRt_mat, 2, function(x) mean(x)), 2))
  rt_low <- sprintf("%.2f", round(apply(estRt_mat, 2, function(x) quantile(x, 0.025)), 2))
  rt_up <- sprintf("%.2f", round(apply(estRt_mat, 2, function(x) quantile(x, 0.975)), 2))
  #
  vioplot(estRt_mat[, 1], estRt_mat[, 2], estRt_mat[, 3], estRt_mat[, 4], estRt_mat[, 5], names = NA, ylim = c(0, panel_B_R_ylim), col = c("#BC3C29FF","#0072B5FF", "#E18727FF", "#7876B1FF", "#FFDC91FF"), ylab = "", xlab = "")
  mtext("Outbreak period (2020)", side = 1, line  = 3, cex = 1.01)
  mtext(expression("R"["0"]), side = 2, line = 3, cex = 1.01)
  axis(1, at = c(1, 1.9, 3, 4, 5), tick = F, labels = c("Jan 1-9", "Jan 11-22", "Jan 23-Feb 1", "Feb 2-16", "Feb 17-Mar 8"))
  abline(h = 1, lwd = 2, lty = 3, col = "red")
  #
  text(1, min(estRt_mat[, 1]) - 0.2, labels = rt_mean[1])
  text(1, min(estRt_mat[, 1]) - 0.45, labels = paste("(", rt_low[1], "-", rt_up[1], ")", sep = ""))
  text(2, min(estRt_mat[, 2]) - 0.2, labels = rt_mean[2])
  text(2, min(estRt_mat[, 2]) - 0.45, labels = paste("(", rt_low[2], "-", rt_up[2], ")", sep = ""))
  text(3, max(estRt_mat[, 3]) + 0.4, labels = rt_mean[3])
  text(3, max(estRt_mat[, 3]) + 0.15, labels = paste("(", rt_low[3], "-", rt_up[3], ")", sep = ""))
  text(4, max(estRt_mat[, 4]) + 0.4, labels = rt_mean[4])
  text(4, max(estRt_mat[, 4]) + 0.15, labels = paste("(", rt_low[4], "-", rt_up[4], ")", sep = ""))
  text(5, max(estRt_mat[, 5]) + 0.4, labels = rt_mean[5])
  text(5, max(estRt_mat[, 5]) + 0.15, labels = paste("(", rt_low[5], "-", rt_up[5], ")", sep = ""))
  text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.12, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "B", xpd = T, cex = 2)
  
  ##########################################   Panel C  ##########################################################
  estN_mat <- apply(pars_estimate, 1, function(x) SEIRsimu(pars = x, init_settings = init_settings, num_periods = 4)[, "Onset_expect"])
  estN_mean <- round(apply(estN_mat, 1, mean), 0)
  estN_up <- round(apply(estN_mat, 1, function(x) quantile(x, 0.975)), 0)
  estN_low <- round(apply(estN_mat, 1, function(x) quantile(x, 0.025)), 0)
  # start C
  plot(ptime, estN_mean, ylim = c(0, max(estN_up, onset_obs_all) * 1.05), xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
  mtext("Onset date (2020)", side = 1, line  = 3, cex = 1.01)
  mtext("No. of ascertained cases", side = 2, line = 3, cex = 1.01)
  axis(1, at = seq(1, 68, 11), labels = mydate[seq(1, 68, 11)])
  #
  abline(v = c(10, 23, 33, 48, 61), lty = 3, lwd = 2, col = "darkgrey")
  text(c(10, 23, 33, 48, 61), par()$usr[4], labels = mydate[c(10, 23, 33, 48, 61)], col = "darkgrey", pos = 3, xpd = T)
  #
  polygon(c(ptime[1:48], rev(ptime[1:48])), c(estN_up[1:48], rev(estN_low[1:48])), col = "#F39B7FB2", border = NA)
  polygon(c(ptime[-c(1:47)], rev(ptime[-c(1:47)])), c(estN_up[-c(1:47)], rev(estN_low[-c(1:47)])), col = "#4DBBD5B2", border = NA)
  #
  points(ptime[1:47], estN_mean[1:47], col = "#BC3C29FF", pch = 16, cex = 0.8)
  points(ptime[-c(1:47)], estN_mean[-c(1:47)], col = "#0072B5FF", pch = 17, cex = 0.8)
  points(ptime, onset_obs_all, col = "black", pch = 4, cex = 0.8)
  #
  legend("topleft", legend = c("Observed", "Fitted",  "Predicted"), col = c("black",  "#BC3C29FF", "#0072B5FF"), pch = c(4, 16, 17), bty = "n")
  text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.12, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "C", xpd = T, cex = 2)
 
  ##########################################   Panel D  ##########################################################
  estN_mat <- apply(pars_estimate, 1, function(x) SEIRsimu(pars = x, init_settings = init_settings, num_periods = 3)[, "Onset_expect"])
  estN_mean <- round(apply(estN_mat, 1, mean), 0)
  estN_up <- round(apply(estN_mat, 1, function(x) quantile(x, 0.975)), 0)
  estN_low <- round(apply(estN_mat, 1, function(x) quantile(x, 0.025)), 0)
  # start D
  plot(ptime, estN_mean, ylim = c(0, max(estN_up, onset_obs_all) * 1.01), xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
  mtext("Onset date (2020)", side = 1, line  = 3, cex = 1.01)
  mtext("No. of ascertained cases", side = 2, line = 3, cex = 1.01)
  axis(1, at = seq(1, 68, 11), labels = mydate[seq(1, 68, 11)])
  #
  abline(v = c(10, 23, 33, 48, 61), lty = 3, lwd = 2, col = "darkgrey")
  text(c(10, 23, 33, 48, 61), par()$usr[4], labels = mydate[c(10, 23, 33, 48, 61)], col = "darkgrey", pos = 3, xpd = T)
  #
  polygon(c(ptime[1:33], rev(ptime[1:33])), c(estN_up[1:33], rev(estN_low[1:33])), col = "#F39B7FB2", border = NA)
  polygon(c(ptime[-c(1:32)], rev(ptime[-c(1:32)])), c(estN_up[-c(1:32)], rev(estN_low[-c(1:32)])), col = "#4DBBD5B2", border = NA)
  #
  points(ptime[1:32], estN_mean[1:32], col = "#BC3C29FF", pch = 16, cex = 0.8)
  points(ptime[-c(1:32)], estN_mean[-c(1:32)], col = "#0072B5FF", pch = 17, cex = 0.8)
  points(ptime, onset_obs_all, col = "black", pch = 4, cex = 0.8)
  #
  legend("topleft", legend = c("Observed", "Fitted",  "Predicted"), col = c("black",  "#BC3C29FF", "#0072B5FF"), pch = c(4, 16, 17), bty = "n")
  text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.12, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "D", xpd = T, cex = 2)
  
  ##########################################   Panel E  ##########################################################
  estN_mat <- apply(pars_estimate, 1, function(x) SEIRsimu(pars = x, init_settings = init_settings, num_periods = 2)[, "Onset_expect"])
  estN_mean <- round(apply(estN_mat, 1, mean), 0)
  estN_up <- round(apply(estN_mat, 1, function(x) quantile(x, 0.975)), 0)
  estN_low <- round(apply(estN_mat, 1, function(x) quantile(x, 0.025)), 0)
  # start E
  plot(ptime, estN_mean, ylim = c(0, max(estN_up, onset_obs_all) * 1.01), xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
  mtext("Onset date (2020)", side = 1, line  = 3, cex = 1.01)
  mtext("No. of ascertained cases", side = 2, line = 3, cex = 1.01)
  axis(1, at = seq(1, 68, 11), labels = mydate[seq(1, 68, 11)])
  #
  abline(v = c(10, 23, 33, 48, 61), lty = 3, lwd = 2, col = "darkgrey")
  text(c(10, 23, 33, 48, 61), par()$usr[4], labels = mydate[c(10, 23, 33, 48, 61)], col = "darkgrey", pos = 3, xpd = T)
  #
  polygon(c(ptime[1:23], rev(ptime[1:23])), c(estN_up[1:23], rev(estN_low[1:23])), col = "#F39B7FB2", border = NA)
  polygon(c(ptime[-c(1:22)], rev(ptime[-c(1:22)])), c(estN_up[-c(1:22)], rev(estN_low[-c(1:22)])), col = "#4DBBD5B2", border = NA)
  #
  points(ptime[1:22], estN_mean[1:22], col = "#BC3C29FF", pch = 16, cex = 0.8)
  points(ptime[-c(1:22)], estN_mean[-c(1:22)], col = "#0072B5FF", pch = 17, cex = 0.8)
  points(ptime, onset_obs_all, col = "black", pch = 4, cex = 0.8)
  #
  legend("topleft", legend = c("Observed", "Fitted",  "Predicted"), col = c("black",  "#BC3C29FF", "#0072B5FF"), pch = c(4, 16, 17), bty = "n")
  text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.12, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "E", xpd = T, cex = 2)
  
  ##########################################   Panel F  ##########################################################
  estAIP_mat <- apply(pars_estimate, 1, function(x) SEIRsimu(pars = x, init_settings = init_settings, num_periods = 5)[, c("I", "A", "P")])
  estI_mat <- estAIP_mat[ptime, ]
  estA_mat <- estAIP_mat[ptime + length(ptime), ]
  estP_mat <- estAIP_mat[ptime + length(ptime) * 2, ]
  estI_mean <- apply(estI_mat, 1, mean)
  estA_mean <- apply(estA_mat, 1, mean)
  estP_mean <- apply(estP_mat, 1, mean)
  estAIP_dat <- rbind(estI_mean, estA_mean, estP_mean)
  barpos <- barplot(estAIP_dat, col = c("#BC3C29FF", "#FFDC91FF", "#0072B5FF"), xlab = "", ylab = "", border = "NA")
  mtext("Date (2020)", side = 1, line  = 3, cex = 1.01)
  mtext("No. of active infectious cases", side = 2, line = 3, cex = 1.01)
  axis(1, at = barpos[seq(1, 68, 11)], labels = mydate[seq(1, 68, 11)])
  legend("topleft", legend = c("Presymptomatic (P)", "Unascertained (A)", "Ascertained (I)"), fill = c("#0072B5FF", "#FFDC91FF", "#BC3C29FF"), bty = "n")
  text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.12, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "F", xpd = T, cex = 2)
  ##figure_F finished
  dev.off()
}