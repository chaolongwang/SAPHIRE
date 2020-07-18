rm(list = ls())
## IMPORTANT: Please set code_root variable properly. 
## code_root should be set to the directory where the repository README file is located. 
## For more information, please read the repository README file
code_root="~/SAPHIRE/"
setwd(paste0(code_root, "scripts_resurgence"))
##
source(paste0(code_root, "R/fun_SEIRresurge.R"))
source(paste0(code_root, "R/init_cond.R"))

##   main models
init_sets_list <- get_init_sets_list(r0 = 0.23)

if (!file.exists("../output/pars_est_run_main_analysis.txt")) {
  stop("Outputs from main analysis cannot be found.\n
      Probably scripts_main/Run_SEIR_main_analysis.R has not yet been run or code_root is not set correctly.")
}

pars_est_dat <- read.table("../output/pars_est_run_main_analysis.txt", header = T)
pars_est_dat <- as.matrix(tail(pars_est_dat, 10000))
##
mydate <- c(paste("Jan", 1:31), paste("Feb", 1:29), paste("Mar", 1:31), paste("Apr", 1:30), paste("May", 1:31))

##
pdf("../output/simulate_resurge_fig3_A.pdf", width = 14, height = 7)
est_mat <- SEIRresurge(pars = apply(pars_est_dat, 2, mean), init_settings = init_sets_list, lift_type = 2, zero_days = 14, outbreak_sign = 3000, return.Mat = T)
est_mat <- na.omit(est_mat)
est_mat <- est_mat[1:152, ]
estIAP_mat <- t(est_mat[, c("I", "A", "P")])
#barht <- apply(estIAP_mat, 2, sum)
barht <- estIAP_mat[1, ]
#c("#F39B7FB2", "#FFDC91FF", "#4DBBD5B2")
par(mfrow = c(1, 1), mar = c(4, 5, 3, 2))
barpos <- barplot(estIAP_mat, ylim = c(0, 60000), col = c("#BC3C29FF", "#FFDC91FF", "#0072B5FF"), xlab = "", ylab = "", border = "NA")
mtext("Date (2020)", side = 1, line  = 3, cex = 1.5)
mtext("No. of active infectious cases", side = 2, line = 3, cex = 1.5)
legend("topleft", legend = c("Presymptomatic (P)", "Unascertained (A)", "Ascertained (I)"), fill = c("#0072B5FF", "#FFDC91FF", "#BC3C29FF"), bty = "n")
axis(1, at = barpos[seq(1, length(barpos), 15)], labels = mydate[seq(1, length(barpos), 15)])
##
Ieq0_index <- head(which(est_mat[, "Onset_expect"] == 0), 1)
arrows(barpos[Ieq0_index], 10000, barpos[Ieq0_index], barht[Ieq0_index], length = 0.1, lwd = 2, col = "red")
text(barpos[Ieq0_index], 10000, labels = "First day of \n I = 0", pos = 3)
arrows(barpos[Ieq0_index + 14], 10000, barpos[Ieq0_index+14], barht[Ieq0_index+14], length = 0.1, lwd = 2, col = "red")
text(barpos[Ieq0_index + 14], 11500, labels = "Lift all controls", pos = 3)
##
AIPgt100_index <- head(which(barht[-(1:Ieq0_index)] >= 100) + Ieq0_index, 1)
arrows(barpos[AIPgt100_index], 10000, barpos[AIPgt100_index], sum(estIAP_mat[, AIPgt100_index]), length = 0.1, lwd = 2, col = "red")
text(barpos[AIPgt100_index-2], 10000, labels = "Resurgence \n defined by I>100", pos = 3)
##
par(fig = c(0.44, 0.965, 0.3, 0.95), new=TRUE)
submydate <- mydate[(Ieq0_index - 10) : min(AIPgt100_index + 5, length(barpos))]
sub_estIAP_dat <- estIAP_mat[, (Ieq0_index - 10) : min(AIPgt100_index + 5, length(barpos))]
for(i in 1:ncol(sub_estIAP_dat)) {
  # i <- 1
  tmp_sum <- sum(sub_estIAP_dat[, i])
  if(tmp_sum >= 500) {
    if(sum(sub_estIAP_dat[1:2, i]) >= 500) {
      sub_estIAP_dat[3, i] <- 0
      sub_estIAP_dat[2, i] <- 500 - sub_estIAP_dat[1, i]
    } else {
      sub_estIAP_dat[3, i] <- 500 - sum(sub_estIAP_dat[1:2, i])
    }
  }
  rm(i, tmp_sum)
}
#
barpos <- barplot(sub_estIAP_dat * 100, col = c("#BC3C29FF", "#FFDC91FF", "#0072B5FF"),
                  xlab = "", ylab = "", border = "NA", yaxt = "n")
axis(2, at = c(0, 10000, 20000, 30000, 40000, 50000), labels = c(0, 10000, 20000, 30000, 40000, 50000) / 100)
abline(h = 10000, col = "red", lwd = 2, lty = 2)
axis(1, at = barpos[seq(1, length(submydate), 15)], labels = submydate[seq(1, length(submydate), 15)])
##
arrows(barpos[11], -10000, barpos[11], 0, length = 0.1, lwd = 2, xpd = T, col = "red")
#text(barpos[11], -10000, labels = "First day of \n I = 0", pos = 1, xpd = T)
arrows(barpos[11+14], -10000, barpos[11+14], 0, length = 0.1, lwd = 2, xpd = T, col = "red")
#text(barpos[11+14], -10000, labels = "Lift all controls", pos = 1, xpd = T)
arrows(barpos[length(barpos) - 5], -10000, barpos[length(barpos) - 5], 0, length = 0.1, lwd = 2, xpd = T, col = "red")
#text(barpos[length(barpos) - 5], -10000, labels = "Resurgence defined \n by I>100", pos = 1, xpd = T)
box(which = "plot", lwd = 1)
dev.off()

