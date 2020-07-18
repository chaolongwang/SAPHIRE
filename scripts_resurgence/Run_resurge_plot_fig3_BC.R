rm(list = ls())
## IMPORTANT: Please set code_root variable properly. 
## code_root should be set to the directory where the repository README file is located. 
## For more information, please read the repository README file
code_root="~/SAPHIRE/"
setwd(paste0(code_root, "scripts_resurgence"))


##
s0_dat1 <- as.matrix(read.table("../output/simu_main_resurgence_date_type_1_lift.txt", header = T))
s0_dat2 <- as.matrix(read.table("../output/simu_main_resurgence_date_type_2_lift.txt", header = T))
s8_dat1 <- as.matrix(read.table("../output/simu_s8_resurgence_date_type_1_lift.txt", header = T))
s8_dat2 <- as.matrix(read.table("../output/simu_s8_resurgence_date_type_2_lift.txt", header = T))

simu_outbreak_date1 <- s0_dat1[, 14]
sum(!is.na(simu_outbreak_date1)) /  length(simu_outbreak_date1)
##
outbreak_prop_mat <- matrix(NA, 30, 4)
outbreak_prop_mat[, 2] <- apply(s0_dat1, 2, function(x) length(na.omit(x)) / length(x))
outbreak_prop_mat[, 1] <- apply(s0_dat2, 2, function(x) length(na.omit(x)) / length(x))
outbreak_prop_mat[, 4] <- apply(s8_dat1, 2, function(x) length(na.omit(x)) / length(x))
outbreak_prop_mat[, 3] <- apply(s8_dat2, 2, function(x) length(na.omit(x)) / length(x))
#outbreak_prop_mat[14, ]
##
cal_date <- function(x) {
  # x <- s0_dat1
  x20 <- apply(x[, 1:20], 2, function(x) mean(na.omit(x)))
  x21p <- mean(na.omit(as.vector(x[, -c(1:20)])))
  return(c(x20, x21p))
}
cal_low_date <- function(x) {
  # x <- s0_dat1
  x20 <- apply(x[, 1:20], 2, function(x) quantile(na.omit(x), 0.025))
  x21p <- quantile(na.omit(as.vector(x[, -c(1:20)])), 0.025)
  return(c(x20, x21p))
}

cal_up_date <- function(x) {
  # x <- s0_dat1
  x20 <- apply(x[, 1:20], 2, function(x) quantile(na.omit(x), 0.975))
  x21p <- quantile(na.omit(as.vector(x[, -c(1:20)])), 0.975)
  return(c(x20, x21p))
}

s0_date1 <- cal_date(s0_dat1)
s0_date2 <- cal_date(s0_dat2)
s8_date1 <- cal_date(s8_dat1)
s8_date2 <- cal_date(s8_dat2)
date_mat <- cbind(s0_date2, s0_date1, s8_date2, s8_date1)
# date_mat[14, ]
cal_low_date(s0_dat2)[14]
cal_up_date(s0_dat2)[14]
cal_low_date(s0_dat1)[14]
cal_up_date(s0_dat1)[14]
#
cal_low_date(s8_dat1)[14]
cal_up_date(s8_dat1)[14]



pdf("../output/simulate_resurge_fig3_BC.pdf", width = 14, height = 7)
par(mfrow = c(1, 2), mar = c(4, 5, 3, 2))
plot(1:nrow(outbreak_prop_mat), outbreak_prop_mat[, 1], type = "l", ylim = c(0, 1), col = "#0072B5FF", lwd = 3, xlab = "", ylab = "", xaxt = "n")
axis(1, at = round(seq(1, 30, length.out = 6), 0), labels = round(seq(1, 30, length.out = 6), 0))
mtext("t (days)", 1, line = 2.5, cex = 1.5)
mtext("Probability of resurgence", 2, line = 2.5, cex = 1.5)
lines(1:nrow(outbreak_prop_mat), outbreak_prop_mat[, 2], col = "#BC3C29FF", lwd = 3)
lines(1:nrow(outbreak_prop_mat), outbreak_prop_mat[, 3], col = "#0072B5FF", lwd = 3, lty = 3)
lines(1:nrow(outbreak_prop_mat), outbreak_prop_mat[, 4], col = "#BC3C29FF", lwd = 3, lty = 3)
points(1:nrow(outbreak_prop_mat), outbreak_prop_mat[, 1], pch = 16, col = "#0072B5FF", cex = 1.5)
points(1:nrow(outbreak_prop_mat), outbreak_prop_mat[, 2], pch = 16, col = "#BC3C29FF", cex = 1.5)
points(1:nrow(outbreak_prop_mat), outbreak_prop_mat[, 3], pch = 18, col = "#0072B5FF", cex = 1.5)
points(1:nrow(outbreak_prop_mat), outbreak_prop_mat[, 4], pch = 18, col = "#BC3C29FF", cex = 1.5)
text(30, 0.9, pos = 2, labels = "Legend in panel C")
text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.1, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.05, labels = "B", xpd = T, cex = 2)
#
plot(1:21, s0_date2, type = "l", ylim = c(23, 45), col = "#0072B5FF", lwd = 3, xlab = "", ylab = "", xaxt = "n")
axis(1, at = round(seq(1, 21, length.out = 6), 0), labels = c("1", "5", "9", "13", "17", "21-30"))
mtext("t (days)", 1, line = 2.5, cex = 1.5)
mtext("Conditional expectation of time to resurgence (days)", 2, line = 2.5, cex = 1.5)
lines(1:21, s0_date1, col = "#BC3C29FF", lwd = 3)
lines(1:21, s8_date2, col = "#0072B5FF", lwd = 3, lty = 3)
lines(1:21, s8_date1, col = "#BC3C29FF", lwd = 3, lty = 3)
points(1:21, s0_date2, pch = 16, col = "#0072B5FF", cex = 1.5)
points(1:21, s0_date1, pch = 16, col = "#BC3C29FF", cex = 1.5)
points(1:21, s8_date2, pch = 18, col = "#0072B5FF", cex = 1.5)
points(1:21, s8_date1, pch = 18, col = "#BC3C29FF", cex = 1.5)
legend("bottomright", title = "When to lift all controls" ,legend = c("t days after first I=0 (M)", "t days of I=0 consecutively (M)", "t days after first I=0 (S8)", "t days of I=0 consecutively (S8)"), 
       col = c("#0072B5FF", "#BC3C29FF"), lty = c(1, 1, 3, 3), lwd = 2, bty = "n", pch = c(16, 16, 18, 18))
text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.1, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.05, labels = "C", xpd = T, cex = 2)
dev.off()
























