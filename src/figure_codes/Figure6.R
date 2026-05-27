library(stringr)
library(survival)

# Codes for Figure 6

# Figure 6a: boxplot comparing carcinoid and atypical SCLC, both harboring chromothripsis
ptdf <- read.csv("ctx.net.asclc.csv")

pdf("boxplot.ki67.grade.pdf", height=4, width=3)
boxplot(ptdf$max_ki67 ~ ptdf$group, ylim = c(0,100), frame = F, las = 1, main = paste0("Ki-67 index (p=", round(t.test(ptdf$max_ki67 ~ ptdf$group)$p.value, 7), ")"), ylab = "Ki-67 index (%)", col = c("#800080", "#FF9912"), names = c(paste0("Low-grade only (n=", nrow(ptdf[ptdf$group == "01_carcinoid",]), ")"), paste0("with high-grade part (n=", nrow(ptdf[ptdf$group == "02_asclc",]), ")")), xlab = "", outline = F)
stripchart(ptdf$max_ki67 ~ ptdf$group, pch = 19, col = rgb(0,0,0,.2), method = "jitter", vertical = T, add = T)
dev.off()


# Figure 6b: survival comparisons
fit <- survfit(Surv(os_months2, os_event2) ~ group, data = ptdf)
time_points <- seq(0, 150, by = 50)
groups <- levels(as.factor(ptdf$group))

## Matrix to store number at risk
n_risk_mat <- matrix(NA, nrow = length(groups), ncol = length(time_points))
rownames(n_risk_mat) <- groups
colnames(n_risk_mat) <- time_points

## Loop over groups
for (i in seq_along(groups)) {
  grp <- groups[i]
  dat_sub <- subset(ptdf, group == grp)
  fit_sub <- survfit(Surv(os_months2, os_event2) ~ 1, data = dat_sub)
  s <- summary(fit_sub, times = time_points)
  
  n_risk_vec <- rep(0, length(time_points))
  matched_times <- match(s$time, time_points)
  n_risk_vec[matched_times[!is.na(matched_times)]] <- s$n.risk
  n_risk_mat[i, ] <- n_risk_vec
}

pdf("km.os.grade.pdf", height=5.5, width=5)
par(mar = c(10.1,4.1,2.1,2.1))
plot(fit, col=c("#800080", "#FF9912"), las=1, xlab = "Time (months)", ylab = "Fraction of survival (%)", main = paste0("Log rank p=", survdiff(Surv(os_months2, os_event2) ~ group, data = ptdf)$p), lwd=1, mark.time = T)
legend("topright", legend = c("Low-grade only", "With high-grade part"), lwd = 1, col=c("#800080", "#FF9912"))

par(xpd = TRUE)
y_base <- -0.35
y_step <- -0.05

for (i in 1:length(groups)) {
  text(time_points, y_base + (i - 1) * y_step, labels = n_risk_mat[i, ])
  text(-10, y_base + (i - 1) * y_step, labels = groups[i], adj = 1)
}

dev.off()


# Figure 6c: swimmers' plot
hx <- read.csv("txhistory.csv")
pdf("swimmer.pdf", height=7, width=7)
plot(c(), xlim = c(0, max(ptdf$os_months2*1.1)), ylim = c(0, nrow(ptdf) + 1), frame = F, xaxt = 'n', yaxt = 'n', xlab = "Time (month)", ylab = "Patient")
axis(1, at = seq(0, 180, 12), cex.axis = 0.6)
axis(2, at = seq(0, nrow(ptdf)+1, 1), las = 1, labels = c("", rev(ptdf$study_id), ""), cex.axis = 0.6)
for (i in 1:nrow(ptdf)){
  if (ptdf$os_event2[i] == 0){
    if (ptdf$initial_ds[i] == "non_metastatic"){
      if (!is.na(ptdf$rec_months[i])){
        segments(0, nrow(ptdf)+1-i, ptdf$rec_months[i], nrow(ptdf)+1-i, lwd = 1, col = "#C6E2FF")
        arrows(ptdf$rec_months[i], nrow(ptdf)+1-i, ptdf$os_months2[i]+0.1, nrow(ptdf)+1-i, lwd = 1, col = "#4682B4", length = 0.05)
      } else {
        arrows(0, nrow(ptdf)+1-i, ptdf$os_months2[i]+0.1, nrow(ptdf)+1-i, lwd = 1, col = "#C6E2FF", length = 0.05)
      }
    } else {
      arrows(0, nrow(ptdf)+1-i, ptdf$os_months2[i]+0.1, nrow(ptdf)+1-i, lwd = 1, col = "#4682B4", length = 0.05)
    }
  } else {
    if (ptdf$initial_ds[i] == "non_metastatic"){
      if (!is.na(ptdf$rec_months[i])){
        segments(0, nrow(ptdf)+1-i, ptdf$rec_months[i], nrow(ptdf)+1-i, lwd = 1, col = "#C6E2FF")
        segments(ptdf$rec_months[i], nrow(ptdf)+1-i, ptdf$os_months2[i]+0.1, nrow(ptdf)+1-i, lwd = 1, col = "#4682B4")
      } else {
        segments(0, nrow(ptdf)+1-i, ptdf$os_months2[i]+0.1, nrow(ptdf)+1-i, lwd = 1, col = "#C6E2FF")
      }
    } else {
      segments(0, nrow(ptdf)+1-i, ptdf$os_months2[i]+0.1, nrow(ptdf)+1-i, lwd = 1, col = "#4682B4")
    }
  }
  if (sa$SA_eligible[sa$infofield == ptdf$infofield[i]] == "yes"){
    rect(as.numeric((sa$SA_start[sa$infofield == ptdf$infofield[i]] - ptdf$day1[i])*12/365), nrow(ptdf)+0.55-i, as.numeric((sa$SA_end[sa$infofield == ptdf$infofield[i]] - ptdf$day1[i])*12/365), nrow(ptdf)+1.45-i, col = "red", border = F)
  }
  if (nrow(hx[hx$infofield == ptdf$infofield[i] & hx$info != "overview",]) != 0){
    tempdf <- hx[hx$infofield == ptdf$infofield[i] & hx$info != "overview",]
    for (j in 1:nrow(tempdf)){
      rect(as.numeric((tempdf$start_date[j] - ptdf$day1[i])*12/365), nrow(ptdf)+0.7-i, as.numeric((tempdf$end_date[j] - ptdf$day1[i])*12/365), nrow(ptdf)+1.3-i, col = txcol[tempdf$txsum[j]], border = F)
    }
    rm(tempdf)
  }
}
legend("topright", legend = c(names(table(hx$txsum)), "Octreotide"), fill = c(txcol, "red"), border = F)
dev.off()

# Figure 6d: table of treatment responses

# Figure 6e: pathology images

# Figure 6f: radiographic images

