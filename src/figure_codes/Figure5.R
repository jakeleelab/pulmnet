library(stringr)
library(survival)

# Codes for Figure 5

# Figure 5a: RFS by ki67
ra <- read.csv("rfs.csv")
fit <- survfit(Surv(rfs_month, recurrence) ~ surg_ki67_group, data = ra)
time_points <- seq(0, 200, by = 40) # Time points for number at risk
groups <- levels(as.factor(ra$surg_ki67_group)) # Unique groups based on ki67 3 and 10

## Matrix to store number at risk
n_risk_mat <- matrix(NA, nrow = length(groups), ncol = length(time_points))
rownames(n_risk_mat) <- groups
colnames(n_risk_mat) <- time_points

## Loop over groups
for (i in seq_along(groups)) {
  grp <- groups[i]
  dat_sub <- subset(ra, surg_ki67_group == grp)
  fit_sub <- survfit(Surv(rfs_month, recurrence) ~ 1, data = dat_sub)
  s <- summary(fit_sub, times = time_points)
  
  n_risk_vec <- rep(0, length(time_points))
  matched_times <- match(s$time, time_points)
  n_risk_vec[matched_times[!is.na(matched_times)]] <- s$n.risk
  n_risk_mat[i, ] <- n_risk_vec
}


pdf("km.rfs.surgicalki67.pdf", height=5.5, width=5)
par(mar = c(10.1,4.1,2.1,2.1))
plot(fit, col=c("darkblue", "darkgreen", "darkred"), las=1, xlab = "Time (months)", ylab = "Fraction of no recurrence (%)", main = paste0("Groups 1 vs 2: p=", survdiff(Surv(rfs_month, recurrence) ~ surg_ki67_group, data = ra[!is.na(ra$surg_ki67) & ra$surg_ki67_group %in% c("01_below3", "02_3to10"),])$p, " 2 vs 3: p=", survdiff(Surv(rfs_month, recurrence) ~ surg_ki67_group, data = ra[!is.na(ra$surg_ki67) & ra$surg_ki67_group %in% c("02_3to10", "03_above10"),])$p), lwd=1, mark.time = T)
legend("topright", legend = c("surg Ki-67 < 3%", "3-10%", "> 10%"), lwd = 1, col=c("darkblue", "darkgreen", "darkred"))

par(xpd = TRUE)
y_base <- -0.35
y_step <- -0.05

for (i in 1:length(groups)) {
  text(time_points, y_base + (i - 1) * y_step, labels = n_risk_mat[i, ])
  text(-10, y_base + (i - 1) * y_step, labels = groups[i], adj = 1)
}

dev.off()

## Coxph statistics
survdiff(Surv(rfs_month, recurrence) ~ surg_ki67_group, data = ra[!is.na(ra$surg_ki67) & ra$surg_ki67_group %in% c("01_below3", "02_3to10"),])
survdiff(Surv(rfs_month, recurrence) ~ surg_ki67_group, data = ra[!is.na(ra$surg_ki67) & ra$surg_ki67_group %in% c("02_3to10", "03_above10"),])


# Figure 5b: RFS by genomic risk
fit <- survfit(Surv(rfs_month, recurrence) ~ risk_group, data = ra)
time_points <- seq(0, 180, by = 60)
groups <- levels(as.factor(ra$risk_group))

## Matrix to store number at risk
n_risk_mat <- matrix(NA, nrow = length(groups), ncol = length(time_points))
rownames(n_risk_mat) <- groups
colnames(n_risk_mat) <- time_points

## Loop over groups
for (i in seq_along(groups)) {
  grp <- groups[i]
  dat_sub <- subset(ra, risk_group == grp)
  fit_sub <- survfit(Surv(rfs_month, recurrence) ~ 1, data = dat_sub)
  s <- summary(fit_sub, times = time_points)
  
  n_risk_vec <- rep(0, length(time_points))
  matched_times <- match(s$time, time_points)
  n_risk_vec[matched_times[!is.na(matched_times)]] <- s$n.risk
  n_risk_mat[i, ] <- n_risk_vec
}

pdf("km.rfs.genomicrisk.pdf", height=3.7, width=3.5)
par(mar = c(6.1,4.1,2.1,2.1))
plot(fit, col=c("#00BFFF", "#FF9912", "#FF1493"), las=1, xlab = "Time (months)", ylab = "Fraction of no recurrence (%)", main = paste0("Log rank p=", survdiff(Surv(rfs_month, recurrence) ~ risk_group, data = ra)$p), lwd = 1, mark.time = T, xaxt = 'n')
axis(1, at = seq(0, 180, by = 60))
tbl <- table(ra$risk_group)

par(xpd = TRUE)  # allow writing outside margins
y_base <- -0.35
y_step <- -0.05

for (i in 1:length(groups)) {
  text(time_points, y_base + (i - 1) * y_step, labels = n_risk_mat[i, ])
  text(-10, y_base + (i - 1) * y_step, labels = groups[i], adj = 1)
}

dev.off()

survdiff(Surv(rfs_month, recurrence) ~ risk_group, data = ra[ra$risk_group %in% c("01_low_risk", "02_intermediate_risk"),])$p
survdiff(Surv(rfs_month, recurrence) ~ risk_group, data = ra[ra$risk_group %in% c("02_intermediate_risk", "03_high_risk"),])$p


# Figure 5c: RFS by minimal risk group
fit <- survfit(Surv(rfs_month, recurrence) ~ ult_group, data = ra)
time_points <- seq(0, 180, by = 60)
groups <- levels(as.factor(ra$ult_group))

## Matrix to store number at risk
n_risk_mat <- matrix(NA, nrow = length(groups), ncol = length(time_points))
rownames(n_risk_mat) <- groups
colnames(n_risk_mat) <- time_points

## Loop over groups
for (i in seq_along(groups)) {
  grp <- groups[i]
  dat_sub <- subset(ra, ult_group == grp)
  fit_sub <- survfit(Surv(rfs_month, recurrence) ~ 1, data = dat_sub)
  s <- summary(fit_sub, times = time_points)
  
  n_risk_vec <- rep(0, length(time_points))
  matched_times <- match(s$time, time_points)
  n_risk_vec[matched_times[!is.na(matched_times)]] <- s$n.risk
  n_risk_mat[i, ] <- n_risk_vec
}


pdf("km.rfs.combinedminimalrisk.pdf", height=3.7, width=3.5)
par(mar = c(6.1,4.1,2.1,2.1))
plot(fit, col=c("#0000FF", "#FF4500"), las=1, xlab = "Time (months)", ylab = "Fraction of no recurrence (%)", main = paste0("Log rank p=", survdiff(Surv(rfs_month, recurrence) ~ ult_group, data = ra)$p), lwd = 1, mark.time = T, xaxt = 'n')
axis(1, at = seq(0, 180, by = 60))
tbl <- table(ra$ult_group)

par(xpd = TRUE)
y_base <- -0.35
y_step <- -0.05

for (i in 1:length(groups)) {
  text(time_points, y_base + (i - 1) * y_step, labels = n_risk_mat[i, ])
  text(-10, y_base + (i - 1) * y_step, labels = groups[i], adj = 1)
}
dev.off()


# Figure 5d: Stage IV OS by ki67
fit <- survfit(Surv(os_months2, os_event2) ~ ki67_group, data = ca[ca$stage_simple == "IV" & ca$ki67_group != "04_unknown",])
time_points <- seq(0, 180, by = 60)
groups <- levels(as.factor(ca$ki67_group[ca$stage_simple == "IV" & ca$ki67_group != "04_unknown"]))

## Matrix to store number at risk
n_risk_mat <- matrix(NA, nrow = length(groups), ncol = length(time_points))
rownames(n_risk_mat) <- groups
colnames(n_risk_mat) <- time_points

## Loop over groups
for (i in seq_along(groups)) {
  grp <- groups[i]
  dat_sub <- subset(ca[ca$stage_simple == "IV" & ca$ki67_group != "04_unknown",], ki67_group == grp)
  fit_sub <- survfit(Surv(os_months2, os_event2) ~ 1, data = dat_sub)
  s <- summary(fit_sub, times = time_points)
  
  n_risk_vec <- rep(0, length(time_points))
  matched_times <- match(s$time, time_points)
  n_risk_vec[matched_times[!is.na(matched_times)]] <- s$n.risk
  n_risk_mat[i, ] <- n_risk_vec
}

pdf("km.os.maxki67.stageIVonly.pdf", height=3.7, width=3.5)
par(mar = c(6.1,4.1,2.1,2.1))
plot(fit, col=c("darkblue", "darkgreen", "darkred"), las=1, xlab = "Time (months)", ylab = "Fraction of survival (%)", main = paste0("Log rank p=", survdiff(Surv(os_months2, os_event2) ~ ki67_group, data = ca[ca$stage_simple == "IV" & ca$ki67_group != "04_unknown",])$p), lwd = 1, mark.time = T, xaxt = 'n')
axis(1, at = seq(0, 180, by = 60))
tbl <- table(ca$ki67_group[ca$stage_simple == "IV" & ca$ki67_group != "04_unknown"])

par(xpd = TRUE)  # allow writing outside margins
y_base <- -0.35
y_step <- -0.05

for (i in 1:length(groups)) {
  text(time_points, y_base + (i - 1) * y_step, labels = n_risk_mat[i, ])
  text(-10, y_base + (i - 1) * y_step, labels = groups[i], adj = 1)
}

dev.off()

survdiff(Surv(os_months2, os_event2) ~ ki67_group, data = ca[ca$ki67_group %in% c("01_below3", "02_3to10") & ca$stage_simple == "IV" & ca$ki67_group != "04_unknown",])
survdiff(Surv(os_months2, os_event2) ~ ki67_group, data = ca[ca$ki67_group %in% c("02_3to10", "03_above10") & ca$stage_simple == "IV" & ca$ki67_group != "04_unknown",])


# Figure 5e: Stage IV OS by genomic risk
fit <- survfit(Surv(os_months2, os_event2) ~ risk_group, data = ca[ca$stage_simple == "IV" & ca$ki67_group != "04_unknown",])
time_points <- seq(0, 180, by = 60)
groups <- levels(as.factor(ca$risk_group[ca$stage_simple == "IV" & ca$ki67_group != "04_unknown"]))

## Matrix to store number at risk
n_risk_mat <- matrix(NA, nrow = length(groups), ncol = length(time_points))
rownames(n_risk_mat) <- groups
colnames(n_risk_mat) <- time_points

## Loop over groups
for (i in seq_along(groups)) {
  grp <- groups[i]
  dat_sub <- subset(ca[ca$stage_simple == "IV" & ca$ki67_group != "04_unknown",], risk_group == grp)
  fit_sub <- survfit(Surv(os_months2, os_event2) ~ 1, data = dat_sub)
  s <- summary(fit_sub, times = time_points)
  
  n_risk_vec <- rep(0, length(time_points))
  matched_times <- match(s$time, time_points)
  n_risk_vec[matched_times[!is.na(matched_times)]] <- s$n.risk
  n_risk_mat[i, ] <- n_risk_vec
}

pdf("km.os.genomicrisk.stageIVonly.pdf", height=3.7, width=3.5)
par(mar = c(6.1,4.1,2.1,2.1))
plot(fit, col=c("#00BFFF", "#FF9912", "#FF1493"), las=1, xlab = "Time (months)", ylab = "Fraction of survival (%)", main = paste0("Log rank p=", survdiff(Surv(os_months2, os_event2) ~ risk_group, data = ca[ca$stage_simple == "IV" & ca$ki67_group != "04_unknown",])$p), lwd = 1, mark.time = T, xaxt = 'n')
axis(1, at = seq(0, 180, by = 60))
tbl <- table(ca$risk_group[ca$stage_simple == "IV" & ca$ki67_group != "04_unknown"])

par(xpd = TRUE)
y_base <- -0.35
y_step <- -0.05

for (i in 1:length(groups)) {
  text(time_points, y_base + (i - 1) * y_step, labels = n_risk_mat[i, ])
  text(-10, y_base + (i - 1) * y_step, labels = groups[i], adj = 1)
}

dev.off()

survdiff(Surv(os_months2, os_event2) ~ risk_group, data = ca[ca$risk_group %in% c("01_low_risk", "02_intermediate_risk") & ca$stage_simple == "IV" & ca$ki67_group != "04_unknown",])
survdiff(Surv(os_months2, os_event2) ~ risk_group, data = ca[ca$risk_group %in% c("02_intermediate_risk", "03_high_risk") & ca$stage_simple == "IV" & ca$ki67_group != "04_unknown",])




