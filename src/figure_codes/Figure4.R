library(stringr)
library(survival)
library(ComplexHeatmap)
library(plotrix)

# Codes for Figure 4

# Figure 4a: conceptual diagram

# Figure 4b: Kaplan-Meier analysis of overall survival by genomic subgroup
pdf("km.os.gen_group.pdf", height=5.5, width=5)
par(mar = c(10.1,4.1,2.1,2.1))
plot(survfit(Surv(os_months2, os_event2) ~ gen_group, data = ca), col=c("#88A0DC", "#800080", "#EE2C2C", "#008000", "#0000FF", "#e78429", "#8B7D7B"), las=1, xlab = "Time (months)", ylab = "Fraction of survival (%)", main = paste0("Log rank p=", survdiff(Surv(os_months2, os_event2) ~ gen_group, data = ca)$p), lwd = 1, mark.time = T)
legend("topright", legend = paste0(names(table(ca$gen_group))), col=c("#88A0DC", "#800080", "#EE2C2C", "#008000", "#0000FF", "#e78429", "#8B7D7B"), lwd=1)

## Number at risk below plot
par(xpd = TRUE)
y_base <- -0.35
y_step <- -0.05

for (i in 1:length(groups)) {
  text(time_points, y_base + (i - 1) * y_step, labels = n_risk_mat[i, ], cex = 0.8)
  text(-10, y_base + (i - 1) * y_step, labels = groups[i], adj = 1)
}

dev.off()

# Figure 4c: bubble plot for stage by genomic subgroup
df <- as.data.frame(table(ca$stage_simple, ca$gen_group))
df$x <- c(rep(2, 4), rep(4, 4), rep(6, 4), rep(8, 4), rep(10, 4), rep(12, 4), rep(14, 4))
df$y <- c(rep(seq(2, 8, 2), 7))

df$Var1 <- as.character(df$Var1)
df$Var2 <- as.character(df$Var2)
df$Freq <- as.numeric(df$Freq)

k <- 7
df$z <- sqrt(df$Freq)/k
df$color <- rep(c("#88A0DC", "#800080", "#EE2C2C", "#008000", "#0000FF", "#e78429", "#8B7D7B"), each = 4)

pdf("bubble.stage.gen_group.pdf", height=6, width=6)
plot(c(), xlim = c(0,15), ylim = c(0,10), xlab = "", ylab = "Stage", las = 1, frame = F, xaxt = "n", yaxt = "n")
axis(1, at = c(2, 4, 6, 8, 10, 12, 14), labels = names(table(ca$gen_group)), tick = T, las = 2)
for (i in 1:nrow(df)){
  symbols(x = df$x[i], y = df$y[i], circles = df$z[i], inches = F, add=T, fg = F, bg = df$color[i])
}
axis(2, at = c(2, 4, 6, 8), labels = c("I", "II", "III", "IV"), tick = T, las=1)
symbols(x = 1, y = 8, circles = sqrt(1)/k, inches = F, add = T, fg = "black")
symbols(x = 1, y = 8, circles = sqrt(10)/k, inches = F, add = T, fg = "black")
symbols(x = 1, y = 8, circles = sqrt(30)/k, inches = F, add = T, fg = "black")
symbols(x = 1, y = 8, circles = sqrt(100)/k, inches = F, add = T, fg = "black")
dev.off()

# Figure 4d: box and whisker plot for ki67
pdf("boxplot.ki67.gen_group.pdf", height=4, width=4)
boxplot(ca$max_ki67 ~ ca$gen_group, ylim = c(0,60), frame = F, las = 2, main = paste0("Ki-67 index (n=", nrow(ca), ")"), ylab = "Ki-67 index (max %)", col=c("#88A0DC", "#800080", "#EE2C2C", "#008000", "#0000FF", "#e78429", "#8B7D7B"), names = c("Mutation negative", "Chromothripsis", "KRAS/NF1", "EIF1AX", "MEN1", "ARID1A-only", "Others"), xlab = "", outline = F, boxlwd = 2/3, medlwd = 2/3)
stripchart(ca$max_ki67 ~ ca$gen_group, pch = 19, col = rgb(0,0,0,.2), method = "jitter", vertical = T, add = T)
dev.off()

t.test(ca$max_ki67[ca$gen_group %in% c("0_mut_neg", "1_chromothripsis")] ~ ca$gen_group[ca$gen_group %in% c("0_mut_neg", "1_chromothripsis")])
t.test(ca$max_ki67[ca$gen_group %in% c("1_chromothripsis", "2_kras_nf1")] ~ ca$gen_group[ca$gen_group %in% c("1_chromothripsis", "2_kras_nf1")])
t.test(ca$max_ki67[ca$gen_group %in% c("1_chromothripsis", "3_eif1ax")] ~ ca$gen_group[ca$gen_group %in% c("1_chromothripsis", "3_eif1ax")])
t.test(ca$max_ki67[ca$gen_group %in% c("1_chromothripsis", "4_men1")] ~ ca$gen_group[ca$gen_group %in% c("1_chromothripsis", "4_men1")])
t.test(ca$max_ki67[ca$gen_group %in% c("1_chromothripsis", "5_arid1a_only")] ~ ca$gen_group[ca$gen_group %in% c("1_chromothripsis", "5_arid1a_only")])

# Figure 4e: box and whisker plot for fga
pdf("boxplot.fga.gen_group.pdf", height=4, width=4)
boxplot(ca$facets_fga[!is.na(ca$facets_fga)] ~ ca$gen_group[!is.na(ca$facets_fga)], ylim = c(0,1), frame = F, las = 2, main = paste0("Fraction of nondiploid genome (n=", nrow(ca[!is.na(ca$facets_fga),]), ")"), ylab = "Fraction of nondiploid genome", col=c("#88A0DC", "#800080", "#EE2C2C", "#008000", "#0000FF", "#e78429", "#8B7D7B"), names = c("Mutation negative", "Chromothripsis", "KRAS/NF1", "EIF1AX", "MEN1", "ARID1A-only", "Others"), xlab = "", outline = F, boxlwd = 2/3, medlwd = 2/3)
stripchart(ca$facets_fga[!is.na(ca$facets_fga)] ~ ca$gen_group[!is.na(ca$facets_fga)], pch = 19, col = rgb(0,0,0,.2), method = "jitter", vertical = T, add = T)
dev.off()

t.test(ca$facets_fga[ca$gen_group %in% c("0_mut_neg", "1_chromothripsis")] ~ ca$gen_group[ca$gen_group %in% c("0_mut_neg", "1_chromothripsis")])
t.test(ca$facets_fga[ca$gen_group %in% c("1_chromothripsis", "3_eif1ax")] ~ ca$gen_group[ca$gen_group %in% c("1_chromothripsis", "3_eif1ax")])
t.test(ca$facets_fga[ca$gen_group %in% c("1_chromothripsis", "5_arid1a_only")] ~ ca$gen_group[ca$gen_group %in% c("1_chromothripsis", "5_arid1a_only")])

# Figure 4f: comprehensive IHC
ndf <- read.csv("comprehensive.ihc.csv")
mutix <- ndf[,c(which(colnames(ndf) == "gen_group"),which(colnames(ndf) == "sex"),which(colnames(ndf) == "histology"),which(colnames(ndf) == "bin_otp"),which(colnames(ndf) == "bin_ascl1"),which(colnames(ndf) == "bin_ttf1"),which(colnames(ndf) == "bin_hnf4a"),which(colnames(ndf) == "bin_dll3"),which(colnames(ndf) == "bin_sez6"),which(colnames(ndf) == "ihc_group"))]
mutix <- as.matrix(mutix)
mutix <- t(mutix)
mutix <- gsub("02_atypical_carcinoid", "atypical_carcinoid", mutix)
mutix <- gsub("01_typical_carcinoid", "typical_carcinoid", mutix)
mutix <- gsub("03_no_resection", "", mutix)
mutix <- gsub("1_chromothripsis", "chromothripsis", mutix)
mutix <- gsub("2_kras_nf1", "kras_nf1", mutix)
mutix <- gsub("3_eif1ax", "eif1ax", mutix)
mutix <- gsub("4_men1", "men1", mutix)
mutix <- gsub("5_arid1a_only", "arid1a_only", mutix)
mutix <- gsub("0_mut_neg", "mut_neg", mutix)
mutix <- gsub("6_others", "others", mutix)

col = c(Unknown = "#4169E1", chromothripsis = "#800080", kras_nf1 = "#EE2C2C", eif1ax = "#008000", men1 = "#0000FF", arid1a_only = "#e78429", mut_neg = "#88A0DC", others = "#8B7D7B", atypical_carcinoid = "#FF9912", typical_carcinoid = "#3D59AB", msi = "#6E8B3D", present = "#191970", focalamp = "#EE2C2C", never = "#B0E2FF", former = "#33A1C9", current = "#191970", Male = "#104E8B", Female = "#CD6889", tx_naive = "#F5DEB3", hormone_only = "#8B7E66", chemo_xrt_rrt = "#292421", A1 = "#FFF68F", A2 = "#9ACD32", B = "#1E90FF", not_available = "#E0E0E0", other = "#708090", positive = "#8B4513", negative = "#C6E2FF")

test.oncoprint = oncoPrint(mutix, get_type = function(x) strsplit(x, ";")[[1]],
                           alter_fun = list(
                             background = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#E0E0E0", col = NA)),
                             Male = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Male"], col = NA)),
                             Female = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Female"], col = NA)),
                             atypical_carcinoid = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["atypical_carcinoid"], col = NA)),
                             typical_carcinoid = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["typical_carcinoid"], col = NA)),
                             chromothripsis = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["chromothripsis"], col = NA)),
                             kras_nf1 = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["kras_nf1"], col = NA)),
                             eif1ax = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["eif1ax"], col = NA)),
                             men1 = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["men1"], col = NA)),
                             arid1a_only = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["arid1a_only"], col = NA)),
                             mut_neg = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["mut_neg"], col = NA)),
                             others = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["others"], col = NA)),
                             Unknown = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Unknown"], col = NA)),
                             A1 = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["A1"], col = NA)),
                             A2 = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["A2"], col = NA)),
                             B = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["B"], col = NA)),
                             other = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["other"], col = NA)),
                             not_available = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["not_available"], col = NA)),
                             positive = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["positive"], col = NA)),
                             negative = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["negative"], col = NA)),
                             
                             chemo_xrt_rrt = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["chemo_xrt_rrt"], col = NA))
                           ), 
                           top_annotation = HeatmapAnnotation(
                             age_dx = as.matrix(ndf[,"age_dx"]),
                             col = list(age_dx = colorRamp2(c(15,95), c("#000000", "#FFFFFF")))
                           ),
                           right_annotation = rowAnnotation(
                             row_barplot = anno_oncoprint_barplot(c("positive", "negative", "A1", "A2", "B", "other"), show_fraction = T)
                           ),
                           col = col, row_order = 1:nrow(mutix), column_order = 1:ncol(mutix), remove_empty_columns = FALSE, column_title = paste0("Pulmonary Carcinoid IMPACT and IHC (# patient = ", nrow(ndf), "; # tumor = ", nrow(ndf), ")")
                           , heatmap_legend_param = list(title = "Alterations", nrow = 2, title_position = "leftcenter"))

pdf(paste0("oncoprint.ihc.pdf"), height=3.2, width=9.5)
draw(test.oncoprint, heatmap_legend_side = "bottom")
dev.off()

# Figure 4g: multivariate Cox proportional hazards model
sdf <- read.csv("multivariate.csv")
pdf("coxph.os.ggforest.pdf", height=10, width=10)
cox <- coxph(Surv(os_months2, os_event2) ~ dic_age + sex + stage_simple + tri_path + dic_ki67 + dic_fga + dic_wgd + gen_group, data = sdf)
ggforest(cox, data = sdf)
dev.off()

# Figure 4h: conceptual summary of genomic subgroups





