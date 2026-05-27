library(stringr)
library(ComplexHeatmap)

# Codes for Figure 3
geneofi <- c("EIF1AX", "ARID1A", "MEN1", "DAXX", "ATRX", "ATM", "SETD2", "KRAS", "PTEN", "TSC2", "CDKN1B", "CDKN2A", "CDKN2B", "CCND1", "MDM2", "CDK4", "MDM4", "MYCL", "RB1", "TP53")

# Figure 3a: oncoprint for neuroendocrine tumors from three major anatomical sites
impact <- read.csv("pancreas.net.csv")
impact <- read.csv("lung.net.csv")
impact <- read.csv("gi.net.csv")

mutix <- impact[,c(which(colnames(impact) %in% c("PATIENT_ID", "SAMPLE_TYPE", "WGD")), c(which(colnames(impact)==geneofi[1]):which(colnames(impact)==geneofi[length(geneofi)])))]
mutix <- mutix[,c(which(colnames(mutix)=="SAMPLE_TYPE"):which(colnames(mutix)==geneofi[length(geneofi)]))]
mutix <- as.matrix(mutix)
mutix <- t(mutix)
mutix <- gsub("Intron", "", mutix)
mutix <- gsub("nonframeshift_deletion", "inframe", mutix)
mutix <- gsub("nonsynonymous_SNV", "missense", mutix)

col = c(missense = "#006400", nonsense = "#000000", splice = "#FF34B3", inframe = "#CD0000", frameshift = "#0000CD", biallelic = "#BABABA", Amp = "#EE2C2C", Del = "#00688B", sv = "#9ACD32", upstrea_reg = "#00BFFF", nonstop = "#800080", Primary = "#DC143C", Metastasis = "#008B45", Unknown = "#4169E1", after = "#FA8072", before = "#63B8FF", wgd = "#000000", chromothripsis = "#800080", atypical_carcinoid = "#FF9912", small_cell_carcinoma = "#FF1493", msi = "#6E8B3D")

test.oncoprint = oncoPrint(mutix, get_type = function(x) strsplit(x, ";")[[1]],
                           alter_fun = list(
                             background = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#E0E0E0", col = NA)),
                             Del = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Del"], col = NA)),
                             Amp = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Amp"], col = NA)),
                             sv = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.7, gp = gpar(fill = col["sv"], col = NA)),
                             missense = function(x, y, w, h) grid.rect(x, y, w*0.7, h*0.35, gp = gpar(fill = col["missense"], col = NA)),
                             nonsense = function(x, y, w, h) grid.rect(x, y, w*0.7, h*0.35, gp = gpar(fill = col["nonsense"], col = NA)),
                             splice = function(x, y, w, h) grid.rect(x, y, w*0.7, h*0.35, gp = gpar(fill = col["splice"], col = NA)),
                             inframe = function(x, y, w, h) grid.rect(x, y, w*0.7, h*0.35, gp = gpar(fill = col["inframe"], col = NA)),
                             frameshift = function(x, y, w, h) grid.rect(x, y, w*0.7, h*0.35, gp = gpar(fill = col["frameshift"], col = NA)),
                             upstrea_reg = function(x, y, w, h) grid.rect(x, y, w*0.7, h*0.35, gp = gpar(fill = col["upstrea_reg"], col = NA)),
                             nonstop = function(x, y, w, h) grid.rect(x, y, w*0.7, h*0.35, gp = gpar(fill = col["nonstop"], col = NA)),
                             wgd = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.7, gp = gpar(fill = col["wgd"], col = NA)),
                             msi = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.7, gp = gpar(fill = col["msi"], col = NA)),
                             chromothripsis = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.7, gp = gpar(fill = col["chromothripsis"], col = NA)),
                             after = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.7, gp = gpar(fill = col["after"], col = NA)),
                             before = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.7, gp = gpar(fill = col["before"], col = NA)),
                             Primary = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Primary"], col = NA)),
                             Metastasis = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Metastasis"], col = NA)),
                             Unknown = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Unknown"], col = NA))
                           ), 
                           top_annotation = HeatmapAnnotation(
                             column_barplot = anno_oncoprint_barplot(c("sv", "Amp", "Del", "missense", "nonsense", "splice", "inframe", "frameshift", "upstrea_reg", "nonstop"), show_fraction = T)
                           ),
                           right_annotation = rowAnnotation(
                             row_barplot = anno_oncoprint_barplot(c("sv", "Amp", "Del", "missense", "nonsense", "splice", "inframe", "frameshift", "upstrea_reg", "nonstop"), show_fraction = T)
                           ),
                           col = col, row_order = 1:nrow(mutix), column_order = 1:ncol(mutix), remove_empty_columns = FALSE, column_title = paste0("IMPACT (# patient = ", nrow(impact), ")")
                           , heatmap_legend_param = list(title = "Alterations", nrow = 2, title_position = "leftcenter"))

pdf(paste0("oncoprint.pancreas.net.pdf"), height=4, width=8)
pdf(paste0("oncoprint.lung.net.pdf"), height=4, width=8)
pdf(paste0("oncoprint.gi.net.pdf"), height=4, width=8)
draw(test.oncoprint, heatmap_legend_side = "bottom")
dev.off()


# Figure 3b: chromoothripsis cases of pancreas and gastrointestinal origin
## Pancreatic case, P01
pdf(paste0("P01.chr11.chromothripsis.pdf"), height = 3.6, width = 3.6)
plot(c(), xlim = c(1, hg19$length[hg19$chr == "11"]), ylim = c(0, 50), frame = F, xlab = "Position", ylab = "Copy number", main = "Chromosome 11 chromothripsis in P01", las = 1)
for (i in c("CCND1")){
  abline(v = genedf$start[genedf$gene == i], col = "pink")
}

load("P01_hisens.Rdata")
df <- out$jointseg
df <- df[df$chrom == 11,]
df$tcn_estimate <- NA
f <- 0.95
p <- 2
for (i in 1:nrow(df)){
  df$tcn_estimate[i] <- ((f*p + (1-f)*2)*2^df$cnlr[i] - (1-f)*2)/f 
}
points(df$maploc, df$tcn_estimate, col = "#9FB6CD", pch=20, cex=0.1)

df <- cncf[cncf$study_id == "P01_hisens" & cncf$chrom == 11,]
df$lcn[is.na(df$lcn)] <- 1
segments(df$loc.start, df$tcn, df$loc.end, df$tcn, lwd = 1, col = "#DC143C")
segments(df$loc.start, df$lcn, df$loc.end, df$lcn, lwd = 1, col = "#27408B")

imp_num <- max(df$tcn)
rect(cytoband$start[cytoband$chr == "11"], -imp_num*0.2, cytoband$end[cytoband$chr == "11"], -imp_num*0.15, col = cytoband$color[cytoband$chr == "11"], xpd = T, lwd = 0.5)
dev.off()

pdf(paste0("P01.chr12.chromothripsis.pdf"), height = 3.6, width = 3.6)
plot(c(), xlim = c(1, hg19$length[hg19$chr == "12"]), ylim = c(0, 50), frame = F, xlab = "Position", ylab = "Copy number", main = "Chromosome 12 chromothripsis in P01", las = 1)
for (i in c("CDK4", "MDM2")){
  abline(v = genedf$start[genedf$gene == i], col = "pink")
}

load("P01_hisens.Rdata")
df <- out$jointseg
df <- df[df$chrom == 12,]
df$tcn_estimate <- NA
f <- 0.95
p <- 2
for (i in 1:nrow(df)){
  df$tcn_estimate[i] <- ((f*p + (1-f)*2)*2^df$cnlr[i] - (1-f)*2)/f 
}
points(df$maploc, df$tcn_estimate, col = "#9FB6CD", pch=20, cex=0.1)

df <- cncf[cncf$study_id == "P01_hisens" & cncf$chrom == 12,]
df$lcn[is.na(df$lcn)] <- 1

segments(df$loc.start, df$tcn, df$loc.end, df$tcn, lwd = 1, col = "#DC143C")
segments(df$loc.start, df$lcn, df$loc.end, df$lcn, lwd = 1, col = "#27408B")

imp_num <- max(df$tcn)
rect(cytoband$start[cytoband$chr == "12"], -imp_num*0.2, cytoband$end[cytoband$chr == "12"], -imp_num*0.15, col = cytoband$color[cytoband$chr == "12"], xpd = T, lwd = 0.5)
dev.off()

## Small intestinal case, S01
pdf(paste0("S01.chr11.chromothripsis.pdf"), height = 3.6, width = 3.5)
plot(c(), xlim = c(1, hg19$length[hg19$chr == "11"]), ylim = c(0, 50), frame = F, xlab = "Position", ylab = "Copy number", main = "Chromosome 11 chromothripsis in S01", las = 1)
for (i in c("CCND1", "PAK1")){
  abline(v = genedf$start[genedf$gene == i], col = "pink")
}

load("S01_hisens.Rdata")
df <- out$jointseg
df <- df[df$chrom == 11,]
df$tcn_estimate <- NA
f <- 0.88
p <- 2
for (i in 1:nrow(df)){
  df$tcn_estimate[i] <- ((f*p + (1-f)*2)*2^df$cnlr[i] - (1-f)*2)/f 
}
points(df$maploc, df$tcn_estimate, col = "#9FB6CD", pch=20, cex=0.1)

df <- cncf[cncf$study_id == "P01_hisens" & cncf$chrom == 11,]
df$lcn[is.na(df$lcn)] <- 1
segments(df$loc.start, df$tcn, df$loc.end, df$tcn, lwd = 1, col = "#DC143C")
segments(df$loc.start, df$lcn, df$loc.end, df$lcn, lwd = 1, col = "#27408B")
imp_num <- max(df$tcn)

rect(cytoband$start[cytoband$chr == "11"], -imp_num*0.2, cytoband$end[cytoband$chr == "11"], -imp_num*0.15, col = cytoband$color[cytoband$chr == "11"], xpd = T, lwd = 0.5)
dev.off()

pdf(paste0("S01.chr12.chromothripsis.pdf"), height = 3.6, width = 3.8)

plot(c(), xlim = c(1, hg19$length[hg19$chr == "12"]), ylim = c(0, 50), frame = F, xlab = "Position", ylab = "Copy number", main = "Chromosome 12 chromothripsis in S01", las = 1)
for (i in c("CDK4", "MDM2")){
  abline(v = genedf$start[genedf$gene == i], col = "pink")
}

load("S01_hisens.Rdata")
df <- out$jointseg
df <- df[df$chrom == 12,]
df$tcn_estimate <- NA
f <- 0.88
p <- 2
for (i in 1:nrow(df)){
  df$tcn_estimate[i] <- ((f*p + (1-f)*2)*2^df$cnlr[i] - (1-f)*2)/f 
}
points(df$maploc, df$tcn_estimate, col = "#9FB6CD", pch=20, cex=0.1)

df <- cncf[cncf$ID == "P01_hisens" & cncf$chrom == 12,]
df$lcn[is.na(df$lcn)] <- 1
segments(df$loc.start, df$tcn, df$loc.end, df$tcn, lwd = 1.5, col = "#DC143C")
segments(df$loc.start, df$lcn, df$loc.end, df$lcn, lwd = 1.5, col = "#27408B")
imp_num <- max(df$tcn)

rect(cytoband$start[cytoband$chr == "12"], -imp_num*0.2, cytoband$end[cytoband$chr == "12"], -imp_num*0.15, col = cytoband$color[cytoband$chr == "12"], xpd = T, lwd = 0.5)
dev.off()




