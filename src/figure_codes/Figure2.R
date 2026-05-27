library(stringr)

# Codes for Figure 2

# Figure 2a: genome-wide copy number of a case harboring chromothripsis
pdf(paste0("C11.genome.chromothripsis.pdf"), height = 5, width = 14)
plot(NULL, xlim = c(hg19$add_values[1], hg19$add_values[nrow(hg19)]+hg19$length[nrow(hg19)]), ylim = c(0, 25), ylab = "Allelic copy number", xlab = "", las=1, xaxt = "n", frame = FALSE, lwd = 2/3)
axis(side = 1, labels = F, at = hg19$add_values, lwd = 2/3)
text(x = hg19$text_location[1:24], y = -2.5, labels = hg19$chr[1:24], xpd=NA)
abline(v = hg19$add_values, col = "gray", lwd = 2/3)

## Loading raw data and transforming into integer copy-number estimates
load("C11_hisens.Rdata")
df <- out$jointseg
df$chr <- df$chrom
df$chr <- as.character(df$chr)
df$chr[df$chr == "23"] <- "X"
df$tcn_estimate <- NA
f <- 0.76 # purity estimate for C11
p <- 1.8 # ploidy estimate for C11
for (i in 1:nrow(df)){
  df$tcn_estimate[i] <- ((f*p + (1-f)*2)*2^df$cnlr[i] - (1-f)*2)/f 
}
df$plot.maploc <- 0
for (j in as.character(c(c(1:23), "X"))){
  df$plot.maploc[df$chr == j] <- df$maploc[df$chr == j] + hg19$add_values[hg19$chr == j]
}
points(df$plot.maploc, df$tcn_estimate, col = "#9FB6CD", pch=20, cex=0.01)

## Copy-number segments
cncf <- read.csv("cncf.csv")
df <- cncf[cncf$study_id == "C11_hisens",]
df$lcn[is.na(df$lcn)] <- 1 # Regions with retention of heterozygosity, can be case specific
df <- df[,c("chrom", "loc.start", "loc.end", "tcn", "lcn")]
colnames(df) <- c("chrom", "start", "end", "tcn", "lcn")
df$chrom <- as.character(df$chrom)
df$chrom[df$chrom == "23"] <- "X"
df$plot.start <- 0
df$plot.end <- 0
for (j in unique(df$chrom)){
  df$plot.start[df$chrom == j] <- df$start[df$chrom == j] + hg19$add_values[hg19$chr == j]
  df$plot.end[df$chrom == j] <- df$end[df$chrom == j] + hg19$add_values[hg19$chr == j]
}
segments(df$plot.start, df$tcn, df$plot.end, df$tcn, lwd = 1, col = "#DC143C")
segments(df$plot.start, df$lcn, df$plot.end, df$lcn, lwd = 1, col = "#27408B")

dev.off()

## Inset 1: chromothripsis in chromosome 2
pdf("C11.chr2.chromothripsis.pdf", height = 2.7, width = 3.5)
plot(NULL, xlim = c(0, hg19$length[hg19$chr == "2"]), ylim = c(0, 4), ylab = "Total copy number", xlab = "", las=1, xaxt = "n", frame = FALSE, main = paste0("C11 (maximal Ki-67 = ", ctx$max_ki67[ctx$study_id == "C11"], "%)"), lwd = 2/3)
load("C11_hisens.Rdata")
df <- out$jointseg
df$chr <- df$chrom
df$chr <- as.character(df$chr)
df$chr[df$chr == "23"] <- "X"
df$tcn_estimate <- NA
f <- 0.76 # purity estimate for C11
p <- 1.8 # ploidy estimate for C11
for (i in 1:nrow(df)){
  df$tcn_estimate[i] <- ((f*p + (1-f)*2)*2^df$cnlr[i] - (1-f)*2)/f 
}
points(df$maploc[df$chrom == "2"], df$tcn_estimate[df$chrom == "2"], col = "#9FB6CD", pch=19, cex=0.1)

df <- cncf[cncf$study_id == "C11_hisens",]
df$lcn[is.na(df$lcn)] <- 1
df <- df[,c("chrom", "loc.start", "loc.end", "tcn", "lcn")]
colnames(df) <- c("chrom", "start", "end", "tcn", "lcn")
df$chrom <- as.character(df$chrom)
df$chrom[df$chrom == "23"] <- "X"

segments(df$start[df$chrom == "2"], df$tcn[df$chrom == "2"], df$end[df$chrom == "2"], df$tcn[df$chrom == "2"], lwd = 1, col = "#DC143C")
segments(df$start[df$chrom == "2"], df$lcn[df$chrom == "2"]-0.01, df$end[df$chrom == "2"], df$lcn[df$chrom == "2"]-0.01, lwd = 1, col = "#27408B")
dev.off()

## Inset 2: chromothripsis in chromosome 11
pdf("C11.chr11.chromothripsis.pdf", height = 2.7, width = 3)
plot(NULL, xlim = c(0, hg19$length[hg19$chr == "11"]), ylim = c(0, 4), ylab = "Total copy number", xlab = "", las=1, xaxt = "n", frame = FALSE, main = paste0("C11 (maximal Ki-67 = ", ctx$max_ki67[ctx$study_id == "C11"], "%)"), lwd = 2/3)
load("C11_hisens.Rdata")
df <- out$jointseg
df$chr <- df$chrom
df$chr <- as.character(df$chr)
df$chr[df$chr == "23"] <- "X"
df$tcn_estimate <- NA
f <- 0.76 # purity estimate for C11
p <- 1.8 # ploidy estimate for C11
for (i in 1:nrow(df)){
  df$tcn_estimate[i] <- ((f*p + (1-f)*2)*2^df$cnlr[i] - (1-f)*2)/f 
}
points(df$maploc[df$chrom == "11"], df$tcn_estimate[df$chrom == "11"], col = "#9FB6CD", pch=19, cex=0.1)

df <- cncf[cncf$study_id == "C11_hisens",]
df$lcn[is.na(df$lcn)] <- 1
df <- df[,c("chrom", "loc.start", "loc.end", "tcn", "lcn")]
colnames(df) <- c("chrom", "start", "end", "tcn", "lcn")
df$chrom <- as.character(df$chrom)
df$chrom[df$chrom == "23"] <- "X"

segments(df$start[df$chrom == "11"], df$tcn[df$chrom == "11"], df$end[df$chrom == "11"], df$tcn[df$chrom == "11"], lwd = 1, col = "#DC143C")
segments(df$start[df$chrom == "11"], df$lcn[df$chrom == "11"]-0.01, df$end[df$chrom == "11"], df$lcn[df$chrom == "11"]-0.01, lwd = 1, col = "#27408B")
dev.off()

# Figure 2b: conceptual illustration

# Figure 2c: chromothripsis diagram
ctx <- read.csv("chromothripsis.details.csv")

pdf(paste0("chromothripsis.cn.heatmap.pdf"), height = 5.8, width = 9)
par(mar = c(2, 6.1, 2, 2.1))
plot(NULL, xlim = c(hg19$add_values[1], hg19$add_values[nrow(hg19)]+hg19$length[nrow(hg19)]), ylim = c(0, nrow(ctx)), ylab = "Sample", xlab = "", las=1, xaxt = "n", yaxt = "n", frame = FALSE, main = paste0("Copy number profile of group 1 tumors with QC pass (n=", nrow(ctx), ")"))
axis(side = 1, labels = F, at = hg19$add_values)
axis(side = 2, at = c(1:nrow(ctx))-0.5, labels = ctx$study_id, las = 1, cex.axis = 1)
abline(v = hg19$add_values, col = "darkgray")
text(x = hg19$text_location[1:24], y = -3, labels = hg19$chr[1:24], xpd=NA)
w <- 0
for (i in 1:nrow(ctx)){
  if (!(ctx$dmp_sample[i] %in% unique(cncf$sample))){
    infoline <- list.files(path = paste0("facets_files/"), pattern = ctx$study_id[i], full.names = T)
    infoline <- infoline[grepl("_hisens", infoline)]
    if (length(infoline) == 1){
      df <- read.csv(infoline, header=T, as.is=T, sep="\t")
    } else {
      if (sum(grepl("refit", infoline)) != 0){
        infoline <- infoline[grepl("refit", infoline)]
      }
      df <- read.csv(infoline[length(infoline)], header=T, as.is=T, sep="\t")
    }
    df$sample <- ctx$dmp_sample[i]
  } else if (nrow(cncf[cncf$sample == ctx$dmp_sample[i],]) != 0) {
    df <- cncf[cncf$ID == paste0(facets$tumor_sample_id[facets$id == ctx$facets_id[i]], "_hisens"),]
  } else {
    print(paste0("Sample ", ctx$dmp_sample[i], " has no copy number information"))
    next
  }
  df <- df[,c("chrom", "loc.start", "loc.end", "tcn.em", "lcn.em")]
  colnames(df) <- c("chrom", "start", "end", "tcn", "lcn")
  df$chrom <- as.character(df$chrom)
  df$chrom[df$chrom == "23"] <- "X"
  df$plot.start <- 0
  df$plot.end <- 0
  for (j in unique(df$chrom)){
    df$plot.start[df$chrom == j] <- df$start[df$chrom == j] + hg19$add_values[hg19$chr == j]
    df$plot.end[df$chrom == j] <- df$end[df$chrom == j] + hg19$add_values[hg19$chr == j]
  }
  df$tcn[df$tcn > 11] <- 11
  for (j in 1:nrow(df)){
    rect(df$plot.start[j], w, df$plot.end[j], w+1, col = modcol[df$tcn[j]+1], border = NA)
  }
  w <- w + 1
}

## Mark chromosomes affected by chromothripsis
w <- 0
for (i in 1:nrow(ctx)){
  for (j in strsplit(ctx$chr_ctx[i], ",", fixed = T)[[1]]){
    if (!(j %in% strsplit(ctx$amp_ctx[i], ",", fixed = T)[[1]])){
      polygon(c(hg19$add_values[hg19$chr == j], hg19$add_values[hg19$chr == j] + hg19$length[hg19$chr == j], hg19$add_values[hg19$chr == j] + hg19$length[hg19$chr == j], hg19$add_values[hg19$chr == j]), c(w,w,w+1,w+1), lwd = 2/3, border = "black")
    }
  }
  w <- w + 1
}

## Mark chromosomes affected by chromothripsis and focal amplification
w <- 0
for (i in 1:nrow(ctx)){
  for (j in strsplit(ctx$chr_ctx[i], ",", fixed = T)[[1]]){
    if (j %in% strsplit(ctx$amp_ctx[i], ",", fixed = T)[[1]]){
      polygon(c(hg19$add_values[hg19$chr == j], hg19$add_values[hg19$chr == j] + hg19$length[hg19$chr == j], hg19$add_values[hg19$chr == j] + hg19$length[hg19$chr == j], hg19$add_values[hg19$chr == j]), c(w,w,w+1,w+1), lwd = 2/3, border = "#FF1493")
    } 
  }
  w <- w + 1
}
legend("topright", legend = c(0:11), fill=modcol)
dev.off()

# Figure 2d: chromothripsis amplifying oncogenes and inactivating tumor suppressor genes
pdf(paste0("C26.chr1.chromothripsis.pdf"), height = 4, width = 6)
plot(c(), xlim = c(1, hg19$length[hg19$chr == "1"]), ylim = c(0, 20), frame = F, xlab = "Position", ylab = "Copy number", main = "Chromosome 1 chromothripsis in C26", las = 1)

load("C26_hisens.Rdata")
df <- out$jointseg
df <- df[df$chrom == 1,]
df$tcn_estimate <- NA
f <- 0.7
p <- 1.8
for (i in 1:nrow(df)){
  df$tcn_estimate[i] <- ((f*p + (1-f)*2)*2^df$cnlr[i] - (1-f)*2)/f 
}
points(df$maploc, df$tcn_estimate, col = "#9FB6CD", pch=20, cex=0.1)

df <- cncf[cncf$study_id == "C26_hisens" & cncf$chrom == 1,]
df$lcn[is.na(df$lcn)] <- 1
segments(df$loc.start, df$tcn, df$loc.end, df$tcn, lwd = 1, col = "#DC143C")
segments(df$loc.start, df$lcn, df$loc.end, df$lcn, lwd = 1, col = "#27408B")
for (i in c("MYCL", "ARID1A")){
  abline(v = genedf$start[genedf$gene == i], col = "pink")
}

rect(cytoband$start[cytoband$chr == "1"], -2, cytoband$end[cytoband$chr == "1"], -1, col = cytoband$color[cytoband$chr == "1"], xpd = T, lwd = 0.5)
dev.off()

# Figure 2e: pathology images

# Figure 2f: pathology images


