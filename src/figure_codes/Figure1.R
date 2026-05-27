library(stringr)
library(ComplexHeatmap)
library(Rediscover)

# Codes for Figure 1

## IMPACT plotting
impact <- read.csv("pulmnet.csv") # File containing genomic and mutational information

## Defining gene sets
putative_drivers <- c("EIF1AX", "ARID1A", "ARID1B", "ARID2", "PBRM1", "MEN1", "KMT2A", "KMT2D", "KDM5C", "ATM", "BAP1", "KRAS", "NF1", "PTEN", "AKT3", "TSC1", "CCND1", "CDK4", "MDM2", "MDM4", "MYCL", "CDKN2A", "RB1", "TP53", "SMARCA4", "SMARCB1", "TERT", "MGA", "SF3B1", "U2AF1", "CDKN1B", "CDKN2B")
geneofi <- c("EIF1AX", "ARID1A", "ARID1B", "ARID2", "PBRM1", "MEN1", "KMT2A", "KMT2D", "KDM5C", "ATM", "BAP1", "KRAS", "NF1", "PTEN", "AKT3", "TSC1", "CCND1", "CDK4", "MDM2", "MDM4", "MYCL", "CDKN2A", "RB1", "TP53") # core gene set for plotting
otherpulm <- c("SMARCA4", "SMARCB1", "TERT", "MGA", "SF3B1", "U2AF1", "CDKN1B", "CDKN2B")

# Figure 1a: Oncoprint
intermed <- impact[,colnames(impact)[!(colnames(impact) %in% otherpulm)]]
list1 <- c(which(colnames(intermed) == geneofi[length(geneofi)]):which(colnames(intermed) == geneofi[1]))
for (i in list1){
  intermed <- intermed[order(intermed[,i], decreasing = T),]
}
intermed <- intermed[order(intermed$chromothripsis, decreasing = T),]

stage_id_levels <- names(table(intermed$stage_simple))
stage_id_colors <- setNames(c(met.brewer(name="Tam", n=4)[c(1:4)]), stage_id_levels)

mutix <- intermed[,c(which(colnames(intermed) == "sex"),which(colnames(intermed) == "smoking"),which(colnames(intermed) == "class_sample"),which(colnames(intermed) == "histology"),which(colnames(intermed) == "necrosis"),which(colnames(intermed) == "ihc_group"),which(colnames(intermed) == "bin_otp"),which(colnames(intermed) == "bin_ascl1"),which(colnames(intermed) == "bin_hnf4a"),which(colnames(intermed) == "WGD"),which(colnames(intermed) == "chromothripsis"),c(which(colnames(intermed)==geneofi[1]):which(colnames(intermed)==geneofi[length(geneofi)])))]
mutix$necrosis[mutix$necrosis != "present"] <- ""
mutix <- as.matrix(mutix)
mutix <- t(mutix)
mutix <- gsub("nonframeshift_deletion", "inframe", mutix)
mutix <- gsub("nonsynonymous_SNV", "missense", mutix)
mutix <- gsub("02_atypical_carcinoid", "atypical_carcinoid", mutix)
mutix <- gsub("01_typical_carcinoid", "typical_carcinoid", mutix)
mutix <- gsub("03_no_resection", "", mutix)

col = c(missense = "#006400", nonsense = "#000000", splice = "#FF34B3", inframe = "#CD0000", frameshift = "#0000CD", germline = "#FF9912", LOH = "#BABABA", Amp = "#EE2C2C", Del = "#00688B", sv = "#9ACD32", upstrea_reg = "#00BFFF", nonstop = "#800080", Unknown = "#4169E1", wgd = "#000000", chromothripsis = "#800080", atypical_carcinoid = "#FF9912", typical_carcinoid = "#3D59AB", msi = "#6E8B3D", present = "#191970", focalamp = "#EE2C2C", never = "#B0E2FF", former = "#33A1C9", current = "#191970", Male = "#104E8B", Female = "#CD6889", tx_naive = "#F5DEB3", hormone_only = "#8B7E66", chemo_xrt_rrt = "#292421", A1 = "#FFF68F", A2 = "#9ACD32", B = "#1E90FF", not_available = "#E0E0E0", other = "#708090", positive = "#8B4513", negative = "#C6E2FF")

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
                             germline = function(x, y, w, h) grid.rect(x, y, w*0.7, h*0.35, gp = gpar(fill = col["germline"], col = NA)),
                             wgd = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.7, gp = gpar(fill = col["wgd"], col = NA)),
                             msi = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.7, gp = gpar(fill = col["msi"], col = NA)),
                             present = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["present"], col = NA)),
                             focalamp = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.7, gp = gpar(fill = col["focalamp"], col = NA)),
                             never = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["never"], col = NA)),
                             former = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["former"], col = NA)),
                             current = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["current"], col = NA)),
                             Male = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Male"], col = NA)),
                             Female = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Female"], col = NA)),
                             atypical_carcinoid = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["atypical_carcinoid"], col = NA)),
                             typical_carcinoid = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["typical_carcinoid"], col = NA)),
                             chromothripsis = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.7, gp = gpar(fill = col["chromothripsis"], col = NA)),
                             Unknown = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Unknown"], col = NA)),
                             A1 = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["A1"], col = NA)),
                             A2 = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["A2"], col = NA)),
                             B = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["B"], col = NA)),
                             other = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["other"], col = NA)),
                             not_available = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["not_available"], col = NA)),
                             positive = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["positive"], col = NA)),
                             negative = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["negative"], col = NA)),
                             tx_naive = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["tx_naive"], col = NA)),
                             hormone_only = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["hormone_only"], col = NA)),
                             chemo_xrt_rrt = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["chemo_xrt_rrt"], col = NA))
                           ), 
                           top_annotation = HeatmapAnnotation(
                             column_barplot = anno_oncoprint_barplot(c("sv", "Amp", "Del", "missense", "nonsense", "splice", "inframe", "frameshift", "upstrea_reg", "nonstop", "germline"), show_fraction = T),
                             age_dx = as.matrix(intermed[,"age_dx"]),
                             smoking_py = as.matrix(intermed[,"py"]),
                             stage = as.matrix(intermed[,"stage_simple"]),
                             mitosis = as.matrix(intermed[,"max_mitosis"]),
                             ki67 = as.matrix(intermed[,"max_ki67"]),
                             tmb = as.matrix(intermed[,"tmb"]),
                             fga = as.matrix(intermed[,"facets_fga"]),
                             ploidy = as.matrix(intermed[,"facets_ploidy"]),
                             purity = as.matrix(intermed[,"facets_purity"]),
                             col = list(age_dx = colorRamp2(c(15,95), c("#000000", "#FFFFFF")), smoking_py = colorRamp2(c(0,70), c("#FFFFFF", "#000080")), stage = stage_id_colors, mitosis = colorRamp2(c(0,50), c("#FFFFFF", "#B0171F")), tmb = colorRamp2(c(0,7.5), c("#FFE1FF", "#4B0082")), ki67 = colorRamp2(c(0,3,15,60), c("#FDE725", "#55C667", "#404788", "#440154")), fga = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), viridis_pal(option = "F", direction = -1)(5)), ploidy = colorRamp2(c(1.5, 2, 4, 5), c("#0000FF", "#FFFFFF", "#CD2626", "#872657")), purity = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), viridis_pal(option = "G", direction = -1)(5)))
                           ),
                           right_annotation = rowAnnotation(
                             row_barplot = anno_oncoprint_barplot(c("sv", "Amp", "Del", "missense", "nonsense", "splice", "inframe", "frameshift", "upstrea_reg", "nonstop", "germline", "tx_naive", "hormone_only", "chemo_xrt_rrt", "A1", "A2", "B", "other"), show_fraction = T)
                           ),
                           col = col, row_order = 1:nrow(mutix), column_order = 1:ncol(mutix), remove_empty_columns = FALSE, column_title = paste0("Pulmonary Carcinoid IMPACT (# patient = ", nrow(impact), "; # tumor = ", nrow(net), ")")
                           , heatmap_legend_param = list(title = "Alterations", nrow = 2, title_position = "leftcenter"))

pdf(paste0("oncoprint.pulmnet.pdf"), height=10, width=14)
draw(test.oncoprint, heatmap_legend_side = "bottom")
dev.off()

# Figure 1b: Lollipop plots made in cBioportal

# Figure 1c: Biallelic oncoprint for tumor suppressor genes
loh <- impact
fg <- read.csv("segmental.cn.csv")

list1 <- rev(which(colnames(loh) %in% c("ATM", "MEN1", "ARID1A", "NF1")))
for (i in list1){
  loh <- loh[order(loh[,i], decreasing = T),]
}

for (i in 1:nrow(loh)){
  for (j in c("ARID1A", "MEN1", "ATM", "NF1")){
    if (nrow(fg[startsWith(fg$sample, loh$dmp_sample[i]),]) != 0 && !is.na(fg$lcn.em[startsWith(fg$sample, loh$dmp_sample[i]) & fg$gene == j][1]) && fg$lcn.em[startsWith(fg$sample, loh$dmp_sample[i]) & fg$gene == j][1] == "0"){
      loh[i,j] <- ifelse(loh[i,j] == "", "LOH", paste0(loh[i,j], ";LOH"))
    }
  }
}

mutix <- loh[,which(colnames(loh) %in% c("ATM", "MEN1", "ARID1A", "NF1"))]
mutix <- as.matrix(mutix)
mutix <- t(mutix)
mutix <- gsub("nonframeshift_deletion", "inframe", mutix)
mutix <- gsub("nonsynonymous_SNV", "missense", mutix)

col = c(missense = "#006400", nonsense = "#000000", splice = "#FF34B3", inframe = "#CD0000", frameshift = "#0000CD", germline = "#FF9912", LOH = "#BABABA", Amp = "#EE2C2C", Del = "#00688B", sv = "#9ACD32", upstrea_reg = "#00BFFF", nonstop = "#800080")

test.oncoprint = oncoPrint(mutix, get_type = function(x) strsplit(x, ";")[[1]],
                           alter_fun = list(
                             background = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#E0E0E0", col = NA)),
                             LOH = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["LOH"], col = NA)),
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
                             germline = function(x, y, w, h) grid.rect(x, y, w*0.7, h*0.35, gp = gpar(fill = col["germline"], col = NA))
                           ), 
                           right_annotation = rowAnnotation(
                             row_barplot = anno_oncoprint_barplot(c("sv", "Amp", "Del", "missense", "nonsense", "splice", "inframe", "frameshift", "upstrea_reg", "nonstop", "germline"), show_fraction = T)
                           ),
                           col = col, row_order = 1:nrow(mutix), column_order = 1:ncol(mutix), remove_empty_columns = FALSE, column_title = paste0("Pulmonary Carcinoid IMPACT (# patient = ", nrow(loh), "; # tumor = ", nrow(net), ")")
                           , heatmap_legend_param = list(title = "Alterations", nrow = 2, title_position = "leftcenter"))

pdf(paste0("oncoprint.TSGs.pdf"), height=3, width=10)
draw(test.oncoprint, heatmap_legend_side = "bottom")
dev.off()

# Figure 1d: Mutual exclusivity and co-occurrence plot
ca <- read.csv("pulmnet.annotated.csv")
df <- ca[,c("dmp_pt", "chromothripsis", "dic_EIF1AX", "dic_MEN1", "dic_ARID1A", "dic_ATM", "dic_KRAS", "dic_NF1", "dic_DN", "focalamp")]
for (i in 2:ncol(df)){
  df[,i] <- grepl("reference", df[,i])
}
df[-1] <- lapply(df[-1], function(x){ifelse(x == TRUE, 0, ifelse(x == FALSE, 1, x))})
rownames(df) <- df[,1]
df <- df[,colnames(df)[colnames(df) != "dmp_pt"]]
mat <- t(as.matrix(df))
pma <- getPM(mat)
getMutex(A=mat, PM=pma, lower.tail = T) # Mutual exclusivity
getMutex(A=mat, PM=pma, lower.tail = F) # Co-occurrence

tempdf <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(tempdf) <- c("gene1", "gene2", "fisher_p", "fisher_or", "c00", "c11", "c10", "c01")
for (i in c("chromothripsis", "dic_EIF1AX", "dic_MEN1", "dic_ARID1A", "dic_ATM", "dic_KRAS", "dic_NF1", "dic_DN", "focalamp")){
  for (j in c("chromothripsis", "dic_EIF1AX", "dic_MEN1", "dic_ARID1A", "dic_ATM", "dic_KRAS", "dic_NF1", "dic_DN", "focalamp")){
    tempdf[nrow(tempdf)+1,] <- NA
    tempdf$gene1[nrow(tempdf)] <- i
    tempdf$gene2[nrow(tempdf)] <- j
    tempdf$fisher_p[nrow(tempdf)] <- fisher.test(table(df[,i], df[,j]))$p
    tempdf$fisher_or[nrow(tempdf)] <- unname(fisher.test(table(df[,i], df[,j]))$estimate)
    tempdf$c00[nrow(tempdf)] <- table(df[,i], df[,j])[1,1]
    tempdf$c11[nrow(tempdf)] <- table(df[,i], df[,j])[2,2]
    tempdf$c10[nrow(tempdf)] <- table(df[,i], df[,j])[2,1]
    tempdf$c01[nrow(tempdf)] <- table(df[,i], df[,j])[1,2]
  }
}
tempdf$bh_q <- p.adjust(tempdf$fisher_p, method = "BH")
tempdf[tempdf$bh_q < 0.05 & tempdf$gene1 != tempdf$gene2,]

varlist <- c("chromothripsis", "dic_EIF1AX", "dic_MEN1", "dic_ARID1A", "dic_ATM", "dic_KRAS", "dic_NF1", "dic_DN", "focalamp")
mat <- matrix(ncol = 9, nrow = 9)
for (k in 1:nrow(tempdf)){
  mat[which(tempdf$gene2[k] == varlist),10-which(tempdf$gene1[k] == varlist)] <- ifelse(tempdf$gene2[k] == tempdf$gene1[k], 0, ifelse(tempdf$fisher_or[k] >= 1, -log10(tempdf$bh_q[k]), log10(tempdf$bh_q[k])))  
}
mat[mat>10] <- 10 # ceiling value of 10.

colnames(mat) <- rev(varlist)
rownames(mat) <- varlist

sigmat <- matrix(ncol = 9, nrow = 9)
for (k in 1:nrow(tempdf)){
  sigmat[which(tempdf$gene2[k] == varlist),10-which(tempdf$gene1[k] == varlist)] <- ifelse(tempdf$gene2[k] == tempdf$gene1[k], 0, ifelse(tempdf$bh_q[k] < 0.05, 1, 0))  
}
colnames(sigmat) <- rev(varlist)
rownames(sigmat) <- varlist

pdf("mutual.exclusivity.heatmap.pdf", height=4, width=5)
Heatmap(mat, row_order = 1:nrow(mat), column_order = 1:nrow(mat), column_title = "Mutual exclusivity and co-occurrence", col = col_fun, name = "-log10(Q-value)", rect_gp = gpar(col = "lightgray", lwd=4/3), row_names_side = "left",
        cell_fun = function(j, i, x, y, width, height, fill){
          if(mat[i,j] > -log10(0.05) || mat[i,j] < log10(0.05))
            grid.text("*", x = x, y = y, gp = gpar(fontsize = 10))
        })
dev.off()
