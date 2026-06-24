## Requirements

Analyses were performed in **R 4.3.3**. The main R packages required to reproduce the analyses are listed below.

### Required R packages

The following packages were attached during the analysis:

| Package | Version |
|---|---:|
| Rediscover | 0.3.2 |
| matrixStats | 1.5.0 |
| ShiftConvolvePoibin | 1.0.0 |
| PoissonBinomial | 1.2.7 |
| Matrix | 1.6-5 |
| GenomicRanges | 1.54.1 |
| GenomeInfoDb | 1.38.8 |
| IRanges | 2.36.0 |
| S4Vectors | 0.40.2 |
| BiocGenerics | 0.48.1 |
| viridis | 0.6.5 |
| MetBrewer | 0.2.0 |
| survminer | 0.5.0 |
| ggpubr | 0.6.1 |
| ggplot2 | 3.5.2 |
| survival | 3.8-3 |
| circlize | 0.4.16 |
| fields | 16.3.1 |
| viridisLite | 0.4.2 |
| spam | 2.11-1 |
| ComplexHeatmap | 2.16.0 |
| stringr | 1.5.1 |
| dplyr | 1.1.4 |

These can be installed using:

```r
install.packages(c(
  "Rediscover",
  "matrixStats",
  "ShiftConvolvePoibin",
  "PoissonBinomial",
  "Matrix",
  "viridis",
  "MetBrewer",
  "survminer",
  "ggpubr",
  "ggplot2",
  "survival",
  "circlize",
  "fields",
  "viridisLite",
  "spam",
  "stringr",
  "dplyr"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "BiocGenerics",
  "S4Vectors",
  "IRanges",
  "GenomeInfoDb",
  "GenomicRanges",
  "ComplexHeatmap"
))

## Session information

The figure scripts were run using the following R environment:

```r
sessionInfo()
```

```text
R version 4.3.3 (2024-02-29)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS 26.5.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib; LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] Rediscover_0.3.2          matrixStats_1.5.0         ShiftConvolvePoibin_1.0.0
[4] PoissonBinomial_1.2.7     Matrix_1.6-5              GenomicRanges_1.54.1     
[7] GenomeInfoDb_1.38.8       IRanges_2.36.0            S4Vectors_0.40.2         
[10] BiocGenerics_0.48.1       viridis_0.6.5             MetBrewer_0.2.0          
[13] survminer_0.5.0           ggpubr_0.6.1              ggplot2_3.5.2            
[16] survival_3.8-3            circlize_0.4.16           fields_16.3.1            
[19] viridisLite_0.4.2         spam_2.11-1               ComplexHeatmap_2.16.0    
[22] stringr_1.5.1             dplyr_1.1.4              

loaded via a namespace (and not attached):
[1] tidyselect_1.2.1        farver_2.1.2            bitops_1.0-9           
[4] RCurl_1.98-1.17         digest_0.6.37           dotCall64_1.2          
[7] lifecycle_1.0.4         cluster_2.1.8.1         magrittr_2.0.3         
[10] compiler_4.3.3          rlang_1.1.6             tools_4.3.3            
[13] data.table_1.17.8       knitr_1.50              ggsignif_0.6.4         
[16] RColorBrewer_1.1-3      abind_1.4-8             withr_3.0.2            
[19] purrr_1.0.4             xtable_1.8-4            colorspace_2.1-1       
[22] scales_1.4.0            iterators_1.0.14        cli_3.6.5              
[25] crayon_1.5.3            generics_0.1.4          rstudioapi_0.17.1      
[28] km.ci_0.5-6             rjson_0.2.23            DNAcopy_1.76.0         
[31] zlibbioc_1.48.2         splines_4.3.3           maps_3.4.3             
[34] parallel_4.3.3          XVector_0.42.0          survMisc_0.5.6         
[37] vctrs_0.6.5             carData_3.0-5           car_3.1-3              
[40] GetoptLong_1.0.5        rstatix_0.7.2           Formula_1.2-5          
[43] clue_0.3-66             foreach_1.5.2           tidyr_1.3.1            
[46] glue_1.8.0              codetools_0.2-20        stringi_1.8.7          
[49] shape_1.4.6.1           gtable_0.3.6            tibble_3.3.0           
[52] pillar_1.11.0           GenomeInfoDbData_1.2.11 R6_2.6.1               
[55] KMsurv_0.1-6            doParallel_1.0.17       evaluate_1.0.4         
[58] lattice_0.22-7          png_0.1-8               backports_1.5.0        
[61] broom_1.0.8             Rcpp_1.1.0              gridExtra_2.3          
[64] xfun_0.52               maftools_2.18.0         zoo_1.8-14             
[67] pkgconfig_2.0.3         GlobalOptions_0.1.2    
```
