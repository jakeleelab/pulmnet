## Requirements

Analyses were performed in **R 4.3.3** on macOS. The main R packages required to reproduce the analyses are listed below.

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
