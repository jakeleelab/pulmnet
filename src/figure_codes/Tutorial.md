### Figure 1

Run the Figure 1 script with:

```r
source("Figure1.R")
```

This script generates the pulmonary carcinoid genomic landscape figures, including:

```text
oncoprint.pulmnet.pdf
oncoprint.TSGs.pdf
mutual.exclusivity.heatmap.pdf
```

Main input files include:

```text
pulmnet.csv
segmental.cn.csv
pulmnet.annotated.csv
```

The script generates oncoprints using `ComplexHeatmap` and performs mutual exclusivity/co-occurrence analysis using `Rediscover`.

Briefly, `Figure1.R` performs the following steps:

1. Loads genomic and mutational data from `pulmnet.csv`.
2. Defines the core pulmonary carcinoid gene set used for plotting.
3. Generates the main genomic landscape oncoprint.
4. Generates a tumor suppressor gene-focused oncoprint incorporating LOH information from `segmental.cn.csv`.
5. Loads annotated patient-level data from `pulmnet.annotated.csv`.
6. Performs mutual exclusivity and co-occurrence analysis across selected genomic features.
7. Saves all outputs as PDF files in the working directory.
