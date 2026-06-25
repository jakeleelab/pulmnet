### Figure 1

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

### Figure 2

This script generates chromothripsis and copy-number visualization panels, including:

```text
C11.genome.chromothripsis.pdf
C11.chr2.chromothripsis.pdf
C11.chr11.chromothripsis.pdf
chromothripsis.cn.heatmap.pdf
C26.chr1.chromothripsis.pdf
```

Main input files and objects include:

```text
cncf.csv
chromothripsis.details.csv
C11_hisens.Rdata
C26_hisens.Rdata
```

The script generates genome-wide and chromosome-level copy-number plots for representative chromothripsis-positive cases. It uses FACETS-derived copy-number output to visualize raw copy-number estimates and segmented total/minor copy-number calls.

Briefly, `Figure2.R` performs the following steps:

1. Loads genome annotation objects, including chromosome size/position information from `hg19`.
2. Generates a genome-wide copy-number plot for case `C11`.
3. Loads FACETS-derived data from `C11_hisens.Rdata`.
4. Converts copy-number log-ratio values into estimated total copy number using sample-specific purity and ploidy values.
5. Overlays segmented total copy number and minor copy number from `cncf.csv`.
6. Generates chromosome-level inset plots for chromosomes 2 and 11 in case `C11`.
7. Loads chromothripsis annotation data from `chromothripsis.details.csv`.
8. Generates a copy-number heatmap across chromothripsis-positive tumors.
9. Marks chromosomes affected by chromothripsis and chromosomes with focal amplification.
10. Generates a chromosome 1 chromothripsis plot for case `C26`, highlighting genes of interest such as `MYCL` and `ARID1A`.
11. Saves all outputs as PDF files in the working directory.

The copy-number transformation used in this script is:

```r
tcn_estimate <- ((f * p + (1 - f) * 2) * 2^cnlr - (1 - f) * 2) / f
```

where:

```text
f    = tumor purity estimate
p    = tumor ploidy estimate
cnlr = copy-number log-ratio
```

Please note that the `*_hisens.Rdata` files and some sample-level FACETS output files are not included in the initial public repository because they may contain germline SNP-informed information. Therefore, some Figure 2 copy-number panels are not fully reproducible from the public repository alone.

### Figure 3

This script generates cross-site neuroendocrine tumor oncoprints and representative chromothripsis copy-number plots, including:

```text
oncoprint.pancreas.net.pdf
oncoprint.lung.net.pdf
oncoprint.gi.net.pdf
P01.chr11.chromothripsis.pdf
P01.chr12.chromothripsis.pdf
S01.chr11.chromothripsis.pdf
S01.chr12.chromothripsis.pdf
```

Main input files and objects include:

```text
pancreas.net.csv
lung.net.csv
gi.net.csv
cncf.csv
P01_hisens.Rdata
S01_hisens.Rdata
```

The script generates oncoprints for neuroendocrine tumors from major anatomical sites and visualizes chromothripsis-associated copy-number alterations in representative pancreatic and gastrointestinal neuroendocrine tumor cases.

Briefly, `Figure3.R` performs the following steps:

1. Defines a neuroendocrine tumor gene set used for cross-site oncoprint visualization.

```r
geneofi <- c(
  "EIF1AX", "ARID1A", "MEN1", "DAXX", "ATRX", "ATM", "SETD2",
  "KRAS", "PTEN", "TSC2", "CDKN1B", "CDKN2A", "CDKN2B",
  "CCND1", "MDM2", "CDK4", "MDM4", "MYCL", "RB1", "TP53"
)
```

2. Loads processed genomic alteration tables for neuroendocrine tumors from different anatomical sites:

```text
pancreas.net.csv
lung.net.csv
gi.net.csv
```

3. Generates site-specific oncoprints using `ComplexHeatmap`.

4. Harmonizes selected mutation annotations before plotting, for example:

```text
nonframeshift_deletion -> inframe
nonsynonymous_SNV      -> missense
```

5. Loads FACETS-derived copy-number data for representative chromothripsis-positive cases:

```text
P01_hisens.Rdata
S01_hisens.Rdata
```

6. Converts copy-number log-ratio values into estimated total copy number using sample-specific purity and ploidy values.

7. Generates chromosome-level chromothripsis plots for pancreatic case `P01`, including:

```text
P01.chr11.chromothripsis.pdf
P01.chr12.chromothripsis.pdf
```

8. Generates chromosome-level chromothripsis plots for small intestinal/gastrointestinal case `S01`, including:

```text
S01.chr11.chromothripsis.pdf
S01.chr12.chromothripsis.pdf
```

9. Overlays segmented total copy number and minor copy number from `cncf.csv`.

10. Marks genes of interest on copy-number plots, including:

```text
CCND1
CDK4
MDM2
PAK1
```

11. Saves all outputs as PDF files in the working directory.

The copy-number transformation used in this script is:

```r
tcn_estimate <- ((f * p + (1 - f) * 2) * 2^cnlr - (1 - f) * 2) / f
```

where:

```text
f    = tumor purity estimate
p    = tumor ploidy estimate
cnlr = copy-number log-ratio
```

Please note that the `*_hisens.Rdata` files and some sample-level FACETS output files are not included in the initial public repository because they may contain germline SNP-informed information. Therefore, the Figure 3 chromothripsis copy-number panels are not fully reproducible from the public repository alone.




