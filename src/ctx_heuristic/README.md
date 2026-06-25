# Chromothripsis Heuristic Screening Code

This repository contains `ctx_heuristic.R`, an R script used to perform an initial heuristic screen for potential chromothripsis cases from FACETS-derived copy-number segmentation data.

## Purpose

The script identifies sample/chromosome pairs with copy-number patterns suggestive of chromothripsis. Specifically, it searches for chromosomes showing:

1. Multiple copy-number segments
2. Multiple total copy-number states
3. Repeated oscillation between a baseline copy-number state and higher copy-number states
4. Sufficient segment complexity to warrant manual review

This script was used as a screening tool to nominate candidate chromothripsis cases for downstream manual review.

## Important note

This code is **not intended to be a definitive chromothripsis caller**. Candidate cases identified by this heuristic should be reviewed manually using copy-number plots, allele-specific copy-number profiles, pathology context, and other available genomic/clinical information.

## Required input files and objects

The script requires FACETS-derived copy-number segment data and sample annotation tables.

### Main input file

```text
impact_facets_annotated.cncf.txt.gz
```

This file should contain FACETS copy-number segment information, including columns such as:

```text
ID
chrom
loc.start
loc.end
tcn
lcn
```

where:

- `tcn` = total copy number
- `lcn` = lesser/minor copy number

## Workflow

Briefly, the script performs the following steps:

1. Loads FACETS copy-number segment data.
2. Extracts sample identifiers from FACETS IDs.
3. Restricts the analysis to study samples.
4. Iterates through each pulmonary neuroendocrine tumor sample.
5. Retrieves FACETS copy-number segments either from the combined copy-number file or from sample-specific FACETS output files.
6. Identifies chromosomes with at least 5 copy-number segments.
7. Merges adjacent segments with identical total and minor copy-number states.
8. Calculates chromosome-level metrics, including:
   - number of copy-number segments
   - number of total copy-number states
   - inferred baseline copy-number state
   - number of oscillating baseline-copy-number blocks
   - fraction of segment length assigned to baseline/low-copy states
9. Flags chromosomes meeting heuristic thresholds for manual review.
10. Exports per-chromosome and per-case candidate tables.

## Heuristic criteria

The main relaxed threshold used for manual review is:

```r
n_oscil >= 3 & n_seg >= 5
```

where:

- `n_oscil` = number of inferred baseline-copy-number oscillation blocks
- `n_seg` = number of merged copy-number segments on a chromosome

A stricter threshold used in exploratory summaries is:

```r
n_oscil >= 3 & n_seg >= 7
```

## Output files

The script writes two tab-delimited output files:

```text
per_chromosome_summary.tsv
per_case_annotation.tsv
```

### `per_chromosome_summary.tsv`

This file contains chromosome-level candidate calls and summary metrics.

Main columns include:

```text
dmp_sample
chrom
n_seg
n_tcn_state
n_oscil
basecn
n_seg_basecn
length_seg_basecn
length_seg_total
```

### `per_case_annotation.tsv`

This file contains case-level annotations for samples with at least one chromosome meeting the relaxed heuristic threshold.
