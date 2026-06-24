# ctx_finder

Chromothripsis-footprint detection on sparse FACETS (targeted-panel) copy-number data.

## Summary

`ctx_finder` scans each chromosome of a tumor sample and nominates the ones that carry a chromothripsis footprint. It works on the coarse, sparse copy-number profiles produced by FACETS on targeted gene panels, where classic WGS-based callers (ShatterSeek and similar) may not apply.

For every (sample, chromosome) it builds copy-number oscillation features (the longest contiguous chain confined to a small number of CN states, allelic LOH-versus-het interleaving, an altered-footprint extent, and a probe-level changepoint confidence. A fixed-threshold classifier combines these into a per-chromosome call, summarized as a sample-level flag (`sample_flag`).

Time-consuming per-sample tasks that compute the changepoints is ran as a Slurm array. This can be adapted to the cluster-setting as needed (or ran sequentially). 

```
ctx_finder/
├── run_ctx_pipeline.py    # manifest + cohort feature builder   (ADD THIS)
├── ctx_changepoint.py     # per-chromosome changepoint features  (ADD THIS)
├── ctx_features.py        # chromosome features + v6 classifier + evaluation
├── ctx_cp_parallel.py     # Slurm-array launcher / worker / collector
├── ctx_label.py           # per-chromosome inspection plots
└── ctx_plots.py           # cohort-level diagnostic plots
```

## Installation (conda)

```bash
conda create -n ctx python=3.10 numpy pandas matplotlib seaborn -c conda-forge
conda activate ctx
pip install ruptures        # changepoint detection (ctx_changepoint falls back to a slow DP solver if missing)
```

The feature builder regenerates probe-level `jointseg` tables from each sample's FACETS `.Rdata`, so **R with FACETS installed and `Rscript` on `PATH`** is required for that step. If your manifest already points at materialized jointseg tables, R is not needed.

```bash
Rscript -e 'library(facets); packageVersion("facets")'   # sanity check
```

## Input data format

You supply one **manifest TSV**. Each row is one sample and must have these columns:

| column            | meaning                                                                 |
|-------------------|-------------------------------------------------------------------------|
| `SAMPLE_ID`       | unique sample identifier                                                |
| `facets_cncf_path`| path to that sample's FACETS output (`.cncf` / `.Rdata`)                |
| `purity`          | FACETS purity estimate                                                  |
| `ploidy`          | FACETS ploidy estimate                                                  |
| `chr_ctx`         | ground-truth positive chromosomes, used only for evaluation (see below) |

`chr_ctx` is the truth label the cohort builder turns into the per-chromosome column `truth_positive_chrom`. If you are running on new, unlabeled data, the column must still exist but can be left empty: the `v6` classifier uses fixed thresholds and does not need truth to make calls. Truth is only consumed by `evaluate_classifier` (to score performance).

Behind each `facets_cncf_path`, the pipeline reads two FACETS tables per sample:

- **Segment table (`seg_df`, cncf):** `chrom`, `loc.start`, `loc.end`, `num.mark`, `cnlr.median`, `mafR`, `tcn.em`, `lcn.em`, `seg`
- **Probe-level table (`jointseg_df`):** `chrom`, `maploc`, `cnlr`, `lorvar`, `vafT`, `het`, `seg` (the per-probe `tcn_estimate` is added internally)

Coordinates are hg19 in the reference build used here. The plotting helper reads centromeres from an hg19 coordinates CSV; edit the default path in `ctx_label.load_centromeres` if you plot.

## Running guide

The run is two phases because the changepoint pass goes through Slurm. Pick one base output directory; everything is written under it.

### Phase 1: cohort features + launch the changepoint array

```python
import os
import pandas as pd
import run_ctx_pipeline as rp
import ctx_features    as cf
import ctx_cp_parallel as cpx

OUT    = "/path/to/results/my_run"
TABLES = os.path.join(OUT, "tables")    # jointseg cache + feature tables
CP_OUT = os.path.join(OUT, "cp_run")    # changepoint outputs
for d in (TABLES, CP_OUT):
    os.makedirs(d, exist_ok=True)

MANIFEST = "/path/to/manifest.tsv"
man_df = pd.read_csv(MANIFEST, sep="\t")

# 1. which samples are processable (have a readable FACETS path, etc.)
proc, ledger = rp.build_processable_sample_ids(man_df, TABLES)

# 2. build per-chromosome cohort features (regenerates jointseg from Rdata; needs Rscript)
cohort, failed = rp.build_cohort_features(man_df, TABLES, proc, use_cache=True)
cohort.to_csv(os.path.join(TABLES, "cohort_chrom_features.tsv"), sep="\t", index=False)

# 3. write the Slurm array script for the changepoint pass (penalty = 8)
cpx.launch(out_dir=CP_OUT, pens=[8], chunk=4, man_path=MANIFEST, tablesDir=TABLES)
```

Then submit the array job:

```bash
sbatch /path/to/results/my_run/cp_run/run_array.sbatch
squeue -u $USER          # wait until done; logs in cp_run/logs/task_*.out
ls /path/to/results/my_run/cp_run/cp_features/ | wc -l   # ~one .tsv per sample, no .FAILED markers
```


### Phase 2: collect changepoint output, classify, evaluate

```python
cohort = pd.read_csv(os.path.join(TABLES, "cohort_chrom_features.tsv"), sep="\t")

# 4. gather per-sample changepoint tsvs, keep penalty = 8
cpx.collect(CP_OUT, pen_const=8.0)                       # writes cp_run/cohort_cp_features.tsv
cp8 = pd.read_csv(os.path.join(CP_OUT, "cohort_cp_features.tsv"), sep="\t")
cp8 = cp8[cp8["pen_const"] == 8.0].copy()

# 5. merge changepoint features onto the cohort (SAMPLE_ID + normalized chrom)
cohort["_c"] = cohort["chrom"].map(cf._norm_chrom)
cp8["_c"]    = cp8["chrom"].map(cf._norm_chrom)
m = cohort.merge(cp8.drop(columns=["chrom"]), on=["SAMPLE_ID", "_c"], how="left")

# 6. apply the classifier
purity_map = dict(zip(man_df["SAMPLE_ID"].astype(str),
                      pd.to_numeric(man_df["purity"], errors="coerce")))
v6 = cf.apply_classifier_v6(m, purity_map)
v6.to_csv(os.path.join(TABLES, "cohort_chrom_features_v6.tsv"), sep="\t", index=False)

# 7. (optional, needs chr_ctx truth) score performance
cf.evaluate_classifier(v6, nominated_col="nominated_v6", label="v6")
```

### Output

`cohort_chrom_features_v6.tsv` has one row per (sample, chromosome). Important columns include:

- `nominated_v6`: the per-chromosome chromothripsis call (boolean)
- `sample_flag_v6`: True if any chromosome in that sample is nominated
- supporting columns: `ctx_suppress` (suppressed mechanistic score), `longest_osc_chain_k3`, `cp_osc_confidence`, `dot_mad_within` (probe noise), `altered_union_frac`, `n_altered_islands`.

