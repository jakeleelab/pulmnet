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
├── ctx_plots.py           # cohort-level diagnostic plots
└── ctx_sibling.py         # sibling-sample corroboration check
```

## Sibling corroboration (`ctx_sibling.py`)

For patients with multiple samples, a chromosome can be nominated in one sample
but missed in another (e.g. a noisier library pushes it below threshold). Real,
shared structural events reproduce FACETS segment breakpoints at (near-)identical
positions across independent samples from the same tumor -- probe positions are
fixed by the panel, so this is a strong, checkable signature -- while unrelated
samples show no such relationship.

`ctx_sibling.sibling_corroboration_report(cohort_df, man_df)` finds every
discordant (patient, chrom) pair in a classifier output table and, for each,
computes what fraction of the flagged sample's breakpoints have a match (within
`tol`, default 200kb) in the "missed" sibling, then ranks that match fraction
against an empirical null built from every other sample in the cohort on that
same chromosome. A high percentile with a nontrivial absolute match fraction is
quantitative support for a real, shared-but-missed event; a percentile near 0
means the "missed" sample most likely genuinely lacks the event (a real
negative, not a detection failure). Validated against carcinoid sibling pairs
already resolved by manual review -- known-corroborating pairs scored 100th
percentile, a known-non-matching pair scored 0th, and a pair the domain expert
called only "supporting evidence" (not fully convincing) scored an intermediate
84.7th percentile with thin absolute overlap, tracking the expert's own
confidence.

Caution: with very few breakpoints on a chromosome, a single coincidental match
can produce a moderate-looking percentile purely because the null distribution
is sparse -- always read `percentile` together with `match_frac` and
`n_bp_flagged`, not `percentile` alone.

```python
import ctx_sibling as csib
report = csib.sibling_corroboration_report(v7, man_df)   # v7 = apply_classifier_v7 output
```

## Thin-chain reliability filter (`ctx_features.apply_thin_chain_filter`)

Manual review of classical-SCLC false positives (2026-08) found a specific
failure mode: the oscillating tcn.em labels driving a `nominated_v7` call were
sometimes carried by a single thin FACETS segment with an unreliable
cellular-fraction estimate (`cf.em`), sandwiched between two ordinary segments
-- not a real recurring structural event. Segment genomic span is unreliable
evidence of this on panel data (a segment can span megabases with almost no
probes in it), and even het-SNP count (`nhet`) alone doesn't discriminate --
real chromothripsis, including the loss-dominant/shard pathway, legitimately
produces small, low-nhet segments too. `cf.em` on the specific segment(s)
forming the chain does discriminate.

**This is v2 (2026-08-22), which supersedes v1.** v1 found the actual
oscillation-chain window (`longest_osc_chain_k3`'s segment run) and vetoed the
**entire chromosome** if any "spike" segment in that window (rounded tcn.em
differs from both neighbors) had `cf.em` below threshold (0.15) -- correct
when the chromosome's only evidence *is* that weak segment, but wrong when a
real, well-supported chain survives around it. GI NET's C016-T02-IM7 chr3
(`2,1,4,1,4,1,2,1,2`) is real chromothripsis per manual review despite one bad
`tcn=4` pair, and v1 wholesale-vetoes it anyway. v2 fixes this: treat the
unreliable spike as a **break** in the segment sequence (not a whole-chromosome
veto) -- excise it, recompute the longest surviving `<=k`-distinct-state chain
on each side, and only drop the chromosome if the **best surviving chain** is
below `min_chain` (default 4, matching `apply_classifier_v7`'s own
`min_osc_chain`). This keeps GI NET chr3 (surviving chain lands at exactly 4)
while still dropping chromosomes whose only evidence is the weak segment
itself (C002-T01-IM5 chr11: surviving chain drops to 3, correctly cut).

Validated on the **carcinoid calibration cohort** against corrected ground
truth (two chromosomes originally believed weak turned out, on closer
domain-expert review with cross-sibling corroboration, to be real): of 73
conventional-pathway nominated rows, v2 at `min_chain=4` drops only
C002-T01-IM5 chr11 (solid confirmed-fake) and correctly restores
C017-T01-IM7 chr11 and C020-T01-IM7 chr16 (both now confirmed real).
A third v1-vetoed row, C017-T02-IM7 chr7, remains genuinely
borderline/unresolved either way. Chromosome-level: carcinoid 89 (raw v7) ->
85 (v1) -> 88 (v2). **Sample-level is identical between v1 and v2 (53/238
either way)** -- v2 is a strict chromosome-level accuracy improvement with no
sample-level cost on the calibration cohort. Also re-run (2026-08-22) on the
other 4 cohorts: classical SCLC 227/679 -> 226/679, pancreatic NET/NEC
81/506 -> 80/506, GI NET 37/245 -> 37/245 (no change), aSCLC 14/25 -> 14/25 (no
change) -- v2's chromosome-level effect there is much smaller than v1's was,
since v2 only removes chromosomes whose *entire* surviving evidence is below
`min_chain`, not merely whose weakest single segment is unreliable.

The loss-dominant/shard pathway (`loss_dominant_shard_hit`) is deliberately
excluded from this filter -- it has its own noise/deficit logic already, and
was empirically untouched by v1 at every threshold tested (carcinoid's own
shard hits sit at `cf.em` >= 0.117), consistent with those calls legitimately
producing small, thin segments this metric would otherwise misread as
unreliable.

```python
import ctx_features as cf
filtered = cf.apply_thin_chain_filter(v7, man_df, tablesDir)   # v7 = apply_classifier_v7 output
# use filtered["nominated_v7_filtered"] / filtered["sample_flag_v7_filtered"]
```

Or via the shared driver script:

```bash
python ctx_run_common/phase5_thin_chain_filter.py --manifest <path> --out <cohort_dir> [--v7-file ...] [--min-chain 4]
```

## Loss-segment filter (`ctx_features.apply_loss_segment_filter`)

Second post-hoc filter, applied **after** the thin-chain filter above.
Manual review of 15 randomly-selected classical-SCLC/pancreatic-NET/GI-NET
chromosomes that survive v7 + the thin-chain filter but carry only 2 "loss
segments" (see below) found them well-supported (no noise concerns) but
qualitatively simple -- almost always a single broad-arm loss, sometimes on
each arm, not the fragmented multi-breakpoint pattern that characterizes real
chromothripsis. Domain-expert judgment: on sparse panel data, where segment
size is an unreliable proxy for genomic extent, it's hard to justify 2 loss
segments as chromothripsis.

A raw FACETS segment counts as a "loss segment" if `round(tcn.em)==1` **or**
`round(lcn.em)==0` -- complete LOH at any total copy number, not just haploid
loss (catches copy-neutral/amplified LOH too). Adjacent raw segments that are
both loss-flagged merge into one block regardless of which condition
triggered each or whether their exact tcn/lcn differ. **Deliberately has no
`cf.em` reliability gate** -- tested filtering unreliable segments out before
merging and rejected it: the "unreliable" divider segments sitting between
confirmed-real loss blocks (C009-T01-IM6 chr1, C017-T01-IM7 chr11,
C020-T01-IM7 chr16) are real, qualitatively distinct amplified/retained
islands (often the exact same `cf.em` value recurring at multiple genomic
locations -- one uncertain clonal cluster, not noise), and filtering them out
before merging would fuse separate loss events together and erase the most
diagnostic feature in all three cases.

`apply_loss_segment_filter(cohort_df, man_df, min_loss_segments=3)` keeps a
currently-nominated, non-shard chromosome only if **either**
`n_loss_segments >= min_loss_segments` **or** the row independently qualifies
for the same WGD-masked allelic rescue `apply_classifier_v7` already uses
(`frac_loh < v7_allelic_alt_frac_loh_max` and
`n_allelic_state_changes >= v7_allelic_alt_lo`) -- a chromosome with real
minor-allele complexity substituting for literal LOH (e.g. C007-T01-IM6
chr11/chr12, chromothripsis arising after a genome-doubling event) should not
be penalized for having zero loss segments by construction. Loss-dominant/
shard-pathway calls are trusted unconditionally, same as the thin-chain
filter. It deliberately does **not** search for chromosomes v7 never
nominated at all -- a cohort-wide scan for such "hidden" chromosomes did
surface some real-looking cases, but also produced far more noise-driven
false leads from poorly-fit samples, so that kind of search is left as a
manual-review exercise rather than baked into this filter.

Validated (2026-08-22) on top of the v2 thin-chain filter, `min_loss_segments=3`:

| Cohort | Chromosome-level | Sample-level |
|---|---|---|
| Carcinoid (calibration) | 88 -> 69 | 53/238 -> 46/238 (22.3% -> 19.3%) |
| Classical SCLC | 376 -> 206 | 226/679 -> 143/679 (33.3% -> 21.1%) |
| Pancreatic NET/NEC | 165 -> 103 | 80/506 -> 61/506 (15.8% -> 12.1%) |
| GI NET | 68 -> 47 | 37/245 -> 32/245 (15.1% -> 13.1%) |
| aSCLC | 25 -> 21 | 14/25 -> 12/25 (56.0% -> 48.0%) |

Every carcinoid chromosome dropped sits at exactly 2 loss segments -- none at
0 or 1. Spot-checked 5 random flipped samples each in classical SCLC,
pancreatic NET/NEC, and GI NET and confirmed convincing: none looked noisy,
all read as simple broad-arm loss rather than chromothripsis. On carcinoid,
two of the flipped samples (C022-T02-IM7, C023-T01-IM7) are known
to have a *different*, never-nominated-by-v7 chromosome that independently
clears the bar on manual review (chr11 and chr8 respectively) -- this filter
cannot find such rows by design (see above), so the true carcinoid count
accounting for those two is 48/238, not 46/238.

```python
import ctx_features as cf
final = cf.apply_loss_segment_filter(filtered, man_df)   # filtered = apply_thin_chain_filter output
# use final["nominated_final"] / final["sample_flag_final"]
```

Or via the shared driver script:

```bash
python ctx_run_common/phase6_loss_segment_filter.py --manifest <path> --out <cohort_dir> [--v7-filtered-file ...] [--min-loss-segments 3]
```

## Shard threshold: baseline-relative + embedded-in-loss (`compute_baseline_shard_features`)

The loss-dominant/shard pathway originally required a surviving fragment
("shard") to clear `4x` the sample's genome-wide **ploidy**. On a low-ploidy
sample, that sets an unreachably high bar for a chromosome whose local
backdrop is already a deep loss well below ploidy -- so the threshold is now
relative to the chromosome's own loss-level `baseline_cn` instead (floored at
1 to avoid a vacuous `0x` threshold when `baseline_cn` is itself 0, e.g. a
homozygous-deletion backdrop): `shard_thr = shard_baseline_mult * max(baseline_cn, 1)`.

Lowering the anchor this way surfaces far more candidates, most of them not
real shards -- e.g. C002-T01-IM5 chr6, a 21kb-scale segment that only
cleared the new threshold because the chromosome's own baseline was low, but
which sits next to an unrelated *gain* segment rather than loss on one side.
A real shard, on manual review, is not just small -- it's genomically
**embedded** in the loss backdrop: a baseline-state segment exists somewhere
upstream of it AND somewhere downstream of it (not necessarily the immediate
neighbor -- that stricter version was tested and rejected: it also vetoed 15
of the then-33 production shard hits, including unambiguous ones like a
tcn=52 shard, purely because the immediate neighboring segment happened to be
a different elevated state rather than exactly the baseline). Requiring
embedding independently recovered every manual by-eye judgment call made
during review, without needing to also pick a specific multiplier -- `4x` and
`5x` land on the identical hit set once embedding is required. `5x` was kept
as the default (a small margin of safety over the original `4x`, since the
two are equivalent on the calibration cohort).

Validated (2026-08-22) against the carcinoid calibration cohort, using the
production ploidy-based-threshold hit set (pre-embedding) as ground truth:

| | Chromosome-level | Sample-level |
|---|---|---|
| Production (ploidy-based threshold, no embedding) | 16 | 14 |
| + embedding required (still ploidy-based threshold) | 14 (-2, both already-positive chromosome-edge cases) | 14 (no change) |
| Baseline-relative threshold + embedding (`4x` or `5x` -- identical) | 23 | 21 |

Net effect on the full carcinoid pipeline (on top of v2 thin-chain + the
loss-segment filter): 46/238 -> 50/238. The 4 newly-nominated samples are
C003-T01-IM5 (chr10), C008-T01-IM6 (chr11), C011-T01-IM6
(chr3), and C022-T02-IM7 (chr11) -- the last of these is the same
chr11 exception noted in the loss-segment filter section above as a known,
previously-manual-only rescue; it is now found automatically and no manual
override is needed for it.

Also re-run (2026-08-22), with the shipped `5x` default, on the other 4
cohorts, but **not yet spot-checked** against real segment data the way the
loss-segment filter was -- treat these numbers as provisional: classical SCLC
143/679 -> 152/679 (21.1% -> 22.4%), pancreatic NET/NEC 61/506 -> 66/506
(12.1% -> 13.0%), GI NET 32/245 -> 33/245 (13.1% -> 13.5%), aSCLC 12/25 ->
12/25 (no change). Unlike carcinoid, `4x` and `5x` are NOT identical in these
other cohorts -- `4x` gives slightly higher counts (167/679, 71/506, 36/245)
-- so the choice of `5x` does matter here and has not been separately
validated the way it was on carcinoid.

```python
import ctx_features as cf
shard_df = cf.compute_baseline_shard_features(seg_df, ploidy)   # per-sample; shard_baseline_mult=5.0 default
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

