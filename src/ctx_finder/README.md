# ctx_finder

Chromothripsis-footprint detection on sparse FACETS (targeted-panel) copy-number data.

**Version 0.8.** Versions: v0.1 is the code deposited with the original manuscript submission; v0.8 is
the version used for the revision (classifier v7 plus the filters described below). Internal names such as
`nominated_v7` refer to the classifier generation, not the package version.

## Summary

`ctx_finder` scans each chromosome of a tumor sample and nominates the ones that carry a chromothripsis footprint. It works on the coarse, sparse copy-number profiles produced by FACETS on targeted gene panels, where classic WGS-based callers (ShatterSeek and similar) may not apply.

For every (sample, chromosome) it builds copy-number oscillation features (the longest contiguous chain confined to a small number of CN states, allelic LOH-versus-het interleaving, an altered-footprint extent, and a probe-level changepoint confidence), plus a separate loss-dominant "shard" description. Fixed-threshold rules combine these into a per-chromosome call (`nominated_final`), summarized as a sample-level flag (`ctx_compatible`).

The per-sample changepoint pass is run as a Slurm array; it can be adapted to other clusters or run sequentially.

```
ctx_finder/
├── run_ctx_pipeline.py    # manifest + cohort feature builder
├── ctx_changepoint.py     # per-chromosome changepoint features
├── ctx_features.py        # chromosome features, classifiers (v6, v7) and post-hoc filters
├── ctx_cp_parallel.py     # Slurm-array launcher / worker / collector
├── ctx_label.py           # per-chromosome inspection plots
├── ctx_plots.py           # cohort-level diagnostic plots
├── ctx_sibling.py         # sibling-sample corroboration check
└── pipeline/              # command-line drivers, phases 1-7 (see "Running guide")
```

## Algorithm, step by step

1. **Input:** per sample, the FACETS segments (cncf) and probe-level data (jointseg), with FACETS purity and ploidy.
2. **Probe copy number:** each probe's log-ratio is converted to an absolute total copy number using purity and ploidy.
3. **Chromosome features:** copy-number states, oscillation chains, allele-specific state changes, altered footprint and probe noise per chromosome.
4. **Changepoint re-segmentation:** probes are re-segmented independently of FACETS (Slurm array) and the oscillation chain is re-scored on those changepoints.
5. **Oscillation pathway (classifier v7):** a chromosome is nominated if its combined oscillation/allelic score is >= 0.35, its oscillation chain has >= 4 segments, and its altered footprint is not weak (footprint < 15% of the chromosome with <= 3 islands is vetoed).
6. **Shard pathway:** a chromosome is also nominated if it is loss-dominant (modal state >= 0.5 copies below ploidy, covering >= 50% of the chromosome) and carries an embedded, non-centromeric fragment at >= 5x that baseline; with two or more such fragments the loss requirement is waived (multi-shard oscillation).
7. **Noise veto:** nominations from either pathway are rejected on chromosomes with probe noise (`dot_mad_within`) > 1.0.
8. **Thin-chain filter:** segments with FACETS cell fraction (`cf.em`) < 0.15 are cut out of the chain; an oscillation call is dropped if the longest surviving chain has < 4 segments.
9. **Loss-segment filter:** an oscillation call needs >= 3 loss segments on the chromosome, unless whole-genome doubling explains their absence.
10. **Within-sample contrast filter:** an oscillation call is dropped if the chromosome does not stand out from the same sample's other chromosomes (contrast score < 2.0).
11. **QC exclusion:** calls are cleared in samples whose FACETS purity was set manually on review (`purity_source == "reviewer_override"`) or that are on an optional manual QC list.
12. **Roll-up:** a sample is positive if any chromosome is called (`ctx_compatible`).

Steps 8-10 apply only to oscillation-pathway calls; shard-pathway calls are exempt. Each step is described in detail in the sections below.

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
python pipeline/phase5_thin_chain_filter.py --manifest <path> --out <cohort_dir> [--v7-file ...] [--min-chain 4]
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
python pipeline/phase6_loss_segment_filter.py --manifest <path> --out <cohort_dir> [--v7-filtered-file ...] [--min-loss-segments 3]
```

## Within-sample contrast filter (`ctx_features.apply_within_sample_contrast_filter`)

Third post-hoc filter, applied after the loss-segment filter, with the same rule in every
cohort. Adopted 2026-10-02.

**Motivation.** Manual review of the pancreatic NET calls against the pathology sign-out
CN plots found many oscillation-pathway calls in noisy samples (CN jitter across the whole
genome) where the called chromosome did not look different from the sample's other
chromosomes. No sample-level noise metric (FACETS `sd_cnlr_residual`, probe-to-probe
jitter, subclonal-segment fraction, FACETS QC flags) separated those from real calls
without also removing reviewed carcinoid true positives.

**Rule.** For each chromosome and each of six complexity features (`cp_osc_chain`,
`cp_n_changepoints`, `n_allelic_state_changes`, `longest_osc_chain_k3`,
`n_altered_islands`, `segments_per_probe`), a robust z-score against the same sample's
other chromosomes: `(x - median(others)) / max(1.4826 * MAD(others), floor)`. The floors
are fixed constants (`CONTRAST_FEATURE_FLOORS`). `contrast_score` is the mean of the
z-scores clipped to [-5, 10]. An oscillation-pathway call is kept only if
`contrast_score >= 2.0`. Shard-pathway calls are exempt, because oscillation features don't
describe them. The median and MAD of the *other* chromosomes stay robust when a minority of
chromosomes are shattered, so multi-chromosome chromothripsis is kept.

**Validation** (`pulmnet/ctx_contrast_20260919/`):
- **Carcinoid:** 48/48 positive samples kept. Only C017-T02 chr2/7/11/19 were dropped,
  and the sample stays positive via its chr12/18 shard calls.
- **Pancreatic NET (reviewed calls):** 21/22 endorsed samples kept. C013-T01 chr3
  (score 1.8) was lost, which was accepted. 7/23 noisy and 3/9 overcalled samples removed.

**Rejected alternatives.**
- Per-sample 2-cluster splitting lost 12/48 carcinoid positives.
- Merging LOH runs across lcn = NA segments (option C) lost 18 carcinoid patients. Short
  retained or gained chromothripsis fragments are often lcn = NA.
- A sample-level `sd_cnlr_residual` cutoff removed reviewed carcinoid true positives.

| Cohort (frozen 2026_09_19) | Patients before -> after |
|---|---|
| Lung carcinoid | 41/236 -> 41/236 (17.4%) |
| Pancreatic NET | 54/450 -> 44/450 (9.8%) |
| GI NET | 29/241 -> 28/241 (11.6%) |
| Classical SCLC | 145/658 -> 123/658 (18.7%) |
| aSCLC | 12/19 -> 12/19 |

Applied automatically by `pipeline/phase6_loss_segment_filter.py`. Use
`--no-contrast-filter` to reproduce the pre-2026-10-02 output, and `--min-contrast-score` to
change the cutoff.

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
override is needed for it. (C011-T01-IM6's chr3 was subsequently found
to be a centromeric-shard artifact -- see next section -- and no longer
contributes to this sample's nomination; the patient stays chromothripsis
positive via a second, unrelated sample. Current carcinoid count: 49/238.)

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

## Centromeric-shard artifact filter (`CENTROMERE_HG19_MB`, `compute_baseline_shard_features`)

Manual review (2026-08-22) of loss-dominant/shard-pathway hits found a
recurring artifact distinct from the small-fragment cases above: a single
"shard" that is actually the one FACETS segment spanning the low-mappability
pericentromeric region -- a probe-count cliff (5-40 markers) against the
confidently-called flanking arms (hundreds), not a real amplified fragment.
Example, C011-T01-IM6 chr3: p-arm 1.45-89.5 Mb (tcn=1, 697 markers),
then one segment 89.5-116.5 Mb (tcn=5, only 39 markers, straddling the
centromere), then q-arm 118.6-197.4 Mb (tcn=1, 673 markers) -- the "shard"
is the low-confidence middle segment, not a real fragment.

A candidate shard is now rejected if its genomic range overlaps
`CENTROMERE_HG19_MB` for its chromosome (an approximate hg19
centromere-plus-immediate-pericentromeric-heterochromatin table). Rejecting
only the centromeric candidate, rather than vetoing the whole chromosome,
correctly leaves a genuine second shard on the same chromosome untouched.

The pattern is rare in carcinoid but common in the other cohorts, and the
same few pericentromeric coordinates recur across unrelated patients (GI
NET: all 4 hits at chr18 ~14-19 Mb; pancreatic NET/NEC: chr10 ~39-42.6 Mb
x5, chr2 ~92.3 Mb x3) -- recurrence at identical coordinates across
unrelated patients is what makes this a panel/assay artifact rather than
scattered noise, unlike the small-fragment false positives above.

Validated (2026-08-22) across all 5 cohorts -- "removed" counts the
`loss_dominant_shard_hit` rows rejected by the centromere check; sample-
and patient-level flips are the net effect after re-running the full
downstream pipeline (a chromosome can still end up `nominated_final` via
an independent pathway or a real shard on the same chromosome, as happened
for pancreatic NET/NEC below):

| Cohort | Loss-dominant hits (before) | Centromeric (removed) | Sample flips | Patient flips | Final count |
|---|---|---|---|---|---|
| Carcinoid | 34 | 3 | 1 | 0 | 50/238 -> 49/238 |
| Classical SCLC | 71 | 7 | 3 | 3 | 152/679 -> 149/679 |
| GI NET | 15 | 4 | 2 | 2 | 33/245 -> 31/245 |
| Pancreatic NET/NEC | 56 | 13 | 0 | 0 | 66/506 -> 66/506 |
| aSCLC | 6 | 0 | 0 | 0 | 12/25 -> 12/25 |

Pancreatic NET/NEC's 10 final-positive artifact rows caused zero flips: 9
were on samples that already had another, unrelated real final chromosome,
and the 10th (C018-T01-IM7 chr18) still nominates on its own merits
through the conventional oscillation pathway once the shard-pathway hit is
removed -- the classifier correctly falls back rather than losing the call.

```python
import ctx_features as cf
cf.CENTROMERE_HG19_MB            # {"1": (121.5, 128.9), ...} -- override via
shard_df = cf.compute_baseline_shard_features(seg_df, ploidy, centromere_mb=...)  # centromere_mb kwarg if needed
```

## Multi-shard oscillation (`compute_baseline_shard_features`, `multi_shard_min`)

Adopted 2026-10-03 for every cohort. When a chromosome carries **two or more embedded
shards**, the alternation baseline -> shard -> baseline -> shard -> baseline is itself an
oscillation between two copy-number states. For such chromosomes the loss-dominance
requirements (`deficit_min`, `baseline_frac_min`) are waived; everything else about a shard
still applies (>= 5x the baseline, <= 15% of the chromosome, embedded in the baseline state,
non-centromeric), as does the v7 noise veto. A chromosome with a single shard keeps the
original loss-dominant rule.

New columns in `cohort_baseline_shard_features.tsv`: `loss_dominant_hit` (original rule),
`multi_shard_hit` (this rule); `loss_dominant_shard_hit` (the shard-pathway flag used
downstream) is now their OR. `multi_shard_min=None` restores the previous behaviour.

Motivation: expert review of pancreatic and GI NET misses where the shattered chromosome sits
on a copy-neutral or gained backdrop with several high-level fragments (C006-T01 chr6,
C019-T02 chr6, C015-T01 chr12; all confirmed on clinical CN plots).

Effect (2026_09_19 runs; patients with `ctx_compatible`):

| Cohort | New shard-pathway chromosomes | Of which called (rest noise-vetoed) | Patients before -> after |
|---|---|---|---|
| Lung carcinoid | 4 | 4 (2 already called otherwise) | 41/236 -> 41/236 |
| Pancreatic NET | 4 | 4 (2 already called) | 44/450 -> 46/450 |
| GI NET | 3 | 1 | 28/241 -> 29/241 |
| Classical SCLC | 10 | 4 (3 already called) | 123/658 -> 123/658 |
| aSCLC | 3 | 2 (both already called) | 12/19 -> 12/19 |

New carcinoid chromosomes (already-positive patients): C017-T02 chr11, C021-T01 chr12.

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

Run the drivers in `pipeline/` in order. Each takes the manifest (`--manifest`) and one output directory
(`--out`); everything is written under `<out>/tables` and `<out>/cp_run`. Run from any directory; the
scripts locate the package themselves. Restrict the manifest to the samples of the cohort you want to call.

```bash
M=/path/to/manifest.tsv
O=/path/to/results/my_cohort

# 1. per-chromosome features + Slurm array script for the changepoint pass (penalty 8)
python pipeline/phase1_build_features.py --manifest $M --out $O \
    --partition <slurm_partition> \
    --env-setup 'source <miniconda>/etc/profile.d/conda.sh; conda activate ctx'
sbatch $O/cp_run/run_array.sbatch       # wait until done; logs in cp_run/logs/

# 2. collect changepoint features, v6 stage
python pipeline/phase2_collect_v6.py --manifest $M --out $O
# 3. shard features + classifier v7 (oscillation and shard pathways, noise veto)
python pipeline/phase3_v7_classify.py --manifest $M --out $O
# 4. (report only) sibling-sample breakpoint corroboration
python pipeline/phase4_sibling_corroboration.py --manifest $M --out $O --roster <sample_table.tsv>
# 5. thin-chain filter
python pipeline/phase5_thin_chain_filter.py --manifest $M --out $O
# 6. loss-segment filter, then within-sample contrast filter
python pipeline/phase6_loss_segment_filter.py --manifest $M --out $O
# 7. reviewer-purity / manual QC exclusion
python pipeline/phase7_reviewer_purity_exclusion.py --manifest $M --out $O [--qc-file <exclusions.tsv>]
```

Phases 2-7 take seconds to minutes; phase 1 and the changepoint array dominate the run time. The library
functions can also be called directly from Python (see the sections above for each filter's API).

### Output

`<out>/tables/cohort_chrom_features_final.qc_adjusted.tsv` has one row per (sample, chromosome):

- `nominated_final`: the final per-chromosome chromothripsis call (boolean)
- `ctx_compatible`: True on every row of a sample with at least one called chromosome
- `loss_dominant_shard_hit`: the chromosome passed the shard criteria (set before the noise veto, so it can be True on an uncalled chromosome)
- intermediate calls from each stage: `nominated_v7`, `nominated_v7_filtered`, plus `thin_chain_veto`, `loss_segment_veto`, `contrast_score` / `contrast_veto`
- supporting features: `longest_osc_chain_k3`, `cp_osc_chain`, `dot_mad_within` (probe noise), `altered_union_frac`, `n_altered_islands`, and others

## Version history

- **v0.8** (revision): classifier v7 (manual-review-calibrated thresholds, alternate allelic path, loss-dominant shard pathway, second-look confidence); sibling corroboration report; thin-chain filter v2; loss-segment filter; baseline-relative, embedded shard rule with centromeric-artifact filter; within-sample contrast filter; multi-shard oscillation rule; command-line pipeline drivers.
- **v0.1** (original submission): feature builder, changepoint pass and classifier v6.
