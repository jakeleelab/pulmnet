# ctx_finder

Detection of chromothripsis from FACETS copy-number profiles of targeted sequencing panels.

`ctx_finder` is a rule-based algorithm that evaluates every chromosome of a tumor sample for two
patterns of chromothripsis: oscillation of copy number between a limited number of states
(**CN oscillation-based candidates**), and focal high-level amplification on a chromosome that is
otherwise lost or that carries several such amplifications (**amplification-based candidates**).
Candidates then pass a series of filters and a final probe-level noise veto, and a sample is called
positive if at least one chromosome remains called. It is designed for the sparse copy-number
profiles that FACETS produces from targeted gene panels (e.g. MSK-IMPACT), where whole-genome
callers such as ShatterSeek are not applicable.

Current version: **v1.0**.

## Requirements

- Python ≥ 3.10 with `numpy`, `pandas` and `ruptures` (changepoint detection).
- R with the `facets` package and `Rscript` on `PATH`, only if probe-level tables must be extracted
  from FACETS `.Rdata` files (they are cached after the first run).
- A Slurm cluster is recommended for the changepoint step; it can also be run sequentially.
- Genome build: hg19 / GRCh37.

## Installation

```bash
git clone https://github.com/jakeleelab/pulmnet.git
conda create -n ctx -c conda-forge python=3.10 numpy pandas ruptures
conda activate ctx
```

The modules are used directly from `src/ctx_finder`; no package installation is required.

## Input

A tab-separated manifest with one row per sample:

| Column | Description |
|---|---|
| `SAMPLE_ID` | unique sample identifier |
| `facets_cncf_path` | path to the sample's FACETS segment file (`*.cncf.txt`); the matching `.Rdata` (probe-level data) is located next to it |
| `purity`, `ploidy` | purity and ploidy of the FACETS fit |
| `purity_source` | optional; calls are cleared in samples marked `reviewer_override` (purity set manually) |

Optionally, a tab-separated file with a `SAMPLE_ID` column lists samples excluded on manual quality
review (`--qc-file`); their calls are cleared as well.

## Usage

```bash
M=manifest.tsv
O=results/my_cohort
PL=src/ctx_finder/pipeline

# Step 1: probe tables, chromosome features, and a Slurm array script for the changepoint step
python $PL/01_prepare.py --manifest $M --out $O \
    --partition <slurm_partition> \
    --env-setup 'source <miniconda>/etc/profile.d/conda.sh; conda activate ctx'
sbatch $O/work/changepoints/run_array.sbatch      # wait until all array tasks have finished

# Step 2: candidates, filters, noise veto, sample QC, output tables
python $PL/02_call.py --manifest $M --out $O [--qc-file exclusions.tsv]
```

Without Slurm, add `--run-changepoints-here` to `01_prepare.py` to run the changepoint step in the
same process (about 10 minutes per sample), then run `02_call.py`. Every parameter below can be
changed on the command line (`python $PL/02_call.py --help`); the values used are written to
`ctx_finder_parameters.tsv`.

## Algorithm

### Input data and preprocessing

For each sample, `ctx_finder` uses the FACETS segment table (total copy number, `tcn.em`;
minor-allele copy number, `lcn.em`), the probe-level log-ratios and heterozygosity of the same run,
and the purity (p) and ploidy (ψ) of the FACETS fit. Autosomes and chromosome X are analyzed. Probe
log-ratios are converted to absolute total copy number:

    CN = [ (p·ψ + 2(1 − p)) · 2^cnlr − 2(1 − p) ] / p

**Probe noise** (σ̂, in copies) is estimated per chromosome from the copy-number differences between
neighboring probes. In a clean profile, neighboring probes have nearly the same copy number, so large
and frequent differences indicate technical scatter. A median-based measure keeps the estimate
unaffected by the few true copy-number changes along the chromosome:

    σ̂ = 1.4826 · median|Δ − median(Δ)| / √2,   Δ = differences between neighboring probes

(1.4826 converts a median absolute deviation to a standard deviation for normally distributed data;
√2 accounts for each difference carrying the noise of two probes.)

### CN oscillation-based candidates

FACETS segments of each chromosome are ordered by position, and each segment's integer total copy
number is taken as its copy-number state. The **oscillating stretch** is the largest number of
consecutive segments that together visit at most three distinct copy-number states. It is long when
a chromosome alternates repeatedly between a few states and short for stepwise changes.

Allelic oscillation is summarized by the number of switches between loss of heterozygosity (LOH)
and heterozygosity, and by the number of changes in minor-allele copy number. The **altered
fraction** is the fraction of the chromosome outside its most frequent copy-number state, and
**altered islands** are the contiguous runs of such segments.

To obtain evidence independent of the FACETS segmentation, probe log-ratios are re-segmented with
the pruned exact linear time (PELT) changepoint algorithm (L2 cost, minimum segment of 5 probes).
The penalty is scaled to each chromosome's noise, β = λ · σ² · ln(n probes), with λ = 8 and σ² the
larger of the median FACETS log-ratio variance and the squared successive-difference noise of the
log-ratios. Re-segmented levels are converted to integer copy-number states, and the
**changepoint confidence** is

    C = (oscillating stretch on the re-segmented states) · min(1, median(|level shift| / σ) / 3)

that is, the oscillating stretch weighted by how clearly the level shifts exceed the noise.

The **oscillation score** is the product of three terms, each scaled to [0, 1]:

    S = F_osc · F_allelic · ramp(C; 0, 3)
    F_osc     = ramp(oscillating stretch; 2, 6)
    F_allelic = 0.25 + 0.75 · max( ramp(LOH↔heterozygosity switches; 1, 6), A_postWGD )
    A_postWGD = ramp(minor-allele changes; 8, 15)   if the LOH fraction is < 0.05, else 0

where ramp(x; a, b) rises linearly from 0 at a to 1 at b. F_allelic is 0.25 when the chromosome has
fewer than 4 heterozygous probes. The A_postWGD term credits repeated minor-allele changes without
LOH, to capture chromothripsis that occurs after whole-genome doubling, where losses to LOH are not
expected.

**Second look.** The default penalty is deliberately strict so that noisy samples do not produce
spurious changepoints. A chromosome with S before changepoint weighting ≥ 0.8 but changepoint
confidence 0 is re-segmented once with λ = 5, and the higher confidence is used.

A chromosome is a CN oscillation-based candidate if S ≥ 0.35, its oscillating stretch spans at least
4 segments, and its altered footprint is not negligible (negligible: altered fraction < 0.15 and at
most 3 altered islands).

### Amplification-based candidates

Chromothripsis can also give rise to high-level, often extrachromosomal, amplification of a few
fragments while much of the chromosome is lost; at panel probe density such events appear as a few
amplified segments rather than as oscillation. For each chromosome, the **baseline** is the
copy-number state covering the largest length. A **focal amplification** is a segment that

- covers at most 15% of the chromosome,
- has total copy number at least 5 × max(baseline, 1),
- has baseline-state segments on both sides, and
- does not overlap the centromere or pericentromeric region, where FACETS frequently places
  sparse-probe segments at recurrent coordinates.

A chromosome is an amplification-based candidate if it shows **chromosome-wide copy-number loss**
(baseline at least 0.5 copies below the sample ploidy and covering at least 50% of the chromosome)
and carries at least one focal amplification, or if it carries two or more focal amplifications
regardless of its baseline.

### Filters and noise veto

Two filters are applied to CN oscillation-based candidates. Amplification-based candidates are
exempt, because their small, high-level and often subclonal segments would otherwise be misread by
these filters.

1. **Loss-segment filter.** Segments with total copy number 1 or minor-allele copy number 0 are
   classified as losses, and adjacent losses are merged into blocks. At least 3 loss blocks are
   required on the chromosome, unless it meets the post-whole-genome-doubling allelic criterion
   (LOH fraction < 0.05 and at least 8 minor-allele changes); chromosomes with one or two losses
   typically reflect simple arm-level events.
2. **Within-sample contrast filter.** To remove calls driven by genome-wide copy-number scatter,
   six complexity features of the chromosome (oscillating stretch, changepoint count and oscillating
   stretch on the re-segmented states, minor-allele changes, altered islands, and segments per
   probe) are compared with the other chromosomes of the same sample:
   z = (x − median(others)) / max(1.4826 · MAD(others), floor), clipped to [−5, 10]. The
   **contrast score** is the mean of the six z-scores, and candidates with a contrast score below
   2.0 are removed.

Finally, the **noise veto** rejects candidates of both types if probe noise exceeds σ̂ = 1.0
copies, unless the median copy-number step between adjacent segments is at least 3 × σ̂ and the
contrast score is at least 1.5, which retains high-amplitude events that stand out in an otherwise
noisy sample. Because each step acts on each chromosome independently, the order of the filters and
the noise veto does not affect the calls.

### Sample calls

Calls are cleared in samples whose FACETS purity was set manually (`purity_source` =
`reviewer_override`), indicating a profile that FACETS could not fit automatically, and in samples
on the optional exclusion list. A sample is called positive if at least one chromosome remains
called.

## Parameters

Tunable parameters (defaults are the values used for the published calls; each is a command-line
option of `02_call.py`, except `changepoint_penalty`, which belongs to `01_prepare.py`):

| Parameter | Default | Meaning |
|---|---|---|
| `changepoint_penalty` | 8 | PELT penalty multiplier λ |
| `score_threshold` | 0.35 | minimum oscillation score S |
| `min_oscillating_stretch` | 4 | minimum oscillating stretch (segments) |
| `footprint_min_fraction` | 0.15 | negligible footprint: altered fraction below this ... |
| `footprint_max_islands` | 3 | ... and at most this many altered islands |
| `post_wgd_min_changes` | 8 | post-WGD allelic credit starts at this many minor-allele changes ... |
| `post_wgd_full_changes` | 15 | ... and is full at this many |
| `post_wgd_max_loh_fraction` | 0.05 | ... only when the LOH fraction is below this |
| `second_look_score` | 0.8 | second look: S before changepoint weighting ≥ this ... |
| `second_look_penalty` | 5 | ... and changepoint confidence 0: re-segment with this λ |
| `focal_amp_multiple` | 5 | focal amplification: copy number ≥ this × max(baseline, 1) |
| `focal_amp_max_fraction` | 0.15 | focal amplification covers at most this fraction of the chromosome |
| `loss_min_deficit` | 0.5 | chromosome-wide loss: baseline ≥ this many copies below ploidy ... |
| `loss_min_baseline_fraction` | 0.5 | ... and covering ≥ this fraction of the chromosome |
| `multi_focal_amp_min` | 2 | this many focal amplifications qualify on any baseline |
| `min_loss_blocks` | 3 | loss-segment filter |
| `min_contrast_score` | 2.0 | within-sample contrast filter |
| `max_probe_noise` | 1.0 | noise veto: maximum probe noise σ̂ (copies) |
| `step_to_noise_min` | 3.0 | noise-veto exception: median copy-number step / σ̂ ≥ this ... |
| `exception_contrast_min` | 1.5 | ... and contrast score ≥ this |

Fixed constants that define the shape of the oscillation score (`SCORE_CONSTANTS` in
`ctx_features.py`): at most 3 states in an oscillating stretch; F_osc ramp 2–6 segments;
LOH↔heterozygosity ramp 1–6 switches; allelic floor 0.25; at least 4 heterozygous probes for allelic
evidence; changepoint-confidence ramp 0–3; PELT minimum segment 5 probes.

## Output

All files are written under `--out`.

**`ctx_finder_calls.tsv`**: one row per sample and chromosome.

| Column | Description |
|---|---|
| `sample_id`, `chrom` | sample and chromosome (1–22, X) |
| `chromothripsis` | final chromosome-level call |
| `candidate_type` | `cn_oscillation`, `amplification`, `both` or `none` (candidates before filters and noise veto) |
| `oscillation_score` | oscillation score S |
| `oscillating_stretch` | oscillating stretch of the FACETS segments |
| `changepoint_confidence` | changepoint confidence C (after the second look) |
| `altered_fraction`, `n_altered_islands` | altered footprint |
| `baseline_cn` | baseline copy-number state |
| `chromosome_wide_loss` | chromosome-wide loss criterion met |
| `n_focal_amplifications` | number of focal amplifications |
| `n_loss_blocks` | number of loss blocks |
| `post_wgd_allelic` | post-whole-genome-doubling allelic criterion met |
| `contrast_score` | within-sample contrast score |
| `probe_noise` | probe noise σ̂ (copies) |
| `cn_step_to_noise` | median copy-number step between adjacent segments / σ̂ (used by the noise-veto exception) |
| `filter` | step that removed a candidate: `loss_segment`, `contrast`, `noise_veto` or `sample_qc`; empty otherwise |

**`ctx_finder_samples.tsv`**: one row per manifest sample.

| Column | Description |
|---|---|
| `sample_id`, `purity`, `ploidy` | from the manifest |
| `qc_status` | `ok`, `manual_purity`, `manual_exclusion`, or `not_analysed` (no probe-level data, purity or ploidy) |
| `chromothripsis` | sample-level call (`True`/`False`; `NA` if not analysed) |
| `n_chromosomes_called`, `called_chromosomes` | number and list of called chromosomes |

**`ctx_finder_parameters.tsv`**: the parameter values used.

`<out>/work/` holds intermediate files (cached probe tables, chromosome features, changepoint
results and logs) and can be deleted after the run, except that keeping `work/probe_tables/` avoids
re-extracting the probe tables in later runs.

## Modules

| File | Purpose |
|---|---|
| `run_ctx_pipeline.py` | manifest handling, probe-table extraction, chromosome feature table |
| `ctx_features.py` | parameters, chromosome features, candidates, filters and noise veto |
| `ctx_changepoint.py`, `ctx_cp_parallel.py` | changepoint re-segmentation and its Slurm array |
| `pipeline/01_prepare.py`, `pipeline/02_call.py` | command-line drivers |

## Version history

- **v1.0**: first public release.

## License

MIT; see the repository's `LICENSE`.
