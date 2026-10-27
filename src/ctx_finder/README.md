# ctx_finder

Detection of chromothripsis from FACETS copy-number profiles of targeted sequencing panels.

`ctx_finder` evaluates every chromosome of a tumor sample and nominates those carrying a
chromothripsis footprint. It is designed for the sparse, coarse copy-number profiles that FACETS
produces from targeted gene panels (e.g. MSK-IMPACT), where whole-genome callers such as
ShatterSeek are not applicable. Calls are made per chromosome and summarized per sample.

Current version: **v0.9**.

## Requirements

- Python ≥ 3.10 with `numpy`, `pandas`, `ruptures` (changepoint detection) and `matplotlib` (plotting
  helpers only). `scikit-learn` is needed only for the optional feature-evaluation helpers.
- R with the `facets` package and `Rscript` on `PATH`, only if probe-level tables must be extracted
  from FACETS `.Rdata` files (they are cached after the first run).
- A Slurm cluster is recommended for the changepoint step; it can also be run sequentially.
- Genome build: hg19 / GRCh37.

## Installation

```bash
git clone https://github.com/jakeleelab/pulmnet.git
conda create -n ctx -c conda-forge python=3.10 numpy pandas matplotlib ruptures
conda activate ctx
```

The modules are used directly from `src/ctx_finder`; no package installation is required.

## Input

A tab-separated manifest with one row per sample:

| Column | Description |
|---|---|
| `SAMPLE_ID` | unique sample identifier |
| `facets_cncf_path` | path to the sample's FACETS segment file (`*.cncf.txt`); the matching `.Rdata` is located next to it |
| `purity`, `ploidy` | FACETS purity and ploidy |
| `chr_ctx` | known chromothriptic chromosomes, used only for evaluation (the column must exist; may be empty) |
| `purity_source` | optional; samples marked `reviewer_override` (manually set purity) are excluded from calling |

## Algorithm

1. **Probe-level copy number.** Each probe's log-ratio is converted to an absolute total copy number
   using the sample's purity and ploidy.
2. **Chromosome features.** For every chromosome, ctx_finder summarizes the copy-number states,
   the longest chain of segments oscillating between a small number of states, allele-specific
   state changes (LOH versus heterozygosity), the extent of the altered footprint, and probe noise.
3. **Changepoint re-segmentation.** Probes are re-segmented independently of FACETS (PELT, with a
   penalty scaled to each chromosome's probe noise), and the oscillation chain is re-scored on the
   inferred changepoints.
4. **Oscillation pathway.** A chromosome is nominated when its combined oscillation and allelic
   score is ≥ 0.35, its oscillation chain spans ≥ 4 segments, and its altered footprint is not
   negligible. A separate allelic criterion captures oscillation after whole-genome doubling, where
   LOH is not reached.
5. **Shard pathway.** A chromosome is also nominated when it is loss-dominant (its modal copy-number
   state lies ≥ 0.5 copies below the sample ploidy and covers ≥ 50% of the chromosome) and retains
   an embedded, non-centromeric fragment at ≥ 5× the modal state. A chromosome with two or more
   such fragments is nominated regardless of its background state, as the alternation between
   background and fragments is itself an oscillation.
6. **Loss-segment filter.** An oscillation call requires ≥ 3 separate loss blocks (total copy
   number 1 or complete LOH) on the chromosome, unless it qualifies through the post-doubling
   allelic criterion.
7. **Within-sample contrast filter.** An oscillation call is kept only if the chromosome stands out
   from the same sample's other chromosomes (mean robust z-score across six complexity features
   ≥ 2.0), which removes calls driven by genome-wide noise.
8. **Noise veto.** Calls from either pathway are rejected when the chromosome's probe-level noise
   exceeds a fixed limit, unless its copy-number steps are at least three times that noise and its
   contrast score is ≥ 1.5 (high-amplitude events that stand out in an otherwise noisy sample).
9. **Sample-level QC.** Calls are cleared in samples with manually set purity or on an optional
   manual exclusion list.
10. **Sample call.** A sample is positive if at least one chromosome is called.

Steps 6–7 apply to oscillation-pathway calls only. Every step acts on each chromosome independently,
so the order of steps 6–8 does not change the calls. Optionally, for patients with several samples,
a sibling-corroboration report tests whether breakpoints of a call recur in the other samples.

## Usage

The command-line drivers in `pipeline/` run the steps in order. All outputs are written under the
directory given by `--out`.

```bash
M=manifest.tsv
O=results/my_cohort
PL=src/ctx_finder/pipeline

# Steps 1-2: features, and a Slurm array script for the changepoint step
python $PL/phase1_build_features.py --manifest $M --out $O \
    --partition <slurm_partition> --env-setup 'conda activate ctx'
sbatch $O/cp_run/run_array.sbatch                 # step 3; wait for completion

python $PL/phase2_collect_v6.py   --manifest $M --out $O    # collect changepoint features
python $PL/phase3_v7_classify.py  --manifest $M --out $O    # steps 4-5
python $PL/phase4_sibling_corroboration.py --manifest $M --out $O --roster <samples.tsv>  # optional report
python $PL/phase5_filters_and_noise_veto.py --manifest $M --out $O   # steps 6-8
python $PL/phase6_sample_qc.py    --manifest $M --out $O [--qc-file <exclusions.tsv>]  # steps 9-10
```

## Output

`<out>/tables/cohort_chrom_features_final.qc_adjusted.tsv`, one row per sample and chromosome:

| Column | Description |
|---|---|
| `nominated_final` | final chromosome-level chromothripsis call |
| `ctx_compatible` | sample-level call (repeated on every row of the sample) |
| `loss_dominant_shard_hit` | chromosome meets the shard-pathway criteria (before the noise veto) |
| `nominated_v7`, `loss_segment_veto`, `contrast_score`, `contrast_veto`, `noise_veto`, `noise_veto_exception`, `step_to_noise` | intermediate calls and filter results |
| `longest_osc_chain_k3`, `cp_osc_chain`, `n_allelic_state_changes`, `dot_mad_within`, ... | supporting features |

## Modules

| File | Purpose |
|---|---|
| `run_ctx_pipeline.py` | manifest handling, probe-table extraction, cohort feature table |
| `ctx_features.py` | chromosome features, classifier and filters |
| `ctx_changepoint.py`, `ctx_cp_parallel.py` | changepoint re-segmentation and its Slurm array |
| `ctx_sibling.py` | sibling-sample breakpoint corroboration |
| `ctx_label.py`, `ctx_plots.py` | plotting helpers |
| `pipeline/` | command-line drivers |

## Version history

- **v0.9**: thin-chain filter removed (redundant with the later filters); noise veto moved to the
  last chromosome-level step with a step-to-noise exception; command-line drivers renumbered.
- **v0.8**: version used for the manuscript revision. Adds the shard pathway, post-doubling allelic
  criterion, thin-chain, loss-segment and within-sample contrast filters, sibling corroboration,
  and the command-line drivers.
- **v0.1**: version deposited with the original manuscript submission.

## License

MIT; see the repository's `LICENSE`.
