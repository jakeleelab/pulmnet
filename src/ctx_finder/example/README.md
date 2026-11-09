# ctx_finder example data

Ten **simulated** tumor samples for trying out `ctx_finder`. No patient data were used: the
copy-number profiles were designed by hand, and the read counts, probe positions and germline
genotypes are synthetic.

The profiles are modeled on the patterns of chromothripsis we observed in lung carcinoids:
chromothripsis of chromosomes 3 and 11 with or without high-level *CCND1* amplification, and
chromothripsis of chromosome 1p with or without *MYCL* amplification. They also include
cases designed to show the other parts of the algorithm: an amplification-based candidate,
chromothripsis after whole-genome doubling, a sample with manually set purity, and negatives.

## How the data were made

1. **Probe layout.** About 72,000 SNP sites on GRCh37, laid out like a targeted capture panel:
   clusters of closely spaced, deeply covered on-target sites plus a sparse, shallow off-target
   backbone. About 12% of sites are heterozygous in the germline, assigned at random.
2. **Designed profiles.** Each sample has a designed total and minor-allele copy number along
   the genome and a tumor purity.
3. **Read counts.** Tumor and normal read counts were simulated from the profiles, with
   realistic depth variation and counting noise.
4. **FACETS.** The counts were processed with FACETS 0.6.2 using the settings of our MSK-IMPACT
   runs (seed 100, `snp.nbhd` 250, `ndepth` 35, purity run `cval` 100, segmentation `cval` 50,
   `min.nhet` 15). Purity, ploidy, segments and integer copy numbers are FACETS' own estimates,
   so they carry ordinary FACETS imperfections.

After FACETS filtering each sample has about 39,000 probes, and probe noise is in the range
seen in real MSK-IMPACT samples with chromothripsis.

## Samples

| Sample | Designed profile | Expected ctx_finder result |
|---|---|---|
| EX01 | chromothripsis of whole chr3 and 11q; *CCND1* amplification (~25 copies) | chr3, chr11 |
| EX02 | chromothripsis of 3p and whole chr11; *CCND1* amplification (~32 copies) | chr3, chr11 |
| EX03 | chromothripsis of whole chr3 and 11q after whole-genome doubling (no LOH); *CCND1* amplification (~45 copies) | chr3, chr11 (through the post-doubling allelic criterion) |
| EX04 | chromothripsis of 1p; *MYCL* amplification (~10 copies) | chr1 |
| EX09 | chromothripsis of whole chr3 and 11q; no amplification | chr3, chr11 |
| EX10 | chromothripsis of 1p; no amplification | chr1 |
| EX05 | chromosome-wide loss of chr12 with focal *CDK4* and *MDM2* amplifications | chr12 (amplification-based candidate) |
| EX06 | arm-level copy-number changes only | negative |
| EX07 | chromothripsis of 3p in a sample whose purity was set manually | negative (cleared by sample QC) |
| EX08 | low-depth, noisy library with arm-level changes only | negative |

## Files

| File | Content |
|---|---|
| `manifest.tsv` | sample sheet: `SAMPLE_ID`, `facets_cncf_path`, `facets_jointseg_path`, `purity`, `ploidy`, `purity_source` (paths relative to this folder) |
| `data/EXnn.cncf.txt` | FACETS segment table (hisens run) |
| `data/EXnn.jointseg.tsv.gz` | FACETS probe-level table (`out$jointseg`) |
| `expected_calls.tsv` | expected sample-level results |
| `run_example.sh` | runs the pipeline and compares the calls with `expected_calls.tsv` |

EX07 is marked `reviewer_override` in `purity_source`, which stands for a purity set manually
during review; its calls are therefore cleared.

## Running the example

From the repository root, with the `ctx` environment active (see the
[ctx_finder README](../README.md) for installation):

```bash
bash src/ctx_finder/example/run_example.sh example_output
```

This runs both steps without Slurm. The changepoint step takes about 10 minutes per sample, so
the whole example takes about 1.5 hours on one core. The script ends by printing each
sample's calls next to the expected results.

To use a Slurm cluster instead, run the two steps yourself:

```bash
E=src/ctx_finder/example
python src/ctx_finder/pipeline/01_prepare.py --manifest $E/manifest.tsv --out example_output \
    --partition <slurm_partition> --env-setup 'source <miniconda>/etc/profile.d/conda.sh; conda activate ctx'
sbatch example_output/work/changepoints/run_array.sbatch     # wait until all tasks have finished
python src/ctx_finder/pipeline/02_call.py --manifest $E/manifest.tsv --out example_output
```

The results are written to `example_output/ctx_finder_calls.tsv` (one row per sample and
chromosome) and `example_output/ctx_finder_samples.tsv` (one row per sample); the
[ctx_finder README](../README.md#output) describes every column.

## Note

These samples are for demonstration and testing only. They are idealized compared with real
tumors and should not be used to estimate how well `ctx_finder` performs.
