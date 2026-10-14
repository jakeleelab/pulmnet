"""
Phase 1 (any cohort): build per-chromosome cohort features + write the Slurm
array script for the changepoint pass. Only --manifest/--out differ between
cohorts.

Usage:
  python pipeline/phase1_build_features.py --manifest <path> --out <dir> [--pen 8] [--chunk 4]
      [--partition <slurm partition>] [--env-setup "source <conda.sh>; conda activate ctx"]
"""
import argparse
import os
import sys

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)  # src/ctx_finder
sys.path.insert(0, REPO)

import pandas as pd
import run_ctx_pipeline as rp
import ctx_cp_parallel as cpx


def _inject_env_setup(sbatch_path, env_setup):
    """cpx.launch()'s generated sbatch script has no env activation -- add the
    given shell lines right before the worker invocation so `sbatch` submissions
    work standalone (compute nodes don't inherit the login shell's environment)."""
    with open(sbatch_path) as fh:
        text = fh.read()
    marker = "mkdir -p"
    idx = text.index(marker)
    line_end = text.index("\n", idx) + 1
    text = text[:line_end] + env_setup.rstrip("\n") + "\n" + text[line_end:]
    with open(sbatch_path, "w") as fh:
        fh.write(text)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pen", type=float, default=8.0)
    ap.add_argument("--chunk", type=int, default=4)
    ap.add_argument("--partition", default="componc_cpu", help="Slurm partition for the changepoint array")
    ap.add_argument("--env-setup", default="",
                    help="shell line(s) that activate the Python environment on compute nodes, "
                         "e.g. 'source <miniconda>/etc/profile.d/conda.sh; conda activate ctx'")
    args = ap.parse_args()

    tables = os.path.join(args.out, "tables")
    cp_out = os.path.join(args.out, "cp_run")
    for d in (tables, cp_out):
        os.makedirs(d, exist_ok=True)

    man_df = pd.read_csv(args.manifest, sep="\t")

    proc, ledger = rp.build_processable_sample_ids(man_df, tables)
    print(f"Processable samples: {len(proc)} / {len(man_df)}")
    ledger.to_csv(os.path.join(tables, "processable_ledger.tsv"), sep="\t", index=False)

    cohort, failed = rp.build_cohort_features(man_df, tables, proc, use_cache=True)
    cohort.to_csv(os.path.join(tables, "cohort_chrom_features.tsv"), sep="\t", index=False)
    print(f"Cohort feature rows: {len(cohort)}, failed samples: {len(failed)}")
    if len(failed):
        print("Failed samples (up to 20):")
        print(failed.head(20).to_string(index=False))

    script = cpx.launch(out_dir=cp_out, pens=[args.pen], chunk=args.chunk, partition=args.partition,
                        man_path=args.manifest, tablesDir=tables)
    if args.env_setup:
        _inject_env_setup(script, args.env_setup)
        print(f"Slurm array script written (environment setup injected): {script}")
    else:
        print(f"Slurm array script written: {script} -- no --env-setup given; make sure the "
              f"array tasks can import pandas/numpy (activate the env in the script if needed)")
    print(f"Submit with: sbatch {script}")


if __name__ == "__main__":
    main()
