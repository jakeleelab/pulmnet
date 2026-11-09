"""
Step 1 of 2: extract probe-level tables, compute chromosome features, and write the
changepoint re-segmentation job (a Slurm array script, or run it sequentially with
--run-changepoints-here).

Usage:
  python pipeline/01_prepare.py --manifest <path> --out <dir>
      [--partition <slurm partition>] [--env-setup "source <conda.sh>; conda activate ctx"]
      [--chunk 4] [--changepoint-penalty 8] [--run-changepoints-here]

Writes <out>/work/ (probe tables, chromosome features, changepoint job). Then submit
<out>/work/changepoints/run_array.sbatch and, when it has finished, run 02_call.py.
"""
import argparse
import os
import sys

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)  # src/ctx_finder
sys.path.insert(0, REPO)

import ctx_cp_parallel as cpx
import ctx_features as cf
import run_ctx_pipeline as rp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--changepoint-penalty", type=float, default=cf.PARAMETERS["changepoint_penalty"])
    ap.add_argument("--chunk", type=int, default=4, help="samples per Slurm array task")
    ap.add_argument("--partition", default="componc_cpu", help="Slurm partition for the changepoint array")
    ap.add_argument("--env-setup", default="",
                    help="shell line(s) that activate the Python environment on compute nodes, "
                         "e.g. 'source <miniconda>/etc/profile.d/conda.sh; conda activate ctx'")
    ap.add_argument("--run-changepoints-here", action="store_true",
                    help="run the changepoint step sequentially in this process instead of writing a Slurm array")
    args = ap.parse_args()

    work = os.path.join(args.out, "work")
    probe_dir = os.path.join(work, "probe_tables")
    cp_dir = os.path.join(work, "changepoints")
    for d in (probe_dir, cp_dir):
        os.makedirs(d, exist_ok=True)

    man_df = rp.load_manifest(args.manifest)
    proc, ledger = rp.build_processable_sample_ids(man_df, probe_dir)
    ledger.to_csv(os.path.join(work, "processable_samples.tsv"), sep="\t", index=False)
    print(f"Processable samples: {len(proc)} / {len(man_df)}")

    cohort, failed = rp.build_cohort_features(man_df, probe_dir, proc)
    cohort.to_csv(os.path.join(work, "chrom_features.tsv"), sep="\t", index=False)
    failed.to_csv(os.path.join(work, "feature_failures.tsv"), sep="\t", index=False)
    print(f"Chromosome features: {len(cohort)} rows, {len(failed)} samples failed")
    if len(failed):
        print(failed.head(20).to_string(index=False))

    ok = sorted(set(proc) - set(failed["SAMPLE_ID"]))
    if not ok:
        sys.exit("No sample could be processed; see the errors above.")
    with open(os.path.join(work, "changepoint_penalty.txt"), "w") as fh:
        fh.write(f"{args.changepoint_penalty}\n")
    if args.run_changepoints_here:
        cpx.run_worker(ok, cp_dir, probe_dir, args.changepoint_penalty, args.manifest)
        print("Changepoint step done; next: 02_call.py")
        return
    script, n_tasks = cpx.launch(cp_dir, ok, args.changepoint_penalty, args.chunk, args.manifest, probe_dir,
                                 partition=args.partition, env_setup=args.env_setup)
    print(f"Changepoint array: {script} ({n_tasks} tasks)")
    if not args.env_setup:
        print("No --env-setup given: make sure the array tasks can import pandas/numpy/ruptures.")
    print(f"Submit with: sbatch {script}; when it has finished, run 02_call.py")


if __name__ == "__main__":
    main()
