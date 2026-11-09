"""
Changepoint re-segmentation of every processable sample, as a Slurm array or sequentially.

  launch : write the sample list and an sbatch --array script (called by pipeline/01_prepare.py)
  worker : changepoint features for a chunk of samples, one TSV per sample (skips samples done)
  collect: concatenate the per-sample TSVs (called by pipeline/02_call.py)

The worker can also be run without Slurm:
  python ctx_cp_parallel.py worker --samples <out>/work/changepoints/samples.txt \
      --out <out>/work/changepoints --probe-dir <out>/work/probe_tables \
      --manifest <manifest.tsv> --pen 8 --chunk 100000 --array-index 0
"""

import os
import argparse

import pandas as pd

import run_ctx_pipeline as rp
import ctx_changepoint as cc


def _norm_chrom(x):
    return str(x).strip().replace("chr", "").replace("CHR", "")


def compute_one_sample(man_df, sample_id, probe_dir, pen):
    _, jointseg_df = rp.prepare_cn_data_for_sample(man_df, sample_id, probe_dir)
    row = man_df.loc[man_df["SAMPLE_ID"].astype(str).eq(str(sample_id))]
    purity = float(row["purity"].iloc[0]) if row.shape[0] and pd.notna(row["purity"].iloc[0]) else None
    ploidy = float(row["ploidy"].iloc[0]) if row.shape[0] and pd.notna(row["ploidy"].iloc[0]) else 2.0
    js = jointseg_df.copy()
    js["_c"] = js["chrom"].map(_norm_chrom)
    rows = []
    for c, sub in js.groupby("_c"):
        f = cc.probe_changepoint_features(sub, purity=purity, ploidy=ploidy, pen_const=pen)
        f["SAMPLE_ID"] = sample_id
        f["chrom"] = c
        rows.append(f)
    return pd.DataFrame(rows)


def run_worker(sample_ids, out_dir, probe_dir, pen, man_path):
    cp_dir = os.path.join(out_dir, "per_sample")
    os.makedirs(cp_dir, exist_ok=True)
    man_df = rp.load_manifest(man_path)
    for sid in sample_ids:
        out_tsv = os.path.join(cp_dir, f"{sid}.tsv")
        if os.path.exists(out_tsv):
            print(f"[skip] {sid} (exists)", flush=True)
            continue
        try:
            df = compute_one_sample(man_df, sid, probe_dir, pen)
            df.to_csv(out_tsv, sep="\t", index=False)
            print(f"[done] {sid} -> {out_tsv} ({df.shape[0]} rows)", flush=True)
        except Exception as e:
            # write a marker so collect can report failures without crashing
            with open(os.path.join(cp_dir, f"{sid}.FAILED"), "w") as fh:
                fh.write(str(e))
            print(f"[FAIL] {sid}: {e}", flush=True)


def launch(out_dir, sample_ids, pen, chunk, man_path, probe_dir, partition="componc_cpu",
           mem="16G", time="04:00:00", env_setup=""):
    os.makedirs(os.path.join(out_dir, "logs"), exist_ok=True)
    samples = os.path.join(out_dir, "samples.txt")
    with open(samples, "w") as fh:
        fh.write("\n".join(sample_ids) + "\n")
    n_tasks = (len(sample_ids) + chunk - 1) // chunk
    script = os.path.join(out_dir, "run_array.sbatch")
    with open(script, "w") as fh:
        fh.write(f"""#!/bin/bash
#SBATCH --job-name=ctx_cp
#SBATCH --partition={partition}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem={mem}
#SBATCH --time={time}
#SBATCH --array=0-{n_tasks - 1}
#SBATCH --output={out_dir}/logs/task_%a.out

{env_setup.rstrip()}
python {os.path.abspath(__file__)} worker \\
    --samples {samples} \\
    --out {out_dir} \\
    --probe-dir {probe_dir} \\
    --manifest {os.path.abspath(man_path)} \\
    --pen {pen} \\
    --chunk {chunk} \\
    --array-index $SLURM_ARRAY_TASK_ID
""")
    return script, n_tasks


def collect(out_dir):
    cp_dir = os.path.join(out_dir, "per_sample")
    tsvs = sorted(os.path.join(cp_dir, f) for f in os.listdir(cp_dir) if f.endswith(".tsv"))
    fails = sorted(f[:-len(".FAILED")] for f in os.listdir(cp_dir) if f.endswith(".FAILED"))
    cp = pd.concat([pd.read_csv(t, sep="\t") for t in tsvs], ignore_index=True)
    return cp, fails


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    w = sub.add_parser("worker")
    w.add_argument("--samples", required=True)
    w.add_argument("--out", required=True)
    w.add_argument("--probe-dir", required=True)
    w.add_argument("--manifest", required=True)
    w.add_argument("--pen", type=float, required=True)
    w.add_argument("--chunk", type=int, default=1)
    w.add_argument("--array-index", type=int, required=True)
    args = ap.parse_args()

    with open(args.samples) as fh:
        all_ids = [x.strip() for x in fh if x.strip()]
    i0 = args.array_index * args.chunk
    chunk_ids = all_ids[i0:i0 + args.chunk]
    print(f"[worker] task {args.array_index}: samples {i0}..{i0 + len(chunk_ids) - 1}", flush=True)
    run_worker(chunk_ids, args.out, args.probe_dir, args.pen, args.manifest)


if __name__ == "__main__":
    main()
