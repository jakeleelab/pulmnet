"""
Phase 2 (any cohort): collect the Slurm array's per-sample changepoint
features, merge onto the cohort feature table, apply the v6 classifier.

Usage:
  python pipeline/phase2_collect_v6.py --manifest <path> --out <dir> [--pen 8]
"""
import argparse
import os
import sys

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)  # src/ctx_finder
sys.path.insert(0, REPO)

import pandas as pd
import ctx_features as cf
import ctx_cp_parallel as cpx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pen", type=float, default=8.0)
    args = ap.parse_args()

    tables = os.path.join(args.out, "tables")
    cp_out = os.path.join(args.out, "cp_run")

    man_df = pd.read_csv(args.manifest, sep="\t")
    cohort = pd.read_csv(os.path.join(tables, "cohort_chrom_features.tsv"), sep="\t")

    cpx.collect(cp_out, pen_const=args.pen)
    cp = pd.read_csv(os.path.join(cp_out, "cohort_cp_features.tsv"), sep="\t")
    cp = cp[cp["pen_const"] == args.pen].copy()

    cohort["_c"] = cohort["chrom"].map(cf._norm_chrom)
    cp["_c"] = cp["chrom"].map(cf._norm_chrom)
    m = cohort.merge(cp.drop(columns=["chrom"]), on=["SAMPLE_ID", "_c"], how="left")

    purity_map = dict(zip(man_df["SAMPLE_ID"].astype(str),
                          pd.to_numeric(man_df["purity"], errors="coerce")))
    v6 = cf.apply_classifier_v6(m, purity_map)
    v6.to_csv(os.path.join(tables, "cohort_chrom_features_v6.tsv"), sep="\t", index=False)
    print(f"Wrote {os.path.join(tables, 'cohort_chrom_features_v6.tsv')} ({len(v6)} rows)")
    print(f"Samples with sample_flag_v6=True: {v6.loc[v6['sample_flag_v6'], 'SAMPLE_ID'].nunique()}")


if __name__ == "__main__":
    main()
