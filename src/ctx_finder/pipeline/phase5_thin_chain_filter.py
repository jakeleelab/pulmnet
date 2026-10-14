"""
Phase 5 (any cohort): thin-chain reliability filter on top of nominated_v7.

v2 (2026-08-22): treats an unreliable oscillation-chain segment (cf.em below
threshold) as a BREAK in the segment sequence and rechecks whether the best
surviving chain still clears min_chain, rather than v1's all-or-nothing veto
of the whole chromosome. See ctx_features.apply_thin_chain_filter docstring
and the ctx_finder README ("Thin-chain reliability filter") for the full
v1-vs-v2 validation narrative -- min_chain=4 (matching apply_classifier_v7's
own min_osc_chain) was checked against manual domain-expert review on the
carcinoid calibration cohort before superseding v1's threshold=0.15 veto.

Usage:
  python pipeline/phase5_thin_chain_filter.py --manifest <path> --out <dir> [--v7-file cohort_chrom_features_v7.tsv] [--min-chain 4]
"""
import argparse
import os
import sys

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)  # src/ctx_finder
sys.path.insert(0, REPO)

import pandas as pd
import ctx_features as cf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--v7-file", default="cohort_chrom_features_v7.tsv",
                    help="filename (within <out>/tables) of the v7 classifier output to filter")
    ap.add_argument("--min-chain", type=int, default=4)
    ap.add_argument("--result-suffix", default="",
                    help="optional suffix for the output filename, e.g. '.roster' -> cohort_chrom_features_v7_filtered.roster.tsv")
    args = ap.parse_args()

    tables = os.path.join(args.out, "tables")
    man_df = pd.read_csv(args.manifest, sep="\t")
    man_df["SAMPLE_ID"] = man_df["SAMPLE_ID"].astype(str)
    v7 = pd.read_csv(os.path.join(tables, args.v7_file), sep="\t")
    v7 = v7.drop(columns=[c for c in v7.columns if c in
                 ["thin_chain_best_chain", "thin_chain_veto", "nominated_v7_filtered",
                  "sample_flag_v7_filtered", "thin_chain_min_chain", "chain_min_cfem_spike",
                  "n_chain_spikes", "thin_chain_threshold"]])

    n_samples = v7["SAMPLE_ID"].nunique()
    n_chrom_before = int(v7["nominated_v7"].sum())
    sample_before = v7.groupby("SAMPLE_ID")["nominated_v7"].any()
    n_samp_before = int(sample_before.sum())

    filtered = cf.apply_thin_chain_filter(v7, man_df, tables, min_chain=args.min_chain)

    out_path = os.path.join(tables, f"cohort_chrom_features_v7_filtered{args.result_suffix}.tsv")
    filtered.to_csv(out_path, sep="\t", index=False)

    n_chrom_after = int(filtered["nominated_v7_filtered"].sum())
    sample_after = filtered.groupby("SAMPLE_ID")["nominated_v7_filtered"].any()
    n_samp_after = int(sample_after.sum())
    flipped = sample_before.index[(sample_before) & (~sample_after)]

    print(f"Wrote {out_path} ({len(filtered)} chromosome rows)")
    print(f"min_chain: {args.min_chain}")
    print(f"Chromosome-level nominated calls: {n_chrom_before} -> {n_chrom_after} "
          f"(-{n_chrom_before - n_chrom_after}, "
          f"{(n_chrom_before - n_chrom_after) / n_chrom_before * 100:.1f}% removed)"
          if n_chrom_before else "no nominated rows")
    print(f"Sample-level flagged: {n_samp_before}/{n_samples} ({n_samp_before / n_samples * 100:.1f}%) "
          f"-> {n_samp_after}/{n_samples} ({n_samp_after / n_samples * 100:.1f}%)")
    print(f"Samples flipped flagged->unflagged: {len(flipped)}")
    if len(flipped):
        print(list(flipped))


if __name__ == "__main__":
    main()
