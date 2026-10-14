"""
Phase 6 (any cohort): loss-segment filter, then within-sample contrast filter, on top of
nominated_v7_filtered (i.e. run AFTER phase5_thin_chain_filter.py).

1. Loss-segment filter (LOCKED 2026-08-22): keep a conventional-pathway (non-shard) call
   only if the chromosome has >=min_loss_segments merged "loss" blocks (tcn.em==1 or
   lcn.em==0, no cf.em reliability gate -- see ctx_features._n_loss_segments docstring for
   why), OR the row independently qualifies for apply_classifier_v7's own WGD-masked allelic
   rescue (frac_loh<v7_allelic_alt_frac_loh_max AND n_allelic_state_changes>=v7_allelic_alt_lo).
   min_loss_segments=3 was checked against manual domain-expert review (carcinoid full cohort
   + random 5-sample spot-checks in classical SCLC, pancreatic NET/NEC, GI NET).
2. Within-sample contrast filter (ADOPTED 2026-10-02, same rule for every cohort): keep an
   oscillation-pathway call only if the chromosome stands out from the same sample's other
   chromosomes (contrast_score >= min_contrast_score). See
   ctx_features.apply_within_sample_contrast_filter and pulmnet/ctx_contrast_20260919/.
Loss-dominant/shard-pathway calls are trusted unconditionally by both filters.

Usage:
  python pipeline/phase6_loss_segment_filter.py --manifest <path> --out <dir> [--v7-filtered-file cohort_chrom_features_v7_filtered.tsv] [--min-loss-segments 3] [--min-contrast-score 2.0] [--no-contrast-filter]
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
    ap.add_argument("--v7-filtered-file", default="cohort_chrom_features_v7_filtered.tsv",
                    help="filename (within <out>/tables) of the thin-chain-filtered output to filter further")
    ap.add_argument("--min-loss-segments", type=int, default=3)
    ap.add_argument("--min-contrast-score", type=float, default=2.0)
    ap.add_argument("--no-contrast-filter", action="store_true",
                    help="skip the within-sample contrast filter (reproduces the pre-2026-10-02 output)")
    ap.add_argument("--result-suffix", default="",
                    help="optional suffix for the output filename, e.g. '.roster' -> cohort_chrom_features_final.roster.tsv")
    args = ap.parse_args()

    tables = os.path.join(args.out, "tables")
    man_df = pd.read_csv(args.manifest, sep="\t")
    man_df["SAMPLE_ID"] = man_df["SAMPLE_ID"].astype(str)
    v7f = pd.read_csv(os.path.join(tables, args.v7_filtered_file), sep="\t")
    v7f = v7f.drop(columns=[c for c in v7f.columns if c in
                   ["n_loss_segments", "wgd_rescue_eligible", "loss_segment_veto",
                    "nominated_final", "sample_flag_final", "loss_segment_min_count",
                    "contrast_score", "contrast_veto", "contrast_min_score"]])

    n_samples = v7f["SAMPLE_ID"].nunique()
    n_chrom_before = int(v7f["nominated_v7_filtered"].sum())
    sample_before = v7f.groupby("SAMPLE_ID")["nominated_v7_filtered"].any()
    n_samp_before = int(sample_before.sum())

    final = cf.apply_loss_segment_filter(v7f, man_df, min_loss_segments=args.min_loss_segments)
    n_chrom_loss = int(final["nominated_final"].sum())
    n_samp_loss = int(final.groupby("SAMPLE_ID")["nominated_final"].any().sum())
    if not args.no_contrast_filter:
        final = cf.apply_within_sample_contrast_filter(final, min_contrast_score=args.min_contrast_score)

    out_path = os.path.join(tables, f"cohort_chrom_features_final{args.result_suffix}.tsv")
    final.to_csv(out_path, sep="\t", index=False)

    n_chrom_after = int(final["nominated_final"].sum())
    sample_after = final.groupby("SAMPLE_ID")["nominated_final"].any()
    n_samp_after = int(sample_after.sum())
    flipped = sample_before.index[(sample_before) & (~sample_after)]

    print(f"Wrote {out_path} ({len(final)} chromosome rows)")
    print(f"min_loss_segments: {args.min_loss_segments} | contrast filter: "
          f"{'off' if args.no_contrast_filter else f'min_contrast_score {args.min_contrast_score}'}")
    print(f"Chromosome-level nominated calls: {n_chrom_before} -> {n_chrom_loss} (loss-segment)"
          f" -> {n_chrom_after} (contrast)")
    print(f"Sample-level flagged: {n_samp_before}/{n_samples} -> {n_samp_loss}/{n_samples} (loss-segment) "
          f"-> {n_samp_after}/{n_samples} ({n_samp_after / n_samples * 100:.1f}%) (contrast)")
    print(f"Samples flipped flagged->unflagged: {len(flipped)}")
    if len(flipped):
        print(list(flipped))


if __name__ == "__main__":
    main()
