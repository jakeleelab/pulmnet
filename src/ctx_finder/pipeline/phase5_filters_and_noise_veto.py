"""
Phase 5 (any cohort): oscillation-call filters, then the probe-level noise veto, on top of the
phase-3 classifier output (nominated_v7, nominated without the noise veto).

1. Loss-segment filter: keep an oscillation-pathway (non-shard) call only if the chromosome has
   >=min_loss_segments merged "loss" blocks (tcn.em==1 or lcn.em==0; see
   ctx_features._n_loss_segments), OR the row qualifies for apply_classifier_v7's post-WGD allelic
   criterion (frac_loh<v7_allelic_alt_frac_loh_max AND n_allelic_state_changes>=v7_allelic_alt_lo).
2. Within-sample contrast filter: keep an oscillation-pathway call only if the chromosome stands
   out from the same sample's other chromosomes (contrast_score >= min_contrast_score). See
   ctx_features.apply_within_sample_contrast_filter.
3. Noise veto (both pathways): remove a call when the chromosome's probe noise exceeds
   noise_max, unless its median copy-number step is >= step_to_noise_min x the noise AND its
   contrast_score >= exception_contrast_min. See ctx_features.apply_noise_veto.
Loss-dominant/shard-pathway calls are exempt from filters 1 and 2.

Usage:
  python pipeline/phase5_filters_and_noise_veto.py --manifest <path> --out <dir> [--v7-file cohort_chrom_features_v7.tsv] [--min-loss-segments 3] [--min-contrast-score 2.0] [--noise-max 1.0] [--step-to-noise-min 3.0] [--exception-contrast-min 1.5] [--no-noise-exception]
"""
import argparse
import os
import sys

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)  # src/ctx_finder
sys.path.insert(0, REPO)

import pandas as pd
import ctx_features as cf

DERIVED = ["n_loss_segments", "wgd_rescue_eligible", "loss_segment_veto", "nominated_final",
           "sample_flag_final", "loss_segment_min_count", "contrast_score", "contrast_veto",
           "contrast_min_score", "cn_step_median", "step_to_noise", "noise_veto",
           "noise_veto_exception", "noise_max", "step_to_noise_min", "exception_contrast_min"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--v7-file", default="cohort_chrom_features_v7.tsv",
                    help="filename (within <out>/tables) of the phase-3 classifier output")
    ap.add_argument("--min-loss-segments", type=int, default=3)
    ap.add_argument("--min-contrast-score", type=float, default=2.0)
    ap.add_argument("--noise-max", type=float, default=1.0)
    ap.add_argument("--step-to-noise-min", type=float, default=3.0)
    ap.add_argument("--exception-contrast-min", type=float, default=1.5)
    ap.add_argument("--no-noise-exception", action="store_true",
                    help="apply the noise veto without the step-to-noise exception")
    ap.add_argument("--result-suffix", default="",
                    help="optional suffix for the output filename, e.g. '.roster' -> cohort_chrom_features_final.roster.tsv")
    args = ap.parse_args()

    tables = os.path.join(args.out, "tables")
    man_df = pd.read_csv(args.manifest, sep="\t")
    man_df["SAMPLE_ID"] = man_df["SAMPLE_ID"].astype(str)
    v7 = pd.read_csv(os.path.join(tables, args.v7_file), sep="\t")
    v7 = v7.drop(columns=[c for c in v7.columns if c in DERIVED])
    if v7["noise_veto_v7"].any():
        sys.exit("phase-3 output already has the noise veto applied; re-run phase 3 (noise_max=None)")

    n_samples = v7["SAMPLE_ID"].nunique()

    def counts(df, col):
        return int(df[col].sum()), int(df.groupby("SAMPLE_ID")[col].any().sum())

    c0, s0 = counts(v7, "nominated_v7")
    final = cf.apply_loss_segment_filter(v7, man_df, min_loss_segments=args.min_loss_segments)
    c1, s1 = counts(final, "nominated_final")
    final = cf.apply_within_sample_contrast_filter(final, min_contrast_score=args.min_contrast_score)
    c2, s2 = counts(final, "nominated_final")
    final = cf.apply_noise_veto(final, man_df, noise_max=args.noise_max,
                                step_to_noise_min=None if args.no_noise_exception else args.step_to_noise_min,
                                exception_contrast_min=args.exception_contrast_min)
    c3, s3 = counts(final, "nominated_final")

    out_path = os.path.join(tables, f"cohort_chrom_features_final{args.result_suffix}.tsv")
    final.to_csv(out_path, sep="\t", index=False)

    print(f"Wrote {out_path} ({len(final)} chromosome rows)")
    print(f"min_loss_segments {args.min_loss_segments} | min_contrast_score {args.min_contrast_score} | "
          f"noise_max {args.noise_max} | step-to-noise exception "
          f"{'off' if args.no_noise_exception else f'>= {args.step_to_noise_min} with contrast >= {args.exception_contrast_min}'}")
    print(f"Chromosome-level calls: {c0} -> {c1} (loss-segment) -> {c2} (contrast) -> {c3} (noise veto; "
          f"{int(final['noise_veto_exception'].sum())} kept by the exception)")
    print(f"Sample-level flagged: {s0} -> {s1} -> {s2} -> {s3} of {n_samples} ({s3 / n_samples * 100:.1f}%)")


if __name__ == "__main__":
    main()
