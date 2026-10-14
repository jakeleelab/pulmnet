"""
Phase 7 (any cohort): exclude samples whose FACETS purity was manually
reviewer-injected -- both the purity-cval run and the hisens run failed to
converge on their own, so a human set purity by hand (often via SNV VAF,
a reasonable but sometimes inaccurate fallback for an overly flat CN
profile). That's a systematic sign of a poor CN profile, so the sample's
chromothripsis call is not trusted even though it's still fully processed
and its raw features stay in the table. Run this AFTER phase6_loss_segment_filter.py.

See the manifest's purity_source column (written by the R cohort-building
scripts): "reviewer_override" is the flag this script acts on.

Usage:
  python pipeline/phase7_reviewer_purity_exclusion.py --manifest <path> --out <dir> [--final-file cohort_chrom_features_final.tsv] [--qc-file <path>] [--result-suffix .qc_adjusted]
"""
import argparse
import os

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--final-file", default="cohort_chrom_features_final.tsv",
                    help="filename (within <out>/tables) of the loss-segment-filtered output to QC-adjust")
    ap.add_argument("--qc-file", default=None,
                    help="optional TSV with a SAMPLE_ID column of hand-curated manual QC exclusions")
    ap.add_argument("--result-suffix", default=".qc_adjusted")
    args = ap.parse_args()

    tables = os.path.join(args.out, "tables")
    man = pd.read_csv(args.manifest, sep="\t")
    man["SAMPLE_ID"] = man["SAMPLE_ID"].astype(str)
    final = pd.read_csv(os.path.join(tables, args.final_file), sep="\t")

    reviewer_override_samples = set(man.loc[man.get("purity_source") == "reviewer_override", "SAMPLE_ID"])

    excluded_samples = set()
    if args.qc_file:
        excl = pd.read_csv(args.qc_file, sep="\t")
        excluded_samples = set(excl["SAMPLE_ID"].astype(str))

    final_adj = final.copy()
    mask = final_adj["SAMPLE_ID"].isin(excluded_samples | reviewer_override_samples)
    final_adj.loc[mask, "nominated_final"] = False
    final_adj.loc[mask, "sample_flag_final"] = False
    final_adj["ctx_compatible"] = final_adj["sample_flag_final"]
    final_adj["PATIENT_ID"] = final_adj["SAMPLE_ID"].str.split("-T").str[0]

    stem, ext = os.path.splitext(args.final_file)
    out_path = os.path.join(tables, f"{stem}{args.result_suffix}{ext}")
    final_adj.to_csv(out_path, sep="\t", index=False)

    n_samples = final_adj["SAMPLE_ID"].nunique()
    n_before = int(final.groupby("SAMPLE_ID")["sample_flag_final"].any().sum())
    n_after = int(final_adj.groupby("SAMPLE_ID")["ctx_compatible"].any().sum())
    n_patients = final_adj["PATIENT_ID"].nunique()
    n_pat_after = int(final_adj.groupby("PATIENT_ID")["ctx_compatible"].any().sum())

    print(f"Wrote {out_path} ({len(final_adj)} chromosome rows)")
    print(f"Manual QC exclusions: {len(excluded_samples)} sample(s) -> {sorted(excluded_samples)}")
    print(f"Reviewer-injected-purity exclusions: {len(reviewer_override_samples)} sample(s) -> {sorted(reviewer_override_samples)}")
    print(f"Sample-level flagged: {n_before}/{n_samples} -> {n_after}/{n_samples}")
    print(f"Patient-level flagged: {n_pat_after}/{n_patients}")


if __name__ == "__main__":
    main()
