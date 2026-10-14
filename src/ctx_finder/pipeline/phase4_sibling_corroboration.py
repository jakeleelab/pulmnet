"""
Phase 4 (any cohort): sibling breakpoint-corroboration report, restricted to
the samples listed in --roster (the processed pool can be broader than the
cohort of interest). Informational only; does not change any call.

Usage:
  python pipeline/phase4_sibling_corroboration.py --manifest <path> --out <dir> --roster <path> [--roster-id-col SAMPLE_ID]
"""
import argparse
import os
import sys

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)  # src/ctx_finder
sys.path.insert(0, REPO)

import pandas as pd
import ctx_sibling as csib


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--roster", required=True)
    ap.add_argument("--roster-id-col", default="SAMPLE_ID")
    args = ap.parse_args()

    tables = os.path.join(args.out, "tables")
    man_df = pd.read_csv(args.manifest, sep="\t")
    v7 = pd.read_csv(os.path.join(tables, "cohort_chrom_features_v7.tsv"), sep="\t")

    roster = pd.read_csv(args.roster, sep="\t")
    roster_ids = set(roster[args.roster_id_col].astype(str))
    if args.roster_id_col == "SAMPLE_ID":
        keep = v7["SAMPLE_ID"].astype(str).isin(roster_ids)
    else:
        v7_patient = v7["SAMPLE_ID"].astype(str).str.split("-T").str[0]
        keep = v7_patient.isin(roster_ids)
    n_before = v7["SAMPLE_ID"].nunique()
    v7 = v7[keep]
    n_after = v7["SAMPLE_ID"].nunique()
    print(f"Roster filter: {n_before} -> {n_after} samples ({n_before - n_after} excluded)")

    report = csib.sibling_corroboration_report(v7, man_df)
    report.to_csv(os.path.join(tables, "sibling_corroboration_report.tsv"), sep="\t", index=False)

    print(f"Discordant sibling pairs checked: {len(report)}")
    if len(report):
        print(report.sort_values("percentile", ascending=False).to_string(index=False))


if __name__ == "__main__":
    main()
