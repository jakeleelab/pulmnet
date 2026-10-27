"""
Phase 3 (any cohort): amplification-based candidate features, gated second-look
confidence, then the v7 classifier (fixed thresholds). The noise veto is deferred to phase 5.

Usage:
  python pipeline/phase3_v7_classify.py --manifest <path> --out <dir>
"""
import argparse
import os
import sys

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)  # src/ctx_finder
sys.path.insert(0, REPO)

import pandas as pd
import ctx_features as cf
import run_ctx_pipeline as rp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    tables = os.path.join(args.out, "tables")
    man_df = pd.read_csv(args.manifest, sep="\t")
    v6 = pd.read_csv(os.path.join(tables, "cohort_chrom_features_v6.tsv"), sep="\t")

    purity_map = dict(zip(man_df["SAMPLE_ID"].astype(str),
                          pd.to_numeric(man_df["purity"], errors="coerce")))
    ploidy_map = dict(zip(man_df["SAMPLE_ID"].astype(str),
                          pd.to_numeric(man_df["ploidy"], errors="coerce")))

    sample_ids = sorted(v6["SAMPLE_ID"].unique())
    shard_frames, failed = [], []
    for sid in sample_ids:
        ploidy = ploidy_map.get(sid)
        try:
            seg_df = rp.load_cncf(man_df, sid)
        except Exception as e:
            failed.append((sid, str(e)))
            continue
        sf = cf.compute_baseline_shard_features(seg_df, ploidy)
        sf.insert(0, "SAMPLE_ID", sid)
        shard_frames.append(sf)

    shard_df = pd.concat(shard_frames, ignore_index=True)
    shard_df.to_csv(os.path.join(tables, "cohort_baseline_shard_features.tsv"), sep="\t", index=False)
    print(f"Amplification-based candidate features: {len(shard_df)} chromosome rows across {len(sample_ids)} samples "
          f"({len(failed)} failed to load)")
    if failed:
        print("Failed samples (up to 20):", failed[:20])

    scored = cf.compute_ctx_score(v6, purity_map, use_quality=False, mech_combine="product")
    second_look_df = cf.compute_second_look_confidence(scored, man_df, tables)
    second_look_df.to_csv(os.path.join(tables, "cohort_second_look_confidence.tsv"), sep="\t", index=False)
    print(f"Second-look confidence: {len(second_look_df)} candidate rows")

    v7 = cf.apply_classifier_v7(v6, purity_map, shard_df, second_look_df=second_look_df,
                            noise_max=None)  # noise veto applied last (phase 5)
    v7.to_csv(os.path.join(tables, "cohort_chrom_features_v7.tsv"), sep="\t", index=False)

    print(f"v7 output: {len(v7)} rows, {v7['SAMPLE_ID'].nunique()} samples")
    print(f"Nominated calls (v7): {int(v7['nominated_v7'].sum())}")
    print(f"Flagged samples (v7): {v7.loc[v7['sample_flag_v7'], 'SAMPLE_ID'].nunique()} / {v7['SAMPLE_ID'].nunique()}")


if __name__ == "__main__":
    main()
