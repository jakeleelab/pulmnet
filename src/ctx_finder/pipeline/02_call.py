"""
Step 2 of 2: nominate CN oscillation-based and amplification-based candidates, apply the
loss-segment filter, the within-sample contrast filter, the probe-level noise veto and
sample QC, and write the results. Run after 01_prepare.py and its changepoint job.

Usage:
  python pipeline/02_call.py --manifest <path> --out <dir> [--qc-file <exclusions.tsv>]
      [--<parameter> <value> ...]   (see --help; defaults are ctx_features.PARAMETERS)

Writes <out>/ctx_finder_calls.tsv (sample x chromosome), <out>/ctx_finder_samples.tsv
(sample) and <out>/ctx_finder_parameters.tsv (the values used).
"""
import argparse
import os
import sys

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)  # src/ctx_finder
sys.path.insert(0, REPO)

import numpy as np
import pandas as pd

import ctx_cp_parallel as cpx
import ctx_features as cf
import run_ctx_pipeline as rp

CALL_COLUMNS = {
    "SAMPLE_ID": "sample_id",
    "chrom": "chrom",
    "nominated_final": "chromothripsis",
    "candidate_type": "candidate_type",
    "ctx_suppress_v7": "oscillation_score",
    "longest_osc_chain_k3": "oscillating_stretch",
    "changepoint_confidence": "changepoint_confidence",
    "altered_union_frac": "altered_fraction",
    "n_altered_islands": "n_altered_islands",
    "baseline_cn": "baseline_cn",
    "chromosome_wide_loss": "chromosome_wide_loss",
    "n_shards": "n_focal_amplifications",
    "n_loss_segments": "n_loss_blocks",
    "post_wgd_allelic": "post_wgd_allelic",
    "contrast_score": "contrast_score",
    "dot_mad_within": "probe_noise",
    "step_to_noise": "cn_step_to_noise",
    "filter": "filter",
}


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--qc-file", default=None,
                    help="optional TSV with a SAMPLE_ID column of samples excluded on manual quality review")
    for name, value in cf.PARAMETERS.items():
        if name != "changepoint_penalty":
            ap.add_argument(f"--{name.replace('_', '-')}", type=type(value), default=value)
    args = ap.parse_args()
    par = {n: getattr(args, n) for n in cf.PARAMETERS if n != "changepoint_penalty"}

    work = os.path.join(args.out, "work")
    probe_dir = os.path.join(work, "probe_tables")
    with open(os.path.join(work, "changepoint_penalty.txt")) as fh:
        par["changepoint_penalty"] = float(fh.read().strip())
    man_df = rp.load_manifest(args.manifest)
    ploidy_map = dict(zip(man_df["SAMPLE_ID"], pd.to_numeric(man_df["ploidy"], errors="coerce")))

    # chromosome features + changepoint features
    feat = pd.read_csv(os.path.join(work, "chrom_features.tsv"), sep="\t")
    feat["chrom"] = feat["chrom"].map(cf._norm_chrom)
    cp, cp_failed = cpx.collect(os.path.join(work, "changepoints"))
    if cp_failed:
        print(f"Changepoint step failed for {len(cp_failed)} samples: {cp_failed[:10]}")
    cp["chrom"] = cp["chrom"].map(cf._norm_chrom)
    d = feat.merge(cp, on=["SAMPLE_ID", "chrom"], how="left")

    # per-chromosome segment features: focal amplifications, loss blocks, copy-number steps
    shard_frames, block_frames = [], []
    for sid in sorted(d["SAMPLE_ID"].unique()):
        seg_df = rp.load_cncf(man_df, sid)
        sf = cf.compute_baseline_shard_features(
            seg_df, ploidy_map.get(sid),
            shard_baseline_mult=par["focal_amp_multiple"], deficit_min=par["loss_min_deficit"],
            baseline_frac_min=par["loss_min_baseline_fraction"], small_frac_max=par["focal_amp_max_fraction"],
            multi_shard_min=par["multi_focal_amp_min"])
        sf.insert(0, "SAMPLE_ID", sid)
        shard_frames.append(sf)
        bf = cf.segment_block_features(seg_df)
        bf.insert(0, "SAMPLE_ID", sid)
        block_frames.append(bf)
    shard_df = pd.concat(shard_frames, ignore_index=True)
    d = d.merge(pd.concat(block_frames, ignore_index=True), on=["SAMPLE_ID", "chrom"], how="left")

    # candidates
    second_look = cf.compute_second_look_confidence(cf.compute_ctx_score(d), man_df, probe_dir,
                                                    score_min=par["second_look_score"],
                                                    pen_fallback=par["second_look_penalty"])
    d = cf.apply_classifier_v7(d, shard_df, second_look_df=second_look,
                               threshold=par["score_threshold"], min_osc_chain=par["min_oscillating_stretch"],
                               footprint_union_min=par["footprint_min_fraction"],
                               footprint_island_max=par["footprint_max_islands"],
                               allelic_alt_lo=par["post_wgd_min_changes"],
                               allelic_alt_hi=par["post_wgd_full_changes"],
                               allelic_alt_frac_loh_max=par["post_wgd_max_loh_fraction"])
    d = d.merge(shard_df.drop(columns=["loss_dominant_shard_hit"]).assign(chrom=lambda x: x["chrom"].astype(str)),
                on=["SAMPLE_ID", "chrom"], how="left")

    # filters (CN oscillation-based candidates), noise veto (all), sample QC
    d = cf.apply_loss_segment_filter(d, min_loss_segments=par["min_loss_blocks"])
    d = cf.apply_within_sample_contrast_filter(d, min_contrast_score=par["min_contrast_score"])
    d = cf.apply_noise_veto(d, noise_max=par["max_probe_noise"], step_to_noise_min=par["step_to_noise_min"],
                            exception_contrast_min=par["exception_contrast_min"])
    qc_status = pd.Series("ok", index=man_df["SAMPLE_ID"])
    if "purity_source" in man_df.columns:
        qc_status[man_df.set_index("SAMPLE_ID")["purity_source"].eq("reviewer_override")] = "manual_purity"
    if args.qc_file:
        excluded = set(pd.read_csv(args.qc_file, sep="\t")["SAMPLE_ID"].astype(str))
        qc_status[qc_status.index.isin(excluded)] = "manual_exclusion"
    qc_cleared = d["SAMPLE_ID"].map(qc_status).ne("ok")
    sample_qc_veto = d["nominated_final"] & qc_cleared
    d["nominated_final"] = d["nominated_final"] & ~qc_cleared

    # outputs
    osc, amp = d["nominated_oscillation"], d["loss_dominant_shard_hit"]
    d["candidate_type"] = np.select([osc & amp, osc, amp], ["both", "cn_oscillation", "amplification"], "none")
    d["filter"] = np.select([d["loss_segment_veto"], d["contrast_veto"], d["noise_veto"], sample_qc_veto],
                            ["loss_segment", "contrast", "noise_veto", "sample_qc"], "")
    d["_order"] = d["chrom"].map(cf._chrom_order)
    d = d.sort_values(["SAMPLE_ID", "_order"])
    d["chrom"] = d["chrom"].replace({"23": "X", "24": "Y"})
    calls = d[list(CALL_COLUMNS)].rename(columns=CALL_COLUMNS).round(4)
    calls.to_csv(os.path.join(args.out, "ctx_finder_calls.tsv"), sep="\t", index=False)

    called = calls[calls["chromothripsis"]].groupby("sample_id")["chrom"].agg(list)
    samples = man_df[["SAMPLE_ID", "purity", "ploidy"]].rename(columns={"SAMPLE_ID": "sample_id"})
    analysed = set(calls["sample_id"])
    samples["qc_status"] = np.where(samples["sample_id"].isin(analysed),
                                    samples["sample_id"].map(qc_status), "not_analysed")
    samples["n_chromosomes_called"] = samples["sample_id"].map(called.str.len()).fillna(0).astype(int)
    samples["called_chromosomes"] = samples["sample_id"].map(called.str.join(",")).fillna("")
    samples["chromothripsis"] = np.where(samples["sample_id"].isin(analysed),
                                         (samples["n_chromosomes_called"] > 0).astype(str), "NA")
    samples = samples[["sample_id", "purity", "ploidy", "qc_status", "chromothripsis",
                       "n_chromosomes_called", "called_chromosomes"]]
    samples.to_csv(os.path.join(args.out, "ctx_finder_samples.tsv"), sep="\t", index=False)
    pd.DataFrame({"parameter": list(par), "value": list(par.values())}).to_csv(
        os.path.join(args.out, "ctx_finder_parameters.tsv"), sep="\t", index=False)

    n_an = len(analysed)
    n_pos = int((samples["chromothripsis"] == "True").sum())
    print(f"Samples analysed: {n_an} of {len(samples)}; chromothripsis-positive: {n_pos}")
    print(f"Chromosome calls: {int(calls['chromothripsis'].sum())}; candidates removed by step: "
          + ", ".join(f"{k} {v}" for k, v in calls["filter"].value_counts().items() if k))
    print(f"Wrote {args.out}/ctx_finder_calls.tsv, ctx_finder_samples.tsv, ctx_finder_parameters.tsv")


if __name__ == "__main__":
    main()
