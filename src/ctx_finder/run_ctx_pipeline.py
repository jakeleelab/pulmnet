"""
Self-contained chromothripsis-footprint pipeline.

A fresh Python session needs only:
    - this file
    - ctx_features.py   (in the same directory or on PYTHONPATH)
    - access to the FACETS manifest + output dirs on the cluster
    - the cached {SAMPLE_ID}.hisens.jointseg.tsv files in tablesDir
      (regenerated from .Rdata via Rscript only if a cache file is missing)

Run:
    python run_ctx_pipeline.py
or import build_everything() / main() from a notebook.

It reproduces the loader logic from explore_cases.ipynb (cells 6, 22) so you do
not need anything else from the notebook to get going again after a restart.
"""

import os
import tempfile
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

import ctx_features as cf


# ======================================================================
# CONFIG  -- edit these paths if they move
# ======================================================================

MAN_PATH = "path/to/samples.pulmnet.facets.path.tsv"
OUT_DIR  = "path/to/out_dir"

# total/minor copy number columns to use from the cncf table
TCN_COL = "tcn.em"
LCN_COL = "lcn.em"

# new candidate features (computed by ctx_features.chrom_features_one_sample)
NEW_FEATURES = [
    "state_changes_per_probe",
    "oscillation_excess",
    "two_state_dominance",
    "n_allelic_state_changes",
    "loh_het_interleave",
    "frac_loh",
    "vaf_imbalance_sd",
    "cnlr_snr",
    "n_resolved_levels",
    "n_singleton_runs",
    "median_run_length",
]


# ======================================================================
# directory setup
# ======================================================================

def setup_dirs(out_dir):
    figuresDir = os.path.join(out_dir, "figures")
    tablesDir  = os.path.join(out_dir, "tables")
    dataDir    = os.path.join(out_dir, "data")
    for d in (figuresDir, tablesDir, dataDir):
        os.makedirs(d, exist_ok=True)
    return figuresDir, tablesDir, dataDir


# ======================================================================
# manifest
# ======================================================================

def load_manifest(man_path=MAN_PATH):
    man_df = pd.read_csv(man_path, on_bad_lines="skip", sep="\t")
    for col in ("SAMPLE_ID", "facets_cncf_path", "chr_ctx", "purity", "ploidy"):
        if col not in man_df.columns:
            raise KeyError(f"manifest missing expected column: {col}")
    return man_df


# ======================================================================
# FACETS loaders  (lifted from notebook cell 6, kept behaviour-identical)
# ======================================================================

def chrom_key(x):
    s = str(x).replace("chr", "").replace("CHR", "")
    if s in ["X", "x", "23"]:
        return 23
    if s in ["Y", "y", "24"]:
        return 24
    if s in ["M", "MT", "m", "mt", "25"]:
        return 25
    return int(s)


def get_sample_row(man_df, sample_id):
    sub = man_df.loc[man_df["SAMPLE_ID"].eq(sample_id)].copy()
    if sub.shape[0] == 0:
        raise ValueError(f"No matching SAMPLE_ID found: {sample_id}")
    if sub.shape[0] > 1:
        raise ValueError(f"Multiple matching SAMPLE_ID rows found: {sample_id}")
    return sub.iloc[0]


def get_facets_paths(man_df, sample_id):
    row = get_sample_row(man_df, sample_id)
    cncf_path = str(row["facets_cncf_path"])
    if pd.isna(cncf_path) or cncf_path == "nan":
        raise ValueError(f"facets_cncf_path is missing for {sample_id}")
    run_dir = Path(cncf_path).parent
    stem = Path(cncf_path).name.replace(".cncf.txt", "")
    rdata_path = run_dir / f"{stem}.Rdata"
    if not rdata_path.exists():
        hits = list(run_dir.glob("*_hisens.Rdata"))
        if len(hits) == 0:
            raise FileNotFoundError(f"No Rdata file found in {run_dir}")
        rdata_path = hits[0]
    return str(cncf_path), str(rdata_path)


def load_cncf(man_df, sample_id):
    cncf_path, _ = get_facets_paths(man_df, sample_id)
    seg_df = pd.read_csv(cncf_path, sep="\t")
    if "ID" in seg_df.columns:
        seg_df = seg_df.drop(columns=["ID"])
    seg_df["_chrom_order"] = seg_df["chrom"].map(chrom_key)
    seg_df = seg_df.sort_values(["_chrom_order", "loc.start", "loc.end"]).reset_index(drop=True)
    return seg_df


def write_jointseg_from_rdata(rdata_path, jointseg_path):
    """Only called when the cached TSV is missing. Requires Rscript."""
    r_code = f"""
load("{rdata_path}")
write.table(out$jointseg, file="{jointseg_path}", sep="\\t", quote=FALSE, row.names=FALSE)
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False) as f:
        f.write(r_code)
        r_script = f.name
    res = subprocess.run(["Rscript", r_script], capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(
            f"Rscript failed extracting {rdata_path}:\n{res.stderr}\n"
            "If R is unavailable here, regenerate the jointseg cache where R exists."
        )
    return jointseg_path


def load_or_write_jointseg(man_df, sample_id, tablesDir, use_cache=True, overwrite=False):
    _, rdata_path = get_facets_paths(man_df, sample_id)
    jointseg_path = os.path.join(tablesDir, f"{sample_id}.hisens.jointseg.tsv")
    if overwrite or (not use_cache) or (not os.path.exists(jointseg_path)):
        write_jointseg_from_rdata(rdata_path, jointseg_path)
    jointseg_df = pd.read_csv(jointseg_path, sep="\t")
    jointseg_df["_chrom_order"] = jointseg_df["chrom"].map(chrom_key)
    jointseg_df = jointseg_df.sort_values(["_chrom_order", "maploc"]).reset_index(drop=True)
    return jointseg_df, jointseg_path


def add_tcn_estimate(jointseg_df, man_df, sample_id, cnlr_col="cnlr", out_col="tcn_estimate"):
    row = get_sample_row(man_df, sample_id)
    purity = float(row["purity"])
    ploidy = float(row["ploidy"])
    d = jointseg_df.copy()
    d[out_col] = ((purity * ploidy + (1 - purity) * 2) * (2 ** d[cnlr_col]) - (1 - purity) * 2) / purity
    return d


def prepare_cn_data_for_sample(man_df, sample_id, tablesDir, use_cache=True, overwrite_jointseg=False):
    """seg_df + jointseg_df (with tcn_estimate) for one sample."""
    seg_df = load_cncf(man_df, sample_id)
    jointseg_df, jointseg_path = load_or_write_jointseg(
        man_df, sample_id, tablesDir, use_cache=use_cache, overwrite=overwrite_jointseg)
    jointseg_df = add_tcn_estimate(jointseg_df, man_df, sample_id)
    return seg_df, jointseg_df, jointseg_path


# ======================================================================
# truth labels
# ======================================================================

def parse_chr_ctx_value(x):
    if pd.isna(x):
        return set()
    s = str(x).strip()
    if s == "":
        return set()
    return set(v.strip().replace("chr", "").replace("CHR", "")
               for v in s.split(",") if v.strip() != "")


def truth_set_for_sample(man_df, sample_id):
    row = get_sample_row(man_df, sample_id)
    return parse_chr_ctx_value(row.get("chr_ctx", np.nan))


# ======================================================================
# processability ledger  (condensed from notebook cell 22)
# ======================================================================

def build_processable_sample_ids(man_df, tablesDir):
    """
    A sample is processable if it has a cncf path that exists AND
    (a cached jointseg TSV exists OR a matching .Rdata exists) AND purity/ploidy.
    Returns (processable_ids, ledger_df).
    """
    rows = []
    for _, row in man_df.iterrows():
        sid = row["SAMPLE_ID"]
        cncf = row.get("facets_cncf_path", np.nan)
        has_path = pd.notna(cncf) and str(cncf).strip() not in ("", "nan", "None")
        cncf_exists = Path(str(cncf)).exists() if has_path else False

        cache_tsv = Path(tablesDir) / f"{sid}.hisens.jointseg.tsv"
        cache_exists = cache_tsv.exists()

        rdata_exists = False
        if has_path:
            run_dir = Path(str(cncf)).parent
            stem = Path(str(cncf)).name.replace(".cncf.txt", "")
            cand = run_dir / f"{stem}.Rdata"
            if cand.exists():
                rdata_exists = True
            elif list(run_dir.glob("*_hisens.Rdata")):
                rdata_exists = True

        has_pp = pd.notna(row.get("purity", np.nan)) and pd.notna(row.get("ploidy", np.nan))
        processable = bool(has_path and cncf_exists and (cache_exists or rdata_exists) and has_pp)
        rows.append(dict(
            SAMPLE_ID=sid, has_path=has_path, cncf_exists=cncf_exists,
            cache_exists=cache_exists, rdata_exists=rdata_exists,
            has_purity_ploidy=has_pp, processable=processable,
            has_chr_ctx=len(parse_chr_ctx_value(row.get("chr_ctx", np.nan))) > 0,
        ))
    ledger = pd.DataFrame(rows)
    processable_ids = ledger.loc[ledger["processable"], "SAMPLE_ID"].astype(str).tolist()
    return processable_ids, ledger


# ======================================================================
# cohort feature table
# ======================================================================

def build_cohort_features(man_df, tablesDir, sample_ids, use_cache=True):
    rows, failed = [], []
    for sid in sample_ids:
        try:
            seg_df, jointseg_df, _ = prepare_cn_data_for_sample(
                man_df, sid, tablesDir, use_cache=use_cache)
            feat = cf.chrom_features_one_sample(
                seg_df, jointseg_df, tcn_col=TCN_COL, lcn_col=LCN_COL)
            feat["SAMPLE_ID"] = sid
            truth = truth_set_for_sample(man_df, sid)
            feat["truth_positive_chrom"] = feat["chrom"].astype(str).isin(truth)
            rows.append(feat)
        except Exception as e:
            failed.append(dict(SAMPLE_ID=sid, error=str(e)))
    cohort = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    failed_df = pd.DataFrame(failed)
    return cohort, failed_df


# ======================================================================
# top-level driver
# ======================================================================

def build_everything(man_path=MAN_PATH, out_dir=OUT_DIR, use_cache=True,
                     n_splits=5, save=True):
    figuresDir, tablesDir, dataDir = setup_dirs(out_dir)
    man_df = load_manifest(man_path)

    processable_ids, ledger = build_processable_sample_ids(man_df, tablesDir)
    print(f"[ledger] {len(processable_ids)} processable / {len(man_df)} manifest rows")
    print(f"[ledger] processable & chr_ctx+: "
          f"{int((ledger['processable'] & ledger['has_chr_ctx']).sum())}")

    cohort, failed_df = build_cohort_features(man_df, tablesDir, processable_ids, use_cache=use_cache)
    if failed_df.shape[0]:
        print(f"[features] {failed_df.shape[0]} samples failed feature extraction "
              f"(see returned failed_df)")
    print(f"[features] cohort table: {cohort.shape[0]} chromosome rows, "
          f"{int(cohort['truth_positive_chrom'].sum())} positive chromosomes")

    # confounder correction + within-sample localization contrast
    if "n_probes" in cohort.columns:
        cohort = cf.add_probe_count_residual(cohort, "state_changes_per_probe", "n_probes")
    present = [c for c in NEW_FEATURES if c in cohort.columns]
    cohort = cf.add_within_sample_contrast(cohort, present)

    # per-feature screen
    screen = cf.feature_screen(cohort, present)
    print("\n[screen] per-feature discrimination vs chr_ctx (sorted by AP):")
    print(screen.to_string(index=False))

    # learned combiner -> out-of-fold per-chromosome probability
    use_cols = present + [c for c in cohort.columns if c.endswith("__ratio_to_sample")]
    pred, coef = cf.learned_combiner_oof(cohort, use_cols, n_splits=n_splits)
    print("\n[combiner] L1-logistic coefficients (full-data fit):")
    print(coef.head(15).to_string(index=False))

    # literature-grounded multiplicative score (oscillation x allelic x quality)
    purity_map = dict(zip(man_df["SAMPLE_ID"].astype(str),
                          pd.to_numeric(man_df["purity"], errors="coerce")))
    scored = cf.compute_ctx_score(cohort, purity_map, cf.ScoreParams())
    print("\n[score] ctx_score components by truth label:")
    print(scored.groupby("truth_positive_chrom")[["F_osc", "F_allelic", "Q_weight", "ctx_score"]]
          .mean().to_string())

    score_screen = cf.feature_screen(
        scored, ["ctx_score", "F_osc", "F_allelic", "mech_evidence",
                 "longest_osc_chain_k3", "longest_osc_chain_k2"])
    print("\n[score] discrimination of score & components vs chr_ctx:")
    print(score_screen.to_string(index=False))

    # sample-level threshold sweep on the multiplicative score (the headline)
    score_sweep = cf.threshold_sweep(scored, score_col="ctx_score",
                                     thresholds=np.linspace(0.05, 0.6, 12))
    print("\n[score] sample-level sweep on ctx_score (compare to your rule's 17 errors):")
    print(score_sweep[["threshold", "pos_recovered", "pos_missed",
                       "neg_correct", "neg_false_positive", "sample_errors"]].to_string(index=False))

    # sample-level threshold sweep on the learned combiner (for comparison)
    sweep = cf.threshold_sweep(pred)
    print("\n[combiner] sample-level sweep on p_ctx_oof:")
    print(sweep[["threshold", "pos_recovered", "pos_missed",
                 "neg_correct", "neg_false_positive", "sample_errors"]].to_string(index=False))

    if save:
        ledger.to_csv(os.path.join(tablesDir, "cohort_facets_processing_ledger.tsv"), sep="\t", index=False)
        cohort.to_csv(os.path.join(tablesDir, "cohort_chrom_features_v2.tsv"), sep="\t", index=False)
        screen.to_csv(os.path.join(tablesDir, "cohort_feature_screen_v2.tsv"), sep="\t", index=False)
        coef.to_csv(os.path.join(tablesDir, "cohort_combiner_coefs_v2.tsv"), sep="\t", index=False)
        pred.to_csv(os.path.join(tablesDir, "cohort_chrom_predictions_v2.tsv"), sep="\t", index=False)
        sweep.to_csv(os.path.join(tablesDir, "cohort_sample_sweep_v2.tsv"), sep="\t", index=False)
        scored.to_csv(os.path.join(tablesDir, "cohort_chrom_ctx_score.tsv"), sep="\t", index=False)
        score_screen.to_csv(os.path.join(tablesDir, "cohort_ctx_score_screen.tsv"), sep="\t", index=False)
        score_sweep.to_csv(os.path.join(tablesDir, "cohort_ctx_score_sweep.tsv"), sep="\t", index=False)
        if failed_df.shape[0]:
            failed_df.to_csv(os.path.join(tablesDir, "cohort_feature_failures_v2.tsv"), sep="\t", index=False)
        print(f"\n[save] wrote outputs to {tablesDir}")

    return dict(man_df=man_df, ledger=ledger, cohort=cohort, failed_df=failed_df,
                screen=screen, pred=pred, coef=coef, sweep=sweep,
                scored=scored, score_screen=score_screen, score_sweep=score_sweep,
                figuresDir=figuresDir, tablesDir=tablesDir, dataDir=dataDir)


def main():
    return build_everything()


if __name__ == "__main__":
    main()
