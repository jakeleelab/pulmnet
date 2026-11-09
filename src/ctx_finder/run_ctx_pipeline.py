"""
Manifest handling, FACETS loaders and the cohort feature table.

The probe-level table (FACETS out$jointseg) of each sample is read from the optional
manifest column facets_jointseg_path (TSV, may be gzipped) when given; otherwise it is
extracted once from the FACETS .Rdata file next to the segment file (requires Rscript) and
cached as {SAMPLE_ID}.hisens.jointseg.tsv in the probe-table directory.

Relative paths in the manifest are resolved from the manifest's own directory.
"""

import os
import tempfile
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

import ctx_features as cf


# total/minor copy number columns to use from the cncf table
TCN_COL = "tcn.em"
LCN_COL = "lcn.em"

REQUIRED_MANIFEST_COLUMNS = ("SAMPLE_ID", "facets_cncf_path", "purity", "ploidy")
PATH_COLUMNS = ("facets_cncf_path", "facets_jointseg_path")   # second one optional


def _given(x):
    return pd.notna(x) and str(x).strip() not in ("", "nan", "None", "NA")


# ======================================================================
# manifest
# ======================================================================

def load_manifest(man_path):
    man_df = pd.read_csv(man_path, on_bad_lines="skip", sep="\t")
    for col in REQUIRED_MANIFEST_COLUMNS:
        if col not in man_df.columns:
            raise KeyError(f"manifest missing expected column: {col}")
    man_df["SAMPLE_ID"] = man_df["SAMPLE_ID"].astype(str)
    base = Path(man_path).resolve().parent
    for col in PATH_COLUMNS:
        if col in man_df.columns:
            man_df[col] = [str(base / str(x)) if _given(x) and not os.path.isabs(str(x)) else x
                           for x in man_df[col]]
    return man_df


# ======================================================================
# FACETS loaders
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


def get_cncf_path(man_df, sample_id):
    cncf_path = get_sample_row(man_df, sample_id)["facets_cncf_path"]
    if not _given(cncf_path):
        raise ValueError(f"facets_cncf_path is missing for {sample_id}")
    return str(cncf_path)


def get_rdata_path(man_df, sample_id):
    """FACETS .Rdata next to the segment file (only needed to extract probe tables)."""
    cncf_path = get_cncf_path(man_df, sample_id)
    run_dir = Path(cncf_path).parent
    stem = Path(cncf_path).name.replace(".cncf.txt", "")
    rdata_path = run_dir / f"{stem}.Rdata"
    if not rdata_path.exists():
        hits = list(run_dir.glob("*_hisens.Rdata"))
        if len(hits) == 0:
            raise FileNotFoundError(f"No probe-level table for {sample_id}: no facets_jointseg_path "
                                    f"in the manifest and no .Rdata file in {run_dir}")
        rdata_path = hits[0]
    return str(rdata_path)


def load_cncf(man_df, sample_id):
    cncf_path = get_cncf_path(man_df, sample_id)
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
            "If R is unavailable here, regenerate the probe-table cache where R exists."
        )
    return jointseg_path


def given_jointseg_path(man_df, sample_id):
    """Probe table named in the manifest (facets_jointseg_path), or None."""
    if "facets_jointseg_path" not in man_df.columns:
        return None
    x = get_sample_row(man_df, sample_id)["facets_jointseg_path"]
    return str(x) if _given(x) else None


def load_or_write_jointseg(man_df, sample_id, probe_dir):
    jointseg_path = given_jointseg_path(man_df, sample_id) or \
        os.path.join(probe_dir, f"{sample_id}.hisens.jointseg.tsv")
    if not os.path.exists(jointseg_path):
        write_jointseg_from_rdata(get_rdata_path(man_df, sample_id), jointseg_path)
    jointseg_df = pd.read_csv(jointseg_path, sep="\t")
    jointseg_df["_chrom_order"] = jointseg_df["chrom"].map(chrom_key)
    jointseg_df = jointseg_df.sort_values(["_chrom_order", "maploc"]).reset_index(drop=True)
    return jointseg_df


def add_tcn_estimate(jointseg_df, man_df, sample_id, cnlr_col="cnlr", out_col="tcn_estimate"):
    """Probe log-ratio -> absolute total copy number, using the sample's purity and ploidy."""
    row = get_sample_row(man_df, sample_id)
    purity = float(row["purity"])
    ploidy = float(row["ploidy"])
    d = jointseg_df.copy()
    d[out_col] = ((purity * ploidy + (1 - purity) * 2) * (2 ** d[cnlr_col]) - (1 - purity) * 2) / purity
    return d


def prepare_cn_data_for_sample(man_df, sample_id, probe_dir):
    """seg_df + jointseg_df (with tcn_estimate) for one sample."""
    seg_df = load_cncf(man_df, sample_id)
    jointseg_df = add_tcn_estimate(load_or_write_jointseg(man_df, sample_id, probe_dir), man_df, sample_id)
    return seg_df, jointseg_df


# ======================================================================
# processable samples
# ======================================================================

def build_processable_sample_ids(man_df, probe_dir):
    """
    A sample is processable if it has a cncf path that exists AND (a probe table named
    in the manifest exists OR a cached probe table exists OR a matching .Rdata exists)
    AND purity/ploidy.
    Returns (processable_ids, ledger_df).
    """
    rows = []
    for _, row in man_df.iterrows():
        sid = row["SAMPLE_ID"]
        cncf = row.get("facets_cncf_path", np.nan)
        has_path = _given(cncf)
        cncf_exists = Path(str(cncf)).exists() if has_path else False

        given = row.get("facets_jointseg_path", np.nan)
        cache_exists = (Path(probe_dir) / f"{sid}.hisens.jointseg.tsv").exists() or \
            (_given(given) and Path(str(given)).exists())

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
        ))
    ledger = pd.DataFrame(rows)
    processable_ids = ledger.loc[ledger["processable"], "SAMPLE_ID"].astype(str).tolist()
    return processable_ids, ledger


# ======================================================================
# cohort feature table
# ======================================================================

def build_cohort_features(man_df, probe_dir, sample_ids):
    rows, failed = [], []
    for sid in sample_ids:
        try:
            seg_df, jointseg_df = prepare_cn_data_for_sample(man_df, sid, probe_dir)
            feat = cf.chrom_features_one_sample(seg_df, jointseg_df, tcn_col=TCN_COL, lcn_col=LCN_COL)
            feat["SAMPLE_ID"] = sid
            rows.append(feat)
        except Exception as e:
            failed.append(dict(SAMPLE_ID=sid, error=str(e)))
    cohort = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    failed_df = pd.DataFrame(failed, columns=["SAMPLE_ID", "error"])
    return cohort, failed_df
