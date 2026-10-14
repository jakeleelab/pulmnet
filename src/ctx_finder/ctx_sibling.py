"""
Sibling-sample corroboration check for multi-sample patients.

When one sample from a patient is nominated for chromothripsis on a chromosome
and another sample from the SAME patient is not, this checks whether the
"missed" sample actually carries the same structural event -- just below the
classifier's threshold (e.g. due to a noisier library) -- rather than a
genuinely different or absent event.

Method: FACETS segment boundaries on a given chromosome are constrained to the
fixed positions of the IMPACT panel's probes, so a real, shared structural
event reproduces breakpoints at (near-)identical positions across independent
samples from the same tumor, while unrelated samples show no such relationship.
For each discordant (patient, chrom, flagged_sample, unflagged_sample) pair,
the check computes what fraction of the flagged sample's breakpoints have a matching
breakpoint (within `tol` bp) in the unflagged sample, then compare that match
fraction against an empirical null built from every OTHER sample in the cohort
on that same chromosome. A high percentile is quantitative support for a real,
shared-but-missed event; a percentile near 0 means the "missed" sample likely
just doesn't carry the event (a genuine negative, not a detection failure).
The report is informational and does not change any call.
"""
import numpy as np
import pandas as pd

import ctx_features as cf
import run_ctx_pipeline as rp


def get_breakpoints(man_df, sample_id, chrom, tcn_col="tcn.em", _cache=None):
    """
    Internal segment-boundary positions (loc.start of every segment after the
    first) for one sample's chromosome. These are the reproducible landmarks
    two samples sharing a real structural event will agree on.
    """
    key = (sample_id, str(chrom))
    if _cache is not None and key in _cache:
        return _cache[key]
    seg = rp.load_cncf(man_df, sample_id)
    seg["_c"] = seg["chrom"].map(cf._norm_chrom)
    g = seg.loc[seg["_c"].eq(cf._norm_chrom(chrom))].sort_values("loc.start")
    bps = g["loc.start"].to_numpy()[1:].astype(float) if len(g) >= 2 else np.array([])
    if _cache is not None:
        _cache[key] = bps
    return bps


def _match_frac(bps_a, bps_b, tol):
    """Fraction of bps_a with a matching breakpoint in bps_b within tol bp."""
    if len(bps_a) == 0:
        return np.nan, 0, 0
    if len(bps_b) == 0:
        return 0.0, 0, len(bps_a)
    matched = sum(1 for bp in bps_a if np.min(np.abs(bps_b - bp)) <= tol)
    return matched / len(bps_a), matched, len(bps_a)


def find_discordant_pairs(cohort_df, sample_col="SAMPLE_ID", chrom_col="chrom",
                          nominated_col="nominated_v7"):
    """
    Every (patient, chrom, flagged_sample, unflagged_sample) combination where
    the patient has >=2 samples and they disagree on that chromosome's call.
    One row per flagged/unflagged pair -- a patient with 2 flagged + 1
    unflagged sample yields 2 rows (the unflagged one vs each flagged one).
    """
    d = cohort_df[[sample_col, chrom_col, nominated_col]].copy()
    d["_patient"] = d[sample_col].astype(str).str.split("-T").str[0]
    d[chrom_col] = d[chrom_col].astype(str)

    pairs = []
    for (pid, chrom), g in d.groupby(["_patient", chrom_col]):
        if g[sample_col].nunique() < 2 or g[nominated_col].nunique() < 2:
            continue
        flagged = g.loc[g[nominated_col], sample_col].tolist()
        unflagged = g.loc[~g[nominated_col], sample_col].tolist()
        for f in flagged:
            for u in unflagged:
                pairs.append((pid, chrom, f, u))
    return pd.DataFrame(pairs, columns=["PATIENT_ID", chrom_col, "flagged_sample", "unflagged_sample"])


def sibling_corroboration_report(cohort_df, man_df, *, tol=200_000, n_null_cap=400,
                                 sample_col="SAMPLE_ID", chrom_col="chrom",
                                 nominated_col="nominated_v7"):
    """
    Run the breakpoint-concordance check on every discordant sibling pair in
    `cohort_df` (output of apply_classifier_v6/v7). Returns one row per pair
    with the true match fraction, an empirical null summary from every other
    cohort sample on that chromosome, and the true pair's percentile rank
    against that null. High percentile (~>=85 in validation) + nontrivial
    absolute match = rescue candidate; percentile near 0 = the unflagged
    sample likely doesn't carry the event at all (a real negative).
    """
    pairs = find_discordant_pairs(cohort_df, sample_col, chrom_col, nominated_col)
    if pairs.empty:
        return pairs.assign(n_bp_flagged=[], matched=[], match_frac=[],
                            null_n=[], null_mean=[], null_median=[], null_max=[],
                            percentile=[])

    all_samples = sorted(cohort_df[sample_col].astype(str).unique())
    pool = all_samples if len(all_samples) <= n_null_cap else all_samples[:n_null_cap]
    bp_cache = {}

    rows = []
    for _, row in pairs.iterrows():
        pid, chrom, flagged, unflagged = row["PATIENT_ID"], row[chrom_col], row["flagged_sample"], row["unflagged_sample"]
        bps_flagged = get_breakpoints(man_df, flagged, chrom, _cache=bp_cache)
        bps_unflagged = get_breakpoints(man_df, unflagged, chrom, _cache=bp_cache)
        true_frac, m, n = _match_frac(bps_flagged, bps_unflagged, tol)

        null_fracs = []
        for other in pool:
            if other in (flagged, unflagged):
                continue
            bps_other = get_breakpoints(man_df, other, chrom, _cache=bp_cache)
            f, _, _ = _match_frac(bps_flagged, bps_other, tol)
            if pd.notna(f):
                null_fracs.append(f)
        null_fracs = np.array(null_fracs)
        pctile = float((null_fracs < true_frac).mean() * 100) if len(null_fracs) and pd.notna(true_frac) else np.nan

        rows.append(dict(
            PATIENT_ID=pid, **{chrom_col: chrom},
            flagged_sample=flagged, unflagged_sample=unflagged,
            n_bp_flagged=n, matched=m,
            match_frac=round(true_frac, 3) if pd.notna(true_frac) else np.nan,
            null_n=len(null_fracs),
            null_mean=round(float(null_fracs.mean()), 3) if len(null_fracs) else np.nan,
            null_median=round(float(np.median(null_fracs)), 3) if len(null_fracs) else np.nan,
            null_max=round(float(null_fracs.max()), 3) if len(null_fracs) else np.nan,
            percentile=round(pctile, 1) if pd.notna(pctile) else np.nan,
        ))
    return pd.DataFrame(rows)
