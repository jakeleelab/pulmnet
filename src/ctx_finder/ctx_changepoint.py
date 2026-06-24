"""
Probe-level changepoint features (Fork 1) for chromothripsis-footprint detection.

Computed from probe-level cnlr + lorvar via PELT (ruptures), INDEPENDENT of the
FACETS segmentation. Replaces the FACETS-segment-derived oscillation chain with
a noise-aware one, and adds within-segment dispersion -- the genuinely new axis
that segment features cannot see.

Per chromosome it returns:
    cp_n_changepoints     # inferred changepoint count (noise-aware fragmentation)
    cp_osc_chain          # longest <=3-state run on INFERRED segment levels
    cp_within_dispersion  # residual probe scatter within inferred segments, in
                          # lorvar (sd) units -- catches spikes/instability FACETS flattened
    cp_osc_confidence     # chain length weighted by how much each level shift beat noise
    cp_n_levels           # number of distinct rounded inferred CN levels

Design / anti-flood guardrails:
  - cnlr is the segmentation signal; per-probe lorvar sets the noise scale.
  - PELT penalty beta = pen_const * sigma2 * log(n_probes), where sigma2 is the
    median lorvar on that chromosome. A level shift must beat the noise model to
    be accepted, so quiet diploids yield ~0 changepoints (unlike the rolling
    median, which had no noise model and flooded).
  - pen_const is the ONE knob; calibrate it on negative chromosomes first
    (see calibrate_penalty) before scoring anything.

Levels are derived on the tcn_estimate scale if available (purity/ploidy aware),
else on cnlr; the chain logic matches longest_le_k_distinct_run in ctx_features.
"""

import numpy as np
import pandas as pd

try:
    import ruptures as rpt
    _HAVE_RUPTURES = True
except Exception:
    _HAVE_RUPTURES = False


def _norm_chrom(x):
    return str(x).strip().replace("chr", "").replace("CHR", "")


def _chrom_order(x):
    s = _norm_chrom(x)
    if s in ("X", "x", "23"): return 23
    if s in ("Y", "y", "24"): return 24
    if s in ("M", "MT", "m", "mt", "25"): return 25
    try: return int(s)
    except ValueError: return 99


def _longest_le_k_distinct_run(states, k):
    """Longest contiguous window with <= k distinct values, on a run-length-
    encoded integer sequence. (Mirror of ctx_features.longest_le_k_distinct_run.)"""
    s = [v for v in states if v is not None and not (isinstance(v, float) and np.isnan(v))]
    if not s:
        return 0
    collapsed = [s[0]]
    for v in s[1:]:
        if v != collapsed[-1]:
            collapsed.append(v)
    best, left, counts = 0, 0, {}
    for right in range(len(collapsed)):
        counts[collapsed[right]] = counts.get(collapsed[right], 0) + 1
        while len(counts) > k:
            counts[collapsed[left]] -= 1
            if counts[collapsed[left]] == 0:
                del counts[collapsed[left]]
            left += 1
        best = max(best, right - left + 1)
    return int(best)


def _pelt_breakpoints(x, sigma2, pen_const, min_size=5, jump=1):
    """
    PELT changepoint detection on 1-D signal x with an L2 (mean-shift) cost and
    a noise-scaled penalty. Returns the list of breakpoint indices (end of each
    segment, ruptures convention; last element == len(x)).
    """
    n = len(x)
    if n < 2 * min_size:
        return [n]
    beta = pen_const * float(sigma2) * np.log(max(n, 2))
    if _HAVE_RUPTURES:
        algo = rpt.Pelt(model="l2", min_size=min_size, jump=jump).fit(x.reshape(-1, 1))
        return algo.predict(pen=beta)
    # ---- fallback exact DP (only used where ruptures is unavailable, e.g. CI) ----
    # O(n^2) penalized least-squares segmentation; fine for test-sized inputs.
    csum = np.concatenate([[0.0], np.cumsum(x)])
    csum2 = np.concatenate([[0.0], np.cumsum(x * x)])
    def seg_cost(a, b):  # cost of x[a:b]
        m = b - a
        s = csum[b] - csum[a]
        s2 = csum2[b] - csum2[a]
        return s2 - (s * s) / m
    F = np.full(n + 1, np.inf); F[0] = -beta
    prev = np.zeros(n + 1, dtype=int)
    for b in range(min_size, n + 1):
        for a in range(0, b - min_size + 1):
            c = F[a] + seg_cost(a, b) + beta
            if c < F[b]:
                F[b] = c; prev[b] = a
    bkps, b = [], n
    while b > 0:
        bkps.append(b); b = prev[b]
    return sorted(bkps)


def _cnlr_to_cn(mean_cnlr, purity, ploidy):
    """Back-transform a segment-mean cnlr to estimated total copy number using
    the same formula as add_tcn_estimate. Returns NaN if purity missing."""
    if purity is None or not np.isfinite(purity) or purity <= 0:
        # no purity: fall back to a monotone proxy (2 * 2**cnlr), good enough for
        # rounding to distinct levels
        return 2.0 * (2.0 ** mean_cnlr)
    return ((purity * ploidy + (1 - purity) * 2) * (2.0 ** mean_cnlr) - (1 - purity) * 2) / purity


def probe_changepoint_features(pts_chr, signal_col="cnlr",
                               cnlr_col="cnlr", lorvar_col="lorvar",
                               purity=None, ploidy=2.0,
                               pen_const=3.0, min_size=5, k_states=3):
    """
    Per-chromosome changepoint features from a probe-level table (one chromosome).

    Segmentation is done in CNLR space (the native space of lorvar). The noise
    scale is the LARGER of median(lorvar) and a robust diff-based estimate, so
    the penalty is never set below the actual probe-to-probe scatter (this is
    what prevents the flood -- lorvar alone understates real scatter ~3x).

    Integer CN levels for the oscillation chain are obtained by back-transforming
    each segment's mean cnlr to copy number via purity/ploidy, then rounding.
    """
    out = dict(cp_n_changepoints=0, cp_osc_chain=0, cp_within_dispersion=0.0,
               cp_osc_confidence=0.0, cp_n_levels=0, cp_n_probes=int(pts_chr.shape[0]))
    d = pts_chr.sort_values("maploc") if "maploc" in pts_chr.columns else pts_chr
    cnlr = pd.to_numeric(d.get(cnlr_col), errors="coerce").to_numpy(dtype=float)
    ok = np.isfinite(cnlr)
    cnlr = cnlr[ok]
    n = len(cnlr)
    if n < 2 * min_size:
        return out

    # noise scale: max of modeled (lorvar) and empirical (robust diff) variance
    if lorvar_col in d.columns and d[lorvar_col].notna().any():
        sigma2_model = float(np.nanmedian(pd.to_numeric(d[lorvar_col], errors="coerce")))
    else:
        sigma2_model = 0.0
    diff_sd = 1.4826 * np.median(np.abs(np.diff(cnlr))) / np.sqrt(2) if n > 1 else 1.0
    sigma2_emp = float(diff_sd ** 2)
    sigma2 = max(sigma2_model, sigma2_emp, 1e-4)

    bkps = _pelt_breakpoints(cnlr, sigma2, pen_const, min_size=min_size)
    starts = [0] + bkps[:-1]
    seg_mean_cnlr, seg_resid_mads = [], []
    for a, b in zip(starts, bkps):
        seg = cnlr[a:b]
        if len(seg) == 0:
            continue
        mu = np.median(seg)
        seg_mean_cnlr.append(mu)
        seg_resid_mads.append(1.4826 * np.median(np.abs(seg - mu)))

    out["cp_n_changepoints"] = max(len(seg_mean_cnlr) - 1, 0)
    # back-transform segment means to CN, round to integer levels
    cn_levels = [int(round(_cnlr_to_cn(m, purity, ploidy))) for m in seg_mean_cnlr]
    cn_levels = [max(c, 0) for c in cn_levels]
    out["cp_n_levels"] = len(set(cn_levels))
    out["cp_osc_chain"] = _longest_le_k_distinct_run(cn_levels, k_states)

    sd_noise = np.sqrt(sigma2)
    if seg_resid_mads:
        out["cp_within_dispersion"] = float(np.median(seg_resid_mads) / sd_noise)

    if len(seg_mean_cnlr) >= 2:
        shifts = np.abs(np.diff(seg_mean_cnlr)) / sd_noise
        sep = float(np.clip(np.median(shifts) / 3.0, 0, 1))
        out["cp_osc_confidence"] = out["cp_osc_chain"] * sep
    return out



def add_changepoint_features(cohort_df, man_df, tablesDir, prepare_fn,
                             pen_const=3.0, signal_col="tcn_estimate",
                             sample_col="SAMPLE_ID"):
    """
    Compute changepoint features for every (sample, chromosome) in cohort_df and
    merge them on [SAMPLE_ID, chrom].

    prepare_fn(man_df, sample_id, tablesDir) -> (seg_df, jointseg_df, path)
    (pass run_ctx_pipeline.prepare_cn_data_for_sample).
    """
    rows = []
    for sid, _ in cohort_df.groupby(sample_col):
        try:
            _, jointseg_df, _ = prepare_fn(man_df, sid, tablesDir)
        except Exception:
            continue
        row = man_df.loc[man_df[sample_col].astype(str).eq(str(sid))]
        purity = float(row["purity"].iloc[0]) if row.shape[0] and pd.notna(row["purity"].iloc[0]) else None
        ploidy = float(row["ploidy"].iloc[0]) if row.shape[0] and pd.notna(row["ploidy"].iloc[0]) else 2.0
        js = jointseg_df.copy()
        js["_c"] = js["chrom"].map(_norm_chrom)
        for c, sub in js.groupby("_c"):
            f = probe_changepoint_features(sub, purity=purity, ploidy=ploidy, pen_const=pen_const)
            f[sample_col] = sid
            f["chrom"] = c
            rows.append(f)
    cp = pd.DataFrame(rows)
    merged = cohort_df.copy()
    merged["_c"] = merged["chrom"].map(_norm_chrom)
    cp["_c"] = cp["chrom"].map(_norm_chrom)
    merged = merged.merge(cp.drop(columns=["chrom"]), on=[sample_col, "_c"], how="left")
    return merged.drop(columns=["_c"])


def calibrate_penalty(cohort_df, man_df, tablesDir, prepare_fn,
                      pen_grid=(1.0, 2.0, 3.0, 5.0, 8.0),
                      truth_col="truth_positive_chrom", sample_col="SAMPLE_ID",
                      n_sample_chroms=400, seed=0):
    """
    ANTI-FLOOD CALIBRATION. For each pen_const, report the median/mean
    cp_n_changepoints on NEGATIVE chromosomes vs POSITIVE chromosomes, on a
    random subset of samples (for speed). A good penalty keeps negatives near
    0-1 changepoints while positives stay elevated. Pick the smallest pen_const
    that keeps negative median <= ~1.
    """
    rng = np.random.default_rng(seed)
    sids = cohort_df[sample_col].drop_duplicates().tolist()
    rng.shuffle(sids)
    rows = []
    cache = {}
    for pen in pen_grid:
        neg_cp, pos_cp = [], []
        seen = 0
        for sid in sids:
            if sid not in cache:
                try:
                    _, js, _ = prepare_fn(man_df, sid, tablesDir)
                except Exception:
                    cache[sid] = None
                    continue
                js = js.copy(); js["_c"] = js["chrom"].map(_norm_chrom)
                cache[sid] = js
            js = cache[sid]
            if js is None:
                continue
            truth_chr = set(cohort_df.loc[(cohort_df[sample_col] == sid) &
                                          cohort_df[truth_col].astype(bool), "chrom"].map(_norm_chrom))
            prow = man_df.loc[man_df[sample_col].astype(str).eq(str(sid))]
            pur = float(prow["purity"].iloc[0]) if prow.shape[0] and pd.notna(prow["purity"].iloc[0]) else None
            plo = float(prow["ploidy"].iloc[0]) if prow.shape[0] and pd.notna(prow["ploidy"].iloc[0]) else 2.0
            for c, sub in js.groupby("_c"):
                f = probe_changepoint_features(sub, purity=pur, ploidy=plo, pen_const=pen)
                (pos_cp if c in truth_chr else neg_cp).append(f["cp_n_changepoints"])
                seen += 1
            if seen >= n_sample_chroms:
                break
        rows.append(dict(pen_const=pen,
                         neg_median=np.median(neg_cp) if neg_cp else np.nan,
                         neg_mean=np.mean(neg_cp) if neg_cp else np.nan,
                         neg_p90=np.percentile(neg_cp, 90) if neg_cp else np.nan,
                         pos_median=np.median(pos_cp) if pos_cp else np.nan,
                         pos_mean=np.mean(pos_cp) if pos_cp else np.nan,
                         n_neg=len(neg_cp), n_pos=len(pos_cp)))
    return pd.DataFrame(rows)
