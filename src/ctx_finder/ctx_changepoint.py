"""
Probe-level changepoint re-segmentation, independent of the FACETS segmentation.

Per chromosome, probe log-ratios (cnlr) are re-segmented with PELT (ruptures; L2
cost) and the oscillating stretch is recomputed on the inferred segment levels:

    cp_n_changepoints     # inferred changepoint count
    cp_osc_chain          # oscillating stretch (longest <=3-state run) on inferred levels
    cp_osc_confidence     # cp_osc_chain weighted by how clearly the level shifts exceed noise

PELT penalty beta = pen_const * sigma2 * log(n_probes), where sigma2 is the larger of
the median FACETS log-ratio variance (lorvar) and the squared robust successive-
difference SD of the log-ratios, so a level shift must beat the probe noise to be
accepted.
"""

import numpy as np
import pandas as pd

try:
    import ruptures as rpt
    _HAVE_RUPTURES = True
except Exception:
    _HAVE_RUPTURES = False


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
    # ---- fallback exact DP (only used where ruptures is unavailable) ----
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
    """Back-transform a segment log-ratio to estimated total copy number using
    purity and ploidy (same formula as run_ctx_pipeline.add_tcn_estimate)."""
    if purity is None or not np.isfinite(purity) or purity <= 0:
        # no purity: fall back to a monotone proxy (2 * 2**cnlr), good enough for
        # rounding to distinct levels
        return 2.0 * (2.0 ** mean_cnlr)
    return ((purity * ploidy + (1 - purity) * 2) * (2.0 ** mean_cnlr) - (1 - purity) * 2) / purity


def probe_changepoint_features(pts_chr, cnlr_col="cnlr", lorvar_col="lorvar",
                               purity=None, ploidy=2.0,
                               pen_const=8.0, min_size=5, k_states=3):
    """
    Changepoint features from a probe-level table of one chromosome.

    Segmentation is done in log-ratio space. Integer copy-number levels for the
    oscillating stretch are obtained by converting each segment's median log-ratio to
    copy number with purity/ploidy and rounding.
    """
    out = dict(cp_n_changepoints=0, cp_osc_chain=0, cp_osc_confidence=0.0)
    d = pts_chr.sort_values("maploc") if "maploc" in pts_chr.columns else pts_chr
    cnlr = pd.to_numeric(d.get(cnlr_col), errors="coerce").to_numpy(dtype=float)
    cnlr = cnlr[np.isfinite(cnlr)]
    n = len(cnlr)
    if n < 2 * min_size:
        return out

    # noise scale: max of modeled (lorvar) and empirical (robust diff) variance
    if lorvar_col in d.columns and d[lorvar_col].notna().any():
        sigma2_model = float(np.nanmedian(pd.to_numeric(d[lorvar_col], errors="coerce")))
    else:
        sigma2_model = 0.0
    diff_sd = 1.4826 * np.median(np.abs(np.diff(cnlr))) / np.sqrt(2) if n > 1 else 1.0
    sigma2 = max(sigma2_model, float(diff_sd ** 2), 1e-4)

    bkps = _pelt_breakpoints(cnlr, sigma2, pen_const, min_size=min_size)
    starts = [0] + bkps[:-1]
    seg_median_cnlr = [np.median(cnlr[a:b]) for a, b in zip(starts, bkps) if b > a]

    out["cp_n_changepoints"] = max(len(seg_median_cnlr) - 1, 0)
    cn_levels = [max(int(round(_cnlr_to_cn(m, purity, ploidy))), 0) for m in seg_median_cnlr]
    out["cp_osc_chain"] = _longest_le_k_distinct_run(cn_levels, k_states)

    if len(seg_median_cnlr) >= 2:
        shifts = np.abs(np.diff(seg_median_cnlr)) / np.sqrt(sigma2)
        sep = float(np.clip(np.median(shifts) / 3.0, 0, 1))
        out["cp_osc_confidence"] = out["cp_osc_chain"] * sep
    return out
