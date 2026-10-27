"""
Chromosome-level features, classifiers and post-hoc filters for chromothripsis
detection on sparse FACETS (targeted-panel) copy-number data.

Inputs per sample:

    seg_df       : FACETS segment table (cncf); columns include chrom, loc.start,
                   loc.end, num.mark, cnlr.median, mafR, tcn.em, lcn.em, cf.em, seg
    jointseg_df  : FACETS probe-level table; columns include chrom, maploc, cnlr,
                   lorvar, vafT, het, seg (tcn_estimate is added internally)
"""

import numpy as np
import pandas as pd


# hg19 centromere + immediate pericentromeric-heterochromatin region, Mb.
# Approximate (UCSC cytoBand/gap track order of magnitude, generously padded
# on the flanks) -- used only to reject candidate "shards" that fall in this
# low-mappability, probe-desert zone, not for any coordinate-precise purpose.
# See compute_baseline_shard_features, "centromeric-shard artifact" filter.
CENTROMERE_HG19_MB = {
    "1": (121.5, 128.9), "2": (92.3, 95.3), "3": (90.5, 93.5), "4": (49.7, 52.7),
    "5": (46.5, 50.7), "6": (58.6, 61.6), "7": (58.1, 61.1), "8": (43.0, 46.8),
    "9": (47.0, 65.0), "10": (39.8, 42.7), "11": (51.6, 55.8), "12": (35.5, 37.8),
    "13": (17.7, 19.0), "14": (17.1, 19.0), "15": (19.0, 20.7), "16": (36.6, 38.9),
    "17": (24.0, 26.5), "18": (15.4, 19.0), "19": (24.4, 28.6), "20": (26.4, 29.7),
    "21": (10.9, 14.3), "22": (12.9, 17.9), "23": (58.6, 61.6),
}
# The chr8 and chr10 windows extend slightly beyond the centromere to cover small,
# sparse-probe segments (chr8:43.09-43.10 Mb, chr10:42.38-42.60 Mb) that recur at
# identical coordinates across unrelated samples, i.e. a panel artifact.


# ----------------------------------------------------------------------
# small helpers
# ----------------------------------------------------------------------

def _norm_chrom(x):
    return str(x).strip().replace("chr", "").replace("CHR", "")


def _chrom_order(x):
    s = _norm_chrom(x)
    if s in ("X", "x", "23"):
        return 23
    if s in ("Y", "y", "24"):
        return 24
    if s in ("M", "MT", "m", "mt", "25"):
        return 25
    try:
        return int(s)
    except ValueError:
        return 99


def _ordered_segments(seg_chr):
    """cncf rows for one chromosome, ordered by genomic start."""
    return seg_chr.sort_values(["loc.start", "loc.end"]).reset_index(drop=True)


def _transitions(state_seq):
    """Number of adjacent changes in a sequence of (hashable) states."""
    if len(state_seq) <= 1:
        return 0
    return int(np.sum([state_seq[i] != state_seq[i - 1] for i in range(1, len(state_seq))]))


# ----------------------------------------------------------------------
# per-(sample, chromosome) feature computation
# ----------------------------------------------------------------------


def probe_reseg_features(tcn_probe, k_states=3, noise=None, min_seg=5):
    """
    FACETS-independent oscillation features from the probe-level estimated total
    copy number of ONE chromosome (ordered by position).

    Robustly self-segments via a rolling median, quantizes to integer CN states,
    run-length-encodes, and measures the longest <=k-state oscillation chain on
    THAT sequence -- so it recovers oscillation that FACETS merged away. Also
    returns a heavy-tail "amplification spike" fraction that captures focal amps
    a coarse 2-segment FACETS call flattens (e.g. C22 chr3).
    """
    x = np.asarray(pd.to_numeric(pd.Series(tcn_probe), errors="coerce").dropna(), dtype=float)
    n = len(x)
    out = dict(probe_osc_chain=0, probe_n_state_changes=0,
               probe_amp_spike_frac=0.0, probe_reseg_n_states=0)
    if n < min_seg:
        return out
    # robust local noise scale if not provided
    if noise is None or not np.isfinite(noise) or noise <= 0:
        noise = 1.4826 * np.median(np.abs(np.diff(x))) / np.sqrt(2) if n > 1 else 1.0
        noise = max(noise, 1e-3)
    # rolling median to denoise (window scaled to chromosome size, odd)
    w = int(max(5, min(51, (n // 30) | 1)))
    s = pd.Series(x)
    rm = s.rolling(w, center=True, min_periods=max(3, w // 3)).median().to_numpy()
    rm = np.where(np.isfinite(rm), rm, x)
    # quantize to integer CN states (clip at 0)
    q = np.clip(np.rint(rm), 0, None).astype(int)
    # run-length encode
    rle = [q[0]]
    for v in q[1:]:
        if v != rle[-1]:
            rle.append(v)
    out["probe_n_state_changes"] = int(len(rle) - 1)
    out["probe_reseg_n_states"] = int(len(set(rle)))
    out["probe_osc_chain"] = int(longest_le_k_distinct_run(rle, k_states))
    # heavy-tail amplification spikes vs robust baseline (focal amps FACETS merged)
    base = np.median(x)
    mad = 1.4826 * np.median(np.abs(x - base))
    mad = max(mad, 1e-3)
    out["probe_amp_spike_frac"] = float(np.mean(x > base + 4.0 * mad))
    return out


def _probe_dispersion(values):
    """Total and within-segment dispersion of position-ordered probe values.

    `values` must already be position-sorted. Returns sd/MAD of the whole cloud
    (signal + noise) and of the lag-1 successive differences (within-segment
    noise, robust to breakpoints). MAD versions are scaled to a Gaussian sd.
    """
    x = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)
    x = x[~np.isnan(x)]
    if x.size < 2:
        return dict(sd_total=np.nan, mad_total=np.nan, sd_within=np.nan, mad_within=np.nan)
    dif = np.diff(x)
    return dict(
        sd_total=float(np.std(x, ddof=1)),
        mad_total=float(1.4826 * np.median(np.abs(x - np.median(x)))),
        sd_within=float(np.std(dif, ddof=1) / np.sqrt(2)),
        mad_within=float(1.4826 * np.median(np.abs(dif - np.median(dif))) / np.sqrt(2)),
    )


def chrom_features_one_sample(seg_df, jointseg_df,
                              tcn_col="tcn.em", lcn_col="lcn.em",
                              point_col="tcn_estimate"):
    """
    Return a DataFrame with one row per chromosome of NEW candidate features.

    These are intended to *complement* your existing cohort_df features
    (n_segments, n_segment_state_changes, point_sd, oscillation_fraction, ...),
    not replace them. Merge on chrom.
    """
    seg = seg_df.copy()
    pts = jointseg_df.copy()
    seg["_c"] = seg["chrom"].map(_norm_chrom)
    pts["_c"] = pts["chrom"].map(_norm_chrom)

    rows = []
    for c, seg_chr in seg.groupby("_c"):
        seg_chr = _ordered_segments(seg_chr)
        pts_chr = pts.loc[pts["_c"].eq(c)].sort_values("maploc")
        _span = focal_span_features(seg_chr, pts_chr, tcn_col=tcn_col, lcn_col=lcn_col, k_states=3)

        n_probes = int(pts_chr.shape[0])
        het = pts_chr["het"] if "het" in pts_chr.columns else pd.Series([], dtype=float)
        n_het = int(np.nansum(het.astype(float))) if len(het) else 0

        # ---- total-CN state structure (T, L) ordered along the chromosome ----
        T = pd.to_numeric(seg_chr.get(tcn_col), errors="coerce")
        L = pd.to_numeric(seg_chr.get(lcn_col), errors="coerce")
        states = list(zip(T.tolist(), L.tolist()))
        A = _transitions(states)                       # total-CN state changes
        D = len({s for s in states if not (pd.isna(s[0]) and pd.isna(s[1]))})

        # oscillation *excess*: transitions beyond a monotone visit of D states.
        # staircase -> ~0 ; revisiting few states -> ->1
        osc_excess = (A - max(D - 1, 0)) / A if A > 0 else 0.0

        # ShatterSeek-style longest contiguous chain confined to <=K CN states
        # (rounded total CN). Discriminates oscillation (long chain) from
        # monotone/BFB staircases (short chain).
        int_states = [int(round(t)) for t in T.tolist() if not pd.isna(t)]
        longest_osc_chain_k2 = longest_le_k_distinct_run(int_states, 2)
        longest_osc_chain_k3 = longest_le_k_distinct_run(int_states, 3)
        longest_osc_chain_k4 = longest_le_k_distinct_run(int_states, 4)
        longest_osc_chain_k5 = longest_le_k_distinct_run(int_states, 5)        
        # CN amplitude spanned by the chromosome's total-CN states (max-min of
        # rounded tcn.em). A 2<->1 allelic-flip wobble has amplitude 1; a real
        # oscillation 2<->4<->2<->5 spans >=2.
        chain_cn_amplitude = (max(int_states) - min(int_states)) if int_states else 0

        # FACETS-independent probe-level re-segmentation (catches merged signal)
        _pn = float(np.sqrt(np.nanmedian(pts_chr["lorvar"]))) if "lorvar" in pts_chr.columns and pts_chr["lorvar"].notna().any() else None
        _pr = probe_reseg_features(pts_chr[point_col] if point_col in pts_chr.columns else [], k_states=3, noise=_pn)
        # within-chromosome probe dispersion: the plotted dots (tcn_estimate) and
        # the stable log-ratio noise scale (cnlr)
        _dsp      = _probe_dispersion(pts_chr[point_col]) if point_col in pts_chr.columns else _probe_dispersion([])
        _dsp_cnlr = _probe_dispersion(pts_chr["cnlr"])    if "cnlr"   in pts_chr.columns else _probe_dispersion([])

        # two-state dominance among transitions (CT oscillates between ~2 states)
        if A > 0:
            pair_counts = {}
            for i in range(1, len(states)):
                if states[i] != states[i - 1]:
                    key = tuple(sorted([states[i], states[i - 1]], key=str))
                    pair_counts[key] = pair_counts.get(key, 0) + 1
            top = max(pair_counts.values()) if pair_counts else 0
            two_state_dom = top / A
        else:
            two_state_dom = 0.0

        # run-length structure of the total-CN state sequence
        run_lengths = []
        if states:
            cur = 1
            for i in range(1, len(states)):
                if states[i] == states[i - 1]:
                    cur += 1
                else:
                    run_lengths.append(cur)
                    cur = 1
            run_lengths.append(cur)
        median_run = float(np.median(run_lengths)) if run_lengths else np.nan
        n_singleton_runs = int(np.sum(np.array(run_lengths) == 1)) if run_lengths else 0

        # ---- allelic / LOH structure ----
        loh = (L == 0).astype(float)                   # 1 where minor CN == 0
        n_allelic_changes = _transitions(L.tolist())
        loh_het_interleave = _transitions(loh.tolist())  # LOH <-> het switches
        frac_loh = float(np.nanmean(loh)) if len(loh) else np.nan

        # probe-level allelic dispersion among het sites (CT scrambles VAF)
        if "vafT" in pts_chr.columns and n_het > 0:
            vaf = pd.to_numeric(pts_chr.loc[het.astype(bool), "vafT"], errors="coerce")
            # fold around 0.5 so 0.2 and 0.8 are equally "imbalanced"
            vaf_imb = np.abs(vaf - 0.5)
            vaf_imb_sd = float(np.nanstd(vaf_imb)) if vaf.notna().sum() > 1 else np.nan
        else:
            vaf_imb_sd = np.nan

        # ---- purity-robust amplitude (raw cnlr SNR) ----
        # spread of segment-level cnlr.median vs typical probe noise (lorvar).
        seg_cnlr = pd.to_numeric(seg_chr.get("cnlr.median"), errors="coerce").dropna()
        if "lorvar" in pts_chr.columns and pts_chr["lorvar"].notna().any():
            noise = float(np.sqrt(np.nanmedian(pts_chr["lorvar"])))
            probe_noise = float(np.nanmedian(pts_chr["lorvar"]))
        else:
            noise = np.nan
            probe_noise = np.nan
        cnlr_spread = float(seg_cnlr.max() - seg_cnlr.min()) if len(seg_cnlr) else np.nan
        cnlr_snr = (cnlr_spread / noise) if (noise and noise > 0) else np.nan

        # number of CN levels that are *resolved* above noise (separated by >2*noise)
        if len(seg_cnlr) and noise and noise > 0:
            levels = np.sort(seg_cnlr.unique())
            resolved = 1
            last = levels[0]
            for v in levels[1:]:
                if v - last > 2 * noise:
                    resolved += 1
                    last = v
            n_resolved_levels = int(resolved)
        else:
            n_resolved_levels = np.nan

        rows.append(dict(
            chrom=c,
            _chrom_order=_chrom_order(c),
            n_probes=n_probes,
            n_het_probes=n_het,
            # fragmentation normalized by probe support
            state_changes_per_probe=(A / n_probes) if n_probes else np.nan,
            segments_per_probe=(seg_chr.shape[0] / n_probes) if n_probes else np.nan,
            # oscillation quality
            n_distinct_cn_states=D,
            oscillation_excess=osc_excess,
            longest_osc_chain_k2=longest_osc_chain_k2,
            longest_osc_chain_k3=longest_osc_chain_k3,
            longest_osc_chain_k4=longest_osc_chain_k4,
            longest_osc_chain_k5=longest_osc_chain_k5,
            dot_sd_total=_dsp["sd_total"],
            dot_mad_total=_dsp["mad_total"],
            dot_sd_within=_dsp["sd_within"],
            dot_mad_within=_dsp["mad_within"],
            cnlr_sd_within=_dsp_cnlr["sd_within"],
            cnlr_mad_within=_dsp_cnlr["mad_within"],            
            chain_cn_amplitude=chain_cn_amplitude,
            probe_osc_chain=_pr["probe_osc_chain"],
            probe_n_state_changes=_pr["probe_n_state_changes"],
            probe_reseg_n_states=_pr["probe_reseg_n_states"],
            probe_amp_spike_frac=_pr["probe_amp_spike_frac"],
            two_state_dominance=two_state_dom,
            median_run_length=median_run,
            n_singleton_runs=n_singleton_runs,
            # allelic / LOH
            n_allelic_state_changes=n_allelic_changes,
            loh_het_interleave=loh_het_interleave,
            frac_loh=frac_loh,
            vaf_imbalance_sd=vaf_imb_sd,
            # purity-robust amplitude
            cnlr_snr=cnlr_snr,
            n_resolved_levels=n_resolved_levels,
            probe_noise=probe_noise,
            **_span,
        ))

    return pd.DataFrame(rows).sort_values("_chrom_order").reset_index(drop=True)


def longest_le_k_distinct_run_bounds(states, k):
    collapsed = []
    for i, v in enumerate(states):
        if v is None or pd.isna(v):
            continue
        v = int(round(v))
        if len(collapsed) == 0 or collapsed[-1]["state"] != v:
            collapsed.append({"state": v, "i0": i, "i1": i})
        else:
            collapsed[-1]["i1"] = i
    if len(collapsed) == 0:
        return dict(chain_len=0, i0=np.nan, i1=np.nan)
    best = (0, 0, 0)
    left = 0
    counts = {}
    for right, r in enumerate(collapsed):
        counts[r["state"]] = counts.get(r["state"], 0) + 1
        while len(counts) > k:
            counts[collapsed[left]["state"]] -= 1
            if counts[collapsed[left]["state"]] == 0:
                del counts[collapsed[left]["state"]]
            left += 1
        cur = right - left + 1
        if cur > best[0]:
            best = (cur, left, right)
    _, l, r = best
    return dict(chain_len=int(best[0]), i0=int(collapsed[l]["i0"]), i1=int(collapsed[r]["i1"]))


def _chain_restricted_min_cfem_spike(seg_chrom, tcn_col="tcn.em", chain_k=3):
    """
    Within the segment run that forms the sample's actual longest_osc_chain_k{chain_k}
    (the same window `longest_le_k_distinct_run_bounds` scores), find any 'spike'
    segment -- one whose rounded tcn.em state differs from BOTH immediate neighbors,
    i.e. an isolated blip that creates a false transition in/out -- and return the
    minimum FACETS `cf.em` (clonal-cluster confidence) among those spikes.

    Restricting to the chain window (not just any spike anywhere on the chromosome)
    matters: a chromosome can carry one weak, low-confidence segment alongside many
    other strong, well-supported oscillating segments elsewhere -- vetoing the whole
    chromosome for its single weakest link anywhere would wrongly kill real, dramatic
    chromothripsis calls. A chromosome whose only evidence is the weak spike is
    still flagged.

    Returns (min_cfem, n_spikes_used); min_cfem is NaN if no scoreable spike exists
    inside the chain window -- callers should treat NaN as "no evidence against
    this row" and not veto on it.
    """
    s = _ordered_segments(seg_chrom)
    T = pd.to_numeric(s.get(tcn_col), errors="coerce")
    valid = T.notna()
    s = s.loc[valid].reset_index(drop=True)
    T = T.loc[valid].reset_index(drop=True)
    if len(s) < 3:
        return np.nan, 0
    int_states = [int(round(t)) for t in T.tolist()]
    bounds = longest_le_k_distinct_run_bounds(int_states, chain_k)
    i0, i1 = bounds["i0"], bounds["i1"]
    if pd.isna(i0):
        return np.nan, 0
    states = np.array(int_states)
    cfems = s["cf.em"].to_numpy() if "cf.em" in s.columns else np.full(len(s), np.nan)
    vals = []
    for i in range(max(1, int(i0)), min(len(s) - 1, int(i1) + 1)):
        if states[i] != states[i - 1] and states[i] != states[i + 1] and pd.notna(cfems[i]):
            vals.append(cfems[i])
    return (min(vals) if vals else np.nan), len(vals)


def focal_span_features(seg_chr, pts_chr, tcn_col="tcn.em", lcn_col="lcn.em", k_states=3):
    seg_chr = _ordered_segments(seg_chr)
    starts = pd.to_numeric(seg_chr["loc.start"], errors="coerce").to_numpy(dtype=float)
    ends = pd.to_numeric(seg_chr["loc.end"], errors="coerce").to_numpy(dtype=float)
    lengths = np.clip(ends - starts, 0, None)
    T = pd.to_numeric(seg_chr.get(tcn_col), errors="coerce")
    L = pd.to_numeric(seg_chr.get(lcn_col), errors="coerce")
    t_states = [int(round(x)) if pd.notna(x) else np.nan for x in T]
    tl_states = [(int(round(t)), int(round(l))) if pd.notna(t) and pd.notna(l) else None for t, l in zip(T, L)]
    chr_start = np.nanmin(starts) if len(starts) else np.nan
    chr_end = np.nanmax(ends) if len(ends) else np.nan
    chr_span = chr_end - chr_start if np.isfinite(chr_start) and np.isfinite(chr_end) else np.nan
    chr_len_sum = np.nansum(lengths)
    b = longest_le_k_distinct_run_bounds(t_states, k_states)
    if np.isfinite(b["i0"]) and np.isfinite(b["i1"]):
        i0, i1 = int(b["i0"]), int(b["i1"])
        chain_start, chain_end = starts[i0], ends[i1]
        chain_span = chain_end - chain_start
        chain_len_sum = np.nansum(lengths[i0:(i1 + 1)])
        chain_n_segments = i1 - i0 + 1
        pts_map = pd.to_numeric(pts_chr["maploc"], errors="coerce") if "maploc" in pts_chr.columns else pd.Series([], dtype=float)
        chain_n_probes = int(((pts_map >= chain_start) & (pts_map <= chain_end)).sum())
    else:
        chain_span = np.nan
        chain_len_sum = np.nan
        chain_n_segments = 0
        chain_n_probes = 0
    valid_t = [x for x in t_states if pd.notna(x)]
    modal_t = pd.Series(valid_t).mode().iloc[0] if len(valid_t) else np.nan
    altered = np.array([pd.notna(x) and x != modal_t for x in t_states], dtype=bool)
    altered_union = float(np.nansum(lengths[altered])) if len(altered) else np.nan
    if altered.any():
        altered_span = float(np.nanmax(ends[altered]) - np.nanmin(starts[altered]))
        runs = []
        cur = 0
        for a, l in zip(altered, lengths):
            if a:
                cur += l
            elif cur > 0:
                runs.append(cur)
                cur = 0
        if cur > 0:
            runs.append(cur)
        n_altered_islands = len(runs)
        largest_altered_island = float(np.nanmax(runs)) if runs else 0.0
    else:
        altered_span = 0.0
        n_altered_islands = 0
        largest_altered_island = 0.0
    change_pos = []
    for i in range(1, len(tl_states)):
        if tl_states[i] != tl_states[i - 1]:
            change_pos.append(starts[i])
    if len(change_pos) >= 2:
        state_change_span = float(np.nanmax(change_pos) - np.nanmin(change_pos))
    else:
        state_change_span = 0.0
    chain_transitions = max(b["chain_len"] - 1, 0)
    chain_span_mb = chain_span / 1e6 if pd.notna(chain_span) else np.nan
    return dict(
        chrom_span_mb=chr_span / 1e6 if pd.notna(chr_span) else np.nan,
        chrom_seg_len_mb=chr_len_sum / 1e6 if pd.notna(chr_len_sum) else np.nan,
        osc_chain_span_mb=chain_span_mb,
        osc_chain_len_sum_mb=chain_len_sum / 1e6 if pd.notna(chain_len_sum) else np.nan,
        osc_chain_frac_chrom=chain_span / chr_span if chr_span and chr_span > 0 else np.nan,
        osc_chain_len_frac_chrom=chain_len_sum / chr_len_sum if chr_len_sum and chr_len_sum > 0 else np.nan,
        osc_chain_n_segments=chain_n_segments,
        osc_chain_frac_segments=chain_n_segments / len(seg_chr) if len(seg_chr) else np.nan,
        osc_chain_n_probes=chain_n_probes,
        osc_chain_frac_probes=chain_n_probes / len(pts_chr) if len(pts_chr) else np.nan,
        osc_chain_density_per_mb=b["chain_len"] / chain_span_mb if chain_span_mb and chain_span_mb > 0 else np.nan,
        osc_transition_density_per_mb=chain_transitions / chain_span_mb if chain_span_mb and chain_span_mb > 0 else np.nan,
        bp_per_osc_transition_mb=chain_span_mb / chain_transitions if chain_transitions > 0 else np.nan,
        altered_union_mb=altered_union / 1e6 if pd.notna(altered_union) else np.nan,
        altered_union_frac=altered_union / chr_len_sum if chr_len_sum and chr_len_sum > 0 else np.nan,
        altered_span_mb=altered_span / 1e6,
        altered_span_frac=altered_span / chr_span if chr_span and chr_span > 0 else np.nan,
        n_altered_islands=n_altered_islands,
        largest_altered_island_mb=largest_altered_island / 1e6,
        largest_altered_island_frac=largest_altered_island / chr_len_sum if chr_len_sum and chr_len_sum > 0 else np.nan,
        state_change_span_mb=state_change_span / 1e6,
        state_change_span_frac=state_change_span / chr_span if chr_span and chr_span > 0 else np.nan,
    )


# ----------------------------------------------------------------------
# cohort-level normalization: localization + confounder residuals
# ----------------------------------------------------------------------

def add_within_sample_contrast(cohort_df, cols, sample_col="SAMPLE_ID"):
    """
    For each feature, add a 'ratio-to-genome' (value / within-sample median) and
    a within-sample rank. Chromothripsis is localized, so a chromosome that is an
    outlier *relative to the rest of its own genome* is more suspicious than one
    that is merely large in absolute terms.
    """
    d = cohort_df.copy()
    for col in cols:
        med = d.groupby(sample_col)[col].transform("median")
        d[f"{col}__ratio_to_sample"] = d[col] / med.replace(0, np.nan)
        d[f"{col}__rank_in_sample"] = d.groupby(sample_col)[col].rank(
            ascending=False, method="min")
    return d


def add_probe_count_residual(cohort_df, frag_col="n_segment_state_changes",
                             probe_col="n_probes"):
    """
    Regress a fragmentation feature on log(probe count) across the whole cohort
    and keep the residual. This strips the 'gene-dense chromosomes look fragmented'
    confounder so the residual reflects fragmentation beyond what probe density alone
    predicts.
    """
    d = cohort_df.copy()
    x = np.log1p(pd.to_numeric(d[probe_col], errors="coerce"))
    y = pd.to_numeric(d[frag_col], errors="coerce")
    ok = x.notna() & y.notna()
    if ok.sum() > 10:
        b1, b0 = np.polyfit(x[ok], y[ok], 1)
        d[f"{frag_col}__resid_vs_probes"] = y - (b0 + b1 * x)
    else:
        d[f"{frag_col}__resid_vs_probes"] = np.nan
    return d


# ----------------------------------------------------------------------
# evaluation
# ----------------------------------------------------------------------

def feature_screen(cohort_df, feature_cols, truth_col="truth_positive_chrom"):
    """
    Per-feature chromosome-level discrimination against the gold label.
    Reports average precision (preferred under heavy imbalance) and AUROC.
    """
    from sklearn.metrics import average_precision_score, roc_auc_score
    y = cohort_df[truth_col].astype(int).values
    out = []
    for col in feature_cols:
        s = pd.to_numeric(cohort_df[col], errors="coerce").values
        ok = np.isfinite(s)
        if y[ok].sum() == 0 or (1 - y[ok]).sum() == 0 or ok.sum() < 20:
            continue
        # try the feature and its sign-flip; report the better-oriented one
        ap_pos = average_precision_score(y[ok], s[ok])
        ap_neg = average_precision_score(y[ok], -s[ok])
        direction = "high" if ap_pos >= ap_neg else "low"
        s_dir = s[ok] if direction == "high" else -s[ok]
        out.append(dict(
            feature=col,
            direction=direction,
            average_precision=max(ap_pos, ap_neg),
            auroc=roc_auc_score(y[ok], s_dir),
            n=int(ok.sum()),
            n_pos=int(y[ok].sum()),
        ))
    return pd.DataFrame(out).sort_values("average_precision", ascending=False).reset_index(drop=True)


def learned_combiner_oof(cohort_df, feature_cols,
                         truth_col="truth_positive_chrom",
                         sample_col="SAMPLE_ID",
                         n_splits=5, C=1.0, seed=0):
    """
    Out-of-fold chromosome-level probabilities from an L1-logistic combiner,
    with folds grouped by SAMPLE_ID (chromosomes within a sample are correlated,
    so a sample's chromosomes must not straddle train/test even if samples are
    treated as independent of each other).

    Returns the cohort_df with an added 'p_ctx_oof' column plus the fitted
    coefficient table from a final full-data fit (for interpretation).
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline
    from sklearn.model_selection import GroupKFold
    from sklearn.impute import SimpleImputer

    d = cohort_df.copy().reset_index(drop=True)
    # drop features with no information cohort-wide (all-NaN or constant)
    feat = []
    for c in feature_cols:
        v = pd.to_numeric(d[c], errors="coerce")
        if v.notna().sum() > 0 and v.nunique(dropna=True) > 1:
            feat.append(c)
    feature_cols = feat

    X = d[feature_cols].apply(pd.to_numeric, errors="coerce").values
    y = d[truth_col].astype(int).values
    groups = d[sample_col].values

    pipe = lambda: make_pipeline(
        # keep_empty_features keeps the column count stable across folds so the
        # coefficient vector always aligns with feature_cols
        SimpleImputer(strategy="median", keep_empty_features=True),
        StandardScaler(),
        LogisticRegression(penalty="l1", solver="liblinear", C=C,
                           class_weight="balanced", max_iter=2000),
    )

    oof = np.full(len(d), np.nan)
    gkf = GroupKFold(n_splits=n_splits)
    for tr, te in gkf.split(X, y, groups):
        model = pipe()
        model.fit(X[tr], y[tr])
        oof[te] = model.predict_proba(X[te])[:, 1]
    d["p_ctx_oof"] = oof

    # final full-data fit for coefficient inspection
    final = pipe()
    final.fit(X, y)
    coefs = final.named_steps["logisticregression"].coef_.ravel()
    coef_tbl = (pd.DataFrame({"feature": feature_cols, "coef": coefs})
                .assign(abs_coef=lambda t: t["coef"].abs())
                .sort_values("abs_coef", ascending=False)
                .reset_index(drop=True))
    return d, coef_tbl


def sample_level_eval(pred_df, score_col="p_ctx_oof", threshold=0.5,
                      truth_col="truth_positive_chrom", sample_col="SAMPLE_ID"):
    """
    Turn chromosome-level scores into the sample-level outcome you actually care
    about: a sample is called positive if any chromosome scores above threshold;
    nominated chromosomes are those above threshold. Then apply your .tex criteria:
      recovered positive  : truth nonempty AND nominated intersects truth
      correct negative     : truth empty AND nothing nominated
    Use OOF scores so the numbers are honest, then sweep threshold for a curve.
    """
    d = pred_df.copy()
    d["_nom"] = d[score_col] >= threshold
    d["_truth"] = d[truth_col].astype(bool)

    rows = []
    for sid, g in d.groupby(sample_col):
        truth_chrs = set(g.loc[g["_truth"], "chrom"])
        nom_chrs = set(g.loc[g["_nom"], "chrom"])
        is_pos = len(truth_chrs) > 0
        if is_pos:
            status = "recovered" if (nom_chrs & truth_chrs) else "missed"
        else:
            status = "correct_negative" if len(nom_chrs) == 0 else "false_positive"
        rows.append(dict(SAMPLE_ID=sid, is_positive=is_pos,
                         n_truth=len(truth_chrs), n_nominated=len(nom_chrs),
                         status=status))
    ev = pd.DataFrame(rows)
    pos = ev[ev["is_positive"]]
    neg = ev[~ev["is_positive"]]
    summary = dict(
        threshold=threshold,
        n_pos=len(pos), n_neg=len(neg),
        pos_recovered=int((pos["status"] == "recovered").sum()),
        pos_missed=int((pos["status"] == "missed").sum()),
        neg_correct=int((neg["status"] == "correct_negative").sum()),
        neg_false_positive=int((neg["status"] == "false_positive").sum()),
    )
    summary["sample_errors"] = summary["pos_missed"] + summary["neg_false_positive"]
    return ev, summary


def threshold_sweep(pred_df, score_col="p_ctx_oof",
                    thresholds=None, **kwargs):
    if thresholds is None:
        thresholds = np.linspace(0.05, 0.95, 19)
    out = [sample_level_eval(pred_df, score_col=score_col, threshold=t, **kwargs)[1]
           for t in thresholds]
    return pd.DataFrame(out)


# ======================================================================
# Literature-grounded chromothripsis score
# ----------------------------------------------------------------------
# Structure (per chromosome):
#
#     score = Q * (w_osc * F_osc + w_allelic * F_allelic)
#
#   F_osc      mechanistic oscillation evidence in [0,1]. Built on the
#              ShatterSeek quantity "longest run of contiguous segments
#              confined to <=K copy-number states". A monotone/BFB staircase
#              visits many distinct states so its <=K window is short ->
#              F_osc ~ 0, which is what gives specificity. (Cortes-Ciriano
#              2020; Rausch 2012; Korbel & Campbell 2013.)
#   F_allelic  single-haplotype LOH footprint in [0,1]: interspersed
#              LOH<->het transitions, the empirically strongest panel signal.
#   Q          quality / confidence weight in (0,1]. DOWNWEIGHTS (does not
#              abstain) for low purity, sparse probes, high probe noise, so
#              weak evidence on a noisy low-purity arm cannot masquerade as a
#              call. Each component has a floor so Q never hits exactly 0.
#
# All thresholds are TUNABLE params with fixed defaults set LOW relative to
# ShatterSeek because panel probe counts are far below WGS resolution and the
# cohort mixes panel versions of differing density/noise.
# ======================================================================

from dataclasses import dataclass, field


def _ramp(x, lo, hi):
    """Clamped linear ramp: 0 below lo, 1 above hi, linear between."""
    x = np.asarray(x, dtype=float)
    if hi <= lo:
        return (x >= hi).astype(float)
    return np.clip((x - lo) / (hi - lo), 0.0, 1.0)


def longest_le_k_distinct_run(states, k):
    """
    Collapse consecutive identical states, then return the length (in collapsed
    units = number of oscillating segments) of the longest contiguous window
    containing at most k DISTINCT copy-number states.

    Oscillation A-B-A-B-A  (k=2) -> collapsed length 5  -> long.
    Staircase   A-B-C-D-E  (k=2) -> any window <=2       -> short (2).
    This is the discrimination that separates chromothripsis from BFB.
    """
    s = [v for v in states if v is not None and not (isinstance(v, float) and np.isnan(v))]
    if len(s) == 0:
        return 0
    # collapse runs of identical adjacent states
    collapsed = [s[0]]
    for v in s[1:]:
        if v != collapsed[-1]:
            collapsed.append(v)
    # sliding window, longest subarray with <= k distinct values
    best = 0
    left = 0
    counts = {}
    for right in range(len(collapsed)):
        counts[collapsed[right]] = counts.get(collapsed[right], 0) + 1
        while len(counts) > k:
            counts[collapsed[left]] -= 1
            if counts[collapsed[left]] == 0:
                del counts[collapsed[left]]
            left += 1
        best = max(best, right - left + 1)
    return int(best)


@dataclass
class ScoreParams:
    # oscillation factor: ramp on longest <=K-state chain length (segments)
    k_states: int = 3          # CN states allowed within an oscillation chain
    osc_chain_lo: float = 2.0  # below this -> no oscillation evidence
    osc_chain_hi: float = 6.0  # at/above this -> full (sparse-panel scaled)
    # allelic factor: ramp on LOH<->het interleave count
    allelic_lo: float = 1.0
    allelic_hi: float = 6.0
    allelic_floor: float = 0.25   # when allelic data is weak/absent, don't zero it
    loh_extent_lo: float = 0.30   # frac_loh below this adds no allelic credit
    loh_extent_hi: float = 0.70   # extensive single-haplotype LOH -> full credit
    min_het_for_allelic: int = 4  # need this many het probes to trust allelic
    # mechanistic mixture weights (need not sum to 1; score is rescaled by sweep)
    w_osc: float = 0.6
    w_allelic: float = 0.4
    # quality weight components (each floored so Q downweights, never abstains)
    purity_lo: float = 0.05
    purity_hi: float = 0.35
    purity_floor: float = 0.15
    probe_lo: float = 4.0
    probe_hi: float = 15.0
    probe_floor: float = 0.15
    noise_lo: float = 0.05     # probe-noise (median lorvar) below this is clean
    noise_hi: float = 0.30     # above this is very noisy
    noise_floor: float = 0.15
    osc_primary_floor: float = 0.5  # for mech_combine='osc_primary'
    use_probe_chain: bool = False   # fold probe-level re-seg chain into F_osc (crude; off)
    use_loh_extent: bool = False    # fold frac_loh into F_allelic (NOT CT-specific; off)
    # --- v2 specificity gates (default OFF -> reproduces v1 exactly) ---
    min_chain_cn_amplitude: int = 0   # require oscillating chain to span >= this many CN units (0 = off)
    max_resolved_levels: int = 0      # down-weight chromosomes spanning > this many CN states (0 = off)
    resolved_levels_penalty: float = 0.0  # multiplicative shrink per extra level above max (0 = off)


def compute_ctx_score(cohort_df, purity_map, params=None, use_quality=True,
                      mech_combine="sum", tcn_int_col="_osc_states_json"):
    """
    Add F_osc, F_allelic, mech, Q, and ctx_score columns to a cohort feature
    table. Requires columns produced by chrom_features_one_sample:
        longest_osc_chain_k2, longest_osc_chain_k3, loh_het_interleave,
        n_het_probes, n_probes, probe_noise
    plus a purity_map: {SAMPLE_ID -> purity}.
    """
    p = params or ScoreParams()
    d = cohort_df.copy()
    # ---- F_osc: oscillation evidence (conjunction lives in the chain length)
    # Use the MAX of the FACETS-segment chain and the FACETS-independent
    # probe-level re-segmentation chain, so oscillation that FACETS merged away
    # (e.g. C11 chr9) is still credited.
    #chain_col = "longest_osc_chain_k3" if p.k_states >= 3 else "longest_osc_chain_k2"
    chain_col = f"longest_osc_chain_k{int(np.clip(p.k_states, 2, 5))}"
    seg_chain = pd.to_numeric(d[chain_col], errors="coerce").fillna(0)
    if p.use_probe_chain and "probe_osc_chain" in d.columns:
        # NOTE: the rolling-median re-segmentation is crude and can manufacture
        # spurious oscillation on noisy negatives -- off by default. A proper
        # changepoint posterior is the right way to use probe-level signal.
        probe_chain = pd.to_numeric(d["probe_osc_chain"], errors="coerce").fillna(0)
        chain = np.maximum(seg_chain, probe_chain)
    else:
        chain = seg_chain
    F_osc = _ramp(chain, p.osc_chain_lo, p.osc_chain_hi)
    # v2 gate A: CN-amplitude. Oscillation must move total CN by a real margin,
    # not just a 2<->1 allelic-flip wobble. Zero F_osc below the amplitude floor.
    if p.min_chain_cn_amplitude > 0 and "chain_cn_amplitude" in d.columns:
        amp = pd.to_numeric(d["chain_cn_amplitude"], errors="coerce").fillna(0)
        F_osc = np.where(amp >= p.min_chain_cn_amplitude, F_osc, 0.0)
    # v2 gate B: state-count penalty. Down-weight chromosomes spanning many CN
    # levels (BFB/high-level amplification), enforcing the "limited 2-3 states"
    # criterion. Multiplicative shrink per distinct level above the cap.
    if p.max_resolved_levels > 0 and p.resolved_levels_penalty > 0 and "n_resolved_levels" in d.columns:
        lev = pd.to_numeric(d["n_resolved_levels"], errors="coerce").fillna(0)
        excess = np.clip(lev - p.max_resolved_levels, 0, None)
        F_osc = F_osc * (1.0 - p.resolved_levels_penalty) ** excess
    # ---- F_allelic: LOH<->het interleave, floored & gated by het support
    interleave_raw = _ramp(pd.to_numeric(d["loh_het_interleave"], errors="coerce").fillna(0),
                           p.allelic_lo, p.allelic_hi)
    if p.use_loh_extent:
        # WARNING: extensive LOH is ubiquitous in cancer and NOT CT-specific;
        # enabling this floods F_allelic across negatives. Off by default.
        loh_extent = _ramp(pd.to_numeric(d["frac_loh"], errors="coerce").fillna(0),
                           p.loh_extent_lo, p.loh_extent_hi)
        allelic_raw = np.maximum(interleave_raw, loh_extent)
    else:
        allelic_raw = interleave_raw
    het_ok = pd.to_numeric(d["n_het_probes"], errors="coerce").fillna(0) >= p.min_het_for_allelic
    F_allelic = np.where(het_ok,
                         p.allelic_floor + (1 - p.allelic_floor) * allelic_raw,
                         p.allelic_floor)
    if mech_combine == "sum":
        mech = p.w_osc * F_osc + p.w_allelic * F_allelic
    elif mech_combine == "product":
        mech = F_osc * F_allelic            # AND-gate: need BOTH
    elif mech_combine == "min":
        mech = np.minimum(F_osc, F_allelic)
    elif mech_combine == "geom":
        mech = (F_osc ** p.w_osc) * (F_allelic ** p.w_allelic)
    elif mech_combine == "osc_primary":
        # oscillation is the primary CT signature; allelic modulates with a high
        # floor so a maxed F_osc with moderate allelic still scores (e.g. C11 chr9)
        mech = F_osc * (p.osc_primary_floor + (1 - p.osc_primary_floor) * F_allelic)
    elif mech_combine == "noisy_and":
        # soft AND: high if both high, but one strong leg partially carries
        mech = 1.0 - (1.0 - F_osc) * (1.0 - F_allelic) * (1.0 - np.minimum(F_osc, F_allelic))
    else:
        raise ValueError(f"unknown mech_combine: {mech_combine}")
    # ---- Q: quality / confidence weight (downweight, never abstain) ----
    purity = d["SAMPLE_ID"].map(purity_map).astype(float)
    q_pur = p.purity_floor + (1 - p.purity_floor) * _ramp(purity, p.purity_lo, p.purity_hi)
    q_prb = p.probe_floor + (1 - p.probe_floor) * _ramp(
        pd.to_numeric(d["n_probes"], errors="coerce").fillna(0), p.probe_lo, p.probe_hi)
    noise = pd.to_numeric(d.get("probe_noise"), errors="coerce").fillna(p.noise_hi)
    q_noise = p.noise_floor + (1 - p.noise_floor) * (1.0 - _ramp(noise, p.noise_lo, p.noise_hi))
    Q = q_pur * q_prb * q_noise
    d["F_osc"] = F_osc
    d["F_allelic"] = F_allelic
    d["mech_evidence"] = mech
    d["Q_weight"] = Q
    d["ctx_score"] = (Q * mech) if use_quality else mech
    return d


def apply_classifier_v3(cohort_df, purity_map, *, params=None, threshold=0.558,
                        suppress=(0.0, 3.0, 0.0),
                        noise_col="dot_mad_within", noise_max=None,
                        score_col="ctx_score", probe_conf_col="cp_osc_confidence",
                        truth_col="truth_positive_chrom", sample_col="SAMPLE_ID"):
    """
    v3 chromothripsis classifier = v2 (Q-off score + probe-confidence suppressor)
    plus a per-chromosome probe-noise veto.

    Per (sample, chromosome):
        ctx_score    = compute_ctx_score(..., use_quality=False)        # = mech_evidence
        ctx_suppress = ctx_score * clip((C-c0)/(c1-c0), 0, 1) floored at m0,
                       C = cp_osc_confidence, (c0,c1,m0) = `suppress`
        nominated_v2 = ctx_suppress >= threshold
        noise_veto   = noise_col > noise_max
        nominated_v3 = nominated_v2 AND NOT noise_veto
    A sample is flagged if ANY of its chromosomes is nominated.

    `cohort_df` must carry `probe_conf_col` and `noise_col`. If noise_max is None
    it is auto-calibrated recall-safe (max noise_col among chromosomes v2 catches),
    so the veto never removes a recovering chromosome.
    """
    if noise_col not in cohort_df.columns:
        raise KeyError(f"noise_col '{noise_col}' missing; build it (chrom_features_one_sample now emits it).")

    c0, c1, m0 = suppress
    d = compute_ctx_score(cohort_df, purity_map, params=params, use_quality=False)

    seg  = pd.to_numeric(d[score_col], errors="coerce").fillna(0.0)
    conf = pd.to_numeric(d[probe_conf_col], errors="coerce").fillna(0.0)
    ramp = np.clip((conf - c0) / max(c1 - c0, 1e-9), 0.0, 1.0)
    d["ctx_suppress"] = seg * (m0 + (1.0 - m0) * ramp)
    d["nominated_v2"] = d["ctx_suppress"] >= threshold

    noise = pd.to_numeric(d[noise_col], errors="coerce")
    if noise_max is None:
        if truth_col not in d.columns:
            raise ValueError("noise_max=None needs truth_col for recall-safe calibration.")
        t = d[truth_col]
        t = t if t.dtype == bool else t.astype(str).str.strip().str.lower().isin(["true", "1", "yes"])
        caught = noise[t & d["nominated_v2"]]
        noise_max = float(caught.max()) if caught.notna().any() else float("inf")

    d["noise_max"]      = noise_max
    d["noise_veto"]     = noise.notna() & (noise > noise_max)
    d["nominated_v3"]   = d["nominated_v2"] & (~d["noise_veto"])
    d["sample_flag_v3"] = d.groupby(sample_col)["nominated_v3"].transform("any")
    return d


def apply_classifier_v5(cohort_df, purity_map, *, params=None,
                        threshold=0.35,
                        suppress=(0.0, 3.0, 0.0),
                        noise_col="dot_mad_within",
                        noise_max=1.135,
                        probe_conf_col="cp_osc_confidence",
                        footprint_span_col="altered_span_frac",
                        footprint_union_col="altered_union_frac",
                        footprint_island_col="n_altered_islands",
                        footprint_span_min=0.25,
                        footprint_union_min=0.15,
                        footprint_island_max=3,
                        truth_col="truth_positive_chrom",
                        sample_col="SAMPLE_ID"):
    """
    v5 chromothripsis classifier.

    v5 = product-combiner mechanistic evidence
         + probe oscillation confidence suppressor
         + fixed probe-noise veto
         + altered-footprint veto.

    Per (sample, chromosome):
        F_osc, F_allelic = compute_ctx_score(...)
        M                = F_osc * F_allelic
        S                = M * clip(C/3, 0, 1), C = cp_osc_confidence
        weak_footprint   = altered_span_frac < 0.25
                           OR (altered_union_frac < 0.15 AND n_altered_islands <= 3)
        nominated_v5     = S >= 0.35 AND dot_mad_within <= 1.135 AND NOT weak_footprint

    A sample is flagged if any chromosome is nominated.
    """
    required = [noise_col, probe_conf_col, footprint_span_col, footprint_union_col, footprint_island_col]
    missing = [c for c in required if c not in cohort_df.columns]
    if missing:
        raise KeyError(f"Missing required v5 columns: {missing}")

    c0, c1, m0 = suppress

    d = compute_ctx_score(
        cohort_df,
        purity_map,
        params=params,
        use_quality=False,
        mech_combine="product",
    )

    mech = pd.to_numeric(d["ctx_score"], errors="coerce").fillna(0.0)
    conf = pd.to_numeric(d[probe_conf_col], errors="coerce").fillna(0.0)

    ramp = np.clip((conf - c0) / max(c1 - c0, 1e-9), 0.0, 1.0)
    d["ctx_suppress"] = mech * (m0 + (1.0 - m0) * ramp)

    noise = pd.to_numeric(d[noise_col], errors="coerce")
    span = pd.to_numeric(d[footprint_span_col], errors="coerce").fillna(-np.inf)
    union = pd.to_numeric(d[footprint_union_col], errors="coerce").fillna(-np.inf)
    islands = pd.to_numeric(d[footprint_island_col], errors="coerce").fillna(0)

    d["nominated_score_v5"] = d["ctx_suppress"] >= threshold
    d["noise_veto_v5"] = noise.notna() & (noise > noise_max)
    d["weak_altered_footprint_v5"] = (
        (span < footprint_span_min) |
        ((union < footprint_union_min) & (islands <= footprint_island_max))
    )

    d["nominated_v5"] = (
        d["nominated_score_v5"] &
        (~d["noise_veto_v5"]) &
        (~d["weak_altered_footprint_v5"])
    )

    d["sample_flag_v5"] = d.groupby(sample_col)["nominated_v5"].transform("any")
    d["v5_threshold"] = threshold
    d["v5_noise_max"] = noise_max
    d["v5_footprint_span_min"] = footprint_span_min
    d["v5_footprint_union_min"] = footprint_union_min
    d["v5_footprint_island_max"] = footprint_island_max

    return d


def apply_classifier_v6(cohort_df, purity_map, *, params=None,
                        threshold=0.35,
                        suppress=(0.0, 3.0, 0.0),
                        noise_col="dot_mad_within",
                        noise_max=1.135,
                        probe_conf_col="cp_osc_confidence",
                        chain_col="longest_osc_chain_k3",
                        min_osc_chain=5,
                        footprint_union_col="altered_union_frac",
                        footprint_island_col="n_altered_islands",
                        footprint_union_min=0.15,
                        footprint_island_max=3,
                        sample_col="SAMPLE_ID"):
    required = [
        noise_col, probe_conf_col, chain_col,
        footprint_union_col, footprint_island_col,
    ]
    missing = [c for c in required if c not in cohort_df.columns]
    if missing:
        raise KeyError(f"Missing required v6 columns: {missing}")

    c0, c1, m0 = suppress

    d = compute_ctx_score(
        cohort_df,
        purity_map,
        params=params,
        use_quality=False,
        mech_combine="product",
    )

    mech = pd.to_numeric(d["ctx_score"], errors="coerce").fillna(0.0)
    conf = pd.to_numeric(d[probe_conf_col], errors="coerce").fillna(0.0)
    ramp = np.clip((conf - c0) / max(c1 - c0, 1e-9), 0.0, 1.0)

    d["ctx_suppress"] = mech * (m0 + (1.0 - m0) * ramp)

    noise = pd.to_numeric(d[noise_col], errors="coerce")
    chain = pd.to_numeric(d[chain_col], errors="coerce").fillna(0)
    union = pd.to_numeric(d[footprint_union_col], errors="coerce").fillna(-np.inf)
    islands = pd.to_numeric(d[footprint_island_col], errors="coerce").fillna(0)

    d["nominated_score_v6"] = d["ctx_suppress"] >= threshold
    d["min_osc_chain_pass_v6"] = chain >= min_osc_chain
    d["noise_veto_v6"] = noise.notna() & (noise > noise_max)

    d["weak_altered_footprint_v6"] = (
        (union < footprint_union_min) &
        (islands <= footprint_island_max)
    )

    d["nominated_v6"] = (
        d["nominated_score_v6"] &
        d["min_osc_chain_pass_v6"] &
        (~d["noise_veto_v6"]) &
        (~d["weak_altered_footprint_v6"])
    )

    d["sample_flag_v6"] = d.groupby(sample_col)["nominated_v6"].transform("any")

    d["v6_threshold"] = threshold
    d["v6_noise_max"] = noise_max
    d["v6_min_osc_chain"] = min_osc_chain
    d["v6_footprint_union_min"] = footprint_union_min
    d["v6_footprint_island_max"] = footprint_island_max

    return d


# ======================================================================
# v7 classifier
# ----------------------------------------------------------------------
# Thresholds calibrated against expert-reviewed chromothripsis calls.
# Changes relative to v6:
#
#   1. min_osc_chain 5 -> 4.
#   2. noise_max 1.135 -> 1.0 (tighter probe-noise ceiling).
#   3. An alternate allelic-evidence path (see F_allelic_v7 below) and a
#      parallel loss-dominant-baseline pathway (see
#      compute_baseline_shard_features) for chromothripsis patterns the v6
#      mechanistic score structurally cannot see: post-WGD chromothripsis
#      that never drops minor-allele copy number to literal zero, and
#      loss-dominant chromothripsis where the shattered/lost state IS the
#      chromosome's own local majority, so nothing reads as "altered"
#      relative to itself.
#   4. A gated "second look" on cp_osc_confidence (see
#      compute_second_look_confidence below). The probe-level confidence
#      check runs at an anti-flood pen_const=8.0, calibrated cohort-wide to
#      stop noisy libraries from generating spurious changepoints. That
#      same stiffness can zero out a chromosome with a maximal FACETS-level
#      mechanistic score (ctx_score>=0.8) whose real, well-separated jumps
#      still don't individually clear the fixed per-breakpoint cost. For
#      that narrow case only (ctx_score >= 0.8 AND cp_osc_confidence == 0),
#      re-check at pen_const=5.0 and use it if higher.
# ======================================================================

def compute_baseline_shard_features(seg_df, ploidy, *,
                                    shard_baseline_mult=5.0,
                                    deficit_min=0.5,
                                    baseline_frac_min=0.5,
                                    small_frac_max=0.15,
                                    multi_shard_min=2,
                                    tcn_col="tcn.em",
                                    centromere_mb=CENTROMERE_HG19_MB):
    """
    Per-chromosome "loss-dominant baseline + retained/amplified shard" features.

    Some chromothripsis events are loss-dominant: most of the chromosome
    collapses to one low copy-number state, so oscillation-chain features
    (built relative to that chromosome's own modal state) see nothing
    "altered" -- the lost state IS the local mode. Comparing the chromosome's
    dominant state against the sample's genome-wide ploidy surfaces these
    cases. Requiring a small, highly-amplified surviving fragment (a "shard")
    keeps the signal specific to real shattering rather than an ordinary
    whole-chromosome deletion, which has no such internal structure.

    The shard threshold is relative to the chromosome's OWN loss-level
    baseline_cn (floored at 1 to avoid a vacuous threshold when baseline_cn
    is itself 0, e.g. a homozygous-deletion backdrop), not the sample's
    genome-wide ploidy -- a low-ploidy sample can otherwise set an
    unreachably high bar for a chromosome whose local backdrop is already
    a deep loss. A shard also must be genomically EMBEDDED in the loss
    backdrop (a baseline-state segment somewhere upstream AND somewhere
    downstream of it), not merely small; a small high segment flanked by an
    unrelated gain rather than the loss backdrop is not a shard.

    A candidate is also rejected if it overlaps CENTROMERE_HG19_MB for its
    chromosome: FACETS often emits a single sparse-probe segment spanning the
    low-mappability pericentromeric region, which recurs at the same
    coordinates across unrelated samples and is a panel artifact rather than
    an amplified fragment. Only the centromeric candidate is removed, so
    genuine shards elsewhere on the same chromosome are unaffected.

    Multi-shard oscillation: when a chromosome
    carries >= multi_shard_min embedded shards, the baseline -> shard -> baseline
    -> shard -> baseline alternation is itself an oscillation between two states,
    so the loss-dominance requirements (deficit_min, baseline_frac_min) are waived.
    Every shard requirement still applies (>= shard_baseline_mult x baseline,
    <= small_frac_max of the chromosome, embedded, non-centromeric), as does the
    v7 noise veto downstream. A single shard keeps the original loss-dominant
    criteria. This captures shattered chromosomes on a copy-neutral or gained
    backdrop that carry several high-level fragments. Pass multi_shard_min=None
    to disable.

    Returns one row per chromosome: baseline_cn (the length-weighted modal
    total-CN state), baseline_frac (how much of the chromosome that state
    covers), deficit (ploidy - baseline_cn), n_shards (embedded shards only;
    see n_shard_candidates for the pre-embedding count), the booleans
    loss_dominant_hit (original single-shard rule) and multi_shard_hit (rule
    above), and loss_dominant_shard_hit = either one (the shard-pathway flag
    used downstream) -- plus extent_start/extent_end/extent_length_mb/
    extent_frac: the genomic SPAN (first-to-last position) of the segments
    at the baseline state, as opposed to baseline_frac's total COVERAGE.
    These differ whenever the baseline state is fragmented (interrupted by
    the shard, or by other segments) rather than one contiguous block -- a
    baseline that's scattered in small patches across the whole chromosome
    has high extent but can have much lower coverage; a baseline confined
    to one corner has both low extent and low coverage. extent_frac divides
    by the chromosome's own observed SPAN (first segment start to last
    segment end, across all states) rather than by summed segment length,
    since the two diverge whenever a chromosome carries a large
    probe-desert gap (e.g. spanning the centromere).
    """
    seg = seg_df.copy()
    seg["_c"] = seg["chrom"].map(_norm_chrom)
    seg["_start"] = pd.to_numeric(seg["loc.start"], errors="coerce")
    seg["_end"] = pd.to_numeric(seg["loc.end"], errors="coerce")
    seg["_len"] = (seg["_end"] - seg["_start"]).clip(lower=0)

    rows = []
    for c, g in seg.groupby("_c"):
        g = g.assign(_tcn=pd.to_numeric(g[tcn_col], errors="coerce")).dropna(subset=["_tcn"])
        chrom_len = g["_len"].sum()
        chrom_span = (g["_end"].max() - g["_start"].min()) if len(g) else np.nan
        if g.empty or chrom_len <= 0 or pd.isna(ploidy):
            rows.append(dict(chrom=c, baseline_cn=np.nan, baseline_frac=np.nan,
                              deficit=np.nan, n_shards=0, n_shard_candidates=0,
                              loss_dominant_hit=False, multi_shard_hit=False,
                              loss_dominant_shard_hit=False,
                              extent_start=np.nan, extent_end=np.nan,
                              extent_length_mb=np.nan, extent_frac=np.nan))
            continue
        by_state = g.groupby("_tcn")["_len"].sum().sort_values(ascending=False)
        baseline_cn = float(by_state.index[0])
        baseline_frac = float(by_state.iloc[0] / chrom_len)
        deficit = float(ploidy - baseline_cn)
        g = g.assign(_len_frac=g["_len"] / chrom_len).sort_values("_start").reset_index(drop=True)
        is_baseline = np.isclose(g["_tcn"].values, baseline_cn)
        shard_thr = shard_baseline_mult * max(baseline_cn, 1.0)
        candidate_mask = (g["_len_frac"] <= small_frac_max) & (g["_tcn"] >= shard_thr)
        candidate_idx = g.index[candidate_mask].tolist()
        n_shard_candidates = len(candidate_idx)

        cen = centromere_mb.get(c)
        if cen is not None:
            cen_start_bp, cen_end_bp = cen[0] * 1e6, cen[1] * 1e6
            non_centromeric_idx = [i for i in candidate_idx
                                    if not (g.at[i, "_start"] < cen_end_bp and g.at[i, "_end"] > cen_start_bp)]
        else:
            non_centromeric_idx = candidate_idx

        embedded_idx = [i for i in non_centromeric_idx
                         if is_baseline[:i].any() and is_baseline[i + 1:].any()]
        shards = g.loc[embedded_idx]
        loss_dominant_hit = (deficit >= deficit_min) and (baseline_frac >= baseline_frac_min) and (len(shards) >= 1)
        multi_shard_hit = multi_shard_min is not None and len(shards) >= multi_shard_min
        hit = loss_dominant_hit or multi_shard_hit

        base_segs = g.loc[g["_tcn"] == baseline_cn]
        extent_start = float(base_segs["_start"].min())
        extent_end = float(base_segs["_end"].max())
        extent_length_mb = (extent_end - extent_start) / 1e6
        extent_frac = float((extent_end - extent_start) / chrom_span) if chrom_span and chrom_span > 0 else np.nan

        rows.append(dict(chrom=c, baseline_cn=baseline_cn, baseline_frac=round(baseline_frac, 4),
                          deficit=round(deficit, 3), n_shards=int(len(shards)),
                          n_shard_candidates=n_shard_candidates,
                          loss_dominant_hit=bool(loss_dominant_hit),
                          multi_shard_hit=bool(multi_shard_hit),
                          loss_dominant_shard_hit=bool(hit),
                          extent_start=extent_start, extent_end=extent_end,
                          extent_length_mb=round(extent_length_mb, 2),
                          extent_frac=round(extent_frac, 4) if pd.notna(extent_frac) else np.nan))
    return pd.DataFrame(rows)


def compute_second_look_confidence(cohort_df, man_df, tablesDir, *,
                                   score_col="ctx_score",
                                   conf_col="cp_osc_confidence",
                                   score_min=0.8,
                                   pen_fallback=5.0,
                                   sample_col="SAMPLE_ID",
                                   chrom_col="chrom"):
    """
    Gated re-check of cp_osc_confidence at a gentler penalty for chromosomes
    the anti-flood pen_const=8.0 probe recheck zeroed out despite a strong
    FACETS-level mechanistic score. See the v7 comment block, change 4.

    `cohort_df` must already carry `score_col` (i.e. call after
    compute_ctx_score / apply_classifier_v6). Only rows matching the gate are
    recomputed (typically very few), re-reading cached jointseg TSVs from
    `tablesDir`.

    Returns a DataFrame with columns [sample_col, chrom_col,
    "cp_osc_confidence_v7"] -- the gentler-penalty confidence, ready to
    `.merge(..., how="left")` and `np.maximum` against the original
    `conf_col` in apply_classifier_v7.
    """
    import os
    import ctx_changepoint as ccp

    gate = (pd.to_numeric(cohort_df[score_col], errors="coerce") >= score_min) & \
           (pd.to_numeric(cohort_df[conf_col], errors="coerce") == 0.0)
    cand = cohort_df.loc[gate, [sample_col, chrom_col]].drop_duplicates()

    rows = []
    for _, row in cand.iterrows():
        sid, chrom = row[sample_col], row[chrom_col]
        jpath = os.path.join(tablesDir, f"{sid}.hisens.jointseg.tsv")
        if not os.path.exists(jpath):
            continue
        j = pd.read_csv(jpath, sep="\t")
        j["_c"] = j["chrom"].map(_norm_chrom)
        sub = j.loc[j["_c"] == _norm_chrom(chrom)].sort_values("maploc")
        prow = man_df.loc[man_df[sample_col].astype(str).eq(str(sid))]
        purity = float(prow["purity"].iloc[0]) if prow.shape[0] and pd.notna(prow["purity"].iloc[0]) else None
        ploidy = float(prow["ploidy"].iloc[0]) if prow.shape[0] and pd.notna(prow["ploidy"].iloc[0]) else 2.0
        f = ccp.probe_changepoint_features(sub, purity=purity, ploidy=ploidy, pen_const=pen_fallback)
        rows.append({sample_col: sid, chrom_col: chrom, "cp_osc_confidence_v7": f["cp_osc_confidence"]})
    return pd.DataFrame(rows, columns=[sample_col, chrom_col, "cp_osc_confidence_v7"])


def apply_classifier_v7(cohort_df, purity_map, shard_df, *, second_look_df=None, params=None,
                        threshold=0.35,
                        suppress=(0.0, 3.0, 0.0),
                        noise_col="dot_mad_within",
                        noise_max=1.0,
                        probe_conf_col="cp_osc_confidence",
                        chain_col="longest_osc_chain_k3",
                        min_osc_chain=4,
                        footprint_union_col="altered_union_frac",
                        footprint_island_col="n_altered_islands",
                        footprint_union_min=0.15,
                        footprint_island_max=3,
                        allelic_alt_lo=8.0,
                        allelic_alt_hi=15.0,
                        allelic_alt_frac_loh_max=0.05,
                        sample_col="SAMPLE_ID",
                        chrom_col="chrom"):
    """
    v6 with four calibrated changes -- see module docstring above. `shard_df`
    is the concatenation of `compute_baseline_shard_features` output across
    the cohort, with a SAMPLE_ID column added per sample, merged on
    (SAMPLE_ID, chrom). `second_look_df`, if given, is
    `compute_second_look_confidence`'s output (SAMPLE_ID, chrom,
    cp_osc_confidence_v7), merged the same way; omit to fall back to plain
    v6-style confidence (no second look).

    `noise_max=None` nominates without the noise veto; the pipeline then applies
    `apply_noise_veto` as its last chromosome-level step. Because every step acts
    on each chromosome independently, this gives the same calls as vetoing here.
    """
    required = [
        noise_col, probe_conf_col, chain_col,
        footprint_union_col, footprint_island_col,
        "loh_het_interleave", "frac_loh", "n_allelic_state_changes",
    ]
    missing = [c for c in required if c not in cohort_df.columns]
    if missing:
        raise KeyError(f"Missing required v7 columns: {missing}")

    c0, c1, m0 = suppress

    d = compute_ctx_score(
        cohort_df,
        purity_map,
        params=params,
        use_quality=False,
        mech_combine="product",
    )

    # ---- alternate allelic-evidence path: WGD-masked chromothripsis ----
    # Fires only when there is NO absolute LOH anywhere on the chromosome
    # (frac_loh ~ 0 -- consistent with a doubled genome retaining a
    # wild-type-like copy) but the minor-allele copy number is nonetheless
    # changing state repeatedly (n_allelic_state_changes), i.e. real allelic
    # complexity that never crosses to literal zero. Gated tightly (ramp
    # starts at 8 transitions, not the 1 used for the ordinary LOH<->het
    # interleave) so it doesn't just reproduce blanket F_allelic removal.
    p = params or ScoreParams()
    interleave_ramp = _ramp(pd.to_numeric(d["loh_het_interleave"], errors="coerce").fillna(0),
                            p.allelic_lo, p.allelic_hi)
    frac_loh = pd.to_numeric(d["frac_loh"], errors="coerce")
    n_changes = pd.to_numeric(d["n_allelic_state_changes"], errors="coerce")
    alt_credit = np.where(frac_loh < allelic_alt_frac_loh_max,
                          _ramp(n_changes, allelic_alt_lo, allelic_alt_hi), 0.0)
    het_ok = pd.to_numeric(d["n_het_probes"], errors="coerce").fillna(0) >= p.min_het_for_allelic
    allelic_raw = np.maximum(interleave_ramp, alt_credit)
    F_allelic_v7 = np.where(het_ok,
                            p.allelic_floor + (1 - p.allelic_floor) * allelic_raw,
                            p.allelic_floor)
    d["F_allelic_v7"] = F_allelic_v7
    mech = pd.to_numeric(d["F_osc"], errors="coerce").fillna(0.0) * F_allelic_v7

    conf = pd.to_numeric(d[probe_conf_col], errors="coerce").fillna(0.0)
    if second_look_df is not None and len(second_look_df):
        sl = second_look_df[[sample_col, chrom_col, "cp_osc_confidence_v7"]].copy()
        sl[chrom_col] = sl[chrom_col].astype(str)
        d[chrom_col] = d[chrom_col].astype(str)
        d = d.merge(sl, on=[sample_col, chrom_col], how="left")
        second_look = pd.to_numeric(d["cp_osc_confidence_v7"], errors="coerce").fillna(0.0)
        conf = np.maximum(conf, second_look)
    ramp = np.clip((conf - c0) / max(c1 - c0, 1e-9), 0.0, 1.0)
    d["ctx_suppress_v7"] = mech * (m0 + (1.0 - m0) * ramp)

    noise = pd.to_numeric(d[noise_col], errors="coerce")
    chain = pd.to_numeric(d[chain_col], errors="coerce").fillna(0)
    union = pd.to_numeric(d[footprint_union_col], errors="coerce").fillna(-np.inf)
    islands = pd.to_numeric(d[footprint_island_col], errors="coerce").fillna(0)

    d["nominated_score_v7"] = d["ctx_suppress_v7"] >= threshold
    d["min_osc_chain_pass_v7"] = chain >= min_osc_chain
    # noise_max=None defers the noise veto to apply_noise_veto, the last step of the pipeline
    d["noise_veto_v7"] = (noise.notna() & (noise > noise_max)) if noise_max is not None else False
    d["weak_altered_footprint_v7"] = (union < footprint_union_min) & (islands <= footprint_island_max)

    nominated_main = (
        d["nominated_score_v7"] &
        d["min_osc_chain_pass_v7"] &
        (~d["noise_veto_v7"]) &
        (~d["weak_altered_footprint_v7"])
    )

    # ---- parallel loss-dominant-baseline pathway ----
    # shard_df hits already carry their own deficit/footprint logic
    # (see compute_baseline_shard_features); they are gated by the SAME noise
    # ceiling as the main pathway, which rejects calls from noisy libraries.
    sd = shard_df[[sample_col, chrom_col, "loss_dominant_shard_hit"]].copy()
    sd[chrom_col] = sd[chrom_col].astype(str)
    d[chrom_col] = d[chrom_col].astype(str)
    d = d.merge(sd, on=[sample_col, chrom_col], how="left")
    d["loss_dominant_shard_hit"] = d["loss_dominant_shard_hit"].fillna(False)
    shard_pathway = d["loss_dominant_shard_hit"] & (~d["noise_veto_v7"])

    d["nominated_v7"] = nominated_main | shard_pathway
    d["sample_flag_v7"] = d.groupby(sample_col)["nominated_v7"].transform("any")

    d["v7_threshold"] = threshold
    d["v7_noise_max"] = noise_max
    d["v7_min_osc_chain"] = min_osc_chain
    d["v7_footprint_union_min"] = footprint_union_min
    d["v7_footprint_island_max"] = footprint_island_max
    d["v7_allelic_alt_lo"] = allelic_alt_lo
    d["v7_allelic_alt_hi"] = allelic_alt_hi
    d["v7_allelic_alt_frac_loh_max"] = allelic_alt_frac_loh_max

    return d


def _n_loss_segments(seg_chrom, tcn_col="tcn.em", lcn_col="lcn.em"):
    """
    Counts merged "loss" blocks on a chromosome, where a raw segment counts as
    loss if round(tcn.em)==1 OR round(lcn.em)==0 -- complete LOH at any total
    copy number, not just haploid loss (catches copy-neutral/amplified LOH
    too). Adjacent raw segments that are both loss-flagged merge into one
    block regardless of which condition triggered each, or whether their
    exact tcn/lcn differ.

    There is deliberately no cf.em reliability gate: segments separating loss
    blocks are often retained or amplified islands with an uncertain clonal
    fraction, and removing them before merging would fuse separate loss
    events.
    """
    s = _ordered_segments(seg_chrom)
    T = pd.to_numeric(s.get(tcn_col), errors="coerce")
    valid = T.notna()
    s = s.loc[valid].reset_index(drop=True)
    T = T.loc[valid].reset_index(drop=True)
    L = pd.to_numeric(s.get(lcn_col), errors="coerce") if lcn_col in s.columns else pd.Series(np.nan, index=s.index)

    n, prev_loss = 0, False
    for t, l in zip(T.tolist(), L.tolist()):
        loss = (int(round(t)) == 1) or (pd.notna(l) and int(round(l)) == 0)
        if loss and not prev_loss:
            n += 1
        prev_loss = loss
    return n


def apply_loss_segment_filter(cohort_df, man_df, *,
                              min_loss_segments=3,
                              nominated_col="nominated_v7",
                              shard_col="loss_dominant_shard_hit",
                              frac_loh_col="frac_loh",
                              n_allelic_col="n_allelic_state_changes",
                              alt_frac_loh_max_col="v7_allelic_alt_frac_loh_max",
                              alt_lo_col="v7_allelic_alt_lo",
                              sample_col="SAMPLE_ID",
                              chrom_col="chrom"):
    """
    First oscillation-call filter, applied to `apply_classifier_v7`'s output.

    Oscillation calls with fewer than three loss blocks (see
    `_n_loss_segments`) typically reflect simple broad-arm losses rather than
    the fragmented, multi-breakpoint pattern of chromothripsis; on sparse
    panel data, segment size is an unreliable proxy for genomic extent.

    For each currently-nominated, non-shard (sample, chrom), keeps the call
    only if EITHER `n_loss_segments >= min_loss_segments` (default 3) OR the
    row independently qualifies for the SAME WGD-masked allelic rescue
    `apply_classifier_v7` already uses (`frac_loh < v7_allelic_alt_frac_loh_max`
    AND `n_allelic_state_changes >= v7_allelic_alt_lo`) -- a chromosome with
    real minor-allele complexity substituting for literal LOH (chromothripsis
    after genome doubling) is not penalized for having no loss segments.
    Loss-dominant/shard-pathway calls are exempt.
    The filter only removes calls; it does not add chromosomes that v7 did not
    nominate.

    Returns `cohort_df` with new columns: `n_loss_segments`,
    `wgd_rescue_eligible`, `loss_segment_veto`, `nominated_final`,
    `sample_flag_final`, `loss_segment_min_count`.
    """
    import run_ctx_pipeline as rp

    d = cohort_df.copy()
    d[chrom_col] = d[chrom_col].astype(str)
    shard = d[shard_col].fillna(False) if shard_col in d.columns else pd.Series(False, index=d.index)
    nom = d[(d[nominated_col] == True) & (~shard)]

    rows = []
    cache = {}
    for sid in nom[sample_col].unique():
        if sid not in cache:
            try:
                seg_df = rp.load_cncf(man_df, sid)
                seg_df["_c"] = seg_df["chrom"].map(_norm_chrom)
                cache[sid] = seg_df
            except Exception:
                cache[sid] = None
        seg_df = cache[sid]
        if seg_df is None:
            continue
        for ch in nom.loc[nom[sample_col] == sid, chrom_col].unique():
            s = seg_df[seg_df["_c"] == _norm_chrom(ch)]
            rows.append({sample_col: sid, chrom_col: ch, "n_loss_segments": _n_loss_segments(s)})

    feat = pd.DataFrame(rows, columns=[sample_col, chrom_col, "n_loss_segments"])
    d = d.merge(feat, on=[sample_col, chrom_col], how="left")

    wgd_rescue = (d[frac_loh_col] < d[alt_frac_loh_max_col]) & (d[n_allelic_col] >= d[alt_lo_col])
    d["wgd_rescue_eligible"] = wgd_rescue.fillna(False)

    d_shard = d[shard_col].fillna(False) if shard_col in d.columns else pd.Series(False, index=d.index)
    veto = (~d_shard) & (d[nominated_col] == True) & (d["n_loss_segments"] < min_loss_segments) & (~d["wgd_rescue_eligible"])
    d["loss_segment_veto"] = veto.fillna(False)
    d["nominated_final"] = d[nominated_col] & (~d["loss_segment_veto"])
    d["sample_flag_final"] = d.groupby(sample_col)["nominated_final"].transform("any")
    d["loss_segment_min_count"] = min_loss_segments

    return d


# Features and fixed MAD floors for the within-sample contrast filter. The floors only stop a
# zero MAD (most chromosomes of a sample sharing the same count) from blowing up the z-score:
# one count step (1.4826 = MAD scale of a 1-unit spread) for the count features, and a
# pooled median per-sample MAD for segments_per_probe. They are fixed constants, so the
# filter does not depend on which samples are analysed together.
CONTRAST_FEATURE_FLOORS = {
    "cp_osc_chain": 1.4826,
    "cp_n_changepoints": 1.4826,
    "n_allelic_state_changes": 1.4826,
    "longest_osc_chain_k3": 1.4826,
    "n_altered_islands": 1.4826,
    "segments_per_probe": 0.0008620492439885324,
}


def compute_within_sample_contrast(cohort_df, *, features=CONTRAST_FEATURE_FLOORS,
                                   sample_col="SAMPLE_ID", min_other_chroms=5,
                                   z_clip=(-5.0, 10.0)):
    """
    Per chromosome, how much more complex it is than the SAME sample's other chromosomes.
    For each feature: z = (x - median(others)) / max(1.4826 * MAD(others), floor), where
    `others` are all of the sample's other chromosomes. The median/MAD of the other
    chromosomes is robust to a minority of shattered chromosomes, so multi-chromosome
    chromothripsis does not mask itself. contrast_score = mean of the clipped z-scores.
    Returns a Series aligned to cohort_df.index (NaN when < min_other_chroms are available).
    """
    d = cohort_df
    zs = []
    for f, floor in features.items():
        x = pd.to_numeric(d[f], errors="coerce").values
        z = np.full(len(d), np.nan)
        for _, idx in d.groupby(sample_col).indices.items():
            v = x[idx]
            for k, i in enumerate(idx):
                if np.isnan(v[k]):
                    continue
                others = np.delete(v, k)
                others = others[~np.isnan(others)]
                if len(others) < min_other_chroms:
                    continue
                med = np.median(others)
                mad = 1.4826 * np.median(np.abs(others - med))
                z[i] = (v[k] - med) / max(mad, floor)
        zs.append(pd.Series(z, index=d.index).clip(*z_clip))
    return pd.concat(zs, axis=1).mean(axis=1)


def apply_within_sample_contrast_filter(cohort_df, *, min_contrast_score=2.0,
                                        nominated_col="nominated_final",
                                        shard_col="loss_dominant_shard_hit",
                                        sample_col="SAMPLE_ID", chrom_col="chrom"):
    """
    Second oscillation-call filter, applied on top of `apply_loss_segment_filter`.

    In samples with genome-wide copy-number jitter, oscillation features can be elevated on
    every chromosome. A call is therefore kept only if the chromosome stands out from the
    other chromosomes of the same sample, rather than relying on a sample-level noise cutoff.

    Keeps an oscillation-pathway call only if contrast_score >= min_contrast_score (see
    `compute_within_sample_contrast`; features cp_osc_chain, cp_n_changepoints,
    n_allelic_state_changes, longest_osc_chain_k3, n_altered_islands, segments_per_probe).
    Loss-dominant/shard-pathway calls are exempt, because oscillation features do not
    describe them.

    Expects every chromosome row of each sample (not only nominated rows). Returns
    `cohort_df` with new columns `contrast_score`, `contrast_veto`, `contrast_min_score`, and
    `nominated_final` / `sample_flag_final` updated.
    """
    d = cohort_df.copy()
    d["contrast_score"] = compute_within_sample_contrast(d, sample_col=sample_col)
    shard = d[shard_col].fillna(False).astype(bool) if shard_col in d.columns else pd.Series(False, index=d.index)
    nominated = d[nominated_col].fillna(False).astype(bool)
    veto = nominated & (~shard) & ~(d["contrast_score"] >= min_contrast_score)
    d["contrast_veto"] = veto
    d["nominated_final"] = nominated & (~veto)
    d["sample_flag_final"] = d.groupby(sample_col)["nominated_final"].transform("any")
    d["contrast_min_score"] = min_contrast_score
    return d


def _median_cn_step(seg_chrom, tcn_col="tcn.em"):
    """Median nonzero |total copy-number difference| between adjacent segments of one
    chromosome (0 when the chromosome has no copy-number change)."""
    t = pd.to_numeric(_ordered_segments(seg_chrom)[tcn_col], errors="coerce").dropna().to_numpy(dtype=float)
    steps = np.abs(np.diff(t))
    steps = steps[steps > 0]
    return float(np.median(steps)) if len(steps) else 0.0


def apply_noise_veto(cohort_df, man_df, *, noise_max=1.0,
                     step_to_noise_min=3.0,
                     exception_contrast_min=1.5,
                     noise_col="dot_mad_within",
                     nominated_col="nominated_final",
                     sample_col="SAMPLE_ID", chrom_col="chrom"):
    """
    Probe-level noise veto, the last chromosome-level step (both pathways).

    A call is removed when the chromosome's within-segment probe noise (`noise_col`,
    absolute copy-number units) exceeds `noise_max`, unless the copy-number changes are
    large relative to that noise and the chromosome stands out from the rest of the sample:
        median nonzero |tcn.em difference| between adjacent segments / noise >= step_to_noise_min
        AND contrast_score >= exception_contrast_min.
    The exception keeps high-amplitude events in samples with genome-wide probe noise; the
    contrast condition rejects samples in which every chromosome shows large jumps.
    `step_to_noise_min=None` disables the exception (plain noise veto). Requires
    `contrast_score` (see `apply_within_sample_contrast_filter`) and `man_df` with
    `facets_cncf_path`.

    Returns `cohort_df` with new columns `cn_step_median`, `step_to_noise`, `noise_veto`,
    `noise_veto_exception`, `noise_max`, `step_to_noise_min`, `exception_contrast_min`, and
    `nominated_final` / `sample_flag_final` updated.
    """
    import run_ctx_pipeline as rp

    d = cohort_df.copy()
    d[chrom_col] = d[chrom_col].astype(str)
    noise = pd.to_numeric(d[noise_col], errors="coerce")
    nominated = d[nominated_col].fillna(False).astype(bool)
    noisy = nominated & noise.notna() & (noise > noise_max)

    d["cn_step_median"] = np.nan
    if step_to_noise_min is not None:
        for sid in d.loc[noisy, sample_col].unique():
            try:
                seg_df = rp.load_cncf(man_df, sid)
            except Exception:
                continue
            seg_df["_c"] = seg_df["chrom"].map(_norm_chrom)
            for i in d.index[noisy & (d[sample_col] == sid)]:
                s = seg_df[seg_df["_c"] == _norm_chrom(d.at[i, chrom_col])]
                d.at[i, "cn_step_median"] = _median_cn_step(s)
    d["step_to_noise"] = d["cn_step_median"] / noise
    if step_to_noise_min is not None:
        exception = noisy & (d["step_to_noise"] >= step_to_noise_min) & \
            (pd.to_numeric(d["contrast_score"], errors="coerce") >= exception_contrast_min)
    else:
        exception = pd.Series(False, index=d.index)
    d["noise_veto_exception"] = exception
    d["noise_veto"] = noisy & ~exception
    d["nominated_final"] = nominated & ~d["noise_veto"]
    d["sample_flag_final"] = d.groupby(sample_col)["nominated_final"].transform("any")
    d["noise_max"] = noise_max
    d["step_to_noise_min"] = step_to_noise_min
    d["exception_contrast_min"] = exception_contrast_min
    return d


def sample_level_counts(df, nominated_col="nominated_v3",
                        truth_col="truth_positive_chrom", sample_col="SAMPLE_ID"):
    """Sample-level triage scorecard from a per-(sample, chromosome) call frame.

    Positive sample = recovered if any of its truth chromosomes is nominated;
    negative sample = false positive if any chromosome is nominated.
    Returns: n_positive / recovered / missed / false_positive / errors (=missed+FP).
    """
    t = df[truth_col]
    t = t if t.dtype == bool else t.astype(str).str.strip().str.lower().isin(["true", "1", "yes"])
    g = (pd.DataFrame({sample_col: df[sample_col], "_pos": t,
                       "_nom": df[nominated_col].astype(bool), "_tp": t & df[nominated_col].astype(bool)})
         .groupby(sample_col).agg(pos=("_pos", "any"), nom=("_nom", "any"), rec=("_tp", "any")))
    P = int(g["pos"].sum())
    rec = int((g["pos"] & g["rec"]).sum())
    fp  = int((~g["pos"] & g["nom"]).sum())
    return pd.Series({"n_positive": P, "recovered": rec, "missed": P - rec,
                      "false_positive": fp, "errors": (P - rec) + fp})


def evaluate_classifier(df, nominated_col="nominated_v3",
                        truth_col="truth_positive_chrom", sample_col="SAMPLE_ID", label=""):
    """Chromosome- and sample-level performance, computed LIVE from `nominated_col`.
    Never reads a stored 'outcome' column. Returns {'chrom':..., 'sample':...}."""
    t = df[truth_col]
    t = t if t.dtype == bool else t.astype(str).str.strip().str.lower().isin(["true", "1", "yes"])
    nom = df[nominated_col].astype(bool)
    chrom = pd.Series({"caught":         int((t & nom).sum()),
                       "MISSED":         int((t & ~nom).sum()),
                       "false_positive": int((~t & nom).sum()),
                       "true_negative":  int((~t & ~nom).sum())})
    g = (pd.DataFrame({sample_col: df[sample_col], "_pos": t, "_nom": nom, "_tp": t & nom})
         .groupby(sample_col).agg(pos=("_pos","any"), nom=("_nom","any"), rec=("_tp","any")))
    P = int(g["pos"].sum()); rec = int((g["pos"] & g["rec"]).sum()); fp = int((~g["pos"] & g["nom"]).sum())
    sample = pd.Series({"n_positive": P, "recovered": rec, "missed": P - rec,
                        "false_positive": fp, "errors": (P - rec) + fp})
    if label: print(f"=== {label} ===")
    print("chromosome-level:\n" + chrom.to_string())
    print("\nsample-level:\n"   + sample.to_string() + "\n")
    return {"chrom": chrom, "sample": sample}


def add_outcome(df, nominated_col="nominated_v3", truth_col="truth_positive_chrom", out_col="outcome"):
    """Add/overwrite a per-chromosome outcome label derived LIVE from `nominated_col`:
    caught / MISSED / false_positive / true_negative."""
    t = df[truth_col]
    t = t if t.dtype == bool else t.astype(str).str.strip().str.lower().isin(["true", "1", "yes"])
    nom = df[nominated_col].astype(bool)
    df[out_col] = np.select([t & nom, t & ~nom, ~t & nom],
                            ["caught", "MISSED", "false_positive"], default="true_negative")
    return df



# ======================================================================
# Bootstrap CIs for sample-level comparison
# ----------------------------------------------------------------------
# With ~23 positive samples, a 1-2 sample difference between methods is noise.
# Resample SAMPLE_IDs with replacement, recompute the sample-level outcome each
# time, and report percentile CIs so method comparisons come with error bars.
# ======================================================================

def bootstrap_sample_eval(pred_df, score_col, threshold,
                          truth_col="truth_positive_chrom",
                          sample_col="SAMPLE_ID",
                          n_boot=2000, seed=0):
    """
    Returns a dict of point estimate + 95% percentile CI for pos_recovered,
    pos_missed, neg_false_positive, and sample_errors, resampling whole samples.
    """
    rng = np.random.default_rng(seed)
    d = pred_df.copy()
    d["_nom"] = pd.to_numeric(d[score_col], errors="coerce") >= threshold
    d["_truth"] = d[truth_col].astype(bool)

    # precompute per-sample truth/nominated chrom sets
    per_sample = {}
    for sid, g in d.groupby(sample_col):
        per_sample[sid] = (set(g.loc[g["_truth"], "chrom"]),
                           set(g.loc[g["_nom"], "chrom"]))
    sids = np.array(list(per_sample.keys()))

    def outcome(sid):
        truth, nom = per_sample[sid]
        if len(truth) > 0:
            return "recovered" if (nom & truth) else "missed"
        return "correct_negative" if len(nom) == 0 else "false_positive"

    base = {sid: outcome(sid) for sid in sids}

    def tally(sample):
        rec = miss = fp = 0
        for sid in sample:
            o = base[sid]
            rec += (o == "recovered"); miss += (o == "missed"); fp += (o == "false_positive")
        return rec, miss, fp, miss + fp

    point = tally(sids)
    boot = np.array([tally(rng.choice(sids, size=len(sids), replace=True))
                     for _ in range(n_boot)])
    names = ["pos_recovered", "pos_missed", "neg_false_positive", "sample_errors"]
    out = {"threshold": threshold, "score_col": score_col}
    for j, nm in enumerate(names):
        out[nm] = int(point[j])
        out[f"{nm}_lo"] = float(np.percentile(boot[:, j], 2.5))
        out[f"{nm}_hi"] = float(np.percentile(boot[:, j], 97.5))
    return out


def best_threshold_for_score(pred_df, score_col, thresholds=None, **kwargs):
    """Pick the threshold minimizing sample_errors (ties -> more recovered)."""
    if thresholds is None:
        thresholds = np.round(np.linspace(0.02, 0.98, 49), 3)
    best = None
    for t in thresholds:
        _, s = sample_level_eval(pred_df, score_col=score_col, threshold=t, **kwargs)
        key = (s["sample_errors"], -s["pos_recovered"])
        if best is None or key < best[0]:
            best = (key, t, s)
    return best[1], best[2]


def compare_scores(pred_df, score_cols, n_boot=2000, seed=0, **kwargs):
    """
    For each score column: find its best threshold, then bootstrap-CI the
    sample-level outcome there. Returns a tidy comparison DataFrame.
    """
    rows = []
    for col in score_cols:
        t, _ = best_threshold_for_score(pred_df, col, **kwargs)
        ci = bootstrap_sample_eval(pred_df, col, t, n_boot=n_boot, seed=seed, **kwargs)
        rows.append(ci)
    return pd.DataFrame(rows)
