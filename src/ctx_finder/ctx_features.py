"""
Chromosome-level features, candidate nomination and filters for chromothripsis
detection on sparse FACETS (targeted-panel) copy-number data.

Inputs per sample:

    seg_df       : FACETS segment table (cncf); columns include chrom, loc.start,
                   loc.end, tcn.em, lcn.em
    jointseg_df  : FACETS probe-level table; columns include chrom, maploc, cnlr,
                   lorvar, het (tcn_estimate is added by run_ctx_pipeline)

All thresholds are collected in PARAMETERS (tunable; each has a command-line option in
pipeline/02_call.py) and SCORE_CONSTANTS (fixed shape of the oscillation score).
"""

import os

import numpy as np
import pandas as pd


# Tunable thresholds and their defaults (the values used for the published calls).
PARAMETERS = {
    # changepoint re-segmentation (pipeline/01_prepare.py)
    "changepoint_penalty": 8.0,         # PELT penalty multiplier (lambda)
    # CN oscillation-based candidates
    "score_threshold": 0.35,            # minimum oscillation score S
    "min_oscillating_stretch": 4,       # minimum oscillating stretch (segments)
    "footprint_min_fraction": 0.15,     # negligible footprint: altered fraction below this ...
    "footprint_max_islands": 3,         # ... and at most this many altered islands
    "post_wgd_min_changes": 8.0,        # post-WGD allelic evidence: minor-allele changes ramp start ...
    "post_wgd_full_changes": 15.0,      # ... and full credit
    "post_wgd_max_loh_fraction": 0.05,  # ... only when the LOH fraction is below this
    "second_look_score": 0.8,           # second look: S before changepoint weighting >= this ...
    "second_look_penalty": 5.0,         # ... and changepoint confidence 0 -> re-segment with this penalty
    # amplification-based candidates
    "focal_amp_multiple": 5.0,          # focal amplification: copy number >= this x max(baseline, 1)
    "focal_amp_max_fraction": 0.15,     # ... covering at most this fraction of the chromosome
    "loss_min_deficit": 0.5,            # chromosome-wide loss: baseline >= this many copies below ploidy
    "loss_min_baseline_fraction": 0.5,  # ... with the baseline covering >= this fraction
    "multi_focal_amp_min": 2,           # this many focal amplifications qualify on any baseline
    # filters of CN oscillation-based candidates
    "min_loss_blocks": 3,               # loss-segment filter
    "min_contrast_score": 2.0,          # within-sample contrast filter
    # probe-level noise veto (all candidates)
    "max_probe_noise": 1.0,             # copies
    "step_to_noise_min": 3.0,           # exception: median CN step / probe noise >= this ...
    "exception_contrast_min": 1.5,      # ... and contrast score >= this
}

# Fixed shape of the oscillation score (Methods; not exposed on the command line).
SCORE_CONSTANTS = {
    "max_states": 3,                    # distinct states allowed in an oscillating stretch
    "stretch_ramp": (2.0, 6.0),         # oscillating stretch -> F_osc
    "loh_switch_ramp": (1.0, 6.0),      # LOH <-> heterozygosity switches -> F_allelic
    "allelic_floor": 0.25,              # F_allelic when allelic evidence is weak or absent
    "min_het_probes": 4,                # heterozygous probes needed to use allelic evidence
    "confidence_ramp": (0.0, 3.0),      # changepoint confidence -> weight on S
    "changepoint_min_size": 5,          # PELT minimum segment size (probes)
}


# hg19 centromere + immediate pericentromeric-heterochromatin region, Mb.
# Approximate (UCSC cytoBand/gap track order of magnitude, generously padded
# on the flanks) -- used only to reject candidate focal amplifications that fall
# in this low-mappability, probe-desert zone, not for any coordinate-precise purpose.
# See compute_baseline_shard_features (centromeric-artifact rule).
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
    containing at most k DISTINCT copy-number states (the oscillating stretch).

    Oscillation A-B-A-B-A  (k=2) -> collapsed length 5  -> long.
    Staircase   A-B-C-D-E  (k=2) -> any window <=2       -> short (2).
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


# ----------------------------------------------------------------------
# per-(sample, chromosome) features from FACETS segments and probes
# ----------------------------------------------------------------------

def _probe_noise(values):
    """Within-chromosome probe noise of position-ordered probe copy-number estimates:
    robust SD of lag-1 successive differences (insensitive to genuine breakpoints)."""
    x = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)
    x = x[~np.isnan(x)]
    if x.size < 2:
        return np.nan
    dif = np.diff(x)
    return float(1.4826 * np.median(np.abs(dif - np.median(dif))) / np.sqrt(2))


def _altered_footprint(seg_chr, tcn_col="tcn.em"):
    """Altered fraction (summed length of segments outside the chromosome's most frequent
    state / total segment length) and number of contiguous altered islands."""
    seg_chr = _ordered_segments(seg_chr)
    starts = pd.to_numeric(seg_chr["loc.start"], errors="coerce").to_numpy(dtype=float)
    ends = pd.to_numeric(seg_chr["loc.end"], errors="coerce").to_numpy(dtype=float)
    lengths = np.clip(ends - starts, 0, None)
    T = pd.to_numeric(seg_chr.get(tcn_col), errors="coerce")
    t_states = [int(round(x)) if pd.notna(x) else np.nan for x in T]
    chr_len_sum = np.nansum(lengths)
    valid_t = [x for x in t_states if pd.notna(x)]
    modal_t = pd.Series(valid_t).mode().iloc[0] if len(valid_t) else np.nan
    altered = np.array([pd.notna(x) and x != modal_t for x in t_states], dtype=bool)
    altered_union = float(np.nansum(lengths[altered])) if len(altered) else np.nan
    if altered.any():
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
    else:
        n_altered_islands = 0
    return dict(
        altered_union_frac=altered_union / chr_len_sum if chr_len_sum and chr_len_sum > 0 else np.nan,
        n_altered_islands=n_altered_islands,
    )


def chrom_features_one_sample(seg_df, jointseg_df,
                              tcn_col="tcn.em", lcn_col="lcn.em",
                              point_col="tcn_estimate"):
    """
    One row per chromosome: oscillating stretch, allelic features, altered footprint,
    probe noise and probe support.
    """
    seg = seg_df.copy()
    pts = jointseg_df.copy()
    seg["_c"] = seg["chrom"].map(_norm_chrom)
    pts["_c"] = pts["chrom"].map(_norm_chrom)

    rows = []
    for c, seg_chr in seg.groupby("_c"):
        seg_chr = _ordered_segments(seg_chr)
        pts_chr = pts.loc[pts["_c"].eq(c)].sort_values("maploc")

        n_probes = int(pts_chr.shape[0])
        het = pts_chr["het"] if "het" in pts_chr.columns else pd.Series([], dtype=float)
        n_het = int(np.nansum(het.astype(float))) if len(het) else 0

        # oscillating stretch: longest run of segments confined to <= max_states integer states
        T = pd.to_numeric(seg_chr.get(tcn_col), errors="coerce")
        L = pd.to_numeric(seg_chr.get(lcn_col), errors="coerce")
        int_states = [int(round(t)) for t in T.tolist() if not pd.isna(t)]

        # allelic / LOH structure
        loh = (L == 0).astype(float)                   # 1 where minor CN == 0
        rows.append(dict(
            chrom=c,
            _chrom_order=_chrom_order(c),
            n_probes=n_probes,
            n_het_probes=n_het,
            segments_per_probe=(seg_chr.shape[0] / n_probes) if n_probes else np.nan,
            longest_osc_chain_k3=longest_le_k_distinct_run(int_states, SCORE_CONSTANTS["max_states"]),
            dot_mad_within=_probe_noise(pts_chr[point_col]) if point_col in pts_chr.columns else np.nan,
            n_allelic_state_changes=_transitions(L.tolist()),
            loh_het_interleave=_transitions(loh.tolist()),  # LOH <-> het switches
            frac_loh=float(np.nanmean(loh)) if len(loh) else np.nan,
            **_altered_footprint(seg_chr, tcn_col=tcn_col),
        ))

    return pd.DataFrame(rows).sort_values("_chrom_order").reset_index(drop=True)


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


def _median_cn_step(seg_chrom, tcn_col="tcn.em"):
    """Median nonzero |total copy-number difference| between adjacent segments of one
    chromosome (0 when the chromosome has no copy-number change)."""
    t = pd.to_numeric(_ordered_segments(seg_chrom)[tcn_col], errors="coerce").dropna().to_numpy(dtype=float)
    steps = np.abs(np.diff(t))
    steps = steps[steps > 0]
    return float(np.median(steps)) if len(steps) else 0.0


def segment_block_features(seg_df):
    """Per chromosome: number of loss blocks (loss-segment filter) and median copy-number
    step between adjacent segments (noise-veto exception)."""
    seg = seg_df.copy()
    seg["_c"] = seg["chrom"].map(_norm_chrom)
    return pd.DataFrame([dict(chrom=c, n_loss_segments=_n_loss_segments(g), cn_step_median=_median_cn_step(g))
                         for c, g in seg.groupby("_c")])


def compute_baseline_shard_features(seg_df, ploidy, *,
                                    shard_baseline_mult=5.0,
                                    deficit_min=0.5,
                                    baseline_frac_min=0.5,
                                    small_frac_max=0.15,
                                    multi_shard_min=2,
                                    tcn_col="tcn.em",
                                    centromere_mb=CENTROMERE_HG19_MB):
    """
    Per-chromosome features for amplification-based candidates: a chromosome-wide
    copy-number loss carrying embedded focal amplifications.

    Chromothripsis can produce high-level, often extrachromosomal, amplification of a
    few fragments while most of the chromosome is lost. On panel data this appears as
    a few amplified segments rather than as oscillation, and oscillation features
    (built relative to the chromosome's own modal state) see nothing "altered" -- the
    lost state IS the local mode. Comparing the chromosome's dominant state against the
    sample's genome-wide ploidy surfaces these cases. Requiring a small, highly
    amplified fragment (a focal amplification) keeps the signal specific to
    chromothripsis rather than an ordinary whole-chromosome deletion, which has no such
    internal structure.

    The focal-amplification threshold is relative to the chromosome's OWN baseline_cn
    (floored at 1 to avoid a vacuous threshold when baseline_cn is itself 0, e.g. a
    homozygous-deletion backdrop), not the sample's genome-wide ploidy -- a low-ploidy
    sample can otherwise set an unreachably high bar for a chromosome whose local
    backdrop is already a deep loss. A focal amplification must also be genomically
    EMBEDDED in the baseline (a baseline-state segment somewhere upstream AND somewhere
    downstream of it), not merely small; a small high segment flanked by an unrelated
    gain rather than the baseline is not counted.

    A candidate is also rejected if it overlaps CENTROMERE_HG19_MB for its
    chromosome: FACETS often emits a single sparse-probe segment spanning the
    low-mappability pericentromeric region, which recurs at the same
    coordinates across unrelated samples and is a panel artifact rather than
    an amplified fragment. Only the centromeric candidate is removed, so
    genuine focal amplifications elsewhere on the same chromosome are unaffected.

    Multiple focal amplifications: when a chromosome carries >= multi_shard_min
    embedded focal amplifications, the repeated alternation between baseline and
    amplified states itself resembles copy-number oscillation, so the chromosome-wide-loss
    requirements (deficit_min, baseline_frac_min) are waived. Every focal-amplification
    requirement still applies (>= shard_baseline_mult x baseline, <= small_frac_max of the
    chromosome, embedded, non-centromeric), as does the probe-level noise veto downstream.
    Pass multi_shard_min=None to disable.

    Returns one row per chromosome: baseline_cn (the length-weighted modal total-CN
    state), baseline_frac (how much of the chromosome that state covers), deficit
    (ploidy - baseline_cn), n_shards (embedded focal amplifications), chromosome_wide_loss
    (deficit and baseline_frac criteria met), and loss_dominant_shard_hit (the
    amplification-based candidate flag used downstream).
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
        if g.empty or chrom_len <= 0 or pd.isna(ploidy):
            rows.append(dict(chrom=c, baseline_cn=np.nan, baseline_frac=np.nan, deficit=np.nan,
                             n_shards=0, chromosome_wide_loss=False, loss_dominant_shard_hit=False))
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
        chromosome_wide_loss = (deficit >= deficit_min) and (baseline_frac >= baseline_frac_min)
        loss_dominant_hit = chromosome_wide_loss and (len(shards) >= 1)
        multi_shard_hit = multi_shard_min is not None and len(shards) >= multi_shard_min

        rows.append(dict(chrom=c, baseline_cn=baseline_cn, baseline_frac=round(baseline_frac, 4),
                         deficit=round(deficit, 3), n_shards=int(len(shards)),
                         chromosome_wide_loss=bool(chromosome_wide_loss),
                         loss_dominant_shard_hit=bool(loss_dominant_hit or multi_shard_hit)))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# oscillation score and candidate nomination
# ----------------------------------------------------------------------

def compute_ctx_score(cohort_df):
    """
    Add F_osc, F_allelic and ctx_score = F_osc * F_allelic (oscillation score before
    changepoint weighting) to a cohort feature table.
        F_osc     = ramp(oscillating stretch; 2, 6)
        F_allelic = floor + (1 - floor) * ramp(LOH <-> het switches; 1, 6),
                    floor alone when fewer than min_het_probes heterozygous probes
    """
    k = SCORE_CONSTANTS
    d = cohort_df.copy()
    chain = pd.to_numeric(d["longest_osc_chain_k3"], errors="coerce").fillna(0)
    F_osc = _ramp(chain, *k["stretch_ramp"])
    interleave = _ramp(pd.to_numeric(d["loh_het_interleave"], errors="coerce").fillna(0), *k["loh_switch_ramp"])
    het_ok = pd.to_numeric(d["n_het_probes"], errors="coerce").fillna(0) >= k["min_het_probes"]
    F_allelic = np.where(het_ok, k["allelic_floor"] + (1 - k["allelic_floor"]) * interleave, k["allelic_floor"])
    d["F_osc"] = F_osc
    d["F_allelic"] = F_allelic
    d["ctx_score"] = F_osc * F_allelic
    return d


def compute_second_look_confidence(cohort_df, man_df, probe_dir, *,
                                   score_col="ctx_score",
                                   conf_col="cp_osc_confidence",
                                   score_min=0.8,
                                   pen_fallback=5.0,
                                   sample_col="SAMPLE_ID",
                                   chrom_col="chrom"):
    """
    Second look: the changepoint re-segmentation runs at a stiff penalty
    (changepoint_penalty) so that noisy libraries do not generate spurious
    changepoints. That stiffness can give zero changepoint confidence to a chromosome
    with a maximal FACETS-level score whose real, well-separated jumps still do not
    individually clear the per-breakpoint cost. For that narrow case only
    (ctx_score >= score_min AND cp_osc_confidence == 0), the chromosome is
    re-segmented at pen_fallback and the higher confidence is used.

    `cohort_df` must already carry `score_col` (call after compute_ctx_score). Reads
    the cached probe tables from `probe_dir`. Returns [sample_col, chrom_col,
    "cp_osc_confidence_v7"].
    """
    import ctx_changepoint as ccp

    gate = (pd.to_numeric(cohort_df[score_col], errors="coerce") >= score_min) & \
           (pd.to_numeric(cohort_df[conf_col], errors="coerce") == 0.0)
    cand = cohort_df.loc[gate, [sample_col, chrom_col]].drop_duplicates()

    rows = []
    for _, row in cand.iterrows():
        sid, chrom = row[sample_col], row[chrom_col]
        jpath = os.path.join(probe_dir, f"{sid}.hisens.jointseg.tsv")
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


def apply_classifier_v7(cohort_df, shard_df, *, second_look_df=None,
                        threshold=0.35,
                        min_osc_chain=4,
                        footprint_union_min=0.15,
                        footprint_island_max=3,
                        allelic_alt_lo=8.0,
                        allelic_alt_hi=15.0,
                        allelic_alt_frac_loh_max=0.05,
                        chain_col="longest_osc_chain_k3",
                        probe_conf_col="cp_osc_confidence",
                        footprint_union_col="altered_union_frac",
                        footprint_island_col="n_altered_islands",
                        sample_col="SAMPLE_ID",
                        chrom_col="chrom"):
    """
    Nominates both candidate types (no filters, no noise veto).

    CN oscillation-based: S = F_osc * F_allelic_v7 * ramp(changepoint confidence; 0, 3)
    >= threshold, oscillating stretch >= min_osc_chain, and not a negligible altered
    footprint (altered fraction < footprint_union_min AND <= footprint_island_max islands).
    F_allelic_v7 additionally credits repeated minor-allele changes without LOH
    (post-WGD allelic evidence). Changepoint confidence is the larger of the
    regular and the second-look confidence.

    Amplification-based: `shard_df` (compute_baseline_shard_features per sample, with
    SAMPLE_ID) loss_dominant_shard_hit.

    Adds ctx_suppress_v7 (= S), changepoint_confidence, nominated_oscillation,
    loss_dominant_shard_hit and nominated_v7 (either type).
    """
    required = [probe_conf_col, chain_col, footprint_union_col, footprint_island_col,
                "loh_het_interleave", "frac_loh", "n_allelic_state_changes"]
    missing = [c for c in required if c not in cohort_df.columns]
    if missing:
        raise KeyError(f"Missing required columns: {missing}")
    k = SCORE_CONSTANTS

    d = compute_ctx_score(cohort_df)

    # ---- post-WGD allelic evidence ----
    # Fires only when there is NO absolute LOH anywhere on the chromosome
    # (frac_loh ~ 0 -- consistent with a doubled genome retaining a
    # wild-type-like copy) but the minor-allele copy number is nonetheless
    # changing state repeatedly (n_allelic_state_changes), i.e. real allelic
    # complexity that never crosses to literal zero. Gated tightly (ramp
    # starts at 8 changes, not the 1 used for the ordinary LOH<->het
    # switches).
    interleave_ramp = _ramp(pd.to_numeric(d["loh_het_interleave"], errors="coerce").fillna(0), *k["loh_switch_ramp"])
    frac_loh = pd.to_numeric(d["frac_loh"], errors="coerce")
    n_changes = pd.to_numeric(d["n_allelic_state_changes"], errors="coerce")
    alt_credit = np.where(frac_loh < allelic_alt_frac_loh_max,
                          _ramp(n_changes, allelic_alt_lo, allelic_alt_hi), 0.0)
    het_ok = pd.to_numeric(d["n_het_probes"], errors="coerce").fillna(0) >= k["min_het_probes"]
    allelic_raw = np.maximum(interleave_ramp, alt_credit)
    F_allelic_v7 = np.where(het_ok, k["allelic_floor"] + (1 - k["allelic_floor"]) * allelic_raw, k["allelic_floor"])
    d["F_allelic_v7"] = F_allelic_v7
    mech = pd.to_numeric(d["F_osc"], errors="coerce").fillna(0.0) * F_allelic_v7

    conf = pd.to_numeric(d[probe_conf_col], errors="coerce").fillna(0.0)
    d[chrom_col] = d[chrom_col].astype(str)
    if second_look_df is not None and len(second_look_df):
        sl = second_look_df[[sample_col, chrom_col, "cp_osc_confidence_v7"]].copy()
        sl[chrom_col] = sl[chrom_col].astype(str)
        d = d.merge(sl, on=[sample_col, chrom_col], how="left")
        second_look = pd.to_numeric(d["cp_osc_confidence_v7"], errors="coerce").fillna(0.0)
        conf = np.maximum(conf, second_look)
    d["changepoint_confidence"] = conf
    c0, c1 = k["confidence_ramp"]
    d["ctx_suppress_v7"] = mech * np.clip((conf - c0) / max(c1 - c0, 1e-9), 0.0, 1.0)

    chain = pd.to_numeric(d[chain_col], errors="coerce").fillna(0)
    union = pd.to_numeric(d[footprint_union_col], errors="coerce").fillna(-np.inf)
    islands = pd.to_numeric(d[footprint_island_col], errors="coerce").fillna(0)
    weak_footprint = (union < footprint_union_min) & (islands <= footprint_island_max)
    d["nominated_oscillation"] = (d["ctx_suppress_v7"] >= threshold) & (chain >= min_osc_chain) & (~weak_footprint)

    # ---- amplification-based candidates ----
    sd = shard_df[[sample_col, chrom_col, "loss_dominant_shard_hit"]].copy()
    sd[chrom_col] = sd[chrom_col].astype(str)
    d = d.merge(sd, on=[sample_col, chrom_col], how="left")
    d["loss_dominant_shard_hit"] = d["loss_dominant_shard_hit"].fillna(False).astype(bool)

    d["nominated_v7"] = d["nominated_oscillation"] | d["loss_dominant_shard_hit"]
    d["post_wgd_allelic"] = (pd.to_numeric(d["frac_loh"], errors="coerce") < allelic_alt_frac_loh_max) & \
        (pd.to_numeric(d["n_allelic_state_changes"], errors="coerce") >= allelic_alt_lo)
    return d


# ----------------------------------------------------------------------
# filters (CN oscillation-based candidates) and probe-level noise veto (all)
# ----------------------------------------------------------------------

def apply_loss_segment_filter(cohort_df, *, min_loss_segments=3,
                              nominated_col="nominated_v7",
                              shard_col="loss_dominant_shard_hit",
                              sample_col="SAMPLE_ID"):
    """
    First filter of CN oscillation-based candidates (`n_loss_segments`, see
    `_n_loss_segments`). Calls with fewer than `min_loss_segments` loss blocks typically
    reflect simple broad-arm losses rather than the fragmented, multi-breakpoint pattern
    of chromothripsis. A call is kept if it has enough loss blocks OR qualifies for the
    post-WGD allelic evidence (`post_wgd_allelic`), where losses to LOH are not expected.
    Amplification-based candidates are exempt. Adds `loss_segment_veto` and
    `nominated_final`.
    """
    d = cohort_df.copy()
    shard = d[shard_col].fillna(False).astype(bool)
    veto = (~shard) & (d[nominated_col] == True) & (d["n_loss_segments"] < min_loss_segments) & \
        (~d["post_wgd_allelic"].astype(bool))
    d["loss_segment_veto"] = veto.fillna(False)
    d["nominated_final"] = d[nominated_col] & (~d["loss_segment_veto"])
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
                                        sample_col="SAMPLE_ID"):
    """
    Second filter of CN oscillation-based candidates. In samples with genome-wide
    copy-number scatter, oscillation features can be elevated on every chromosome, so a
    call is kept only if the chromosome stands out from the same sample's other
    chromosomes (contrast_score >= min_contrast_score; see compute_within_sample_contrast).
    Amplification-based candidates are exempt. Expects every chromosome row of each
    sample. Adds `contrast_score` and `contrast_veto`; updates `nominated_final`.
    """
    d = cohort_df.copy()
    d["contrast_score"] = compute_within_sample_contrast(d, sample_col=sample_col)
    shard = d[shard_col].fillna(False).astype(bool)
    nominated = d[nominated_col].fillna(False).astype(bool)
    veto = nominated & (~shard) & ~(d["contrast_score"] >= min_contrast_score)
    d["contrast_veto"] = veto
    d["nominated_final"] = nominated & (~veto)
    return d


def apply_noise_veto(cohort_df, *, noise_max=1.0,
                     step_to_noise_min=3.0,
                     exception_contrast_min=1.5,
                     noise_col="dot_mad_within",
                     nominated_col="nominated_final"):
    """
    Probe-level noise veto, the last chromosome-level step (all candidates).

    A call is removed when the chromosome's probe noise (`noise_col`, copy-number units)
    exceeds `noise_max`, unless the copy-number changes are large relative to that noise
    and the chromosome stands out from the rest of the sample:
        cn_step_median / noise >= step_to_noise_min AND contrast_score >= exception_contrast_min.
    The exception keeps high-amplitude events in samples with genome-wide probe noise; the
    contrast condition rejects samples in which every chromosome shows large jumps.
    `step_to_noise_min=None` disables the exception. Requires `cn_step_median` and
    `contrast_score`. Adds `step_to_noise`, `noise_veto`, `noise_veto_exception`; updates
    `nominated_final`.
    """
    d = cohort_df.copy()
    noise = pd.to_numeric(d[noise_col], errors="coerce")
    nominated = d[nominated_col].fillna(False).astype(bool)
    noisy = nominated & noise.notna() & (noise > noise_max)
    d["step_to_noise"] = pd.to_numeric(d["cn_step_median"], errors="coerce") / noise
    if step_to_noise_min is not None:
        exception = noisy & (d["step_to_noise"] >= step_to_noise_min) & \
            (pd.to_numeric(d["contrast_score"], errors="coerce") >= exception_contrast_min)
    else:
        exception = pd.Series(False, index=d.index)
    d["noise_veto_exception"] = exception
    d["noise_veto"] = noisy & ~exception
    d["nominated_final"] = nominated & ~d["noise_veto"]
    return d
