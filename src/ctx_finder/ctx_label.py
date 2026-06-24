"""
Labeling and inspection helpers for the chr_ctx-positive cases.

Purpose: you have ~43 positive (sample, chromosome) pairs but have never looked
at them individually. This module gives them short, safe labels and lets you
inspect the actual copy-number / allelic profile so you can SEE what the
oscillation and LOH<->het interleaving look like in the real positives, and
judge where a sensible boundary is.

  make_sample_codes(scored)        -> {SAMPLE_ID -> 'C01', ...} for ctx-pos samples
  positive_register(scored, man_df)-> tidy table of the positive chromosomes,
                                       labeled and sorted by score (caught vs missed)
  plot_interleave(seg_df, jointseg_df, chrom, ...) -> per-chromosome inspection plot
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt



def load_centromeres(path="path/to/hg19.coordinates.csv"):
    cc = pd.read_csv(path)
    cc["_c"] = cc["chr"].astype(str).str.strip().str.replace("chr","",regex=False).str.replace("CHR","",regex=False)
    return {r["_c"]: (float(r["centro_start"]), float(r["centro_end"])) for _, r in cc.iterrows()}


def make_sample_codes(scored, sample_col="SAMPLE_ID",
                      truth_col="truth_positive_chrom"):
    """Map each chr_ctx-positive SAMPLE_ID to a short stable code C01..Cnn."""
    pos_samples = (scored.loc[scored[truth_col].astype(bool), sample_col]
                   .astype(str).drop_duplicates().sort_values().tolist())
    width = max(2, len(str(len(pos_samples))))
    return {sid: f"C{str(i + 1).zfill(width)}" for i, sid in enumerate(pos_samples)}


def positive_register(scored, man_df, threshold=None, score_col="ctx_score",
                      sample_col="SAMPLE_ID", truth_col="truth_positive_chrom",
                      extra_cols=()):
    """
    One row per chr_ctx-positive chromosome, with a safe code, the key features,
    the score, and (optionally) whether it clears `threshold`. Sorted by score
    descending so caught vs missed positives are easy to read off.
    """
    code_map = make_sample_codes(scored, sample_col, truth_col)
    d = scored.loc[scored[truth_col].astype(bool)].copy()
    pm = dict(zip(man_df[sample_col].astype(str),
                  pd.to_numeric(man_df["purity"], errors="coerce")))
    d["code"] = d[sample_col].astype(str).map(code_map)
    d["purity"] = d[sample_col].astype(str).map(pm)

    base_cols = [
        "code", "chrom", "purity", "n_probes", "n_het_probes",
        "longest_osc_chain_k3", "longest_osc_chain_k2", "probe_osc_chain",
        "probe_amp_spike_frac",
        "loh_het_interleave", "n_allelic_state_changes", "frac_loh",
        "oscillation_excess", "n_resolved_levels", "probe_noise",
        "F_osc", "F_allelic", "mech_evidence", "Q_weight", score_col,
    ]
    cols = [c for c in base_cols if c in d.columns] + list(extra_cols)
    out = d[cols].copy()
    if threshold is not None:
        out["nominated"] = out[score_col] >= threshold
        out["status"] = np.where(out["nominated"], "caught", "MISSED")
    return out.sort_values(score_col, ascending=False).reset_index(drop=True)


def _chrom_norm(x):
    return str(x).strip().replace("chr", "").replace("CHR", "")


def plot_interleave(seg_df, jointseg_df, chrom, code=None, centromeres=None,
                    point_col="tcn_estimate", tcn_col="tcn.em", lcn_col="lcn.em",
                    save_path=None, figsize=(13, 4)):
    """
    Inspect one chromosome: probe-level estimated total CN (dots) with segment
    total CN (red) and minor CN (blue) overlaid, in rank-x. LOH stretches
    (lcn==0) are shaded so you can SEE the LOH<->het interleaving that drives
    the allelic feature. x is probe rank within the chromosome. If `centromeres`
    is given, the centromere is shown as a shaded band when it contains probes,
    else a dashed line at its position; in-centromere probes are tinted gray.
    """
    c = _chrom_norm(chrom)
    seg = seg_df.copy(); pts = jointseg_df.copy()
    seg = seg.loc[seg["chrom"].map(_chrom_norm).eq(c)].sort_values("loc.start").reset_index(drop=True)
    pts = pts.loc[pts["chrom"].map(_chrom_norm).eq(c)].sort_values("maploc").reset_index(drop=True)
    if pts.shape[0] == 0:
        raise ValueError(f"no probes for chrom {chrom}")

    pts["_x"] = np.arange(pts.shape[0])
    def span(row):
        lo = np.searchsorted(pts["maploc"].values, row["loc.start"])
        hi = np.searchsorted(pts["maploc"].values, row["loc.end"])
        return max(lo, 0), min(hi, pts.shape[0] - 1)

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(pts["_x"], pts[point_col], s=8, c="black", alpha=0.35,
               linewidths=0, label=point_col)

    # centromere: tint in-centromere probes, shade band if any probes there, else mark with a dashed line
    has_centro = bool(centromeres and c in centromeres)
    centro_as_band = False
    if has_centro:
        cs, ce = centromeres[c]
        in_cen = pts["maploc"].between(cs, ce)
        if in_cen.any():
            ax.scatter(pts.loc[in_cen, "_x"], pts.loc[in_cen, point_col],
                       s=12, c="gray", alpha=0.7, zorder=3)
        cx0 = int(np.searchsorted(pts["maploc"].values, cs))
        cx1 = int(np.searchsorted(pts["maploc"].values, ce))
        if cx1 > cx0:
            ax.axvspan(cx0, cx1, color="gray", alpha=0.18, linewidth=0, zorder=0)
            centro_as_band = True
        else:
            ax.axvline(cx0, color="gray", ls="--", alpha=0.5, linewidth=1, zorder=0)

    for _, row in seg.iterrows():
        x0, x1 = span(row)
        if pd.notna(row.get(tcn_col)):
            ax.hlines(row[tcn_col], x0, x1, color="red", linewidth=2.2)
        if pd.notna(row.get(lcn_col)):
            ax.hlines(row[lcn_col], x0, x1, color="blue", linewidth=2.0)
            if row[lcn_col] == 0:    # LOH stretch
                ax.axvspan(x0, x1, color="orange", alpha=0.12, linewidth=0)

    # legend proxies
    from matplotlib.lines import Line2D
    import matplotlib.patches as mpatches
    handles = [
        Line2D([0], [0], color="red", lw=2.2, label=f"{tcn_col} (total CN)"),
        Line2D([0], [0], color="blue", lw=2.0, label=f"{lcn_col} (minor CN)"),
        Line2D([0], [0], marker="o", color="black", lw=0, alpha=0.5, label=point_col),
        Line2D([0], [0], color="orange", lw=8, alpha=0.3, label="LOH (minor=0)"),
    ]
    if has_centro:
        handles.append(
            mpatches.Patch(facecolor="gray", alpha=0.4, label="centromere")
            if centro_as_band else
            Line2D([0], [0], color="gray", ls="--", alpha=0.6, lw=1, label="centromere")
        )
    ax.legend(handles=handles, frameon=False, loc="upper right", fontsize=8)

    title = f"chr{c}" + (f"  [{code}]" if code else "")
    ax.set_title(f"{title}  | probes={pts.shape[0]}, segments={seg.shape[0]}")
    ax.set_xlabel("probe rank within chromosome")
    ax.set_ylabel("copy number")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
    return fig, ax


def false_positive_register(scored, man_df, threshold, score_col="ctx_score",
                            sample_col="SAMPLE_ID", truth_col="truth_positive_chrom",
                            extra_cols=()):
    """
    One row per FALSE-POSITIVE chromosome: a chr_ctx-negative-sample chromosome
    (or a non-truth chromosome in a positive sample) that scores >= threshold.
    Mirror of positive_register, for diagnosing where FPs come from.

    'fp_type' distinguishes:
      'neg_sample'  -> nomination in a sample with empty chr_ctx (true FP that
                       breaks specificity at the sample level)
      'pos_sample_wrong_chrom' -> nomination of a non-truth chromosome in a
                       chr_ctx-positive sample (does not change sample-level
                       recovery, but is an extra/incorrect chromosome call)
    """
    d = scored.copy()
    pm = dict(zip(man_df[sample_col].astype(str),
                  pd.to_numeric(man_df["purity"], errors="coerce")))
    sample_has_truth = d.groupby(sample_col)[truth_col].transform("any")
    d["_truth"] = d[truth_col].astype(bool)
    d["_nom"] = pd.to_numeric(d[score_col], errors="coerce") >= threshold
    fp = d.loc[d["_nom"] & ~d["_truth"]].copy()
    fp["fp_type"] = np.where(sample_has_truth.loc[fp.index], "pos_sample_wrong_chrom", "neg_sample")
    fp["purity"] = fp[sample_col].astype(str).map(pm)

    base_cols = [
        sample_col, "chrom", "fp_type", "purity", "n_probes", "n_het_probes",
        "longest_osc_chain_k3", "longest_osc_chain_k2", "probe_osc_chain",
        "loh_het_interleave", "n_allelic_state_changes", "frac_loh",
        "oscillation_excess", "n_resolved_levels", "probe_amp_spike_frac",
        "probe_noise", "F_osc", "F_allelic", "ctx_score",
    ]
    cols = [c for c in base_cols if c in fp.columns] + list(extra_cols)
    return fp[cols].sort_values(score_col, ascending=False).reset_index(drop=True)
