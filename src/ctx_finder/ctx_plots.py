"""
Diagnostic scatter plots for the chromothripsis-footprint analysis.

Two entry points:

  plot_quality_landscape(scored, man_df, ...)
      Where do the chr_ctx-positive chromosomes sit in purity x quality space,
      relative to all other chromosomes? Directly visualizes the confound that
      positives are enriched in lower-quality / lower-purity samples.

  plot_score_diagnostics(scored, man_df, score_cols, thresholds, ...)
      For each score: score (y) vs sample purity (x), truth highlighted,
      threshold line drawn. Shows separation AND whether missed positives are
      the low-purity ones.

Both save a PNG and return (fig, axes). Agg backend; safe on a headless node.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _purity_map(man_df):
    return dict(zip(man_df["SAMPLE_ID"].astype(str),
                    pd.to_numeric(man_df["purity"], errors="coerce")))


def _attach_sample_cols(scored, man_df):
    d = scored.copy()
    pm = _purity_map(man_df)
    d["purity"] = d["SAMPLE_ID"].astype(str).map(pm)
    # sample is chr_ctx-positive if it has any truth-positive chromosome
    sample_pos = d.groupby("SAMPLE_ID")["truth_positive_chrom"].transform("any")
    d["sample_is_positive"] = sample_pos
    return d


def plot_quality_landscape(scored, man_df, save_path=None,
                           y_metric="probe_noise", ylog=False,
                           figsize=(13, 5)):
    """
    Two panels:
      (A) per-chromosome: x = sample purity, y = quality metric (default probe
          noise; try 'n_probes' or 'Q_weight'). Gray = not flagged,
          red = chr_ctx truth-positive chromosome.
      (B) per-sample: purity distribution, chr_ctx-positive vs negative samples
          (strip plot with jitter), to show if positives skew low-purity.
    """
    d = _attach_sample_cols(scored, man_df)

    fig, axes = plt.subplots(1, 2, figsize=figsize,
                             gridspec_kw={"width_ratios": [2.2, 1]})

    # ---- panel A: purity x quality, truth highlighted ----
    ax = axes[0]
    bg = d.loc[~d["truth_positive_chrom"].astype(bool)]
    fg = d.loc[d["truth_positive_chrom"].astype(bool)]
    ax.scatter(bg["purity"], pd.to_numeric(bg[y_metric], errors="coerce"),
               s=14, c="lightgray", alpha=0.45, linewidths=0, label="not flagged")
    ax.scatter(fg["purity"], pd.to_numeric(fg[y_metric], errors="coerce"),
               s=46, c="crimson", alpha=0.9, edgecolors="black", linewidths=0.4,
               label="chr_ctx truth")
    if ylog:
        ax.set_yscale("log")
    ax.set_xlabel("sample purity")
    ax.set_ylabel(y_metric)
    ax.set_title(f"Ground-truth chromosomes in purity x {y_metric} space")
    ax.legend(frameon=False, loc="best")
    ax.spines[["top", "right"]].set_visible(False)

    # ---- panel B: purity by sample status (one point per sample) ----
    ax = axes[1]
    per_sample = (d.groupby("SAMPLE_ID")
                    .agg(purity=("purity", "first"),
                         is_pos=("sample_is_positive", "first"))
                    .reset_index())
    rng = np.random.default_rng(0)
    for cat, color, label in [(False, "0.6", "chr_ctx-neg sample"),
                              (True, "crimson", "chr_ctx-pos sample")]:
        sub = per_sample.loc[per_sample["is_pos"].eq(cat)]
        x = (1 if cat else 0) + rng.uniform(-0.12, 0.12, size=len(sub))
        ax.scatter(x, sub["purity"], s=30, c=color, alpha=0.75,
                   linewidths=0, label=label)
        if len(sub):
            ax.hlines(np.nanmedian(sub["purity"]),
                      (1 if cat else 0) - 0.2, (1 if cat else 0) + 0.2,
                      color="black", linewidth=1.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["neg", "pos"])
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylabel("sample purity")
    ax.set_title("Purity by sample status\n(black = median)")
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
    return fig, axes


def plot_score_diagnostics(scored, man_df, score_cols, thresholds=None,
                           x_metric="purity", save_path=None, ncols=2,
                           figsize=None):
    """
    One subplot per score: y = score, x = sample purity (default).
    Gray = not flagged, red = chr_ctx truth. Dashed line at the score's
    threshold (pass a dict {score_col: threshold}). Lets you see whether the
    positives separate, and whether any missed positive is simply low-purity.
    """
    d = _attach_sample_cols(scored, man_df)
    thresholds = thresholds or {}
    n = len(score_cols)
    nrows = int(np.ceil(n / ncols))
    figsize = figsize or (6.5 * ncols, 4.2 * nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    axes = axes.ravel()

    for ax, col in zip(axes, score_cols):
        x = pd.to_numeric(d[x_metric], errors="coerce")
        y = pd.to_numeric(d[col], errors="coerce")
        truth = d["truth_positive_chrom"].astype(bool)
        ax.scatter(x[~truth], y[~truth], s=14, c="lightgray", alpha=0.45,
                   linewidths=0, label="not flagged")
        ax.scatter(x[truth], y[truth], s=46, c="crimson", alpha=0.9,
                   edgecolors="black", linewidths=0.4, label="chr_ctx truth")
        if col in thresholds:
            ax.axhline(thresholds[col], color="navy", linestyle="--",
                       linewidth=1.2, label=f"threshold={thresholds[col]:.2f}")
        ax.set_xlabel(x_metric)
        ax.set_ylabel(col)
        ax.set_title(col)
        ax.legend(frameon=False, fontsize=8, loc="best")
        ax.spines[["top", "right"]].set_visible(False)

    for ax in axes[n:]:
        ax.set_visible(False)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
    return fig, axes


def plot_calls_scatter(scored, man_df, threshold, score_col="ctx_score",
                       x_metric="longest_osc_chain_k3", y_metric="loh_het_interleave",
                       save_path=None, figsize=(10, 8), annotate_missed=True,
                       sample_col="SAMPLE_ID", truth_col="truth_positive_chrom"):
    """
    Chromosome-level call scatter in feature space, colored by the 4 outcomes
    given a score threshold:
      gray   = true negative (not flagged, scored below threshold)
      blue   = false positive (not flagged, scored >= threshold)
      green  = caught positive (flagged, scored >= threshold)
      red    = MISSED positive (flagged, scored < threshold)
    Missed positives are annotated with their sample code + chrom.
    """
    import numpy as np
    d = scored.copy()
    pm = dict(zip(man_df[sample_col].astype(str),
                  pd.to_numeric(man_df["purity"], errors="coerce")))
    code_map = None
    try:
        from ctx_label import make_sample_codes
        code_map = make_sample_codes(d, sample_col, truth_col)
    except Exception:
        pass
    d["_truth"] = d[truth_col].astype(bool)
    d["_call"] = pd.to_numeric(d[score_col], errors="coerce") >= threshold

    def cat(r):
        if r["_truth"] and r["_call"]:   return "caught"
        if r["_truth"] and not r["_call"]: return "MISSED"
        if (not r["_truth"]) and r["_call"]: return "false_pos"
        return "true_neg"
    d["_cat"] = d.apply(cat, axis=1)

    style = {"true_neg": ("lightgray", 14, 0.4, "true negative"),
             "false_pos": ("tab:blue", 48, 0.85, "false positive"),
             "caught": ("seagreen", 48, 0.85, "caught positive"),
             "MISSED": ("crimson", 90, 0.95, "MISSED positive")}
    fig, ax = plt.subplots(figsize=figsize)
    jit = np.random.default_rng(0)
    for c in ["true_neg", "false_pos", "caught", "MISSED"]:
        sub = d.loc[d["_cat"].eq(c)]
        if sub.shape[0] == 0:
            continue
        col, s, a, lab = style[c]
        xj = pd.to_numeric(sub[x_metric], errors="coerce") + jit.uniform(-0.15, 0.15, len(sub))
        yj = pd.to_numeric(sub[y_metric], errors="coerce") + jit.uniform(-0.15, 0.15, len(sub))
        ax.scatter(xj, yj, s=s, c=col, alpha=a, linewidths=0.4 if c != "true_neg" else 0,
                   edgecolors="black" if c != "true_neg" else "none",
                   label=f"{lab} (n={sub.shape[0]})")
    if annotate_missed:
        for _, r in d.loc[d["_cat"].eq("MISSED")].iterrows():
            code = code_map.get(str(r[sample_col]), str(r[sample_col])[:6]) if code_map else str(r[sample_col])[:6]
            ax.annotate(f"{code}:chr{r['chrom']}",
                        (pd.to_numeric(r[x_metric], errors="coerce"),
                         pd.to_numeric(r[y_metric], errors="coerce")),
                        fontsize=8, xytext=(5, 5), textcoords="offset points", color="crimson")
    ax.set_xlabel(x_metric); ax.set_ylabel(y_metric)
    ax.set_title(f"Calls in feature space (threshold={threshold}, score={score_col})")
    ax.legend(frameon=False, loc="best")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight"); plt.close(fig)
    return fig, ax


def plot_score_vs_rank(scored, threshold, score_col="ctx_score",
                       truth_col="truth_positive_chrom", save_path=None, figsize=(12, 5)):
    """
    Every chromosome ranked by score (x = rank, y = score). Truth in red,
    negatives gray, threshold line drawn. Shows the global separation and
    exactly where positives fall relative to the cut.
    """
    import numpy as np
    d = scored.copy()
    d["_s"] = pd.to_numeric(d[score_col], errors="coerce")
    d = d.sort_values("_s", ascending=False).reset_index(drop=True)
    d["_rank"] = np.arange(len(d))
    d["_truth"] = d[truth_col].astype(bool)
    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(d.loc[~d["_truth"], "_rank"], d.loc[~d["_truth"], "_s"],
               s=10, c="lightgray", alpha=0.5, linewidths=0, label="not flagged")
    ax.scatter(d.loc[d["_truth"], "_rank"], d.loc[d["_truth"], "_s"],
               s=40, c="crimson", alpha=0.9, edgecolors="black", linewidths=0.4, label="chr_ctx truth")
    ax.axhline(threshold, color="navy", linestyle="--", linewidth=1.2, label=f"threshold={threshold}")
    ax.set_xlabel("chromosome rank by score"); ax.set_ylabel(score_col)
    ax.set_title("Score-ranked chromosomes, truth highlighted")
    ax.legend(frameon=False, loc="best")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight"); plt.close(fig)
    return fig, ax
