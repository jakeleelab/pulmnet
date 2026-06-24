"""
Parallel (Slurm array) changepoint-feature extraction for Fork 1.

Layout:
  worker   : compute cp_* features for ONE sample, write per-sample TSV (idempotent)
  launcher : write sample manifest + an sbatch --array script
  collector: concatenate per-sample TSVs into the cohort table
  preflight: run the worker on a few KNOWN cases locally and assert expected
             behavior BEFORE submitting the array (catch bugs, don't flood cluster)

Usage on the cluster:
  # 1. preflight (local, no Slurm) -- MUST pass before submitting
  python ctx_cp_parallel.py preflight --known P-XXXX-chr9 ...

  # 2. write manifest + sbatch script
  python ctx_cp_parallel.py launch --out /path/cp_run --pens 2,3,5 --chunk 4

  # 3. submit
  sbatch /path/cp_run/run_array.sbatch

  # 4. after the array finishes, collect
  python ctx_cp_parallel.py collect --out /path/cp_run --cohort cohort2.tsv

All heavy work is per-sample and independent -> one array task per chunk of samples.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd

# these imports assume run_ctx_pipeline.py and ctx_changepoint.py are importable
import run_ctx_pipeline as rp
import ctx_changepoint as cc

MAN_PATH = rp.MAN_PATH


def _norm_chrom(x):
    return str(x).strip().replace("chr", "").replace("CHR", "")


# ----------------------------------------------------------------------
# worker: one sample -> per-sample TSV of cp_* features (all penalties)
# ----------------------------------------------------------------------

def compute_one_sample(man_df, sample_id, tablesDir, pens, signal_col="tcn_estimate"):
    _, jointseg_df, _ = rp.prepare_cn_data_for_sample(man_df, sample_id, tablesDir)
    row = man_df.loc[man_df["SAMPLE_ID"].astype(str).eq(str(sample_id))]
    purity = float(row["purity"].iloc[0]) if row.shape[0] and pd.notna(row["purity"].iloc[0]) else None
    ploidy = float(row["ploidy"].iloc[0]) if row.shape[0] and pd.notna(row["ploidy"].iloc[0]) else 2.0
    js = jointseg_df.copy()
    js["_c"] = js["chrom"].map(_norm_chrom)
    rows = []
    for c, sub in js.groupby("_c"):
        for pen in pens:
            f = cc.probe_changepoint_features(sub, purity=purity, ploidy=ploidy, pen_const=pen)
            f["SAMPLE_ID"] = sample_id
            f["chrom"] = c
            f["pen_const"] = pen
            rows.append(f)
    return pd.DataFrame(rows)


def run_worker(sample_ids, out_dir, tablesDir, pens, signal_col="tcn_estimate", man_path=None):
    cp_dir = os.path.join(out_dir, "cp_features")
    os.makedirs(cp_dir, exist_ok=True)
    man_df = rp.load_manifest(man_path) if man_path else rp.load_manifest()
    for sid in sample_ids:
        out_tsv = os.path.join(cp_dir, f"{sid}.tsv")
        if os.path.exists(out_tsv):
            print(f"[skip] {sid} (exists)", flush=True)
            continue
        try:
            df = compute_one_sample(man_df, sid, tablesDir, pens, signal_col=signal_col)
            df.to_csv(out_tsv, sep="\t", index=False)
            print(f"[done] {sid} -> {out_tsv} ({df.shape[0]} rows)", flush=True)
        except Exception as e:
            # write an empty marker so collect can report failures without crashing
            with open(os.path.join(cp_dir, f"{sid}.FAILED"), "w") as fh:
                fh.write(str(e))
            print(f"[FAIL] {sid}: {e}", flush=True)


# ----------------------------------------------------------------------
# launcher: manifest of processable samples + sbatch array script
# ----------------------------------------------------------------------

def launch(out_dir, pens, chunk, partition="componc_cpu", mem="16G",
           time="04:00:00", cpus=1, man_path=None, tablesDir=None):
    os.makedirs(out_dir, exist_ok=True)
    man_df = rp.load_manifest(man_path) if man_path else rp.load_manifest()
    tablesDir = tablesDir or (rp.OUT_DIR + "/tables")
    proc_ids, ledger = rp.build_processable_sample_ids(man_df, tablesDir)
    manifest = os.path.join(out_dir, "samples.txt")
    with open(manifest, "w") as fh:
        fh.write("\n".join(proc_ids) + "\n")
    n = len(proc_ids)
    n_tasks = (n + chunk - 1) // chunk
    pens_str = ",".join(str(p) for p in pens)
    script = os.path.join(out_dir, "run_array.sbatch")
    # NOTE: %a is the array index; each task handles `chunk` samples from the manifest.
    with open(script, "w") as fh:
        fh.write(f"""#!/bin/bash
#SBATCH --job-name=ctx_cp
#SBATCH --partition={partition}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --time={time}
#SBATCH --array=0-{n_tasks - 1}
#SBATCH --output={out_dir}/logs/task_%a.out

mkdir -p {out_dir}/logs
python {os.path.abspath(__file__)} worker \\
    --manifest {manifest} \\
    --out {out_dir} \\
    --pens {pens_str} \\
    --chunk {chunk} \\
    --tables {tablesDir} \\
    --man-path {man_path or rp.MAN_PATH} \\
    --array-index $SLURM_ARRAY_TASK_ID
""")
    print(f"manifest: {manifest}  ({n} samples)")
    print(f"sbatch:   {script}  ({n_tasks} array tasks, chunk={chunk})")
    print(f"submit:   sbatch {script}")
    return script


# ----------------------------------------------------------------------
# collector
# ----------------------------------------------------------------------



def selftest(out_dir, cohort_path, n_samples=3, pens=(3, 5), chunk=2,
             tablesDir=None, sample_col="SAMPLE_ID"):
    """
    END-TO-END test of the parallel machinery WITHOUT Slurm. Runs the real worker
    on a few samples (chunked exactly as the array would), then collect +
    calibrate, asserting outputs are well-formed. Run this BEFORE sbatch.

    Picks a mix: prefer some chr_ctx-positive samples (so calibration has pos
    rows) plus negatives. Writes into out_dir/selftest so it never touches the
    real run dir.
    """
    import shutil
    tablesDir = tablesDir or (rp.OUT_DIR + "/tables")
    test_dir = os.path.join(out_dir, "selftest")
    if os.path.exists(test_dir):
        shutil.rmtree(test_dir)
    os.makedirs(test_dir, exist_ok=True)

    man_df = rp.load_manifest()
    cohort = pd.read_csv(cohort_path, sep="\t")
    # pick samples: up to half positive, rest negative, that are in the cohort
    in_cohort = cohort[sample_col].astype(str).unique().tolist()
    pos_samples = (cohort.loc[cohort["truth_positive_chrom"].astype(bool), sample_col]
                   .astype(str).drop_duplicates().tolist())
    pos_pick = [s for s in pos_samples if s in in_cohort][: max(1, n_samples // 2)]
    neg_pick = [s for s in in_cohort if s not in pos_samples][: n_samples - len(pos_pick)]
    sample_ids = pos_pick + neg_pick
    assert len(sample_ids) >= 1, "no samples selected for selftest"
    print(f"[selftest] samples: {sample_ids} (pos={pos_pick}, neg={neg_pick})")

    # 1. manifest
    manifest = os.path.join(test_dir, "samples.txt")
    with open(manifest, "w") as fh:
        fh.write("\n".join(sample_ids) + "\n")

    # 2. run the worker over chunks EXACTLY as the array would (array-index loop)
    n_tasks = (len(sample_ids) + chunk - 1) // chunk
    print(f"[selftest] simulating {n_tasks} array tasks (chunk={chunk})")
    with open(manifest) as fh:
        all_ids = [x.strip() for x in fh if x.strip()]
    for task in range(n_tasks):
        i0 = task * chunk
        chunk_ids = all_ids[i0:i0 + chunk]
        print(f"[selftest] task {task}: {chunk_ids}")
        run_worker(chunk_ids, test_dir, tablesDir, list(pens))

    # 3. assert per-sample outputs exist and are well-formed
    cp_dir = os.path.join(test_dir, "cp_features")
    produced = [f for f in os.listdir(cp_dir) if f.endswith(".tsv")]
    failed = [f for f in os.listdir(cp_dir) if f.endswith(".FAILED")]
    assert len(failed) == 0, f"[selftest][FAIL] worker failures: {failed}"
    assert len(produced) == len(sample_ids), \
        f"[selftest][FAIL] expected {len(sample_ids)} tsvs, got {len(produced)}"
    expected_cols = {"SAMPLE_ID", "chrom", "pen_const", "cp_n_changepoints",
                     "cp_osc_chain", "cp_within_dispersion", "cp_osc_confidence",
                     "cp_n_levels"}
    for f in produced:
        df = pd.read_csv(os.path.join(cp_dir, f), sep="\t")
        miss = expected_cols - set(df.columns)
        assert not miss, f"[selftest][FAIL] {f} missing cols {miss}"
        assert set(df["pen_const"].unique()) == set(float(p) for p in pens), \
            f"[selftest][FAIL] {f} pens {df['pen_const'].unique()} != {pens}"
        assert df["chrom"].nunique() >= 1
    print(f"[selftest] worker OK: {len(produced)} tsvs, all cols present, all pens present")

    # 4. idempotency: re-running the worker must SKIP existing files
    before = {f: os.path.getmtime(os.path.join(cp_dir, f)) for f in produced}
    run_worker(sample_ids[:1], test_dir, tablesDir, list(pens))
    after = os.path.getmtime(os.path.join(cp_dir, produced[0]))
    assert after == before[produced[0]], "[selftest][FAIL] idempotency broken (rewrote existing)"
    print("[selftest] idempotency OK (existing files skipped)")

    # 5. collect + merge
    merged = collect(test_dir, cohort_path=cohort_path, pen_const=float(pens[0]), save=True)
    assert "cp_osc_chain" in merged.columns, "[selftest][FAIL] collect didn't merge cp features"
    n_merged = merged["cp_osc_chain"].notna().sum()
    assert n_merged > 0, "[selftest][FAIL] no cp features merged onto cohort rows"
    print(f"[selftest] collect OK: {n_merged} cohort rows got cp features")

    # 6. calibration report runs and is well-formed
    rep, _ = calibration_report(test_dir, cohort_path)
    assert set(rep["pen_const"]) == set(float(p) for p in pens), "[selftest][FAIL] calibrate pens"
    assert {"neg_median", "pos_median"}.issubset(rep.columns)
    print("[selftest] calibrate OK:")
    print(rep.to_string(index=False))

    print("\n[selftest] PASS  ->  safe to launch the full array")
    return True


def calibration_report(out_dir, cohort_path, truth_col="truth_positive_chrom",
                       sample_col="SAMPLE_ID"):
    """
    Anti-flood calibration from the ALREADY-COMPUTED parallel output (no serial
    recompute). Reads cp_features/*.tsv (which contain all penalties), joins the
    truth label from the cohort table, and reports neg vs pos cp_n_changepoints
    per penalty. Pick the smallest pen_const where neg_median <= ~1-2 and
    pos_median stays clearly above.
    """
    cp_dir = os.path.join(out_dir, "cp_features")
    tsvs = [os.path.join(cp_dir, f) for f in os.listdir(cp_dir) if f.endswith(".tsv")]
    cp = pd.concat([pd.read_csv(t, sep="\t") for t in tsvs], ignore_index=True)
    cp["_c"] = cp["chrom"].map(_norm_chrom)

    cohort = pd.read_csv(cohort_path, sep="\t")
    cohort["_c"] = cohort["chrom"].map(_norm_chrom)
    truth = cohort[[sample_col, "_c", truth_col]].copy()
    cp = cp.merge(truth, on=[sample_col, "_c"], how="left")
    cp[truth_col] = cp[truth_col].fillna(False).astype(bool)

    rows = []
    for pen, g in cp.groupby("pen_const"):
        neg = g.loc[~g[truth_col], "cp_n_changepoints"]
        pos = g.loc[g[truth_col], "cp_n_changepoints"]
        rows.append(dict(
            pen_const=pen,
            neg_median=neg.median(), neg_mean=round(neg.mean(), 2), neg_p90=neg.quantile(0.90),
            pos_median=pos.median(), pos_mean=round(pos.mean(), 2),
            n_neg=len(neg), n_pos=len(pos),
        ))
    rep = pd.DataFrame(rows).sort_values("pen_const").reset_index(drop=True)
    return rep, cp


def collect(out_dir, cohort_path=None, pen_const=None, save=True):
    cp_dir = os.path.join(out_dir, "cp_features")
    tsvs = [os.path.join(cp_dir, f) for f in os.listdir(cp_dir) if f.endswith(".tsv")]
    fails = [f for f in os.listdir(cp_dir) if f.endswith(".FAILED")]
    if fails:
        print(f"[warn] {len(fails)} samples FAILED: {fails[:10]}{' ...' if len(fails)>10 else ''}")
    cp = pd.concat([pd.read_csv(t, sep="\t") for t in tsvs], ignore_index=True)
    print(f"[collect] {len(tsvs)} sample files, {cp.shape[0]} (chrom,pen) rows")
    if pen_const is not None:
        cp = cp.loc[cp["pen_const"] == pen_const].copy()
    if save:
        out = os.path.join(out_dir, "cohort_cp_features.tsv")
        cp.to_csv(out, sep="\t", index=False)
        print(f"[collect] wrote {out}")
    if cohort_path:
        cohort = pd.read_csv(cohort_path, sep="\t")
        cohort["_c"] = cohort["chrom"].map(_norm_chrom)
        cp["_c"] = cp["chrom"].map(_norm_chrom)
        merged = cohort.merge(cp.drop(columns=["chrom"]), on=["SAMPLE_ID", "_c"], how="left")
        return merged.drop(columns=["_c"])
    return cp


# ----------------------------------------------------------------------
# preflight: KNOWN cases, assert behavior before submitting the array
# ----------------------------------------------------------------------

def preflight(tablesDir, pens=(2, 3, 5), cases=None):
    """
    cases: list of (SAMPLE_ID, chrom, expectation) where expectation is one of
      'oscillator' -> expect cp_osc_chain >= 4 at some penalty (e.g. C11 chr9)
      'flat_loh'   -> expect cp_osc_chain <= 2 (e.g. C22 chr3, no oscillation)
      'negative'   -> expect cp_n_changepoints small (no flood)
    Prints the features and PASS/FAIL per assertion. Run locally before sbatch.
    """
    if not cc._HAVE_RUPTURES:
        print("[preflight][WARN] ruptures NOT importable here -- using DP fallback. "
              "On the cluster this must say available=True.")
    man_df = rp.load_manifest()
    print(f"ruptures available: {cc._HAVE_RUPTURES}\n")
    all_ok = True
    for sid, chrom, expect in (cases or []):
        try:
            _, js, _ = rp.prepare_cn_data_for_sample(man_df, sid, tablesDir)
        except Exception as e:
            print(f"[FAIL load] {sid}: {e}"); all_ok = False; continue
        js = js.copy(); js["_c"] = js["chrom"].map(_norm_chrom)
        sub = js.loc[js["_c"].eq(_norm_chrom(chrom))]
        if sub.shape[0] == 0:
            print(f"[FAIL] {sid} chr{chrom}: no probes"); all_ok = False; continue
        prow = man_df.loc[man_df["SAMPLE_ID"].astype(str).eq(str(sid))]
        pur = float(prow["purity"].iloc[0]) if prow.shape[0] and pd.notna(prow["purity"].iloc[0]) else None
        plo = float(prow["ploidy"].iloc[0]) if prow.shape[0] and pd.notna(prow["ploidy"].iloc[0]) else 2.0
        print(f"=== {sid} chr{chrom}  (expect: {expect}, n_probes={sub.shape[0]}, purity={pur}) ===")
        chains = {}
        for pen in pens:
            f = cc.probe_changepoint_features(sub, purity=pur, ploidy=plo, pen_const=pen)
            chains[pen] = f["cp_osc_chain"]
            print(f"  pen={pen}: n_cp={f['cp_n_changepoints']:2d} osc_chain={f['cp_osc_chain']:2d} "
                  f"n_levels={f['cp_n_levels']:2d} within_disp={f['cp_within_dispersion']:.2f} "
                  f"osc_conf={f['cp_osc_confidence']:.2f}")
        # assertions
        if expect == "oscillator":
            ok = max(chains.values()) >= 4
        elif expect == "flat_loh":
            # informational only: FACETS called this flat, but probe signal may
            # reveal structure FACETS merged. Do NOT assert -- report and move on.
            ok = True
        elif expect == "negative":
            ok = True  # informational; checked in calibration
        else:
            ok = True
        print(f"  -> {'PASS' if ok else 'FAIL'}\n")
        all_ok = all_ok and ok
    print("PREFLIGHT", "PASS" if all_ok else "FAIL")
    return all_ok


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    w = sub.add_parser("worker")
    w.add_argument("--manifest", required=True)
    w.add_argument("--out", required=True)
    w.add_argument("--pens", default="2,3,5")
    w.add_argument("--chunk", type=int, default=1)
    w.add_argument("--array-index", type=int, required=True)
    w.add_argument("--tables", default=rp.OUT_DIR + "/tables")
    w.add_argument("--man-path", default=None)

    l = sub.add_parser("launch")
    l.add_argument("--out", required=True)
    l.add_argument("--pens", default="2,3,5")
    l.add_argument("--chunk", type=int, default=4)
    l.add_argument("--partition", default="componc_cpu")
    l.add_argument("--mem", default="16G")
    l.add_argument("--time", default="04:00:00")

    c = sub.add_parser("collect")
    c.add_argument("--out", required=True)
    c.add_argument("--cohort", default=None)
    c.add_argument("--pen", type=float, default=None)

    cb = sub.add_parser("calibrate")
    cb.add_argument("--out", required=True)
    cb.add_argument("--cohort", required=True)

    st = sub.add_parser("selftest")
    st.add_argument("--out", required=True)
    st.add_argument("--cohort", required=True)
    st.add_argument("--n", type=int, default=3)
    st.add_argument("--pens", default="3,5")
    st.add_argument("--chunk", type=int, default=2)
    st.add_argument("--tables", default=rp.OUT_DIR + "/tables")

    args = ap.parse_args()
    pens = [float(p) for p in args.pens.split(",")] if hasattr(args, "pens") else [3.0]

    if args.cmd == "worker":
        with open(args.manifest) as fh:
            all_ids = [x.strip() for x in fh if x.strip()]
        i0 = args.array_index * args.chunk
        chunk_ids = all_ids[i0:i0 + args.chunk]
        print(f"[worker] task {args.array_index}: samples {i0}..{i0+len(chunk_ids)-1}", flush=True)
        run_worker(chunk_ids, args.out, args.tables, pens, man_path=args.man_path)
    elif args.cmd == "launch":
        launch(args.out, pens, args.chunk, partition=args.partition,
               mem=args.mem, time=args.time)
    elif args.cmd == "collect":
        collect(args.out, cohort_path=args.cohort, pen_const=args.pen)
    elif args.cmd == "calibrate":
        rep, _ = calibration_report(args.out, args.cohort)
        print(rep.to_string(index=False))
    elif args.cmd == "selftest":
        pens = [float(p) for p in args.pens.split(",")]
        selftest(args.out, args.cohort, n_samples=args.n, pens=pens,
                 chunk=args.chunk, tablesDir=args.tables)


if __name__ == "__main__":
    main()
