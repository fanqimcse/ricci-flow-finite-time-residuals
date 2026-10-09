# -*- coding: utf-8 -*-
"""Experiment D analysis: transitivity, manipulated, degree sequence fixed.

The pooled audit nominated transitivity but could not separate it from size,
family, or degree sequence. Experiment D manipulates it directly: the same
configuration-model base graph is rewired by degree-preserving swaps towards
a target transitivity, so the degree sequence is invariant and only the
wiring pattern, and whatever it drags along, changes.

Three questions, in order of importance:

1. Causal: does the plateau rise when transitivity rises, within every base?
2. Attribution: the rewiring also moves assortativity. Experiment B already
   calibrated the assortativity effect under base fixed effects
   (beta = +0.341 in log10 per unit r). How much of the observed rise can
   that calibration account for? Whatever survives is the transitivity route.
3. Magnitude: does pushing transitivity to the real-network median (0.463)
   lift the plateau from the synthetic 0.002 towards the real 0.197?

Caveat carried through: high-transitivity runs do not saturate inside the
step budget, so the reported fold is an upper bound until the long run lands.
"""

import glob
import os

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "synthetic_results")

REAL_MEDIAN_PLATEAU = 0.19668
REAL_MEDIAN_TRANS = 0.4626
ASSORT_BETA_LOG10 = 0.341      # experiment B, base fixed effects


def main():
    frames = []
    for path in sorted(glob.glob(os.path.join(RESULTS, "clust_*.csv"))):
        try:
            frames.append(pd.read_csv(path))
        except Exception:
            pass
    if not frames:
        raise SystemExit("no clustering results yet")
    d = pd.concat(frames, ignore_index=True)
    # keep the longest run for each (base, target)
    d = d.sort_values("steps").drop_duplicates(
        subset=["base", "target_T"], keep="last")
    d = d.sort_values(["base", "target_T"]).reset_index(drop=True)
    d.to_csv(os.path.join(RESULTS, "clustering_merged.csv"), index=False)

    L = []
    A = L.append
    A("EXPERIMENT D: TRANSITIVITY, MANIPULATED (DEGREE SEQUENCE FIXED)")
    A("=" * 70)
    A("cells: {}".format(len(d)))
    A("status: {}".format(dict(d["status"].value_counts())))
    A("")
    A("degree sequence held within base? {}".format(
        d.groupby("base")["degree_seq_hash"].nunique().max() == 1))
    A("edges identical within base? {}".format(
        d.groupby("base")["edges"].nunique().max() == 1))
    A("")
    A(d[["base", "target_T", "achieved_T", "transitivity", "final_r",
         "nodes", "edges", "var_ratio", "steps", "status"]].to_string(
             index=False, float_format=lambda v: "{:+.5f}".format(v)))
    A("")

    A("--- 1. within-base effect ---")
    folds = []
    for base, sub in d.groupby("base"):
        sub = sub.sort_values("transitivity")
        if len(sub) < 2 or sub["transitivity"].nunique() < 2:
            A("  base {}: only one level, fold not estimable".format(base))
            continue
        lo, hi = sub.iloc[0], sub.iloc[-1]
        fold = float(hi["var_ratio"]) / max(float(lo["var_ratio"]), 1e-12)
        folds.append(fold)
        A("  base {}: T {:.3f}->{:.3f}  plateau {:.5f} -> {:.5f}  "
          "fold = {:.1f}x".format(
              base, float(lo["transitivity"]), float(hi["transitivity"]),
              float(lo["var_ratio"]), float(hi["var_ratio"]), fold))
        A("           assortativity {:.3f} -> {:.3f} (drift {:+.3f})".format(
            float(lo["final_r"]), float(hi["final_r"]),
            float(hi["final_r"]) - float(lo["final_r"])))
    if folds:
        A("")
        A("  median fold across bases: {:.1f}x  (n={})".format(
            float(np.median(folds)), len(folds)))
        A("  bases with fold > 1: {}/{}".format(
            sum(1 for f in folds if f > 1), len(folds)))
    A("")

    A("--- 2. attribution: how much is the assortativity drift? ---")
    A("  experiment B calibrated d log10(plateau)/d r = {:+.3f} "
      "(base fixed effects)".format(ASSORT_BETA_LOG10))
    for base, sub in d.groupby("base"):
        sub = sub.sort_values("transitivity")
        if len(sub) < 2:
            continue
        lo, hi = sub.iloc[0], sub.iloc[-1]
        dr = float(hi["final_r"]) - float(lo["final_r"])
        predicted_log10 = ASSORT_BETA_LOG10 * dr
        predicted_fold = 10 ** predicted_log10
        observed_log10 = np.log10(float(hi["var_ratio"]) /
                                  max(float(lo["var_ratio"]), 1e-12))
        share = predicted_log10 / observed_log10 if observed_log10 else np.nan
        A("  base {}: dr={:+.3f} -> assortativity predicts {:.2f}x ; "
          "observed {:.1f}x ; share explained = {:.1%}".format(
              base, dr, predicted_fold, 10 ** observed_log10, share))
    A("")

    A("--- 3. magnitude: does this close the synthetic-real gap? ---")
    lo_cells = d[d["transitivity"] <= 0.20]
    hi_cells = d[d["transitivity"] >= 0.35]
    if len(lo_cells) and len(hi_cells):
        lo_med = float(lo_cells["var_ratio"].median())
        hi_med = float(hi_cells["var_ratio"].median())
        A("  low transitivity  (T<=0.20): n={} median plateau={:.5f}".format(
            len(lo_cells), lo_med))
        A("  high transitivity (T>=0.35): n={} median plateau={:.5f} "
          "(median T={:.3f})".format(
              len(hi_cells), hi_med, float(hi_cells["transitivity"].median())))
        A("  rise = {:.1f}x".format(hi_med / lo_med if lo_med else np.nan))
        A("  real-network median plateau = {:.5f} (median T={:.3f})".format(
            REAL_MEDIAN_PLATEAU, REAL_MEDIAN_TRANS))
        A("  remaining gap after matching transitivity = {:.1f}x".format(
            REAL_MEDIAN_PLATEAU / hi_med if hi_med else np.nan))
        A("  share of the original ~90x gap that transitivity accounts for:")
        total = np.log10(REAL_MEDIAN_PLATEAU / lo_med)
        done = np.log10(hi_med / lo_med)
        A("    log10 total = {:.2f}, log10 explained = {:.2f} -> {:.0%}".format(
            total, done, done / total if total else np.nan))
    A("")

    if len(d) >= 4:
        A("--- 4. pooled rank correlation within the experiment ---")
        for cov in ("transitivity", "final_r", "edges"):
            x = d[cov].astype(float).values
            y = np.log10(d["var_ratio"].astype(float).values)
            m = np.isfinite(x) & np.isfinite(y)
            if m.sum() >= 4 and np.std(x[m]) > 0:
                rho, p = stats.spearmanr(x[m], y[m])
                A("  Spearman(log10 plateau, {:14s}) rho={:+.3f} p={:.3g} "
                  "n={}".format(cov, rho, p, int(m.sum())))
    A("")

    A("--- caveat ---")
    A("High-transitivity cells do not reach 'saturated' within the step")
    A("budget, so their plateau is still falling and the fold reported here")
    A("is an UPPER bound. The 400-step run is the arbiter; until it lands,")
    A("quote the fold as 'at a matched horizon of 170 steps'.")

    text = "\n".join(L)
    with open(os.path.join(RESULTS, "CLUSTERING_REPORT.txt"), "w",
              encoding="utf-8") as fh:
        fh.write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
