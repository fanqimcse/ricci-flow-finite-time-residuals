# -*- coding: utf-8 -*-
"""Robustness of the one within-family result that survives.

Predictive audit 2 found the first structural predictor that works INSIDE a
family: within ENZYMES (n=11), Spearman(transitivity, log10 plateau) = +0.782
(p=0.0045), while assortativity (-0.345) and modularity (-0.364) both fail.
Everything else in the predictive programme is either an ecological
correlation or a failure out of sample, so this single result carries a lot of
weight and has to be stress-tested before it goes in the paper.

Tests:
  1. leave-one-out jackknife on the 11 ENZYMES networks
  2. partial Spearman controlling for size (log10 nodes), density and
     assortativity -- transitivity, size and density move together
  3. multiple-comparison correction: audit 2 screened 7 covariates
  4. is it an artefact of one influential network?
"""

import os

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.join(HERE, "ricci_flow_longrun_results", "analysis")
OUT = os.path.join(ANALYSIS, "WITHIN_FAMILY_ROBUSTNESS.txt")


def spear(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 5 or np.std(x[m]) == 0 or np.std(y[m]) == 0:
        return np.nan, np.nan, int(m.sum())
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


def partial_spearman(x, y, controls):
    """Rank-transform then residualise both sides on the controls."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    Z = np.column_stack([np.asarray(c, dtype=float) for c in controls])
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(Z).all(axis=1)
    if m.sum() < 6:
        return np.nan, np.nan, int(m.sum())
    xr = stats.rankdata(x[m])
    yr = stats.rankdata(y[m])
    Zr = np.column_stack([stats.rankdata(Z[m][:, j])
                          for j in range(Z.shape[1])])
    Z1 = np.column_stack([np.ones(len(Zr)), Zr])
    bx, *_ = np.linalg.lstsq(Z1, xr, rcond=None)
    by, *_ = np.linalg.lstsq(Z1, yr, rcond=None)
    rx = xr - Z1 @ bx
    ry = yr - Z1 @ by
    r, p = stats.pearsonr(rx, ry)
    return float(r), float(p), int(m.sum())


def main():
    meta = pd.read_csv(os.path.join(ANALYSIS, "07_merged_with_structure.csv"))
    keep = meta[(meta["converged"] == True)                        # noqa: E712
                & (meta["initial_variance"] > 1e-20)
                & (meta["nodes"] >= 20) & (meta["edges"] >= 20)
                & np.isfinite(meta["var_ratio_final"])
                & (meta["var_ratio_final"] > 0)].copy()
    keep["y"] = np.log10(keep["var_ratio_final"].astype(float))
    enz = keep[keep["family"] == "ENZYMES"].reset_index(drop=True)

    L = []
    A = L.append
    A("WITHIN-FAMILY ROBUSTNESS: TRANSITIVITY INSIDE ENZYMES")
    A("=" * 66)
    A("networks: {}".format(len(enz)))
    A("")
    A(enz[["dataset", "nodes", "edges", "transitivity", "assortativity",
           "modularity", "density", "var_ratio_final"]].to_string(
               index=False, float_format=lambda v: "{:+.4f}".format(v)))
    A("")

    A("--- 1. leave-one-out jackknife, Spearman(transitivity, log10 plateau) ---")
    full_rho, full_p, _ = spear(enz["transitivity"], enz["y"])
    A("  full sample rho={:+.3f} p={:.3g}".format(full_rho, full_p))
    vals = []
    for i in range(len(enz)):
        sub = enz.drop(index=i)
        r, p, _ = spear(sub["transitivity"], sub["y"])
        vals.append(r)
        A("  drop {:24s} rho={:+.3f} (p={:.3g})".format(
            str(enz.loc[i, "dataset"]), r, p))
    A("  range {:+.3f} .. {:+.3f}; all positive: {}".format(
        min(vals), max(vals), all(v > 0 for v in vals)))
    A("")

    A("--- 2. partial correlation, controlling for confounds ---")
    for label, ctrl in (
            ("size (log10 nodes)", [enz["log10_nodes"]]),
            ("density", [enz["density"]]),
            ("assortativity", [enz["assortativity"]]),
            ("size + density", [enz["log10_nodes"], enz["density"]]),
            ("size + density + assortativity",
             [enz["log10_nodes"], enz["density"], enz["assortativity"]]),
    ):
        r, p, n = partial_spearman(enz["transitivity"], enz["y"], ctrl)
        A("  control {:34s} partial rho={:+.3f} p={:.3g} n={}".format(
            label, r, p, n))
    A("")

    A("--- 3. multiple comparisons ---")
    A("  audit 2 screened 7 covariates within this family.")
    A("  Bonferroni-corrected p for transitivity = {:.4f}".format(
        min(1.0, full_p * 7)))
    A("  BH threshold at q=0.05, rank 1 of 7 = {:.4f}".format(0.05 / 7))
    A("")

    A("--- 4. is one network driving it? ---")
    lo, hi = min(vals), max(vals)
    A("  jackknife range width = {:.3f} (a single influential point would "
      "produce a wide swing)".format(hi - lo))
    A("  worst-case rho after dropping any single network = {:+.3f}".format(lo))
    A("")

    A("--- cross-check: the same test on the other covariates ---")
    for cov in ("assortativity", "modularity", "degree_cv", "log10_nodes",
                "density", "avg_path_len"):
        r, p, n = spear(enz[cov], enz["y"])
        vals2 = []
        for i in range(len(enz)):
            sub = enz.drop(index=i)
            rr, _, _ = spear(sub[cov], sub["y"])
            vals2.append(rr)
        A("  {:14s} rho={:+.3f} p={:.3g} | jackknife {:+.3f}..{:+.3f} "
          "sign-stable={}".format(cov, r, p, min(vals2), max(vals2),
                                  all(v > 0 for v in vals2)
                                  or all(v < 0 for v in vals2)))
    A("")

    A("--- reading ---")
    A("Transitivity is promoted if the jackknife is sign-stable, the partial")
    A("correlation survives size/density control, and the Bonferroni p stays")
    A("below 0.05. If the partial correlation collapses under size control,")
    A("report transitivity and size as not separable inside this family.")

    text = "\n".join(L)
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
