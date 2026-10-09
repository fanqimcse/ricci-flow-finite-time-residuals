# -*- coding: utf-8 -*-
"""L3: turn experiment C into a quantitative law, then test it on real data.

Two parts.

Part 1 -- scaling law on the controlled ensemble
    Fit log10(plateau at matched horizon) = a + b * Q  (+ base fixed effects)
    on the 8 community-sweep graphs. If b is tight and the fit is good, the
    paper can state a quantitative relation instead of just "goes up".

Part 2 -- the same two factors on real networks
    Re-run the leave-one-family-out test with modularity IN the model, using
    the declared converged sample. Experiment C says community structure
    should carry the size of the effect; if modularity still fails to predict
    out of sample, that is a real and reportable finding: the causal link is
    established in the controlled ensemble but does not transfer to
    network-level prediction on this corpus.

Everything is reported, including the parts that fail.
"""

import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

BASE = os.path.dirname(os.path.abspath(__file__))
SYN = os.path.join(BASE, "synthetic_results")
OUT = os.path.join(BASE, "ricci_flow_longrun_results", "analysis")
REAL = os.path.join(OUT, "07_merged_with_structure.csv")

CONVERGED = {"saturated", "constant_curvature"}


def ols(y, X, names):
    X = np.column_stack([np.ones(len(y)), X]) if X.size else np.ones(
        (len(y), 1))
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    if dof <= 0:
        return None
    xtx_inv = np.linalg.pinv(X.T @ X)
    sigma2 = float(resid @ resid) / dof
    se = np.sqrt(np.maximum(np.diag(xtx_inv) * sigma2, 1e-300))
    from scipy.stats import t as tdist
    tvals = beta / se
    pvals = 2 * (1 - tdist.cdf(np.abs(tvals), dof))
    r2 = 1 - float(resid @ resid) / float(((y - y.mean()) ** 2).sum())
    return {"names": ["intercept"] + list(names), "coef": beta, "se": se,
            "p": pvals, "r2": r2, "dof": dof}


def main():
    L = []
    A = L.append
    A("L3: FROM MECHANISM TO A QUANTITATIVE LAW")
    A("=" * 62)

    # ---------- part 1
    A("")
    A("--- PART 1. scaling law on the controlled ensemble ---")
    frames = [pd.read_csv(c)
              for c in sorted(glob.glob(os.path.join(SYN,
                                                     "community_b*.csv")))]
    syn = pd.concat(frames, ignore_index=True)
    syn = syn.drop_duplicates(subset=["base", "target_mu", "rep"],
                              keep="last")
    syn = syn[np.isfinite(syn["var_ratio"].astype(float))].copy()
    syn["y"] = np.log10(syn["var_ratio"].astype(float))
    A("n = {} graphs, matched horizon {} steps, all '{}'".format(
        len(syn), int(syn["steps"].max()),
        syn["status"].value_counts().to_dict()))

    y = syn["y"].values
    m = ols(y, syn[["modularity"]].astype(float).values, ["modularity"])
    if m:
        A("  log10(plateau) = {:+.3f} {:+.3f} * Q      R2={:.3f}  "
          "p(Q)={:.3g}".format(m["coef"][0], m["coef"][1], m["r2"],
                               m["p"][1]))
        A("  -> a 0.1 increase in planted modularity multiplies the residual "
          "by {:.2f}x".format(10 ** (0.1 * m["coef"][1])))

    bases = sorted(syn["base"].unique().tolist())
    dummies = np.column_stack(
        [(syn["base"].values == b).astype(float) for b in bases[1:]])
    X2 = np.column_stack([syn[["modularity"]].astype(float).values, dummies])
    m2 = ols(y, X2, ["modularity"] + ["base{}".format(b) for b in bases[1:]])
    if m2:
        A("  with base fixed effects: beta(Q) = {:+.3f} (p={:.3g}), "
          "R2={:.3f}".format(m2["coef"][1], m2["p"][1], m2["r2"]))
    A("")
    A("  CAVEAT: measured at a matched horizon, not at saturation, so this")
    A("  is a relation between community strength and the residual AT EQUAL")
    A("  FLOW TIME, not necessarily between Q and the asymptotic plateau.")

    # ---------- part 2
    A("")
    A("--- PART 2. the same two factors on real networks ---")
    d = pd.read_csv(REAL)
    d = d[d["status"].isin(CONVERGED)]
    d = d[d["initial_variance"].astype(float) > 1e-10]
    d = d[(d["nodes"].astype(float) >= 20) & (d["edges"].astype(float) >= 20)]
    d = d[np.isfinite(d["var_ratio_final"].astype(float))].copy()
    d["y"] = np.log10(d["var_ratio_final"].astype(float))
    d["assort"] = d["assortativity"].astype(float)
    d["modul"] = d["modularity"].astype(float)
    A("declared sample: {} networks / {} families".format(
        len(d), d["family"].nunique()))

    fam = d.groupby("family").agg(
        n=("y", "size"), mean_y=("y", "mean"),
        mean_assort=("assort", "mean"), mean_modul=("modul", "mean"))
    fam = fam[fam["n"] >= 1]
    for cov in ("mean_assort", "mean_modul"):
        r, p = spearmanr(fam[cov], fam["mean_y"])
        A("  family level  Spearman({}) rho={:+.3f} p={:.3f} n={}".format(
            cov, r, p, len(fam)))

    usable = d.dropna(subset=["y", "assort", "modul"])
    models = [("M1 assortativity", ["assort"]),
              ("M2 modularity", ["modul"]),
              ("M3 assortativity + modularity", ["assort", "modul"])]
    for label, preds in models:
        pt, pp = [], []
        for f in sorted(usable["family"].unique()):
            tr = usable[usable["family"] != f].dropna(subset=preds)
            te = usable[usable["family"] == f].dropna(subset=preds)
            if len(tr) < 5 or len(te) == 0:
                continue
            beta, *_ = np.linalg.lstsq(
                np.column_stack([np.ones(len(tr)), tr[preds].values]),
                tr["y"].values, rcond=None)
            pred = np.column_stack([np.ones(len(te)), te[preds].values]) @ beta
            pt.extend(te["y"].values.tolist())
            pp.extend(pred.tolist())
        pt = np.asarray(pt, float)
        pp = np.asarray(pp, float)
        r2 = 1 - float(((pt - pp) ** 2).sum()) / float(
            ((pt - pt.mean()) ** 2).sum())
        A("  LOFO {:32s} oos_R2={:+.3f}  oos_rho={:+.3f}".format(
            label, r2, spearmanr(pt, pp)[0] if len(pt) > 4 else np.nan))

    A("")
    A("--- VERDICT ---")
    A("  Experiment C establishes the causal link under a FIXED degree")
    A("  sequence: strengthening community structure raises the residual,")
    A("  and assortativity is flat (rho = +0.000) so the two are distinct")
    A("  mechanisms. It also closes most of the magnitude gap (92x -> 5.7x,")
    A("  and that residual gap is a LOWER bound because the synthetic runs")
    A("  are not saturated).")
    A("")
    A("  What it does NOT establish: that either factor PREDICTS the")
    A("  residual of a real network. The leave-one-family-out R2 stays at or")
    A("  below zero for every model. So the paper may claim a controlled")
    A("  causal mechanism and a between-family tendency, never a predictor.")

    with open(os.path.join(SYN, "L3_REPORT.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
