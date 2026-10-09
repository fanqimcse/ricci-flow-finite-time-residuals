# -*- coding: utf-8 -*-
"""L3 predictive audit, second pass: transitivity on the same declared sample.

The pooled real+synthetic audit (n=97) put transitivity ahead of every other
candidate, and it was the only predictor that was strong BOTH pooled
(rho=+0.821) and inside the real networks alone (rho=+0.621, p=9e-04).
Assortativity was the previous best and模块化 was a failure out of sample.
So transitivity has to go through the same gauntlet, on the same declared
sample, before it can be promoted.

Gauntlet (identical to FINAL_PREDICTIVE_AUDIT):
  Q1  within-family  -- the only design free of family confounding
  Q2  family level   -- n=7, with jackknife
  Q3  leave-one-family-out -- the honest out-of-sample test
Plus: the two-factor model (assortativity + transitivity) and its LOFO skill.

Baseline for Q3 is the intercept-only model (predict the held-out family with
the training mean). A model that cannot beat that has no out-of-sample skill,
whatever its in-sample fit looks like.
"""

import os

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.join(HERE, "ricci_flow_longrun_results", "analysis")
OUT = os.path.join(ANALYSIS, "PREDICTIVE_AUDIT2.txt")


def spear(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 5 or np.std(x[m]) == 0 or np.std(y[m]) == 0:
        return np.nan, np.nan, int(m.sum())
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


def ols_fit(y, X):
    X1 = np.column_stack([np.ones(len(y)), X]) if X.shape[1] else \
        np.ones((len(y), 1))
    beta, *_ = np.linalg.lstsq(X1, y, rcond=None)
    return beta, X1.shape[1]


def ols_predict(beta, X):
    X1 = np.column_stack([np.ones(len(X)), X]) if X.shape[1] else \
        np.ones((len(X), 1))
    return X1 @ beta


def main():
    meta = pd.read_csv(os.path.join(ANALYSIS, "07_merged_with_structure.csv"))
    keep = meta[(meta["converged"] == True)                        # noqa: E712
                & (meta["initial_variance"] > 1e-20)
                & (meta["nodes"] >= 20) & (meta["edges"] >= 20)
                & np.isfinite(meta["var_ratio_final"])
                & (meta["var_ratio_final"] > 0)].copy()
    keep["y"] = np.log10(keep["var_ratio_final"].astype(float))
    keep = keep.sort_values("dataset").reset_index(drop=True)

    L = []
    A = L.append
    A("PREDICTIVE AUDIT 2 - TRANSITIVITY ON THE DECLARED SAMPLE")
    A("=" * 66)
    A("sample: {} networks, {} families".format(
        len(keep), keep["family"].nunique()))
    A("")
    A(keep["family"].value_counts().to_string())
    A("")

    covs = ["transitivity", "assortativity", "modularity", "degree_cv",
            "log10_nodes", "density", "avg_path_len"]

    A("--- Q1. WITHIN-FAMILY Spearman(covariate, log10 plateau) ---")
    for cov in covs:
        parts = []
        for fam, sub in keep.groupby("family"):
            if len(sub) < 5:
                parts.append("{} n={} (n<5)".format(fam, len(sub)))
                continue
            r, p, n = spear(sub[cov], sub["y"])
            parts.append("{} n={} rho={:+.3f} (p={:.3g})".format(
                fam, n, r, p))
        A("  {:14s} {}".format(cov, " | ".join(parts)))
    A("")

    A("--- Q2. FAMILY LEVEL (n={}) ---".format(keep["family"].nunique()))
    fam = keep.groupby("family").agg(
        n=("y", "size"), y=("y", "mean"),
        **{c: (c, "mean") for c in covs})
    A(fam.to_string(float_format=lambda v: "{:+.4f}".format(v)))
    A("")
    for cov in covs:
        r, p, n = spear(fam[cov], fam["y"])
        A("  Spearman(family {:14s}, family plateau) rho={:+.3f} p={:.3g}"
          .format(cov, r, p))
    A("")
    for cov in ("transitivity", "assortativity"):
        vals = []
        for drop in fam.index:
            sub = fam.drop(index=drop)
            r, _, _ = spear(sub[cov], sub["y"])
            vals.append(r)
        A("  jackknife {}: full={:+.3f} range {:+.3f}..{:+.3f}".format(
            cov, spear(fam[cov], fam["y"])[0], min(vals), max(vals)))
    A("")

    A("--- Q3. LEAVE-ONE-FAMILY-OUT ---")
    models = {
        "M0 intercept only": [],
        "M1 transitivity": ["transitivity"],
        "M2 assortativity": ["assortativity"],
        "M3 transitivity + assortativity": ["transitivity", "assortativity"],
        "M4 T + assort + modul": ["transitivity", "assortativity",
                                  "modularity"],
        "M5 T + assort + modul + logN + degCV": [
            "transitivity", "assortativity", "modularity", "log10_nodes",
            "degree_cv"],
    }
    for label, cols in models.items():
        preds = np.empty(len(keep))
        preds[:] = np.nan
        for fam in keep["family"].unique():
            test = keep["family"] == fam
            train = keep[~test]
            if len(train) < 4:
                continue
            beta, _ = ols_fit(train["y"].values,
                              train[cols].astype(float).values
                              if cols else np.zeros((len(train), 0)))
            Xt = keep.loc[test, cols].astype(float).values if cols else \
                np.zeros((int(test.sum()), 0))
            preds[test.values] = ols_predict(beta, Xt)
        m = np.isfinite(preds)
        y = keep["y"].values[m]
        p_ = preds[m]
        ss_res = float(((y - p_) ** 2).sum())
        ss_tot = float(((y - y.mean()) ** 2).sum())
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
        rho, pv, _ = spear(p_, y)
        A("  {:42s} n={} oos_R2={:+.3f} oos_rho={:+.3f} (p={:.3g})".format(
            label, int(m.sum()), r2, rho, pv))
    A("")

    A("--- Q4. in-sample fit (for contrast only, NOT evidence) ---")
    for label, cols in models.items():
        if not cols:
            continue
        X = keep[cols].astype(float).values
        beta, _ = ols_fit(keep["y"].values, X)
        pred = ols_predict(beta, X)
        ss_res = float(((keep["y"].values - pred) ** 2).sum())
        ss_tot = float(((keep["y"].values - keep["y"].values.mean()) ** 2)
                       .sum())
        A("  {:42s} in-sample R2={:.3f}".format(label, 1 - ss_res / ss_tot))
    A("")

    A("--- reading ---")
    A("A predictor is promoted only if it survives Q1 or Q2 AND beats the")
    A("intercept baseline in Q3. In-sample R2 is reported only to show how")
    A("misleading it is here.")

    text = "\n".join(L)
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
