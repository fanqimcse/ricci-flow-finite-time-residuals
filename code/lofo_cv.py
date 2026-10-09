# -*- coding: utf-8 -*-
"""Leave-one-family-out prediction test for the residual plateau.

Why this matters
----------------
A correlation says "these two move together". A referee who wants to raise
the paper's tier asks a harder question:

    "Can you PREDICT the residual of a network family you did not fit on?"

That turns a covariate into a quantitative law. Two honesty constraints:

  * the held-out unit is the SOURCE FAMILY, never the single network
    (ENZYMES_g1..g113 are draws from one generator; scc_rt_* one corpus),
  * every fit uses the remaining families only,
  * fit quality is out-of-sample R^2 and Spearman on the held-out family.

The model list is fixed in advance - nothing is selected afterwards on p.

Models
------
  M0  intercept only                                   (null baseline)
  M1  assortativity
  M2  modularity
  M3  assortativity + modularity
  M4  assortativity + modularity + transitivity + log10(nodes)

M2 vs M1 is the decisive comparison for the 92x magnitude gap: if modularity
predicts better than assortativity on real networks, community structure --
not degree mixing -- sets the size of the residual.
"""

import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

BASE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE, "ricci_flow_longrun_results", "analysis")
SRC = os.path.join(OUT_DIR, "07_merged_with_structure.csv")

MIN_VARIANCE = 1e-10
MIN_NODES = 20
USABLE_STATUS = {"saturated", "constant_curvature", "horizon_reached"}


def ols_fit(y, X):
    X = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return beta


def ols_predict(beta, X):
    X = np.column_stack([np.ones(len(X)), X])
    return X @ beta


def main():
    frame = pd.read_csv(SRC)
    lines = []
    lines.append("LEAVE-ONE-FAMILY-OUT PREDICTION TEST")
    lines.append("=" * 62)
    lines.append("source: 07_merged_with_structure.csv")
    lines.append("")

    n0 = len(frame)
    steps = []
    if "initial_variance" in frame.columns:
        before = len(frame)
        frame = frame[frame["initial_variance"].astype(float) > MIN_VARIANCE]
        steps.append(("E2 degenerate initial variance", before, len(frame)))
    if "nodes" in frame.columns:
        before = len(frame)
        frame = frame[frame["nodes"].astype(float) >= MIN_NODES]
        steps.append(("E3 nodes >= 20", before, len(frame)))
    if "status" in frame.columns:
        before = len(frame)
        frame = frame[frame["status"].isin(USABLE_STATUS)]
        steps.append(("E1 usable trajectory", before, len(frame)))
    frame = frame[np.isfinite(frame["var_ratio_final"].astype(float))]
    steps.append(("E4 finite ratio", None, len(frame)))

    lines.append("exclusion chain: {} -> {}".format(n0, len(frame)))
    for label, before, after in steps:
        if before is None:
            lines.append("  {:32s} -> {}".format(label, after))
        else:
            lines.append("  {:32s} {} -> {} (-{})".format(
                label, before, after, before - after))
    lines.append("")
    lines.append("families retained:")
    lines.append(frame.groupby("family").size().to_string())
    lines.append("")

    frame = frame.copy()
    frame["y"] = np.log10(frame["var_ratio_final"].astype(float).values)
    frame["assort"] = frame["assortativity"].astype(float)
    frame["modul"] = frame["modularity"].astype(float)
    frame["trans"] = frame["transitivity"].astype(float)
    frame["logn"] = np.log10(frame["nodes"].astype(float).clip(lower=1.0))
    frame["dcv"] = frame["degree_cv"].astype(float)

    models = [
        ("M0 intercept only", []),
        ("M1 assortativity", ["assort"]),
        ("M2 modularity", ["modul"]),
        ("M3 assortativity + modularity", ["assort", "modul"]),
        ("M4 assort + modul + trans + logN",
         ["assort", "modul", "trans", "logn"]),
    ]

    core = ["y", "assort", "modul", "trans", "logn"]
    usable = frame.dropna(subset=core)
    lines.append("complete-case sample: {} networks / {} families".format(
        len(usable), usable["family"].nunique()))
    lines.append("")

    families = sorted(usable["family"].unique().tolist())
    results = {}

    for label, predictors in models:
        pooled_true, pooled_pred = [], []
        per_family = []
        for fam in families:
            train = usable[usable["family"] != fam]
            test = usable[usable["family"] == fam]
            if predictors:
                train = train.dropna(subset=predictors)
                test = test.dropna(subset=predictors)
            if len(train) < 8 or len(test) == 0:
                continue
            beta = ols_fit(train["y"].values, train[predictors].values)
            pred = ols_predict(beta, test[predictors].values)
            truth = test["y"].values
            ss_res = float(((truth - pred) ** 2).sum())
            ss_tot = float(((truth - truth.mean()) ** 2).sum())
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
            if len(truth) >= 4 and np.std(truth) > 0 and np.std(pred) > 0:
                rho, p = spearmanr(truth, pred)
            else:
                rho, p = np.nan, np.nan
            per_family.append({"family": fam, "n_test": len(truth),
                               "oos_r2": r2, "oos_rho": rho,
                               "mae_log10": float(np.mean(np.abs(truth - pred)))})
            pooled_true.extend(truth.tolist())
            pooled_pred.extend(pred.tolist())

        pt = np.asarray(pooled_true, dtype=float)
        pp = np.asarray(pooled_pred, dtype=float)
        ss_res = float(((pt - pp) ** 2).sum())
        ss_tot = float(((pt - pt.mean()) ** 2).sum())
        pooled_r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
        if len(pt) > 4 and np.std(pt) > 0 and np.std(pp) > 0:
            pooled_rho, pooled_p = spearmanr(pt, pp)
        else:
            pooled_rho, pooled_p = np.nan, np.nan
        results[label] = {
            "per_family": pd.DataFrame(per_family),
            "pooled_r2": pooled_r2,
            "pooled_rho": pooled_rho,
            "pooled_p": pooled_p,
            "n_pred": len(pt),
        }

        lines.append("--- {} ---".format(label))
        lines.append(pd.DataFrame(per_family).to_string(
            index=False, float_format=lambda v: "{:+.3f}".format(v)))
        lines.append("  POOLED out-of-sample: n={} R2={:+.3f} rho={:+.3f} "
                     "p={:.3g}".format(len(pt), pooled_r2, pooled_rho, pooled_p))
        lines.append("")

    lines.append("--- model comparison (pooled out-of-sample) ---")
    comp = pd.DataFrame([
        {"model": k, "n_pred": v["n_pred"], "oos_R2": v["pooled_r2"],
         "oos_rho": v["pooled_rho"], "p": v["pooled_p"]}
        for k, v in results.items()])
    lines.append(comp.to_string(index=False,
                                float_format=lambda v: "{:+.3f}".format(v)))
    lines.append("")
    lines.append("Reading rule: a model earns the word 'predictive' only if")
    lines.append("pooled out-of-sample R2 > 0 AND it beats M0. R2 < 0 means the")
    lines.append("held-out families are predicted worse than their own mean.")

    os.makedirs(OUT_DIR, exist_ok=True)
    with open(os.path.join(OUT_DIR, "LOFO_REPORT.txt"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    comp.to_csv(os.path.join(OUT_DIR, "16_lofo_comparison.csv"), index=False)
    for label, res in results.items():
        safe = (label.replace(" ", "_").replace("+", "plus")
                     .replace("(", "").replace(")", ""))
        res["per_family"].to_csv(
            os.path.join(OUT_DIR, "16_lofo_{}.csv".format(safe)), index=False)
    print("written LOFO_REPORT.txt")


if __name__ == "__main__":
    main()
