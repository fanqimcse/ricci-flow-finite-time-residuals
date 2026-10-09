# -*- coding: utf-8 -*-
"""One declared sample, three questions.

Prior analyses in this project used three different samples:
  * converged subset            -> assortativity rho = +0.821 (n=7 families)
  * all networks with a ratio   -> assortativity rho = +0.291 (n=13 families)
  * LOFO with a status filter   -> out-of-sample R2 = +0.04 (ns)
Sample switching is how a weak effect gets mistaken for a strong one. This
script fixes ONE sample up front and answers three questions on it:

  Q1  within-family   : among sibling networks, does higher assortativity
                        go with a higher residual plateau?
  Q2  between-family  : do families with higher mean assortativity retain
                        more? And is modularity better?
  Q3  out-of-sample   : leave-one-family-out R2 for a fixed model list.

Declared sample
---------------
  status in {saturated, constant_curvature}   (the flow actually finished)
  initial_variance > 1e-10                    (ratio is not 0/0 noise)
  nodes >= 20 and edges >= 20                 (covariates estimable)
  var_ratio_final finite
"""

import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

BASE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE, "ricci_flow_longrun_results", "analysis")
SRC = os.path.join(OUT_DIR, "07_merged_with_structure.csv")

CONVERGED = {"saturated", "constant_curvature"}
MIN_VARIANCE = 1e-10
MIN_NODES = 20
MIN_EDGES = 20


def rho_p(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 5 or np.std(x[m]) == 0 or np.std(y[m]) == 0:
        return np.nan, np.nan, int(m.sum())
    r, p = spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


def ols_fit(y, X):
    X = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return beta


def ols_predict(beta, X):
    X = np.column_stack([np.ones(len(X)), X])
    return X @ beta


def main():
    raw = pd.read_csv(SRC)
    lines = []
    lines.append("FINAL PREDICTIVE AUDIT  (one declared sample)")
    lines.append("=" * 62)

    frame = raw.copy()
    chain = []
    before = len(frame)
    frame = frame[frame["status"].isin(CONVERGED)]
    chain.append(("S1 flow converged", before, len(frame)))
    before = len(frame)
    frame = frame[frame["initial_variance"].astype(float) > MIN_VARIANCE]
    chain.append(("S2 non-degenerate variance", before, len(frame)))
    before = len(frame)
    frame = frame[(frame["nodes"].astype(float) >= MIN_NODES) &
                  (frame["edges"].astype(float) >= MIN_EDGES)]
    chain.append(("S3 size >= 20/20", before, len(frame)))
    before = len(frame)
    frame = frame[np.isfinite(frame["var_ratio_final"].astype(float))]
    chain.append(("S4 finite ratio", before, len(frame)))

    lines.append("declared sample: {} -> {} networks, {} families".format(
        len(raw), len(frame), frame["family"].nunique()))
    for label, b, a in chain:
        lines.append("  {:28s} {:3d} -> {:3d}".format(label, b, a))
    lines.append("")
    lines.append(frame.groupby("family").size().to_string())
    lines.append("")

    frame = frame.copy()
    frame["y"] = np.log10(frame["var_ratio_final"].astype(float))
    frame["assort"] = frame["assortativity"].astype(float)
    frame["modul"] = frame["modularity"].astype(float)
    frame["trans"] = frame["transitivity"].astype(float)
    frame["logn"] = np.log10(frame["nodes"].astype(float).clip(lower=1.0))

    # ---- Q1 within-family
    lines.append("--- Q1. WITHIN-FAMILY Spearman(covariate, log10 plateau) ---")
    wf_rows = []
    for fam, sub in frame.groupby("family"):
        sub = sub.dropna(subset=["assort", "modul", "y"])
        if len(sub) < 5:
            lines.append("  {:10s} n={:<3d} (n<5, not estimable)".format(
                fam, len(sub)))
            continue
        ra, pa, _ = rho_p(sub["assort"], sub["y"])
        rm, pm, _ = rho_p(sub["modul"], sub["y"])
        wf_rows.append({"family": fam, "n": len(sub),
                        "rho_assort": ra, "p_assort": pa,
                        "rho_modul": rm, "p_modul": pm})
        lines.append("  {:10s} n={:<3d} assort rho={:+.3f} (p={:.3f}) | "
                     "modul rho={:+.3f} (p={:.3f})".format(
                         fam, len(sub), ra, pa, rm, pm))
    wf = pd.DataFrame(wf_rows)
    if len(wf):
        lines.append("  median within-family: assort {:+.3f} ({}/{} positive) | "
                     "modul {:+.3f} ({}/{} positive)".format(
                         float(np.nanmedian(wf["rho_assort"])),
                         int((wf["rho_assort"] > 0).sum()), len(wf),
                         float(np.nanmedian(wf["rho_modul"])),
                         int((wf["rho_modul"] > 0).sum()), len(wf)))
    lines.append("")

    # ---- Q2 between-family
    fam_means = frame.groupby("family").agg(
        n=("y", "size"),
        mean_y=("y", "mean"),
        mean_assort=("assort", "mean"),
        mean_modul=("modul", "mean"),
        mean_trans=("trans", "mean"),
        mean_logn=("logn", "mean"),
    ).reset_index()
    lines.append("--- Q2. BETWEEN-FAMILY (family means, n={}) ---".format(
        len(fam_means)))
    lines.append(fam_means.to_string(index=False,
                                     float_format=lambda v: "{:+.3f}".format(v)))
    for cov in ("mean_assort", "mean_modul", "mean_trans", "mean_logn"):
        r, p, n = rho_p(fam_means[cov], fam_means["mean_y"])
        lines.append("  Spearman({:12s}, mean plateau) rho={:+.3f} p={:.3f}".format(
            cov, r, p))
    lines.append("")
    lines.append("  jackknife (drop one family), assortativity:")
    base_r, _, _ = rho_p(fam_means["mean_assort"], fam_means["mean_y"])
    jr = []
    for i in range(len(fam_means)):
        sub = fam_means.drop(fam_means.index[i])
        r, _, _ = rho_p(sub["mean_assort"], sub["mean_y"])
        jr.append(r)
        lines.append("    drop {:10s} -> {:+.3f}".format(
            str(fam_means.iloc[i]["family"]), r))
    lines.append("    full {:+.3f} | range {:+.3f}..{:+.3f}".format(
        base_r, float(np.nanmin(jr)), float(np.nanmax(jr))))
    lines.append("")

    # ---- Q3 leave-one-family-out
    lines.append("--- Q3. LEAVE-ONE-FAMILY-OUT ---")
    models = [
        ("M0 intercept only", []),
        ("M1 assortativity", ["assort"]),
        ("M2 modularity", ["modul"]),
        ("M3 assortativity + modularity", ["assort", "modul"]),
        ("M4 assort + modul + trans + logN",
         ["assort", "modul", "trans", "logn"]),
    ]
    usable = frame.dropna(subset=["y", "assort", "modul", "trans", "logn"])
    families = sorted(usable["family"].unique().tolist())
    comp_rows = []
    for label, preds in models:
        pt, pp = [], []
        for fam in families:
            train = usable[usable["family"] != fam]
            test = usable[usable["family"] == fam]
            if preds:
                train = train.dropna(subset=preds)
                test = test.dropna(subset=preds)
            if len(train) < 5 or len(test) == 0:
                continue
            beta = ols_fit(train["y"].values, train[preds].values)
            pred = ols_predict(beta, test[preds].values)
            pt.extend(test["y"].values.tolist())
            pp.extend(pred.tolist())
        pt = np.asarray(pt, float)
        pp = np.asarray(pp, float)
        ss_res = float(((pt - pp) ** 2).sum())
        ss_tot = float(((pt - pt.mean()) ** 2).sum())
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
        rho = spearmanr(pt, pp)[0] if len(pt) > 4 else np.nan
        comp_rows.append({"model": label, "n_pred": len(pt),
                          "oos_R2": r2, "oos_rho": rho})
        lines.append("  {:34s} n={:<3d} oos_R2={:+.3f} oos_rho={:+.3f}".format(
            label, len(pt), r2, rho))
    lines.append("")

    lines.append("--- VERDICT ---")
    best = max(comp_rows[1:], key=lambda r: r["oos_R2"])
    lines.append("  best non-null model: {} (oos_R2={:+.3f})".format(
        best["model"], best["oos_R2"]))
    if best["oos_R2"] <= 0:
        lines.append("  => NO model beats the held-out family mean.")
        lines.append("     Do not call any covariate a predictor of the")
        lines.append("     residual. Describe direction and significance,")
        lines.append("     and report the null out-of-sample result.")
    else:
        lines.append("  => out-of-sample skill is positive but must be")
        lines.append("     reported with the per-family spread and n.")

    with open(os.path.join(OUT_DIR, "FINAL_PREDICTIVE_AUDIT.txt"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    pd.DataFrame(comp_rows).to_csv(
        os.path.join(OUT_DIR, "18_final_lofo.csv"), index=False)
    wf.to_csv(os.path.join(OUT_DIR, "18_within_family.csv"), index=False)
    fam_means.to_csv(os.path.join(OUT_DIR, "18_family_means.csv"), index=False)
    print("ok")


if __name__ == "__main__":
    main()
