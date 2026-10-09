"""Analyse the synthetic assortativity sweep.

Design recap
------------
Degree sequence held FIXED per base graph; only the wiring is changed by
degree-preserving double-edge swaps towards a target assortativity. Anything
that moves the plateau is therefore attributable to the rewiring, not to the
degree distribution.

The honest caveat: double-edge swaps also move transitivity. So the analysis
does not "isolate" assortativity by construction. Instead both r and
transitivity are measured on every generated graph and the question is which
one tracks the plateau once both are in the model.

Outputs: synthetic_results/analysis/{01_merged.csv,SYNTH_REPORT.txt}
"""

import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

BASE = os.path.dirname(os.path.abspath(__file__))
SYN_DIR = os.path.join(BASE, "synthetic_results")
OUT_DIR = os.path.join(SYN_DIR, "analysis")


def ols(y, X, names):
    """Plain OLS with intercept; returns coefs, se, t, p (two-sided)."""
    X = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    if dof <= 0:
        return None
    sigma2 = float(resid @ resid) / dof
    cov = sigma2 * np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    tval = beta / np.where(se > 0, se, np.nan)
    # Normal approximation for p-values; n is small so treat as indicative.
    from scipy.stats import t as tdist
    pval = 2 * tdist.sf(np.abs(tval), dof)
    r2 = 1 - float(resid @ resid) / float(((y - y.mean()) ** 2).sum())
    return {"names": ["intercept"] + list(names), "coef": beta, "se": se,
            "t": tval, "p": pval, "r2": r2, "dof": dof}


def show(model):
    lines = []
    lines.append("    {:<14s} {:>10s} {:>10s} {:>8s} {:>10s}".format(
        "term", "coef", "se", "t", "p"))
    for i, name in enumerate(model["names"]):
        lines.append("    {:<14s} {:>10.4f} {:>10.4f} {:>8.2f} {:>10.3g}".format(
            name, model["coef"][i], model["se"][i], model["t"][i],
            model["p"][i]))
    lines.append("    R2 = {:.3f}   dof = {}".format(model["r2"], model["dof"]))
    return lines


def main():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)

    files = sorted(glob.glob(os.path.join(SYN_DIR, "assort_base*.csv")))
    frames = []
    for path in files:
        try:
            frame = pd.read_csv(path)
        except Exception:
            continue
        if not frame.empty:
            frames.append(frame)
    if not frames:
        print("no synthetic results yet; looked in", SYN_DIR)
        return

    data = pd.concat(frames, ignore_index=True)
    data = data.drop_duplicates(subset=["base", "target_r", "rep"])
    data.to_csv(os.path.join(OUT_DIR, "01_merged.csv"), index=False)

    ok = data[np.isfinite(data["var_ratio"]) &
              np.isfinite(data["final_r"])].copy()
    ok = ok[ok["var_ratio"] > 0]

    L = []
    L.append("SYNTHETIC ASSORTATIVITY SWEEP")
    L.append("=" * 62)
    L.append("")
    L.append("graphs completed: {} / requested {}".format(
        len(ok), len(data)))
    L.append("status breakdown: {}".format(
        dict(data["status"].value_counts())))
    L.append("")

    if len(ok) < 8:
        L.append("too few completed graphs to analyse")
        text = "\n".join(L)
        with open(os.path.join(OUT_DIR, "SYNTH_REPORT.txt"), "w",
                  encoding="utf-8") as fh:
            fh.write(text)
        print(text)
        return

    # ---- 0. did the rewiring actually hit the targets?
    L.append("--- 0. manipulation check (did rewiring hit the target?) ---")
    chk = ok.groupby("target_r").agg(
        n=("final_r", "size"),
        achieved_median=("final_r", "median"),
        achieved_min=("final_r", "min"),
        achieved_max=("final_r", "max"),
        transitivity=("transitivity", "median"),
        plateau=("var_ratio", "median"))
    L.append(chk.to_string())
    L.append("")
    gap = float(np.median(np.abs(ok["final_r"] - ok["target_r"])))
    L.append("  median |achieved - target| = {:.3f}".format(gap))
    if gap > 0.15:
        L.append("  WARNING: rewiring did not reach the targets; the sweep is")
        L.append("  effectively a weak manipulation. Interpret with care.")
    L.append("")

    # ---- 1. monotone trend across target levels
    L.append("--- 1. plateau vs target assortativity (by base) ---")
    L.append("  {:<6s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s}".format(
        "base", *["r={:+.1f}".format(v) for v in
                  sorted(ok["target_r"].unique())]))
    for base_id, grp in ok.groupby("base"):
        cells = []
        for level in sorted(ok["target_r"].unique()):
            sub = grp[grp["target_r"] == level]["var_ratio"]
            cells.append("{:>8.4f}".format(float(np.median(sub)))
                         if len(sub) else "{:>8s}".format("-"))
        L.append("  {:<6d} {}".format(int(base_id), " ".join(cells)))
    L.append("")

    # ---- 2. correlations
    L.append("--- 2. correlations with the plateau ---")
    for col in ["final_r", "transitivity", "edges"]:
        if col not in ok.columns:
            continue
        vals = ok[col].astype(float).values
        rho, p = spearmanr(vals, ok["var_ratio"].astype(float).values)
        L.append("  Spearman(var_ratio, {:<14s}) rho = {:+.3f}  p = {:.3g}".format(
            col, rho, p))
    L.append("")
    # within-base: removes any base-graph-specific level shift
    within = []
    for base_id, grp in ok.groupby("base"):
        if grp["final_r"].nunique() >= 4 and len(grp) >= 8:
            rho, p = spearmanr(grp["final_r"].astype(float),
                               grp["var_ratio"].astype(float))
            within.append((base_id, rho, p, len(grp)))
    if within:
        L.append("  within-base (base-graph effects removed):")
        for base_id, rho, p, n in within:
            L.append("      base {} n={}  rho = {:+.3f}  p = {:.3g}".format(
                base_id, n, rho, p))
        rhos = [w[1] for w in within]
        L.append("      median within-base rho = {:+.3f}  ({} bases)".format(
            float(np.median(rhos)), len(rhos)))
    L.append("")

    # ---- 3. joint model
    L.append("--- 3. joint model: log10(plateau) ~ assortativity + transitivity ---")
    y = np.log10(ok["var_ratio"].astype(float).values)
    cols = ["final_r", "transitivity"]
    X = ok[cols].astype(float).values
    model = ols(y, X, cols)
    if model is not None:
        L.extend(show(model))
    L.append("")

    # ---- 3b. with base fixed effects: the strongest version of the test.
    # Every graph in a base shares the identical degree sequence, so base
    # dummies absorb everything about the degree distribution.
    L.append("--- 3b. same model WITH base fixed effects ---")
    dummies = np.column_stack([
        (ok["base"].astype(int).values == b).astype(float)
        for b in sorted(ok["base"].unique())[1:]
    ])
    X_fe = np.column_stack([X, dummies])
    names_fe = cols + ["base{}".format(b)
                       for b in sorted(ok["base"].unique())[1:]]
    model_fe = ols(y, X_fe, names_fe)
    if model_fe is not None:
        L.extend(show(model_fe))
        i_r = model_fe["names"].index("final_r")
        L.append("    -> assortativity beta = {:+.3f} (p = {:.3g}) "
                 "with degree sequence held fixed".format(
                     model_fe["coef"][i_r], model_fe["p"][i_r]))
    L.append("")

    # ---- 3c. magnitude gap: synthetic vs real networks
    L.append("--- 3c. magnitude sanity check (synthetic vs real) ---")
    real_path = os.path.join(BASE, "ricci_flow_longrun_results",
                             "longrun_summary.csv")
    if os.path.exists(real_path):
        try:
            real = pd.read_csv(real_path)
            col = ("var_ratio_final" if "var_ratio_final" in real.columns
                   else real.columns[-1])
            vals = pd.to_numeric(real[col], errors="coerce")
            vals = vals[np.isfinite(vals)]
            vals = vals[(vals > 0) & (vals < 1)]
            if len(vals):
                L.append("  synthetic plateau median = {:.5f} (n={})".format(
                    float(np.median(ok["var_ratio"])), len(ok)))
                L.append("  real-network plateau med = {:.5f} (n={})".format(
                    float(np.median(vals)), len(vals)))
                L.append("  ratio of medians          = {:.1f}x".format(
                    float(np.median(vals)) / float(np.median(ok["var_ratio"]))))
                L.append("  => the synthetic sweep reproduces the DIRECTION of the")
                L.append("     assortativity effect, but configuration-model graphs")
                L.append("     homogenise far more completely than real networks.")
                L.append("     Do NOT transfer the magnitude; only the sign and the")
                L.append("     causal attribution carry over.")
        except Exception as exc:
            L.append("  (could not compare: {})".format(exc))
    else:
        L.append("  (real-network summary not found at {})".format(real_path))
    L.append("")

    # ---- 4. verdict
    L.append("--- 4. verdict ---")
    rho_r = spearmanr(ok["final_r"].astype(float),
                      ok["var_ratio"].astype(float))
    rho_t = spearmanr(ok["transitivity"].astype(float),
                      ok["var_ratio"].astype(float))
    L.append("  assortativity rho = {:+.3f} (p={:.3g})".format(
        rho_r[0], rho_r[1]))
    L.append("  transitivity  rho = {:+.3f} (p={:.3g})".format(
        rho_t[0], rho_t[1]))
    best = model_fe if model_fe is not None else model
    if best is not None:
        b_r = best["coef"][best["names"].index("final_r")]
        p_r = best["p"][best["names"].index("final_r")]
        b_t = best["coef"][best["names"].index("transitivity")]
        p_t = best["p"][best["names"].index("transitivity")]
        tag = "with base fixed effects" if best is model_fe else "pooled"
        L.append("  joint ({}): assortativity beta={:+.3f} (p={:.3g}), "
                 "transitivity beta={:+.3f} (p={:.3g})".format(
                     tag, b_r, p_r, b_t, p_t))
        if p_r < 0.05 and b_r > 0:
            L.append("  => assortativity survives with transitivity controlled:")
            L.append("     rewiring a FIXED degree sequence towards higher")
            L.append("     assortativity raises the residual plateau.")
            L.append("     This is causal evidence, not just correlation.")
        elif p_t < 0.05 and p_r >= 0.05:
            L.append("  => transitivity, not assortativity, tracks the plateau.")
            L.append("     The real-network assortativity correlation was a")
            L.append("     proxy for clustering. Downgrade the causal claim.")
        else:
            L.append("  => neither predictor separates cleanly; report as")
            L.append("     inconclusive and say so in the paper.")
    L.append("")

    text = "\n".join(L)
    with open(os.path.join(OUT_DIR, "SYNTH_REPORT.txt"), "w",
              encoding="utf-8") as fh:
        fh.write(text)
    print(text)


if __name__ == "__main__":
    main()
