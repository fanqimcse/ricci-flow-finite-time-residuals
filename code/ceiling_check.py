"""Is the `initial_kappa -> plateau ratio` link real, or a ceiling artefact?

Why this matters
----------------
In the converged subset the initial mean curvature kappa0 is the strongest
predictor of the plateau ratio var_final / var_initial (rho = +0.512, and
+0.786 at family level). But kappa is bounded above by 1. If high-kappa
networks start with a *mechanically* small curvature variance (everything
pressed against the ceiling), then a large ratio could simply reflect a small
denominator rather than a failure to homogenise.

This script separates the three competing readings:

  R1  ratio artefact      kappa0 predicts var_initial, and the ratio inherits it
  R2  genuine residual    kappa0 predicts the *absolute* residual variance too
  R3  collinear proxy     kappa0 is only a proxy for transitivity / assortativity

Outputs: 11_ceiling_diagnostics.csv, 12_nested_models.csv, CEILING_REPORT.txt
"""

import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

BASE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE, "ricci_flow_longrun_results", "analysis")


def load(name):
    return pd.read_csv(os.path.join(OUT_DIR, name))


def rho(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 6:
        return np.nan, np.nan, int(mask.sum())
    if np.std(x[mask]) == 0 or np.std(y[mask]) == 0:
        return np.nan, np.nan, int(mask.sum())
    r, p = spearmanr(x[mask], y[mask])
    return float(r), float(p), int(mask.sum())


def ols(y, X, names):
    y = np.asarray(y, dtype=float)
    X = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    xtx_inv = np.linalg.pinv(X.T @ X)
    sigma2 = float(resid @ resid) / dof
    se = np.sqrt(np.maximum(np.diag(xtx_inv) * sigma2, 1e-300))
    from scipy.stats import t as tdist
    tvals = beta / se
    pvals = 2.0 * (1.0 - tdist.cdf(np.abs(tvals), dof))
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - float(resid @ resid) / ss_tot if ss_tot > 0 else np.nan
    return pd.DataFrame({
        "term": ["intercept"] + list(names),
        "coef": beta, "se": se, "t": tvals, "p": pvals,
    }), float(r2), float(1 - (1 - r2) * (len(y) - 1) / dof), int(len(y))


def main():
    merged = load("07_merged_with_structure.csv")
    converged = merged[merged["converged"] == True].copy()  # noqa: E712
    # same exclusion chain as finalize_structure.py (degenerate 0/0 ratios and
    # graphs too small for topological covariates)
    converged = converged[converged["initial_variance"] >= 1e-10]
    converged = converged[(converged["nodes"] >= 20) & (converged["edges"] >= 20)]
    converged = converged[np.isfinite(converged["var_ratio_final"])]

    # absolute residual variance: reconstruct from ratio * initial variance
    converged["final_variance_abs"] = (converged["var_ratio_final"]
                                       * converged["initial_variance"])
    converged["log_final_variance"] = np.log10(
        converged["final_variance_abs"].clip(lower=1e-12))

    lines = []
    lines.append("CEILING / RE-SCALING DIAGNOSTICS FOR THE kappa0 RESULT")
    lines.append("=" * 64)
    lines.append("")
    lines.append("converged datasets: {}".format(len(converged)))
    lines.append("")

    # ---------- 1. is kappa0 correlated with the denominator? ----------
    lines.append("--- R1: does kappa0 mechanically shrink var_initial? ---")
    r, p, n = rho(converged["initial_kappa"], converged["initial_variance"])
    lines.append("  spearman(kappa0, var_initial)      = {:+.3f}  p={:.3g}  n={}".format(r, p, n))
    r, p, n = rho(converged["initial_variance"], converged["var_ratio_final"])
    lines.append("  spearman(var_initial, ratio)       = {:+.3f}  p={:.3g}  n={}".format(r, p, n))
    hi = converged["initial_kappa"] > converged["initial_kappa"].median()
    lines.append("  median var_initial, high kappa0    = {:.4f}".format(
        converged.loc[hi, "initial_variance"].median()))
    lines.append("  median var_initial, low  kappa0    = {:.4f}".format(
        converged.loc[~hi, "initial_variance"].median()))
    lines.append("  max kappa0 observed                = {:.4f}  (theoretical ceiling = 1)".format(
        converged["initial_kappa"].max()))
    near = (converged["initial_kappa"] > 0.8).sum()
    lines.append("  datasets with kappa0 > 0.8         = {}  (ceiling risk)".format(int(near)))
    lines.append("")

    # ---------- 2. absolute residual ----------
    lines.append("--- R2: does kappa0 also predict the ABSOLUTE residual variance? ---")
    rows = []
    for target in ["var_ratio_final", "final_variance_abs", "log_final_variance"]:
        for cov_name in ["initial_kappa", "transitivity", "assortativity",
                         "degree_cv", "log10_nodes", "initial_variance"]:
            r, p, n = rho(converged[cov_name], converged[target])
            rows.append({"target": target, "covariate": cov_name,
                         "spearman": r, "p_value": p, "n": n})
    diag = pd.DataFrame(rows)
    diag.to_csv(os.path.join(OUT_DIR, "11_ceiling_diagnostics.csv"), index=False)
    for target in ["var_ratio_final", "final_variance_abs", "log_final_variance"]:
        lines.append("  target = {}".format(target))
        sub = diag[diag["target"] == target]
        for _, row in sub.iterrows():
            if row["spearman"] != row["spearman"]:
                continue
            lines.append("    {:<20s} rho={:+.3f}  p={:.3g}".format(
                row["covariate"], row["spearman"], row["p_value"]))
        lines.append("")

    # ---------- 3. collinearity ----------
    lines.append("--- R3: is kappa0 just a proxy for a topological covariate? ---")
    pairs = [("initial_kappa", "transitivity"), ("initial_kappa", "assortativity"),
             ("initial_kappa", "degree_cv"), ("initial_kappa", "log10_nodes"),
             ("initial_kappa", "modularity"), ("initial_kappa", "density"),
             ("transitivity", "assortativity"), ("transitivity", "degree_cv")]
    for a, b in pairs:
        if a in converged.columns and b in converged.columns:
            r, p, n = rho(converged[a], converged[b])
            if r == r:
                lines.append("  spearman({:<14s}, {:<14s}) = {:+.3f}  p={:.3g}".format(
                    a, b, r, p))
    lines.append("")

    # ---------- 4. nested models ----------
    lines.append("--- nested models on the converged subset (target = ratio) ---")
    specs = [
        ("M0 kappa0 only", ["initial_kappa"]),
        ("M1 transitivity only", ["transitivity"]),
        ("M2 kappa0 + transitivity", ["initial_kappa", "transitivity"]),
        ("M3 kappa0 + topology", ["initial_kappa", "transitivity",
                                  "assortativity", "degree_cv"]),
        ("M4 topology only", ["transitivity", "assortativity", "degree_cv"]),
        ("M5 kappa0 + size", ["initial_kappa", "log10_nodes"]),
    ]
    model_rows = []
    for label, terms in specs:
        sub = converged.dropna(subset=terms + ["var_ratio_final"])
        if len(sub) < len(terms) + 6:
            lines.append("  {:<28s} skipped (n={})".format(label, len(sub)))
            continue
        X = np.column_stack([
            (sub[t].astype(float).values - sub[t].astype(float).mean())
            / (sub[t].astype(float).std(ddof=0) + 1e-300) for t in terms])
        res, r2, adj, n = ols(sub["var_ratio_final"].values, X, terms)
        lines.append("  {:<28s} n={:<3d} R2={:.3f} adjR2={:.3f}".format(
            label, n, r2, adj))
        for _, row in res.iterrows():
            if row["term"] == "intercept":
                continue
            lines.append("      {:<22s} beta={:+.3f} p={:.3g}".format(
                row["term"], row["coef"], row["p"]))
            model_rows.append({"model": label, "n": n, "r2": r2, "adj_r2": adj,
                               "term": row["term"], "coef": row["coef"],
                               "p": row["p"]})
    pd.DataFrame(model_rows).to_csv(
        os.path.join(OUT_DIR, "12_nested_models.csv"), index=False)
    lines.append("")

    # ---------- verdict ----------
    lines.append("--- verdict ---")
    r_abs, p_abs, _ = rho(converged["initial_kappa"],
                          converged["final_variance_abs"])
    r_den, p_den, _ = rho(converged["initial_kappa"], converged["initial_variance"])
    lines.append("  kappa0 vs absolute residual : rho={:+.3f} (p={:.3g})".format(
        r_abs if r_abs == r_abs else float("nan"), p_abs))
    lines.append("  kappa0 vs initial variance  : rho={:+.3f} (p={:.3g})".format(
        r_den if r_den == r_den else float("nan"), p_den))
    if r_abs == r_abs and abs(r_abs) > 0.25 and p_abs < 0.10:
        lines.append("  => R2 supported: the link survives on the absolute scale.")
    else:
        lines.append("  => ratio-scale result NOT confirmed on the absolute scale;")
        lines.append("     report kappa0 as a ratio-scale predictor only, with the")
        lines.append("     ceiling caveat stated explicitly.")
    if r_den == r_den and abs(r_den) > 0.4:
        lines.append("  => WARNING: kappa0 is entangled with the denominator.")

    text = "\n".join(lines)
    with open(os.path.join(OUT_DIR, "CEILING_REPORT.txt"), "w", encoding="utf-8") as fh:
        fh.write(text)
    print(text)


if __name__ == "__main__":
    main()
