"""Final structural analysis with a declared, auditable exclusion chain.

What went wrong before
----------------------
The first pass filtered only `initial_variance > 0`. Several datasets carry
initial curvature variances of 1e-32 / 1e-33 - floating point noise, not a
variance - so their plateau "ratio" is a 0/0 artefact pinned at exactly 1.0.
Those datasets are also tiny (N=4, N=5) with transitivity = 1.0 and a high
kappa0, which is precisely the profile that manufactures the kappa0 and
transitivity signals. They must be removed *and the removal must be reported*.

Exclusion chain (declared before looking at the outcome):
  E1  status not in {saturated, constant_curvature, horizon_reached,
      budget_exceeded}   -> no usable trajectory
  E2  initial_variance < 1e-10  -> degenerate, ratio is 0/0 noise
  E3  nodes < 20 or edges < 20  -> topological covariates not estimable
  E4  var_ratio_final not finite

Every step reports how many datasets it drops, and the headline correlations
are recomputed both with and without the chain so the sensitivity is visible.

Outputs: 13_exclusion_audit.csv, 14_final_correlations.csv,
         15_final_models.csv, FINAL_STRUCTURE_REPORT.txt
"""

import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

BASE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE, "ricci_flow_longrun_results", "analysis")

MIN_VARIANCE = 1e-10
MIN_NODES = 20
MIN_EDGES = 20


def rho(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 6 or np.std(x[m]) == 0 or np.std(y[m]) == 0:
        return np.nan, np.nan, int(m.sum())
    r, p = spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


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
    adj = 1.0 - (1.0 - r2) * (len(y) - 1) / dof
    return pd.DataFrame({
        "term": ["intercept"] + list(names), "coef": beta, "se": se,
        "t": tvals, "p": pvals}), float(r2), float(adj), int(len(y))


def rank_ols(x, y):
    """R2 in rank space - robust to the heavy right tail of the ratio."""
    m = np.isfinite(np.asarray(x, float)) & np.isfinite(np.asarray(y, float))
    if m.sum() < 8:
        return np.nan
    rx = rankdata(np.asarray(x, float)[m])
    ry = rankdata(np.asarray(y, float)[m])
    X = np.column_stack([np.ones(len(rx)), rx])
    b, *_ = np.linalg.lstsq(X, ry, rcond=None)
    res = ry - X @ b
    return float(1.0 - res @ res / ((ry - ry.mean()) ** 2).sum())


def bh_fdr(pvals):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(pvals, dtype=float)
    ok = np.isfinite(p)
    out = np.full(len(p), np.nan)
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order]
    m = len(ranked)
    adj = ranked * m / (np.arange(1, m + 1))
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    res = np.empty(m)
    res[order] = np.minimum(adj, 1.0)
    out[ok] = res
    return out


COVARIATES = ["initial_kappa", "transitivity", "assortativity", "degree_cv",
              "modularity", "density", "log10_nodes", "mean_degree",
              "algebraic_connectivity", "lcc_fraction", "n_components"]


def main():
    full = pd.read_csv(os.path.join(OUT_DIR, "07_merged_with_structure.csv"))

    lines = []
    lines.append("FINAL STRUCTURAL ANALYSIS OF THE RESIDUAL CURVATURE VARIANCE")
    lines.append("=" * 66)
    lines.append("")

    # ---------------- exclusion chain ----------------
    steps = []
    current = full.copy()

    def drop(mask, code, reason):
        nonlocal current
        n_before = len(current)
        removed = current[mask]
        current = current[~mask].copy()
        steps.append({"step": code, "reason": reason,
                      "n_before": n_before, "n_removed": int(mask.sum()),
                      "n_after": len(current),
                      "examples": ", ".join(removed["dataset"].head(4).tolist())})

    usable_status = ["saturated", "constant_curvature", "horizon_reached",
                     "budget_exceeded"]
    drop(~current["status"].isin(usable_status), "E1",
         "no usable trajectory (failed / diverged / skipped / insufficient)")
    drop(current["initial_variance"] < MIN_VARIANCE, "E2",
         "initial curvature variance < 1e-10 -> ratio is 0/0 noise")
    drop((current["nodes"] < MIN_NODES) | (current["edges"] < MIN_EDGES), "E3",
         "graph too small for topological covariates (N<20 or E<20)")
    drop(~np.isfinite(current["var_ratio_final"]), "E4",
         "plateau ratio not finite")

    audit = pd.DataFrame(steps)
    audit.to_csv(os.path.join(OUT_DIR, "13_exclusion_audit.csv"), index=False)
    lines.append("--- exclusion chain ---")
    for _, row in audit.iterrows():
        lines.append("  {} {:<52s} -{} -> {}".format(
            row["step"], row["reason"][:52], row["n_removed"], row["n_after"]))
        if row["n_removed"]:
            lines.append("       e.g. {}".format(row["examples"]))
    lines.append("")

    clean = current.copy()
    conv = clean[clean["converged"] == True].copy()  # noqa: E712
    lines.append("clean sample: {} (converged {})".format(len(clean), len(conv)))
    lines.append("families represented: {}".format(clean["family"].nunique()))
    lines.append("")

    # ---------------- correlations ----------------
    rows = []
    for label, frame in [("clean_all", clean), ("clean_converged", conv)]:
        for cov_name in COVARIATES:
            if cov_name not in frame.columns:
                continue
            r, p, n = rho(frame[cov_name], frame["var_ratio_final"])
            r_abs, p_abs, _ = rho(
                frame[cov_name], frame["var_ratio_final"] * frame["initial_variance"])
            rows.append({"stratum": label, "covariate": cov_name, "n": n,
                         "spearman_ratio": r, "p_ratio": p,
                         "rank_r2": rank_ols(frame[cov_name], frame["var_ratio_final"]),
                         "spearman_abs_residual": r_abs, "p_abs": p_abs})
    corr = pd.DataFrame(rows)
    # multiple-comparison control within each stratum
    for label in corr["stratum"].unique():
        m = corr["stratum"] == label
        corr.loc[m, "p_ratio_bh"] = bh_fdr(corr.loc[m, "p_ratio"].values)
    corr.to_csv(os.path.join(OUT_DIR, "14_final_correlations.csv"), index=False)

    lines.append("--- Spearman vs plateau ratio (higher = LESS homogenised) ---")
    for label in ["clean_all", "clean_converged"]:
        sub = corr[corr["stratum"] == label].dropna(subset=["spearman_ratio"])
        sub = sub.reindex(sub["spearman_ratio"].abs().sort_values(ascending=False).index)
        lines.append("")
        lines.append("  [{}]  n={}".format(label, int(sub["n"].max())))
        lines.append("  {:<24s} {:>7s} {:>9s} {:>9s} {:>8s} {:>8s}".format(
            "covariate", "rho", "p", "p_BH", "rankR2", "rho_abs"))
        for _, row in sub.iterrows():
            lines.append("  {:<24s} {:>+7.3f} {:>9.3g} {:>9.3g} {:>8.3f} {:>+8.3f}".format(
                row["covariate"], row["spearman_ratio"], row["p_ratio"],
                row["p_ratio_bh"],
                row["rank_r2"] if row["rank_r2"] == row["rank_r2"] else float("nan"),
                row["spearman_abs_residual"]
                if row["spearman_abs_residual"] == row["spearman_abs_residual"]
                else float("nan")))
    lines.append("")

    # ---------------- nested models ----------------
    lines.append("--- nested models (converged, clean; target = plateau ratio) ---")
    specs = [
        ("M0 kappa0", ["initial_kappa"]),
        ("M1 transitivity", ["transitivity"]),
        ("M2 kappa0+transitivity", ["initial_kappa", "transitivity"]),
        ("M3 topology", ["transitivity", "assortativity", "degree_cv"]),
        ("M4 kappa0+topology", ["initial_kappa", "transitivity",
                                "assortativity", "degree_cv"]),
        ("M5 kappa0+size", ["initial_kappa", "log10_nodes"]),
    ]
    model_rows = []
    for label, terms in specs:
        sub = conv.dropna(subset=terms + ["var_ratio_final"])
        if len(sub) < len(terms) + 6:
            lines.append("  {:<24s} skipped (n={})".format(label, len(sub)))
            continue
        X = np.column_stack([
            (sub[t].astype(float).values - sub[t].astype(float).mean())
            / (sub[t].astype(float).std(ddof=0) + 1e-300) for t in terms])
        res, r2, adj, n = ols(sub["var_ratio_final"].values, X, terms)
        lines.append("  {:<24s} n={:<3d} R2={:.3f} adjR2={:.3f}".format(
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
        os.path.join(OUT_DIR, "15_final_models.csv"), index=False)
    lines.append("")

    # ---------------- family level ----------------
    fam = conv.groupby("family").agg(
        n=("dataset", "count"),
        ratio=("var_ratio_final", "median"),
        initial_kappa=("initial_kappa", "median"),
        transitivity=("transitivity", "median"),
        assortativity=("assortativity", "median"),
        degree_cv=("degree_cv", "median"),
        log10_nodes=("log10_nodes", "median"),
    ).reset_index()
    fam = fam[fam["n"] >= 2]
    fam.to_csv(os.path.join(OUT_DIR, "16_family_level_final.csv"), index=False)
    lines.append("--- family level (one point per source family) ---")
    lines.append("  families with >=2 converged datasets: {}".format(len(fam)))
    for cov_name in ["initial_kappa", "transitivity", "assortativity",
                     "degree_cv", "log10_nodes"]:
        r, p, n = rho(fam[cov_name], fam["ratio"])
        if r == r:
            lines.append("  {:<20s} rho={:+.3f}  p={:.3g}  (n={})".format(
                cov_name, r, p, n))
    lines.append("")

    # ---------------- sensitivity vs the unfiltered version ----------------
    lines.append("--- sensitivity: what the degenerate samples were doing ---")
    old = full[(full["converged"] == True)]  # noqa: E712
    for cov_name in ["initial_kappa", "transitivity"]:
        r_old, p_old, n_old = rho(old[cov_name], old["var_ratio_final"])
        r_new, p_new, n_new = rho(conv[cov_name], conv["var_ratio_final"])
        lines.append("  {:<16s} unfiltered rho={:+.3f} (n={}) -> clean rho={:+.3f} (n={})".format(
            cov_name, r_old, n_old, r_new, n_new))
    lines.append("")

    # ---------------- verdict ----------------
    lines.append("--- verdict ---")
    passed = corr[(corr["stratum"] == "clean_converged")].dropna(
        subset=["spearman_ratio"])
    passed = passed[passed["p_ratio_bh"] < 0.10]
    sub = passed.reindex(
        passed["spearman_ratio"].abs().sort_values(ascending=False).index)
    if len(sub):
        lines.append("  predictors surviving BH(10%) in the converged subset:")
        for _, row in sub.iterrows():
            lines.append("    {:<24s} rho={:+.3f}  p_BH={:.3g}".format(
                row["covariate"], row["spearman_ratio"], row["p_ratio_bh"]))
    else:
        lines.append("  NO covariate survives BH(10%) in the converged subset.")
        lines.append("  => the 'residual is determined by structure' claim is NOT")
        lines.append("     supported at this sample size; report it as exploratory.")
    ranked = corr[(corr["stratum"] == "clean_converged")].dropna(
        subset=["spearman_ratio"])
    top = ranked.reindex(
        ranked["spearman_ratio"].abs().sort_values(ascending=False).index)
    if len(top):
        best = top.iloc[0]
        lines.append("  strongest (uncorrected): {} rho={:+.3f} p={:.3g}".format(
            best["covariate"], best["spearman_ratio"], best["p_ratio"]))

    text = "\n".join(lines)
    with open(os.path.join(OUT_DIR, "FINAL_STRUCTURE_REPORT.txt"), "w",
              encoding="utf-8") as fh:
        fh.write(text)
    print(text)


if __name__ == "__main__":
    main()
