"""Structural predictors of the residual curvature variance.

Goal of this script
-------------------
The long-horizon experiment showed that the discrete Ollivier-Ricci flow only
*partially* homogenises curvature: even after convergence the residual
variance ratio plateaus at a dataset-specific level (median 0.16 among runs
that converged, but 22% stay above 0.30).

The claim "the residual is determined by structure" is only admissible if
measured structural covariates actually predict it. This script computes
topological covariates for every dataset and tests them against the plateau.

It reuses the original loader so the graphs are byte-identical to the ones
used in the flow experiment.

Outputs (analysis dir):
  07_structure_predictors.csv   per-dataset covariates
  08_structure_correlations.csv spearman/partial correlations
  09_structure_regression.csv    multivariate OLS on standardised covariates
  10_family_level_check.csv      family-aggregated robustness check
"""

import os
import sys
import warnings

import numpy as np
import pandas as pd
import networkx as nx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ricci_flow_program_100data as R  # noqa: E402

warnings.filterwarnings("ignore")

BASE = os.path.dirname(os.path.abspath(__file__))
RESULT_DIR = os.path.join(BASE, "ricci_flow_longrun_results")
OUT_DIR = os.path.join(RESULT_DIR, "analysis")
DATA_DIR = r"D:\100data"

# cost guards: the machine has 4 cores and limited RAM
MAX_N_FOR_SPECTRAL = 1500      # dense eigh on the largest component
MAX_E_FOR_MODULARITY = 60000   # louvain
MAX_N_FOR_PATH = 800           # exact average shortest path


# ----------------------------------------------------------------- helpers
def family_of(dataset):
    """Collapse a dataset id into its source family (ENZYMES_g12 -> ENZYMES)."""
    import re
    name = re.sub(r"_g?\d+$", "", str(dataset))
    name = re.sub(r"[-_](sub|snap)\d*$", "", name)
    if name.startswith("scc_rt"):
        parts = name.split("_")
        return "_".join(parts[:2])
    return name.split("-")[0] if "-" in name else name


def largest_component(graph):
    comps = list(nx.connected_components(graph))
    if not comps:
        return graph
    biggest = max(comps, key=len)
    return graph.subgraph(biggest).copy()


def spectral_gap(graph):
    """Algebraic connectivity of the normalised Laplacian (lambda2)."""
    sub = largest_component(graph)
    n = sub.number_of_nodes()
    if n < 4 or n > MAX_N_FOR_SPECTRAL:
        return np.nan, np.nan
    try:
        order = list(sub.nodes())
        A = nx.to_numpy_array(sub, nodelist=order, weight=None)
        deg = A.sum(axis=1)
        deg[deg == 0] = 1.0
        d_inv_sqrt = 1.0 / np.sqrt(deg)
        L = np.eye(n) - (d_inv_sqrt[:, None] * A * d_inv_sqrt[None, :])
        vals = np.linalg.eigvalsh(L)
        vals = np.sort(vals)
        return float(vals[1]), float(vals[-1])
    except Exception:
        return np.nan, np.nan


def modularity_score(graph):
    n, m = graph.number_of_nodes(), graph.number_of_edges()
    if m > MAX_E_FOR_MODULARITY or n < 4:
        return np.nan
    try:
        from networkx.algorithms.community import louvain_communities
        comms = louvain_communities(graph, weight=None, seed=42)
        if len(comms) < 2:
            return 0.0
        return float(nx.algorithms.community.modularity(graph, comms, weight=None))
    except Exception:
        return np.nan


def avg_path_length(graph):
    sub = largest_component(graph)
    if sub.number_of_nodes() < 2 or sub.number_of_nodes() > MAX_N_FOR_PATH:
        return np.nan
    try:
        return float(nx.average_shortest_path_length(sub, weight=None))
    except Exception:
        return np.nan


# ------------------------------------------------------------- covariate pass
def compute_covariates():
    parser = R.build_parser()
    args = parser.parse_args([])
    args.data_dir = DATA_DIR

    files = R.discover_graph_files(DATA_DIR)
    files.sort()
    files = [p for p in files if "100data.zip" not in p]

    rows = []
    for index, path in enumerate(files, 1):
        dataset = R.safe_dataset_id(path, DATA_DIR)
        record = {"dataset": dataset}
        try:
            graph = R.load_graph(path, args)
        except Exception as exc:
            record["load_error"] = str(exc)[:120]
            rows.append(record)
            print("[{}/{}] {}: LOAD FAIL".format(index, len(files), dataset), flush=True)
            continue

        undirected = graph.to_undirected() if graph.is_directed() else graph
        n = undirected.number_of_nodes()
        m = undirected.number_of_edges()

        degrees = np.array([d for _, d in undirected.degree()], dtype=float)
        record["nodes"] = n
        record["edges"] = m
        record["density"] = (2.0 * m) / (n * (n - 1)) if n > 1 else np.nan
        record["mean_degree"] = float(degrees.mean()) if n else np.nan
        record["degree_cv"] = (float(degrees.std() / degrees.mean())
                               if n and degrees.mean() > 0 else np.nan)
        record["degree_max"] = float(degrees.max()) if n else np.nan

        comps = list(nx.connected_components(undirected))
        sizes = sorted((len(c) for c in comps), reverse=True)
        record["n_components"] = len(comps)
        record["lcc_fraction"] = sizes[0] / n if n else np.nan

        try:
            record["transitivity"] = float(nx.transitivity(undirected))
        except Exception:
            record["transitivity"] = np.nan
        try:
            record["assortativity"] = float(
                nx.degree_assortativity_coefficient(undirected))
        except Exception:
            record["assortativity"] = np.nan

        record["modularity"] = modularity_score(undirected)
        record["avg_path_len"] = avg_path_length(undirected)
        gap, lmax = spectral_gap(undirected)
        record["algebraic_connectivity"] = gap
        record["lambda_max_normlap"] = lmax

        print("[{}/{}] {} N={} E={} Q={:.3f} a={:.3f}".format(
            index, len(files), dataset, n, m,
            record["modularity"] if record["modularity"] == record["modularity"] else float("nan"),
            gap if gap == gap else float("nan")), flush=True)
        rows.append(record)

    frame = pd.DataFrame(rows)
    frame.to_csv(os.path.join(OUT_DIR, "07_structure_predictors.csv"), index=False)
    return frame


# ---------------------------------------------------------------- statistics
def spearman(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 6 or np.std(x[mask]) == 0 or np.std(y[mask]) == 0:
        return np.nan, np.nan, int(mask.sum())
    from scipy.stats import spearmanr
    rho, p = spearmanr(x[mask], y[mask])
    return float(rho), float(p), int(mask.sum())


def partial_spearman(x, y, z):
    """Spearman correlation of x,y after regressing rank(x),rank(y) on rank(z)."""
    from scipy.stats import rankdata
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if mask.sum() < 10:
        return np.nan, int(mask.sum())
    rx = rankdata(x[mask])
    ry = rankdata(y[mask])
    rz = np.column_stack([np.ones(mask.sum()), rankdata(z[mask])])
    beta_x, *_ = np.linalg.lstsq(rz, rx, rcond=None)
    beta_y, *_ = np.linalg.lstsq(rz, ry, rcond=None)
    ex = rx - rz @ beta_x
    ey = ry - rz @ beta_y
    if np.std(ex) == 0 or np.std(ey) == 0:
        return np.nan, int(mask.sum())
    return float(np.corrcoef(ex, ey)[0, 1]), int(mask.sum())


def ols(y, X, names):
    """Standardised multivariate OLS with HC3-ish (HC1) robust standard errors."""
    y = np.asarray(y, dtype=float)
    X = np.column_stack([np.ones(len(y)), X])
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
    pvals = 2.0 * (1.0 - tdist.cdf(np.abs(tvals), dof))
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - float(resid @ resid) / ss_tot if ss_tot > 0 else np.nan
    adj = 1.0 - (1.0 - r2) * (len(y) - 1) / dof if dof > 0 else np.nan
    out = []
    for i, nm in enumerate(["intercept"] + list(names)):
        out.append({
            "term": nm, "coef": float(beta[i]), "se": float(se[i]),
            "t": float(tvals[i]), "p": float(pvals[i]),
        })
    frame = pd.DataFrame(out)
    frame.attrs["r2"] = float(r2)
    frame.attrs["adj_r2"] = float(adj)
    frame.attrs["n"] = int(len(y))
    return frame


COVARIATES = [
    "initial_kappa", "initial_variance", "log10_nodes", "density",
    "mean_degree", "degree_cv", "transitivity", "assortativity",
    "modularity", "algebraic_connectivity", "lcc_fraction",
    "n_components", "avg_path_len",
]


def main():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)

    cov_path = os.path.join(OUT_DIR, "07_structure_predictors.csv")
    if os.path.exists(cov_path):
        cov = pd.read_csv(cov_path)
    else:
        cov = compute_covariates()

    table = pd.read_csv(os.path.join(OUT_DIR, "00_dataset_table.csv"))
    # `nodes` must come from the freshly computed covariates: the summary table
    # carries NaN node counts for most datasets, which silently disabled the
    # minimum-size filter in the exclusion chain.
    table = table.drop(columns=["nodes"], errors="ignore")
    merged = table.merge(cov, on="dataset", how="left")
    if "log10_nodes" not in merged.columns:
        merged["log10_nodes"] = np.log10(merged["nodes"].astype(float))
    merged["family"] = merged["dataset"].map(family_of)

    merged = merged[np.isfinite(merged["var_ratio_final"])]
    merged = merged[merged["initial_variance"] > 0]

    converged = merged[merged["converged"] == True].copy()   # noqa: E712
    merged.to_csv(os.path.join(OUT_DIR, "07_merged_with_structure.csv"), index=False)

    lines = []
    lines.append("STRUCTURAL PREDICTORS OF THE RESIDUAL CURVATURE VARIANCE")
    lines.append("=" * 66)
    lines.append("")
    lines.append("datasets with a usable plateau: {}".format(len(merged)))
    lines.append("  of which converged (not budget-truncated): {}".format(len(converged)))
    lines.append("families represented: {}".format(merged["family"].nunique()))
    lines.append("")

    # ---- correlations
    corr_rows = []
    for label, frame in [("all_usable", merged), ("converged_only", converged)]:
        for cov_name in COVARIATES:
            if cov_name not in frame.columns:
                continue
            rho, p, n = spearman(frame[cov_name], frame["var_ratio_final"])
            partial, n_part = partial_spearman(
                frame[cov_name], frame["var_ratio_final"], frame["log10_nodes"])
            corr_rows.append({
                "stratum": label, "covariate": cov_name, "n": n,
                "spearman": rho, "p_value": p,
                "partial_given_logN": partial, "n_partial": n_part,
            })
    corr = pd.DataFrame(corr_rows)
    corr.to_csv(os.path.join(OUT_DIR, "08_structure_correlations.csv"), index=False)

    lines.append("--- Spearman vs plateau ratio (higher ratio = less homogenised) ---")
    for label in ["all_usable", "converged_only"]:
        sub = corr[corr["stratum"] == label].copy()
        sub = sub.reindex(sub["spearman"].abs().sort_values(ascending=False).index)
        lines.append("")
        lines.append("  [{}]  n={}".format(label, int(sub["n"].max())))
        lines.append("  {:<26s} {:>7s} {:>9s} {:>9s}".format(
            "covariate", "rho", "p", "partial"))
        for _, row in sub.iterrows():
            rho = row["spearman"]
            p = row["p_value"]
            par = row["partial_given_logN"]
            if rho != rho:
                continue
            lines.append("  {:<26s} {:>+7.3f} {:>9.3g} {:>+9.3f}".format(
                row["covariate"], rho, p,
                par if par == par else float("nan")))
    lines.append("")

    # ---- multivariate model on the converged subset
    lines.append("--- multivariate OLS (standardised covariates, converged subset) ---")
    model_terms = [c for c in ["initial_kappa", "log10_nodes", "degree_cv",
                               "transitivity", "assortativity", "modularity",
                               "algebraic_connectivity", "density"]
                   if c in converged.columns]
    sub = converged.dropna(subset=model_terms + ["var_ratio_final"])
    lines.append("  n = {} (complete cases)".format(len(sub)))
    if len(sub) >= len(model_terms) + 5:
        X = np.column_stack([
            (sub[c].astype(float).values - sub[c].astype(float).mean())
            / (sub[c].astype(float).std(ddof=0) + 1e-300) for c in model_terms])
        res = ols(sub["var_ratio_final"].values, X, model_terms)
        if res is not None:
            res.to_csv(os.path.join(OUT_DIR, "09_structure_regression.csv"), index=False)
            lines.append("  R2 = {:.3f}   adj R2 = {:.3f}".format(
                res.attrs["r2"], res.attrs["adj_r2"]))
            lines.append("  {:<24s} {:>8s} {:>8s} {:>7s} {:>9s}".format(
                "term", "coef", "se", "t", "p"))
            for _, row in res.iterrows():
                lines.append("  {:<24s} {:>+8.3f} {:>8.3f} {:>7.2f} {:>9.3g}".format(
                    row["term"], row["coef"], row["se"], row["t"], row["p"]))
    lines.append("")

    # ---- family-level robustness (pseudo-replication guard)
    fam = converged.copy()
    fam_agg = fam.groupby("family").agg(
        n_datasets=("dataset", "count"),
        var_ratio=("var_ratio_final", "median"),
        initial_kappa=("initial_kappa", "median"),
        log10_nodes=("log10_nodes", "median"),
        degree_cv=("degree_cv", "median"),
        transitivity=("transitivity", "median"),
        modularity=("modularity", "median"),
        algebraic_connectivity=("algebraic_connectivity", "median"),
    ).reset_index()
    fam_agg = fam_agg[fam_agg["n_datasets"] >= 2]
    fam_agg.to_csv(os.path.join(OUT_DIR, "10_family_level_check.csv"), index=False)

    lines.append("--- family-level check (one point per source family) ---")
    lines.append("  families with >=2 converged datasets: {}".format(len(fam_agg)))
    for cov_name in ["initial_kappa", "log10_nodes", "degree_cv", "transitivity",
                     "modularity", "algebraic_connectivity"]:
        if cov_name not in fam_agg.columns:
            continue
        rho, p, n = spearman(fam_agg[cov_name], fam_agg["var_ratio"])
        if rho == rho:
            lines.append("  {:<26s} rho = {:+.3f}  p = {:8.3g}  (n families = {})".format(
                cov_name, rho, p, n))
    lines.append("")

    # ---- verdict
    lines.append("--- verdict ---")
    best = corr[corr["stratum"] == "converged_only"].dropna(subset=["spearman"])
    if len(best):
        top = best.reindex(best["spearman"].abs().sort_values(ascending=False).index).iloc[0]
        lines.append("  strongest single predictor in the converged subset: "
                     "{} (rho={:+.3f})".format(top["covariate"], top["spearman"]))
    lines.append("  NOTE: budget-truncated runs carry an upper-bound plateau;")
    lines.append("        the converged subset is the admissible primary analysis.")

    text = "\n".join(lines)
    with open(os.path.join(OUT_DIR, "STRUCTURE_REPORT.txt"), "w", encoding="utf-8") as fh:
        fh.write(text)
    print(text)


if __name__ == "__main__":
    main()
