# -*- coding: utf-8 -*-
"""
ricci_flow_six_analysis_final.py

Final statistical analysis for the Ricci-flow universal-network experiment.

INPUT
-----
raw_results.csv:
dataset, step, mean_curvature, curvature_variance, lambda1,
nodes, edges, time, distance_time, sinkhorn_time, lambda_time

dataset_summary.csv:
domain, dataset, file, status, message, nodes, edges, components,
weight_mode, weight_reason, readme, readme_direction, duplicates,
invalid_lines, transport_cells_per_step, shortest_path_work, slope, intercept

This program:
1. Reconstructs dataset-level endpoint metrics from raw trajectories.
2. Performs paired non-parametric tests and effect sizes.
3. Bootstrap confidence intervals.
4. Leave-one-dataset-out robustness analysis.
5. Cross-domain validation using the supplied domain field.
6. Hierarchical/domain-aware analysis where available.
7. Network-size dependence analysis.
8. Permutation/null-model tests.
9. Repairs spectral analysis by explicitly distinguishing valid positive
   spectral gaps from zero/non-identifiable cases.
10. Produces publication-ready CSV tables and figures.

IMPORTANT:
- Dataset is the independent experimental unit.
- The 10 steps within a dataset are NOT treated as independent replicates.
- No TensorFlow is used.
- No sampling of graph trajectories is performed; this script analyzes all
  records present in raw_results.csv.
"""

from __future__ import annotations

import os
import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy import stats

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------
# USER CONFIGURATION
# ---------------------------------------------------------------------

RAW_RESULTS = Path(
    r"D:\2026codes\2026821100data\results_correct_fast\raw_results.csv"
)
DATASET_SUMMARY = Path(
    r"D:\2026codes\2026821100data\results_correct_fast\dataset_summary.csv"
)

OUTPUT_DIR = RAW_RESULTS.parent / "six_analysis_final"

N_BOOT = 10000
N_PERM = 20000
RANDOM_SEED = 20260829
ALPHA = 0.05

# Minimum number of datasets in a domain for domain-level inferential
# summaries. Smaller domains are still retained descriptively.
MIN_DOMAIN_N = 2

# ---------------------------------------------------------------------
# GENERAL UTILITIES
# ---------------------------------------------------------------------

RNG = np.random.default_rng(RANDOM_SEED)


def log(msg=""):
    print(msg, flush=True)


def ensure_dirs():
    for name in [
        "tables",
        "figures",
        "robustness",
        "domain",
        "null_model",
        "spectral",
        "logs",
    ]:
        (OUTPUT_DIR / name).mkdir(parents=True, exist_ok=True)


def clean_numeric(s):
    return pd.to_numeric(s, errors="coerce").replace(
        [np.inf, -np.inf], np.nan
    )


def mean_ci(values, confidence=0.95):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    n = len(x)
    if n < 2:
        return np.nan, np.nan, np.nan
    m = float(np.mean(x))
    se = float(stats.sem(x))
    tcrit = float(stats.t.ppf((1 + confidence) / 2, n - 1))
    return m, m - tcrit * se, m + tcrit * se


def bootstrap_stat(values, func=np.median, n_boot=N_BOOT, seed=RANDOM_SEED):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    n = len(x)
    if n == 0:
        return np.nan, np.nan, np.nan
    if n == 1:
        v = float(func(x))
        return v, v, v

    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_boot, n))
    samples = x[idx]
    boot = np.asarray([func(row) for row in samples], dtype=float)
    point = float(func(x))
    lo = float(np.quantile(boot, 0.025))
    hi = float(np.quantile(boot, 0.975))
    return point, lo, hi


def rank_biserial_from_wilcoxon(x):
    """
    Rank-biserial effect size for a paired Wilcoxon test.

    r_rb = (W_plus - W_minus) / (W_plus + W_minus)

    Positive means x tends to be positive.
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    x = x[x != 0]
    if len(x) == 0:
        return np.nan
    absx = np.abs(x)
    ranks = stats.rankdata(absx, method="average")
    w_plus = float(np.sum(ranks[x > 0]))
    w_minus = float(np.sum(ranks[x < 0]))
    denom = w_plus + w_minus
    return (w_plus - w_minus) / denom if denom else np.nan


def paired_test(x, alternative="two-sided"):
    """
    Paired Wilcoxon signed-rank test against zero.
    Returns statistic, p, rank-biserial effect, n_nonzero.
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    nonzero = x[x != 0]
    n = len(nonzero)

    if n == 0:
        return np.nan, np.nan, np.nan, 0

    try:
        res = stats.wilcoxon(
            nonzero,
            alternative=alternative,
            zero_method="wilcox",
            method="auto",
        )
        p = float(res.pvalue)
        stat = float(res.statistic)
    except Exception:
        # Exact calculation can fail for some tied/large inputs.
        try:
            res = stats.wilcoxon(
                nonzero,
                alternative=alternative,
                zero_method="wilcox",
                method="approx",
            )
            p = float(res.pvalue)
            stat = float(res.statistic)
        except Exception:
            p = np.nan
            stat = np.nan

    return stat, p, rank_biserial_from_wilcoxon(nonzero), n


def paired_permutation_p(x, n_perm=N_PERM, seed=RANDOM_SEED):
    """
    Sign-flip permutation test for a paired endpoint change.

    Null: the sign of each dataset-level change is exchangeable.
    Two-sided p-value based on the absolute mean statistic.
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = len(x)
    if n == 0:
        return np.nan, np.nan

    observed = float(np.mean(x))
    rng = np.random.default_rng(seed)
    signs = rng.choice(np.array([-1.0, 1.0]), size=(n_perm, n))
    perm_stats = np.mean(signs * x[None, :], axis=1)
    p = (1.0 + np.sum(np.abs(perm_stats) >= abs(observed))) / (n_perm + 1.0)
    return observed, float(p)


def spearman_with_ci(x, y, n_boot=N_BOOT, seed=RANDOM_SEED):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = len(x)
    if n < 3:
        return np.nan, np.nan, np.nan, np.nan, n

    rho, p = stats.spearmanr(x, y)

    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(x[idx])) < 2 or len(np.unique(y[idx])) < 2:
            continue
        vals.append(stats.spearmanr(x[idx], y[idx]).statistic)

    if len(vals) >= 100:
        lo, hi = np.quantile(vals, [0.025, 0.975])
    else:
        lo, hi = np.nan, np.nan

    return float(rho), float(p), float(lo), float(hi), n


def benjamini_hochberg(pvalues):
    p = np.asarray(pvalues, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    mask = np.isfinite(p)
    if not np.any(mask):
        return out

    pv = p[mask]
    m = len(pv)
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)

    tmp = np.empty(m)
    tmp[order] = q
    out[mask] = tmp
    return out


def standardized_name(x):
    return str(x).strip()


# ---------------------------------------------------------------------
# LOAD AND RECONSTRUCT DATA
# ---------------------------------------------------------------------

def load_inputs():
    if not RAW_RESULTS.exists():
        raise FileNotFoundError(f"Raw results not found: {RAW_RESULTS}")
    if not DATASET_SUMMARY.exists():
        raise FileNotFoundError(
            f"Dataset summary not found: {DATASET_SUMMARY}"
        )

    log(f"Loading raw results: {RAW_RESULTS}")
    raw = pd.read_csv(RAW_RESULTS)

    log(f"Loading dataset summary: {DATASET_SUMMARY}")
    summary = pd.read_csv(DATASET_SUMMARY)

    required_raw = {
        "dataset",
        "step",
        "mean_curvature",
        "curvature_variance",
        "lambda1",
        "nodes",
        "edges",
    }
    required_summary = {"domain", "dataset"}

    missing_raw = sorted(required_raw - set(raw.columns))
    missing_summary = sorted(required_summary - set(summary.columns))

    if missing_raw:
        raise ValueError(f"raw_results.csv missing columns: {missing_raw}")
    if missing_summary:
        raise ValueError(
            f"dataset_summary.csv missing columns: {missing_summary}"
        )

    for col in [
        "step",
        "mean_curvature",
        "curvature_variance",
        "lambda1",
        "nodes",
        "edges",
    ]:
        raw[col] = clean_numeric(raw[col])

    raw["dataset"] = raw["dataset"].astype(str).str.strip()
    summary["dataset"] = summary["dataset"].astype(str).str.strip()
    summary["domain"] = summary["domain"].fillna("Unknown").astype(str).str.strip()

    # There may be status rows in summary without valid raw trajectories.
    # Keep only the best/first metadata row per dataset.
    summary = (
        summary.sort_values(["dataset"])
        .drop_duplicates(subset=["dataset"], keep="first")
        .copy()
    )

    raw = raw.sort_values(["dataset", "step"]).reset_index(drop=True)

    log(f"Raw records: {len(raw)}")
    log(f"Raw datasets: {raw['dataset'].nunique()}")
    log(f"Summary datasets: {summary['dataset'].nunique()}")

    return raw, summary


def build_dataset_metrics(raw, summary):
    records = []

    for dataset, sub in raw.groupby("dataset", sort=True):
        sub = sub.sort_values("step").copy()

        # A valid trajectory requires at least two distinct steps.
        sub = sub.drop_duplicates(subset=["step"], keep="last")

        if len(sub) < 2:
            continue

        first = sub.iloc[0]
        last = sub.iloc[-1]

        k0 = float(first["mean_curvature"])
        kT = float(last["mean_curvature"])
        v0 = float(first["curvature_variance"])
        vT = float(last["curvature_variance"])
        l0 = float(first["lambda1"]) if np.isfinite(first["lambda1"]) else np.nan
        lT = float(last["lambda1"]) if np.isfinite(last["lambda1"]) else np.nan

        delta_k = kT - k0
        var_red = np.nan
        if np.isfinite(v0) and v0 != 0 and np.isfinite(vT):
            var_red = (v0 - vT) / abs(v0)

        lambda_imp = np.nan
        if np.isfinite(l0) and np.isfinite(lT):
            # Spectral improvement is only identifiable when initial
            # lambda1 is strictly positive. A zero gap is not an
            # "improvement denominator"; it is a non-identifiable case.
            if l0 > 0:
                lambda_imp = (lT - l0) / abs(l0)

        normalized_rate = np.nan
        if np.isfinite(delta_k):
            normalized_rate = delta_k / max(1, len(sub) - 1)

        # Robust convergence score:
        # fraction of step-to-step variance changes that are non-increasing.
        variance_steps = sub["curvature_variance"].values.astype(float)
        finite_v = np.isfinite(variance_steps)
        vv = variance_steps[finite_v]
        if len(vv) >= 2:
            dv = np.diff(vv)
            convergence_score = float(np.mean(dv <= 1e-12))
            variance_monotonic = bool(np.all(dv <= 1e-12))
        else:
            convergence_score = np.nan
            variance_monotonic = False

        # Linear slope of curvature versus step, useful but NOT a primary
        # inferential endpoint.
        x = sub["step"].values.astype(float)
        y = sub["mean_curvature"].values.astype(float)
        mask = np.isfinite(x) & np.isfinite(y)
        if np.sum(mask) >= 2:
            slope, intercept, _, _, _ = stats.linregress(x[mask], y[mask])
        else:
            slope, intercept = np.nan, np.nan

        nodes = first["nodes"]
        edges = first["edges"]

        records.append(
            {
                "dataset": dataset,
                "nodes": float(nodes) if np.isfinite(nodes) else np.nan,
                "edges": float(edges) if np.isfinite(edges) else np.nan,
                "n_steps": int(len(sub)),
                "initial_kappa": k0,
                "final_kappa": kT,
                "delta_kappa": delta_k,
                "initial_variance": v0,
                "final_variance": vT,
                "variance_reduction": var_red,
                "lambda_initial": l0,
                "lambda_final": lT,
                "lambda_improvement": lambda_imp,
                "normalized_rate": normalized_rate,
                "convergence_score": convergence_score,
                "variance_monotonic": variance_monotonic,
                "curvature_slope": slope,
                "curvature_intercept": intercept,
            }
        )

    metrics = pd.DataFrame(records)

    # Merge supplied domain metadata.
    domain_cols = ["dataset", "domain"]
    for c in summary.columns:
        if c not in domain_cols:
            domain_cols.append(c)

    merged = metrics.merge(
        summary[domain_cols],
        on="dataset",
        how="left",
        suffixes=("", "_summary"),
    )

    merged["domain"] = merged["domain"].fillna("Unknown")

    return merged


# ---------------------------------------------------------------------
# ANALYSIS 1: UNIVERSAL EVOLUTION / PAIRED TESTS
# ---------------------------------------------------------------------

PRIMARY_ENDPOINTS = {
    "delta_kappa": "delta_kappa",
    "variance_reduction": "variance_reduction",
    "lambda_improvement": "lambda_improvement",
    "normalized_rate": "normalized_rate",
    "convergence_score": "convergence_score",
}


def analyze_universal(metrics):
    rows = []

    for endpoint, col in PRIMARY_ENDPOINTS.items():
        x = clean_numeric(metrics[col]).dropna().values

        if endpoint == "convergence_score":
            # Treat as a descriptive convergence proportion; also test
            # whether the per-dataset score differs from 0.5.
            diff = x - 0.5
        else:
            diff = x

        med = float(np.median(x)) if len(x) else np.nan
        mean = float(np.mean(x)) if len(x) else np.nan
        sd = float(np.std(x, ddof=1)) if len(x) > 1 else np.nan

        _, p, effect, n_nonzero = paired_test(diff)

        boot_median, boot_lo, boot_hi = bootstrap_stat(
            x, func=np.median, seed=RANDOM_SEED + hash(endpoint) % 10000
        )
        boot_mean, boot_mean_lo, boot_mean_hi = bootstrap_stat(
            x, func=np.mean, seed=RANDOM_SEED + 10000 + hash(endpoint) % 10000
        )

        perm_mean, perm_p = paired_permutation_p(
            diff,
            n_perm=N_PERM,
            seed=RANDOM_SEED + 20000 + hash(endpoint) % 10000,
        )

        mean_m, mean_lo, mean_hi = mean_ci(x)

        rows.append(
            {
                "endpoint": endpoint,
                "n": len(x),
                "mean": mean,
                "sd": sd,
                "median": med,
                "mean_ci95_low_t": mean_lo,
                "mean_ci95_high_t": mean_hi,
                "bootstrap_mean": boot_mean,
                "bootstrap_mean_ci95_low": boot_mean_lo,
                "bootstrap_mean_ci95_high": boot_mean_hi,
                "bootstrap_median": boot_median,
                "bootstrap_median_ci95_low": boot_lo,
                "bootstrap_median_ci95_high": boot_hi,
                "wilcoxon_statistic": _,
                "wilcoxon_p": p,
                "rank_biserial": effect,
                "wilcoxon_nonzero_n": n_nonzero,
                "sign_flip_mean": perm_mean,
                "sign_flip_p": perm_p,
            }
        )

    out = pd.DataFrame(rows)
    out["wilcoxon_q_BH"] = benjamini_hochberg(out["wilcoxon_p"].values)
    out["sign_flip_q_BH"] = benjamini_hochberg(out["sign_flip_p"].values)

    out.to_csv(
        OUTPUT_DIR / "tables" / "01_universal_endpoint_statistics.csv",
        index=False,
    )

    return out


# ---------------------------------------------------------------------
# ANALYSIS 2: ROBUSTNESS / LEAVE-ONE-OUT
# ---------------------------------------------------------------------

def leave_one_out(metrics, endpoint):
    x = metrics[["dataset", endpoint]].copy()
    x[endpoint] = clean_numeric(x[endpoint])
    x = x.dropna()

    if len(x) < 3:
        return pd.DataFrame()

    full_mean = float(x[endpoint].mean())
    full_median = float(x[endpoint].median())

    rows = []
    for i, row in x.iterrows():
        y = x.drop(index=i)[endpoint].values
        rows.append(
            {
                "excluded_dataset": row["dataset"],
                "n_remaining": len(y),
                "mean": float(np.mean(y)),
                "median": float(np.median(y)),
                "mean_change_from_full": float(np.mean(y) - full_mean),
                "median_change_from_full": float(np.median(y) - full_median),
            }
        )

    out = pd.DataFrame(rows)
    out.to_csv(
        OUTPUT_DIR / "robustness" / f"loo_{endpoint}.csv",
        index=False,
    )

    return out


def analyze_robustness(metrics):
    summary_rows = []

    for endpoint in [
        "delta_kappa",
        "variance_reduction",
        "lambda_improvement",
        "normalized_rate",
    ]:
        loo = leave_one_out(metrics, endpoint)
        if loo.empty:
            continue

        x = clean_numeric(metrics[endpoint]).dropna().values
        full_mean = float(np.mean(x))
        full_median = float(np.median(x))

        max_abs_mean_shift = float(
            np.max(np.abs(loo["mean_change_from_full"]))
        )
        max_abs_median_shift = float(
            np.max(np.abs(loo["median_change_from_full"]))
        )

        max_mean_row = loo.iloc[
            int(np.argmax(np.abs(loo["mean_change_from_full"].values)))
        ]
        max_median_row = loo.iloc[
            int(np.argmax(np.abs(loo["median_change_from_full"].values)))
        ]

        summary_rows.append(
            {
                "endpoint": endpoint,
                "n": len(x),
                "full_mean": full_mean,
                "full_median": full_median,
                "max_abs_loo_mean_shift": max_abs_mean_shift,
                "dataset_causing_max_mean_shift": max_mean_row[
                    "excluded_dataset"
                ],
                "max_abs_loo_median_shift": max_abs_median_shift,
                "dataset_causing_max_median_shift": max_median_row[
                    "excluded_dataset"
                ],
            }
        )

    out = pd.DataFrame(summary_rows)
    out.to_csv(
        OUTPUT_DIR / "tables" / "02_leave_one_out_summary.csv",
        index=False,
    )

    # Top influential datasets by absolute standardized endpoint response.
    rows = []
    for endpoint in [
        "delta_kappa",
        "variance_reduction",
        "lambda_improvement",
        "normalized_rate",
    ]:
        x = clean_numeric(metrics[endpoint])
        med = float(x.median())
        mad = float(np.median(np.abs(x.dropna() - med)))
        denom = mad if mad > 0 else float(x.std())
        denom = denom if np.isfinite(denom) and denom > 0 else 1.0

        tmp = metrics[["dataset", endpoint]].copy()
        tmp["robust_z"] = (tmp[endpoint] - med) / denom
        tmp["endpoint"] = endpoint
        rows.append(tmp[["endpoint", "dataset", endpoint, "robust_z"]])

    influence = pd.concat(rows, ignore_index=True)
    influence["abs_robust_z"] = influence["robust_z"].abs()
    influence.sort_values(
        ["endpoint", "abs_robust_z"], ascending=[True, False]
    ).to_csv(
        OUTPUT_DIR / "robustness" / "dataset_influence_scores.csv",
        index=False,
    )

    return out


# ---------------------------------------------------------------------
# ANALYSIS 3: CROSS-DOMAIN VALIDATION
# ---------------------------------------------------------------------

def domain_summary(metrics):
    rows = []

    for domain, sub in metrics.groupby("domain", sort=True):
        n = len(sub)

        row = {
            "domain": domain,
            "n_datasets": n,
            "median_nodes": float(sub["nodes"].median()),
            "median_edges": float(sub["edges"].median()),
        }

        for endpoint in [
            "delta_kappa",
            "variance_reduction",
            "lambda_improvement",
            "normalized_rate",
            "convergence_score",
        ]:
            x = clean_numeric(sub[endpoint]).dropna().values
            row[f"{endpoint}_n"] = len(x)
            row[f"{endpoint}_mean"] = (
                float(np.mean(x)) if len(x) else np.nan
            )
            row[f"{endpoint}_median"] = (
                float(np.median(x)) if len(x) else np.nan
            )
            if len(x) >= 2:
                _, lo, hi = bootstrap_stat(
                    x,
                    func=np.median,
                    n_boot=N_BOOT,
                    seed=RANDOM_SEED + len(rows) * 100 + len(endpoint),
                )
                row[f"{endpoint}_median_ci95_low"] = lo
                row[f"{endpoint}_median_ci95_high"] = hi
            else:
                row[f"{endpoint}_median_ci95_low"] = np.nan
                row[f"{endpoint}_median_ci95_high"] = np.nan

            if endpoint != "convergence_score" and len(x) >= MIN_DOMAIN_N:
                _, p, effect, nz = paired_test(x)
                row[f"{endpoint}_wilcoxon_p"] = p
                row[f"{endpoint}_rank_biserial"] = effect
                row[f"{endpoint}_test_n"] = nz
            else:
                row[f"{endpoint}_wilcoxon_p"] = np.nan
                row[f"{endpoint}_rank_biserial"] = np.nan
                row[f"{endpoint}_test_n"] = np.nan

        rows.append(row)

    out = pd.DataFrame(rows)

    for endpoint in [
        "delta_kappa",
        "variance_reduction",
        "lambda_improvement",
        "normalized_rate",
    ]:
        col = f"{endpoint}_wilcoxon_p"
        out[f"{endpoint}_wilcoxon_q_BH"] = benjamini_hochberg(
            out[col].values
        )

    out.to_csv(
        OUTPUT_DIR / "domain" / "03_cross_domain_summary.csv",
        index=False,
    )

    return out


def domain_effect_heterogeneity(metrics, endpoint):
    """
    Cochran-like heterogeneity diagnostic using inverse-variance weights
    at the domain level. This is a descriptive meta-analytic diagnostic,
    not a claim that datasets are iid.
    """
    rows = []

    for domain, sub in metrics.groupby("domain", sort=True):
        x = clean_numeric(sub[endpoint]).dropna().values
        if len(x) < 2:
            continue

        mean = float(np.mean(x))
        var = float(np.var(x, ddof=1))
        se2 = var / len(x) if var > 0 else np.nan

        rows.append(
            {
                "domain": domain,
                "n": len(x),
                "mean": mean,
                "variance": var,
                "se2": se2,
            }
        )

    d = pd.DataFrame(rows)

    if len(d) < 2:
        return d, np.nan, np.nan, np.nan

    valid = np.isfinite(d["se2"]) & (d["se2"] > 0)
    d = d.loc[valid].copy()

    if len(d) < 2:
        return d, np.nan, np.nan, np.nan

    w = 1.0 / d["se2"].values
    y = d["mean"].values
    pooled = np.sum(w * y) / np.sum(w)
    Q = float(np.sum(w * (y - pooled) ** 2))
    df = len(y) - 1
    p = float(stats.chi2.sf(Q, df))
    I2 = max(0.0, (Q - df) / Q) if Q > 0 else 0.0

    return d, Q, p, I2


def analyze_cross_domain(metrics):
    out = domain_summary(metrics)

    hetero_rows = []
    for endpoint in [
        "delta_kappa",
        "variance_reduction",
        "lambda_improvement",
        "normalized_rate",
    ]:
        d, Q, p, I2 = domain_effect_heterogeneity(metrics, endpoint)
        hetero_rows.append(
            {
                "endpoint": endpoint,
                "domains_in_meta_diagnostic": len(d),
                "Q": Q,
                "Q_p": p,
                "I2": I2,
            }
        )

    hetero = pd.DataFrame(hetero_rows)
    hetero.to_csv(
        OUTPUT_DIR / "domain" / "03b_domain_heterogeneity.csv",
        index=False,
    )

    return out, hetero


# ---------------------------------------------------------------------
# ANALYSIS 4: HIERARCHICAL / DOMAIN-AWARE TESTS
# ---------------------------------------------------------------------

def try_mixed_effects(metrics):
    """
    Fit a random-intercept model with domain as a grouping factor if
    statsmodels is available.

    Model:
        endpoint ~ 1 + log10(nodes) + (1 | domain)

    This is intentionally implemented as an optional analysis. The script
    remains runnable without statsmodels.
    """
    try:
        import statsmodels.formula.api as smf
    except Exception as exc:
        pd.DataFrame(
            [{"status": "statsmodels_unavailable", "message": repr(exc)}]
        ).to_csv(
            OUTPUT_DIR / "tables" / "04_mixed_effects_status.csv",
            index=False,
        )
        return

    rows = []

    df = metrics.copy()
    df["log10_nodes"] = np.log10(df["nodes"].clip(lower=1))

    endpoints = [
        "variance_reduction",
        "delta_kappa",
        "normalized_rate",
        "lambda_improvement",
    ]

    for endpoint in endpoints:
        sub = df[
            ["domain", "log10_nodes", endpoint]
        ].dropna().copy()

        if len(sub) < 10 or sub["domain"].nunique() < 2:
            rows.append(
                {
                    "endpoint": endpoint,
                    "status": "insufficient_data",
                }
            )
            continue

        try:
            model = smf.mixedlm(
                f"{endpoint} ~ log10_nodes",
                sub,
                groups=sub["domain"],
            )
            fit = model.fit(
                reml=True,
                method="lbfgs",
                maxiter=2000,
                disp=False,
            )

            rows.append(
                {
                    "endpoint": endpoint,
                    "status": "ok",
                    "n": len(sub),
                    "domains": sub["domain"].nunique(),
                    "intercept": float(fit.params.get("Intercept", np.nan)),
                    "log10_nodes_beta": float(
                        fit.params.get("log10_nodes", np.nan)
                    ),
                    "log10_nodes_p": float(
                        fit.pvalues.get("log10_nodes", np.nan)
                    ),
                    "domain_random_var": float(
                        fit.cov_re.iloc[0, 0]
                        if fit.cov_re.shape[0]
                        else np.nan
                    ),
                    "aic": float(fit.aic) if np.isfinite(fit.aic) else np.nan,
                    "bic": float(fit.bic) if np.isfinite(fit.bic) else np.nan,
                }
            )
        except Exception as exc:
            rows.append(
                {
                    "endpoint": endpoint,
                    "status": "fit_failed",
                    "message": repr(exc),
                }
            )

    pd.DataFrame(rows).to_csv(
        OUTPUT_DIR / "tables" / "04_mixed_effects_models.csv",
        index=False,
    )


# ---------------------------------------------------------------------
# ANALYSIS 5: NETWORK-SIZE DEPENDENCE
# ---------------------------------------------------------------------

def analyze_size_dependence(metrics):
    rows = []

    for endpoint in [
        "delta_kappa",
        "variance_reduction",
        "lambda_improvement",
        "normalized_rate",
        "convergence_score",
    ]:
        sub = metrics[["nodes", "edges", endpoint]].copy()
        sub = sub.dropna()
        if len(sub) < 3:
            continue

        log_nodes = np.log10(sub["nodes"].clip(lower=1))
        log_edges = np.log10(sub["edges"].clip(lower=1))

        rho_n, p_n, lo_n, hi_n, n_n = spearman_with_ci(
            log_nodes.values,
            sub[endpoint].values,
            seed=RANDOM_SEED + 30000,
        )
        rho_e, p_e, lo_e, hi_e, n_e = spearman_with_ci(
            log_edges.values,
            sub[endpoint].values,
            seed=RANDOM_SEED + 31000,
        )

        rows.append(
            {
                "endpoint": endpoint,
                "n_nodes_analysis": n_n,
                "spearman_rho_log10_nodes": rho_n,
                "p_log10_nodes": p_n,
                "rho_ci95_low_nodes": lo_n,
                "rho_ci95_high_nodes": hi_n,
                "n_edges_analysis": n_e,
                "spearman_rho_log10_edges": rho_e,
                "p_log10_edges": p_e,
                "rho_ci95_low_edges": lo_e,
                "rho_ci95_high_edges": hi_e,
            }
        )

    out = pd.DataFrame(rows)

    if len(out):
        out["q_log10_nodes_BH"] = benjamini_hochberg(
            out["p_log10_nodes"].values
        )
        out["q_log10_edges_BH"] = benjamini_hochberg(
            out["p_log10_edges"].values
        )

    out.to_csv(
        OUTPUT_DIR / "tables" / "05_network_size_dependence.csv",
        index=False,
    )

    return out


# ---------------------------------------------------------------------
# ANALYSIS 6: NULL / PERMUTATION MODEL
# ---------------------------------------------------------------------

def permutation_null_endpoint(values, statistic=np.median, n_perm=N_PERM,
                              seed=RANDOM_SEED):
    """
    Null model for endpoint sign.

    For each dataset, randomly flip the sign of its endpoint. This preserves
    the empirical magnitude distribution and destroys a systematic directional
    effect. The statistic is the median of the sign-flipped sample.

    This is appropriate for paired endpoint changes because the null is
    symmetry around zero.
    """
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    n = len(x)

    if n == 0:
        return np.nan, np.nan, np.nan, np.nan

    observed = float(statistic(x))
    rng = np.random.default_rng(seed)

    batch = 1000
    count = 0
    total = 0

    for start in range(0, n_perm, batch):
        b = min(batch, n_perm - start)
        signs = rng.choice([-1.0, 1.0], size=(b, n))
        vals = signs * x[None, :]
        statv = np.median(vals, axis=1) if statistic is np.median else np.mean(
            vals, axis=1
        )
        count += int(np.sum(np.abs(statv) >= abs(observed)))
        total += b

    p = (1 + count) / (1 + total)
    return observed, float(p), count, total


def analyze_null_models(metrics):
    rows = []

    for endpoint in [
        "delta_kappa",
        "variance_reduction",
        "lambda_improvement",
        "normalized_rate",
    ]:
        x = clean_numeric(metrics[endpoint]).dropna().values
        if len(x) < 3:
            continue

        obs_med, p_med, _, _ = permutation_null_endpoint(
            x,
            statistic=np.median,
            n_perm=N_PERM,
            seed=RANDOM_SEED + 40000 + len(rows),
        )
        obs_mean, p_mean, _, _ = permutation_null_endpoint(
            x,
            statistic=np.mean,
            n_perm=N_PERM,
            seed=RANDOM_SEED + 41000 + len(rows),
        )

        rows.append(
            {
                "endpoint": endpoint,
                "n": len(x),
                "observed_median": obs_med,
                "median_sign_flip_p": p_med,
                "observed_mean": obs_mean,
                "mean_sign_flip_p": p_mean,
            }
        )

    out = pd.DataFrame(rows)
    if len(out):
        out["median_q_BH"] = benjamini_hochberg(
            out["median_sign_flip_p"].values
        )
        out["mean_q_BH"] = benjamini_hochberg(
            out["mean_sign_flip_p"].values
        )

    out.to_csv(
        OUTPUT_DIR / "null_model" / "06_sign_flip_null_model.csv",
        index=False,
    )

    return out


# ---------------------------------------------------------------------
# SPECTRAL REPAIR / ANALYSIS
# ---------------------------------------------------------------------

def analyze_spectral(metrics):
    """
    Spectral analysis is separated into:
    A. valid positive lambda1 cases;
    B. lambda1 == 0 / non-positive cases, which are reported as
       non-identifiable rather than forced into a percentage improvement.

    This prevents the previous bio-CE-GT-type zero-gap issue from being
    silently interpreted as an ordinary spectral improvement.
    """
    df = metrics.copy()

    df["lambda_initial_status"] = np.select(
        [
            df["lambda_initial"].isna(),
            df["lambda_initial"] <= 0,
        ],
        [
            "missing",
            "zero_or_nonpositive",
        ],
        default="positive",
    )

    valid = df[
        (df["lambda_initial"] > 0)
        & np.isfinite(df["lambda_initial"])
        & np.isfinite(df["lambda_final"])
    ].copy()

    invalid = df[df["lambda_initial_status"] != "positive"].copy()

    valid_rows = []
    if len(valid):
        x = valid["lambda_improvement"].values
        stat, p, effect, n = paired_test(x)
        med, lo, hi = bootstrap_stat(
            x,
            func=np.median,
            seed=RANDOM_SEED + 50000,
        )
        mean, mean_lo, mean_hi = bootstrap_stat(
            x,
            func=np.mean,
            seed=RANDOM_SEED + 50001,
        )

        valid_rows.append(
            {
                "status": "positive_initial_gap",
                "n": len(x),
                "mean_improvement": float(np.mean(x)),
                "median_improvement": float(np.median(x)),
                "bootstrap_median": med,
                "bootstrap_median_ci95_low": lo,
                "bootstrap_median_ci95_high": hi,
                "bootstrap_mean": mean,
                "bootstrap_mean_ci95_low": mean_lo,
                "bootstrap_mean_ci95_high": mean_hi,
                "wilcoxon_statistic": stat,
                "wilcoxon_p": p,
                "rank_biserial": effect,
                "wilcoxon_nonzero_n": n,
            }
        )

    invalid_summary = (
        invalid.groupby("lambda_initial_status")
        .size()
        .reset_index(name="n")
        if len(invalid)
        else pd.DataFrame(columns=["lambda_initial_status", "n"])
    )

    spectral_summary = pd.DataFrame(valid_rows)
    spectral_summary.to_csv(
        OUTPUT_DIR / "spectral" / "07_spectral_valid_statistics.csv",
        index=False,
    )
    invalid_summary.to_csv(
        OUTPUT_DIR / "spectral" / "07_spectral_nonidentifiable_cases.csv",
        index=False,
    )

    # Save the dataset-level spectral audit.
    audit_cols = [
        "dataset",
        "domain",
        "lambda_initial",
        "lambda_final",
        "lambda_improvement",
        "lambda_initial_status",
    ]
    df[audit_cols].sort_values(
        ["lambda_initial_status", "dataset"]
    ).to_csv(
        OUTPUT_DIR / "spectral" / "07_spectral_dataset_audit.csv",
        index=False,
    )

    return spectral_summary, invalid_summary


# ---------------------------------------------------------------------
# FIGURES
# ---------------------------------------------------------------------

def savefig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_endpoint_distribution(metrics, endpoint, title, filename,
                               xlabel=None):
    x = clean_numeric(metrics[endpoint]).dropna().values
    if len(x) == 0:
        return

    plt.figure(figsize=(8, 5))
    plt.hist(x, bins=25)
    plt.axvline(0, linestyle="--", linewidth=1)
    plt.xlabel(xlabel or endpoint)
    plt.ylabel("Number of datasets")
    plt.title(title)
    savefig(OUTPUT_DIR / "figures" / filename)


def plot_domain_forest(domain_df, endpoint, title, filename):
    med_col = f"{endpoint}_median"
    lo_col = f"{endpoint}_median_ci95_low"
    hi_col = f"{endpoint}_median_ci95_high"

    d = domain_df[
        ["domain", "n_datasets", med_col, lo_col, hi_col]
    ].dropna()

    if d.empty:
        return

    d = d.sort_values(med_col)
    y = np.arange(len(d))

    plt.figure(figsize=(10, max(5, 0.42 * len(d))))
    x = d[med_col].values
    lo = x - d[lo_col].values
    hi = d[hi_col].values - x

    plt.errorbar(
        x,
        y,
        xerr=np.vstack([lo, hi]),
        fmt="o",
        capsize=3,
    )
    plt.axvline(0, linestyle="--", linewidth=1)
    plt.yticks(y, d["domain"].values)
    plt.xlabel(endpoint)
    plt.title(title)
    savefig(OUTPUT_DIR / "figures" / filename)


def plot_size_relationship(metrics, endpoint, filename, title):
    d = metrics[["nodes", endpoint]].dropna()
    if len(d) < 3:
        return

    x = np.log10(d["nodes"].clip(lower=1).values)
    y = d[endpoint].values

    plt.figure(figsize=(7, 5))
    plt.scatter(x, y, alpha=0.75)
    if len(x) >= 3:
        slope, intercept, _, _, _ = stats.linregress(x, y)
        xx = np.linspace(np.min(x), np.max(x), 100)
        plt.plot(xx, slope * xx + intercept)
    plt.xlabel("log10(number of nodes)")
    plt.ylabel(endpoint)
    plt.title(title)
    savefig(OUTPUT_DIR / "figures" / filename)


def plot_robustness(metrics, endpoint, filename):
    loo = leave_one_out(metrics, endpoint)
    if loo.empty:
        return

    x = np.arange(len(loo))
    y = loo["mean"].values

    plt.figure(figsize=(12, 5))
    plt.plot(x, y, linewidth=1)
    plt.axhline(float(metrics[endpoint].mean()), linestyle="--", linewidth=1)
    plt.xlabel("Leave-one-dataset-out iteration")
    plt.ylabel(f"Mean {endpoint}")
    plt.title(f"Leave-one-dataset-out robustness: {endpoint}")
    savefig(OUTPUT_DIR / "figures" / filename)


def plot_trajectory_envelope(raw, filename="fig_universal_trajectory_envelope.png"):
    """
    Plot all dataset trajectories after per-dataset normalization.

    This does NOT treat steps as independent replicates. It is a visualization
    of complete trajectories.
    """
    curves = []
    for dataset, sub in raw.groupby("dataset"):
        sub = sub.sort_values("step")
        if len(sub) < 2:
            continue
        y = sub["mean_curvature"].values.astype(float)
        if not np.all(np.isfinite(y)):
            continue
        scale = max(abs(y[0]), 1e-12)
        curves.append(y / scale)

    if not curves:
        return

    arr = np.asarray(curves)
    steps = np.arange(arr.shape[1])

    q25, med, q75 = np.quantile(arr, [0.25, 0.5, 0.75], axis=0)

    plt.figure(figsize=(8, 5))
    plt.plot(steps, med, linewidth=2, label="Median")
    plt.fill_between(steps, q25, q75, alpha=0.2, label="IQR")
    plt.axhline(1, linestyle="--", linewidth=1)
    plt.xlabel("Ricci-flow step")
    plt.ylabel("Mean curvature / initial curvature")
    plt.title("Normalized curvature trajectories across datasets")
    plt.legend()
    savefig(OUTPUT_DIR / "figures" / filename)


def plot_variance_trajectories(raw):
    curves = []
    for dataset, sub in raw.groupby("dataset"):
        sub = sub.sort_values("step")
        y = sub["curvature_variance"].values.astype(float)
        if len(y) == 10 and np.all(np.isfinite(y)) and y[0] > 0:
            curves.append(y / y[0])

    if not curves:
        return

    arr = np.asarray(curves)
    steps = np.arange(arr.shape[1])
    q10, q50, q90 = np.quantile(arr, [0.10, 0.50, 0.90], axis=0)

    plt.figure(figsize=(8, 5))
    plt.plot(steps, q50, linewidth=2, label="Median")
    plt.fill_between(steps, q10, q90, alpha=0.2, label="10–90%")
    plt.axhline(1, linestyle="--", linewidth=1)
    plt.xlabel("Ricci-flow step")
    plt.ylabel("Curvature variance / initial variance")
    plt.title("Curvature variance evolution")
    plt.legend()
    savefig(OUTPUT_DIR / "figures" / "fig_variance_evolution.png")


def plot_spectral_trajectories(raw):
    curves = []
    for dataset, sub in raw.groupby("dataset"):
        sub = sub.sort_values("step")
        y = sub["lambda1"].values.astype(float)
        if len(y) == 10 and np.all(np.isfinite(y)) and y[0] > 0:
            curves.append(y / y[0])

    if not curves:
        return

    arr = np.asarray(curves)
    steps = np.arange(arr.shape[1])
    q25, q50, q75 = np.quantile(arr, [0.25, 0.5, 0.75], axis=0)

    plt.figure(figsize=(8, 5))
    plt.plot(steps, q50, linewidth=2, label="Median")
    plt.fill_between(steps, q25, q75, alpha=0.2, label="IQR")
    plt.axhline(1, linestyle="--", linewidth=1)
    plt.xlabel("Ricci-flow step")
    plt.ylabel("λ1 / initial λ1")
    plt.title("Spectral-gap evolution: positive initial λ1 only")
    plt.legend()
    savefig(OUTPUT_DIR / "figures" / "fig_spectral_evolution_repaired.png")


def plot_correlation_matrix(metrics):
    cols = [
        "nodes",
        "edges",
        "initial_kappa",
        "delta_kappa",
        "variance_reduction",
        "lambda_improvement",
        "normalized_rate",
        "convergence_score",
    ]
    d = metrics[cols].copy()
    d["nodes"] = np.log10(d["nodes"].clip(lower=1))
    d["edges"] = np.log10(d["edges"].clip(lower=1))
    corr = d.corr(method="spearman")

    plt.figure(figsize=(9, 7))
    im = plt.imshow(corr.values, aspect="auto", vmin=-1, vmax=1)
    plt.colorbar(im, label="Spearman ρ")
    plt.xticks(range(len(corr.columns)), corr.columns, rotation=70)
    plt.yticks(range(len(corr.index)), corr.index)

    for i in range(len(corr.index)):
        for j in range(len(corr.columns)):
            val = corr.iloc[i, j]
            if np.isfinite(val):
                plt.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7)

    plt.title("Spearman correlation matrix")
    savefig(OUTPUT_DIR / "figures" / "fig_correlation_matrix.png")

    corr.to_csv(
        OUTPUT_DIR / "tables" / "08_spearman_correlation_matrix.csv"
    )


# ---------------------------------------------------------------------
# REPORT / AUDIT
# ---------------------------------------------------------------------

def create_dataset_statistics_table(metrics, summary):
    # Table 1: dataset statistics.
    cols = [
        "dataset",
        "domain",
        "nodes",
        "edges",
        "n_steps",
        "initial_kappa",
        "final_kappa",
        "delta_kappa",
        "initial_variance",
        "final_variance",
        "variance_reduction",
        "lambda_initial",
        "lambda_final",
        "lambda_improvement",
        "normalized_rate",
        "convergence_score",
        "variance_monotonic",
    ]
    table = metrics[cols].sort_values(["domain", "dataset"])
    table.to_csv(
        OUTPUT_DIR / "tables" / "Table_1_dataset_statistics.csv",
        index=False,
    )
    return table


def create_universal_table(metrics):
    """
    Table 2: universal Ricci-flow evolution.
    """
    rows = []

    for endpoint in [
        "delta_kappa",
        "variance_reduction",
        "lambda_improvement",
        "normalized_rate",
    ]:
        x = clean_numeric(metrics[endpoint]).dropna().values
        if len(x) == 0:
            continue

        stat, p, effect, n = paired_test(x)
        med, lo, hi = bootstrap_stat(
            x,
            func=np.median,
            seed=RANDOM_SEED + 60000 + len(rows),
        )
        mean = float(np.mean(x))

        rows.append(
            {
                "endpoint": endpoint,
                "n_datasets": len(x),
                "mean": mean,
                "median": float(np.median(x)),
                "median_bootstrap_ci95_low": lo,
                "median_bootstrap_ci95_high": hi,
                "wilcoxon_statistic": stat,
                "wilcoxon_p": p,
                "rank_biserial_effect": effect,
                "nonzero_n": n,
            }
        )

    out = pd.DataFrame(rows)
    out["wilcoxon_q_BH"] = benjamini_hochberg(out["wilcoxon_p"].values)

    out.to_csv(
        OUTPUT_DIR / "tables" / "Table_2_universal_ricci_flow_evolution.csv",
        index=False,
    )
    return out


def create_trajectory_completeness(raw, metrics):
    counts = raw.groupby("dataset")["step"].nunique()

    completeness = pd.DataFrame(
        {
            "n_datasets": [len(counts)],
            "datasets_with_10_steps": [int(np.sum(counts == 10))],
            "fraction_with_10_steps": [
                float(np.mean(counts == 10))
            ],
            "min_steps": [int(counts.min())],
            "max_steps": [int(counts.max())],
        }
    )
    completeness.to_csv(
        OUTPUT_DIR / "tables" / "09_trajectory_completeness.csv",
        index=False,
    )


def create_reproducibility_manifest(raw, summary, metrics):
    manifest = pd.DataFrame(
        [
            {
                "raw_results_path": str(RAW_RESULTS),
                "dataset_summary_path": str(DATASET_SUMMARY),
                "raw_records": len(raw),
                "raw_datasets": raw["dataset"].nunique(),
                "summary_datasets": summary["dataset"].nunique(),
                "metric_datasets": metrics["dataset"].nunique(),
                "domains": metrics["domain"].nunique(),
                "bootstrap_replicates": N_BOOT,
                "permutation_replicates": N_PERM,
                "random_seed": RANDOM_SEED,
                "independent_unit": "dataset",
                "trajectory_steps_as_independent_replicates": False,
                "tensorflow_used": False,
            }
        ]
    )
    manifest.to_csv(
        OUTPUT_DIR / "tables" / "00_reproducibility_manifest.csv",
        index=False,
    )


def create_quality_audit(raw, summary, metrics):
    rows = []

    for dataset, sub in raw.groupby("dataset"):
        steps = sorted(pd.unique(sub["step"].dropna()))
        row = {
            "dataset": dataset,
            "raw_records": len(sub),
            "unique_steps": len(steps),
            "has_step_0": 0 in steps,
            "has_step_9": 9 in steps,
            "exactly_10_steps": len(steps) == 10,
            "duplicate_step_records": int(len(sub) - len(steps)),
            "missing_curvature": int(sub["mean_curvature"].isna().sum()),
            "missing_variance": int(sub["curvature_variance"].isna().sum()),
            "missing_lambda1": int(sub["lambda1"].isna().sum()),
        }
        rows.append(row)

    audit = pd.DataFrame(rows)
    audit.to_csv(
        OUTPUT_DIR / "tables" / "10_raw_data_quality_audit.csv",
        index=False,
    )

    summary_audit = summary.copy()
    if "status" in summary_audit.columns:
        summary_audit["status"] = summary_audit["status"].astype(str)

    summary_audit.to_csv(
        OUTPUT_DIR / "tables" / "11_dataset_summary_audit_copy.csv",
        index=False,
    )


def write_readme(universal, domain_df, size_df, null_df):
    lines = [
        "RICCI FLOW SIX-ANALYSIS FINAL REPORT",
        "=" * 72,
        "",
        f"Raw results: {RAW_RESULTS}",
        f"Dataset summary: {DATASET_SUMMARY}",
        f"Output directory: {OUTPUT_DIR}",
        "",
        "INDEPENDENT UNIT",
        "Dataset is the independent experimental unit.",
        "The 10 flow steps are trajectory observations, not independent replicates.",
        "",
        "STATISTICAL FRAMEWORK",
        "- Wilcoxon signed-rank tests at dataset level",
        "- rank-biserial effect size",
        "- 10,000-replicate bootstrap confidence intervals",
        "- sign-flip permutation null models",
        "- leave-one-dataset-out robustness",
        "- cross-domain validation using dataset_summary.csv domain",
        "- optional domain random-intercept mixed-effects models",
        "- Spearman network-size dependence with bootstrap CI",
        "- repaired spectral analysis excluding non-positive initial λ1 from",
        "  percentage-improvement denominators",
        "",
        "INTERPRETATION RULE",
        "A positive mean alone is not treated as evidence of universality.",
        "Median, confidence interval, permutation/Wilcoxon evidence,",
        "cross-domain consistency, and leave-one-out robustness are examined.",
        "",
        "UNIVERSAL ENDPOINTS",
    ]

    if universal is not None and len(universal):
        for _, r in universal.iterrows():
            lines.append(
                f"{r['endpoint']}: n={int(r['n'])}, "
                f"mean={r['mean']:.8g}, median={r['median']:.8g}, "
                f"Wilcoxon p={r['wilcoxon_p']:.8g}, "
                f"BH-q={r['wilcoxon_q_BH']:.8g}"
            )

    lines.extend(["", "OUTPUTS"])
    for p in sorted(OUTPUT_DIR.rglob("*")):
        if p.is_file():
            lines.append(str(p.relative_to(OUTPUT_DIR)))

    (OUTPUT_DIR / "FINAL_REPORT_README.txt").write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


# ---------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------

def main():
    ensure_dirs()

    log("=" * 72)
    log("RICCI FLOW FINAL SIX-ANALYSIS")
    log("=" * 72)

    raw, summary = load_inputs()
    metrics = build_dataset_metrics(raw, summary)

    if metrics.empty:
        raise RuntimeError("No valid dataset-level trajectories found.")

    log("")
    log(
        f"Valid metric datasets: {metrics['dataset'].nunique()} | "
        f"Domains: {metrics['domain'].nunique()}"
    )

    # Save reconstructed endpoint data first.
    metrics.to_csv(
        OUTPUT_DIR / "tables" / "00_dataset_level_metrics_reconstructed.csv",
        index=False,
    )

    create_dataset_statistics_table(metrics, summary)
    create_universal_table(metrics)
    create_trajectory_completeness(raw, metrics)
    create_reproducibility_manifest(raw, summary, metrics)
    create_quality_audit(raw, summary, metrics)

    log("")
    log("Analysis 1/6: universal endpoint statistics")
    universal = analyze_universal(metrics)

    log("Analysis 2/6: leave-one-dataset-out robustness")
    robustness = analyze_robustness(metrics)

    log("Analysis 3/6: cross-domain validation")
    domain_df, heterogeneity = analyze_cross_domain(metrics)

    log("Analysis 4/6: domain-aware mixed-effects analysis")
    try_mixed_effects(metrics)

    log("Analysis 5/6: network-size dependence")
    size_df = analyze_size_dependence(metrics)

    log("Analysis 6/6: sign-flip permutation null models")
    null_df = analyze_null_models(metrics)

    log("Spectral repair / audit")
    spectral_df, spectral_invalid = analyze_spectral(metrics)

    # Figures.
    log("Generating figures...")
    plot_endpoint_distribution(
        metrics,
        "delta_kappa",
        "Distribution of dataset-level curvature change",
        "fig_delta_kappa_distribution.png",
        "Final κ − Initial κ",
    )
    plot_endpoint_distribution(
        metrics,
        "variance_reduction",
        "Distribution of curvature variance reduction",
        "fig_variance_reduction_distribution.png",
        "Variance reduction",
    )
    plot_endpoint_distribution(
        metrics,
        "lambda_improvement",
        "Distribution of spectral-gap improvement (valid λ1 only)",
        "fig_lambda_improvement_distribution.png",
        "Relative λ1 improvement",
    )

    plot_trajectory_envelope(raw)
    plot_variance_trajectories(raw)
    plot_spectral_trajectories(raw)

    plot_domain_forest(
        domain_df,
        "variance_reduction",
        "Cross-domain validation: curvature variance reduction",
        "fig_cross_domain_variance_reduction_forest.png",
    )
    plot_domain_forest(
        domain_df,
        "delta_kappa",
        "Cross-domain validation: curvature change",
        "fig_cross_domain_delta_kappa_forest.png",
    )

    plot_size_relationship(
        metrics,
        "variance_reduction",
        "fig_size_dependence_variance_reduction.png",
        "Network size vs curvature variance reduction",
    )
    plot_size_relationship(
        metrics,
        "delta_kappa",
        "fig_size_dependence_delta_kappa.png",
        "Network size vs curvature change",
    )
    plot_size_relationship(
        metrics,
        "normalized_rate",
        "fig_size_dependence_normalized_rate.png",
        "Network size vs normalized curvature evolution rate",
    )

    for endpoint in [
        "delta_kappa",
        "variance_reduction",
        "lambda_improvement",
        "normalized_rate",
    ]:
        plot_robustness(
            metrics,
            endpoint,
            f"fig_loo_{endpoint}.png",
        )

    plot_correlation_matrix(metrics)

    # A compact text report with conclusions driven by actual computed values.
    with open(
        OUTPUT_DIR / "FINAL_STATISTICAL_REPORT.txt",
        "w",
        encoding="utf-8",
    ) as f:
        f.write("RICCI FLOW FINAL STATISTICAL REPORT\n")
        f.write("=" * 72 + "\n\n")
        f.write(
            f"Dataset-level n = {metrics['dataset'].nunique()}\n"
        )
        f.write(
            f"Domains = {metrics['domain'].nunique()}\n"
        )
        f.write(
            f"Raw records = {len(raw)}\n\n"
        )

        f.write("PRIMARY ENDPOINTS\n")
        f.write("-" * 72 + "\n")
        for _, r in universal.iterrows():
            f.write(
                f"{r['endpoint']}: n={int(r['n'])}; "
                f"mean={r['mean']:.8g}; median={r['median']:.8g}; "
                f"bootstrap median 95% CI="
                f"[{r['bootstrap_median_ci95_low']:.8g}, "
                f"{r['bootstrap_median_ci95_high']:.8g}]; "
                f"Wilcoxon p={r['wilcoxon_p']:.8g}; "
                f"BH q={r['wilcoxon_q_BH']:.8g}; "
                f"rank-biserial={r['rank_biserial']:.8g}\n"
            )

        f.write("\nCROSS-DOMAIN HETEROGENEITY\n")
        f.write("-" * 72 + "\n")
        f.write(heterogeneity.to_string(index=False))
        f.write("\n\nNETWORK SIZE DEPENDENCE\n")
        f.write("-" * 72 + "\n")
        f.write(size_df.to_string(index=False))
        f.write("\n\nNULL MODEL\n")
        f.write("-" * 72 + "\n")
        f.write(null_df.to_string(index=False))
        f.write("\n\nSPECTRAL AUDIT\n")
        f.write("-" * 72 + "\n")
        f.write(spectral_df.to_string(index=False))
        f.write("\n\nNON-IDENTIFIABLE SPECTRAL CASES\n")
        f.write("-" * 72 + "\n")
        f.write(spectral_invalid.to_string(index=False))

    write_readme(universal, domain_df, size_df, null_df)

    log("")
    log("=" * 72)
    log("ALL SIX ANALYSES COMPLETED")
    log(f"Output: {OUTPUT_DIR}")
    log("=" * 72)


if __name__ == "__main__":
    main()
