# -*- coding: utf-8 -*-
"""
ricci_flow_6analysis.py

Statistical analysis of universal Ricci-flow evolution across complex networks.

INPUT:
    raw_results.csv
    dataset_summary.csv

This program:
    1. Reads ONLY the two CSV files.
    2. Does NOT read original graph files.
    3. Does NOT run Ricci-flow calculations.
    4. Does NOT use TensorFlow.
    5. Does NOT use load_graph(), analyze_graph(), or run().
    6. Treats dataset as the independent experimental unit.
    7. Uses step=0...9 trajectories for evolution/convergence analysis.
    8. Uses dataset_summary.domain for cross-domain validation.

ANALYSIS:
    A. Universal Ricci Flow Evolution
    B. Statistical Significance
    C. Correlation / Mechanism
    D. Convergence Analysis
    E. Robustness Analysis
    F. Cross-domain Validation

OUTPUT:
    analysis_6results/
        tables/
        statistics/
        figures/
        reports/

Required packages:
    numpy
    pandas
    scipy
    matplotlib

No TensorFlow.
"""

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


# ============================================================
# 0. CONFIGURATION
# ============================================================

RAW_RESULTS_PATH = (
    r"D:\2026codes\2026821100data\results_correct_fast\raw_results.csv"
)

DATASET_SUMMARY_PATH = (
    r"D:\2026codes\2026821100data\results_correct_fast\dataset_summary.csv"
)

OUTPUT_DIR = (
    r"D:\2026codes\2026821100data\analysis_6results"
)

STEPS_EXPECTED = 10
BOOTSTRAP_ITERATIONS = 5000
RANDOM_SEED = 20260829

np.random.seed(RANDOM_SEED)


# ============================================================
# 1. DIRECTORY SETUP
# ============================================================

OUTPUT_DIR = Path(OUTPUT_DIR)

TABLE_DIR = OUTPUT_DIR / "tables"
STAT_DIR = OUTPUT_DIR / "statistics"
FIG_DIR = OUTPUT_DIR / "figures"
REPORT_DIR = OUTPUT_DIR / "reports"

for d in [OUTPUT_DIR, TABLE_DIR, STAT_DIR, FIG_DIR, REPORT_DIR]:
    d.mkdir(parents=True, exist_ok=True)


# ============================================================
# 2. REQUIRED COLUMNS
# ============================================================

RAW_REQUIRED = [
    "dataset",
    "step",
    "mean_curvature",
    "curvature_variance",
    "lambda1",
    "nodes",
    "edges",
    "time",
    "distance_time",
    "sinkhorn_time",
    "lambda_time",
]

SUMMARY_REQUIRED = [
    "domain",
    "dataset",
    "file",
    "status",
    "message",
    "nodes",
    "edges",
    "components",
    "weight_mode",
    "weight_reason",
    "readme",
    "readme_direction",
    "duplicates",
    "invalid_lines",
    "transport_cells_per_step",
    "shortest_path_work",
    "slope",
    "intercept",
]


# ============================================================
# 3. UTILITY FUNCTIONS
# ============================================================

def safe_numeric(df, columns):
    """Convert selected columns to numeric."""
    for c in columns:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def safe_divide(a, b):
    """Element-wise safe division."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    result = np.full_like(a, np.nan, dtype=float)

    mask = np.isfinite(a) & np.isfinite(b) & (np.abs(b) > 1e-15)

    result[mask] = a[mask] / b[mask]

    return result


def benjamini_hochberg(pvalues):
    """
    Benjamini-Hochberg FDR correction.

    Returns adjusted p-values.
    """
    pvalues = np.asarray(pvalues, dtype=float)

    result = np.full(len(pvalues), np.nan)

    valid = np.isfinite(pvalues)

    if valid.sum() == 0:
        return result

    pv = pvalues[valid]
    n = len(pv)

    order = np.argsort(pv)
    ranked = pv[order]

    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)

    restored = np.empty(n)
    restored[order] = adjusted

    result[valid] = restored

    return result


def bootstrap_mean_ci(values, iterations=BOOTSTRAP_ITERATIONS,
                       confidence=0.95, seed=RANDOM_SEED):
    """
    Bootstrap confidence interval for the mean.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    if len(values) < 2:
        return np.nan, np.nan

    rng = np.random.default_rng(seed)

    n = len(values)

    samples = rng.choice(
        values,
        size=(iterations, n),
        replace=True
    )

    means = np.mean(samples, axis=1)

    alpha = 1.0 - confidence

    low = np.quantile(means, alpha / 2)
    high = np.quantile(means, 1 - alpha / 2)

    return float(low), float(high)


def bootstrap_median_ci(values, iterations=BOOTSTRAP_ITERATIONS,
                         confidence=0.95, seed=RANDOM_SEED):
    """
    Bootstrap confidence interval for the median.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    if len(values) < 2:
        return np.nan, np.nan

    rng = np.random.default_rng(seed)

    n = len(values)

    samples = rng.choice(
        values,
        size=(iterations, n),
        replace=True
    )

    medians = np.median(samples, axis=1)

    alpha = 1.0 - confidence

    low = np.quantile(medians, alpha / 2)
    high = np.quantile(medians, 1 - alpha / 2)

    return float(low), float(high)


def rank_biserial_from_wilcoxon(values):
    """
    Rank-biserial effect size for paired Wilcoxon-type analysis.

    For a one-sample signed-rank test:
        r_rb = (W_plus - W_minus) / (W_plus + W_minus)
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    values = values[values != 0]

    if len(values) == 0:
        return np.nan

    abs_values = np.abs(values)

    ranks = stats.rankdata(abs_values)

    w_plus = np.sum(ranks[values > 0])
    w_minus = np.sum(ranks[values < 0])

    denominator = w_plus + w_minus

    if denominator == 0:
        return np.nan

    return float((w_plus - w_minus) / denominator)


def cohens_d_one_sample(values):
    """
    One-sample Cohen's d.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    if len(values) < 2:
        return np.nan

    sd = np.std(values, ddof=1)

    if sd <= 1e-15:
        return np.nan

    return float(np.mean(values) / sd)


def safe_pearson(x, y):
    mask = np.isfinite(x) & np.isfinite(y)

    if mask.sum() < 3:
        return np.nan, np.nan, int(mask.sum())

    r, p = stats.pearsonr(x[mask], y[mask])

    return float(r), float(p), int(mask.sum())


def safe_spearman(x, y):
    mask = np.isfinite(x) & np.isfinite(y)

    if mask.sum() < 3:
        return np.nan, np.nan, int(mask.sum())

    r, p = stats.spearmanr(x[mask], y[mask])

    return float(r), float(p), int(mask.sum())


def regression_slope(x, y):
    mask = np.isfinite(x) & np.isfinite(y)

    if mask.sum() < 3:
        return np.nan, np.nan, np.nan

    result = stats.linregress(x[mask], y[mask])

    return (
        float(result.slope),
        float(result.intercept),
        float(result.rvalue),
    )


def fmt(x, digits=6):
    if x is None:
        return "NA"

    try:
        if not np.isfinite(float(x)):
            return "NA"
        return f"{float(x):.{digits}g}"
    except Exception:
        return "NA"


# ============================================================
# 4. LOAD DATA
# ============================================================

def load_input_data():
    print("=" * 72)
    print("RICCI FLOW — SIX-ANALYSIS STATISTICAL FRAMEWORK")
    print("=" * 72)

    print("\nInput raw results:")
    print(RAW_RESULTS_PATH)

    print("\nInput dataset summary:")
    print(DATASET_SUMMARY_PATH)

    if not os.path.exists(RAW_RESULTS_PATH):
        raise FileNotFoundError(
            f"raw_results.csv not found:\n{RAW_RESULTS_PATH}"
        )

    if not os.path.exists(DATASET_SUMMARY_PATH):
        raise FileNotFoundError(
            f"dataset_summary.csv not found:\n{DATASET_SUMMARY_PATH}"
        )

    raw = pd.read_csv(RAW_RESULTS_PATH)
    summary = pd.read_csv(DATASET_SUMMARY_PATH)

    print("\nRaw columns:")
    print(list(raw.columns))

    print("\nSummary columns:")
    print(list(summary.columns))

    missing_raw = [
        c for c in RAW_REQUIRED
        if c not in raw.columns
    ]

    missing_summary = [
        c for c in SUMMARY_REQUIRED
        if c not in summary.columns
    ]

    if missing_raw:
        raise ValueError(
            "raw_results.csv is missing required columns:\n"
            + "\n".join(missing_raw)
        )

    if missing_summary:
        raise ValueError(
            "dataset_summary.csv is missing required columns:\n"
            + "\n".join(missing_summary)
        )

    numeric_raw = [
        "step",
        "mean_curvature",
        "curvature_variance",
        "lambda1",
        "nodes",
        "edges",
        "time",
        "distance_time",
        "sinkhorn_time",
        "lambda_time",
    ]

    raw = safe_numeric(raw, numeric_raw)

    numeric_summary = [
        "nodes",
        "edges",
        "components",
        "duplicates",
        "invalid_lines",
        "transport_cells_per_step",
        "shortest_path_work",
        "slope",
        "intercept",
    ]

    summary = safe_numeric(summary, numeric_summary)

    raw["dataset"] = raw["dataset"].astype(str).str.strip()
    summary["dataset"] = summary["dataset"].astype(str).str.strip()

    summary["domain"] = (
        summary["domain"]
        .astype(str)
        .str.strip()
    )

    print(f"\nRaw records: {len(raw)}")
    print(f"Summary records: {len(summary)}")
    print(
        f"Raw datasets: {raw['dataset'].nunique()}"
    )
    print(
        f"Summary datasets: {summary['dataset'].nunique()}"
    )

    return raw, summary


# ============================================================
# 5. BUILD DATASET-LEVEL METRICS
# ============================================================

def calculate_dataset_metrics(raw, summary):

    print("\n" + "=" * 72)
    print("BUILDING DATASET-LEVEL METRICS")
    print("=" * 72)

    rows = []

    valid_datasets = []

    for dataset, sub in raw.groupby("dataset", sort=True):

        sub = sub.sort_values("step").copy()

        if len(sub) < 2:
            continue

        sub = sub.dropna(
            subset=[
                "step",
                "mean_curvature",
                "curvature_variance",
            ]
        )

        if len(sub) < 2:
            continue

        # ----------------------------------------------------
        # Initial and final records
        # ----------------------------------------------------

        first = sub.iloc[0]
        last = sub.iloc[-1]

        k0 = float(first["mean_curvature"])
        kT = float(last["mean_curvature"])

        v0 = float(first["curvature_variance"])
        vT = float(last["curvature_variance"])

        l0 = float(first["lambda1"])
        lT = float(last["lambda1"])

        step0 = float(first["step"])
        stepT = float(last["step"])

        T = stepT - step0

        # ----------------------------------------------------
        # Main metrics
        # ----------------------------------------------------

        delta_kappa = kT - k0

        if np.isfinite(v0) and abs(v0) > 1e-15:
            variance_reduction = (v0 - vT) / abs(v0)
        else:
            variance_reduction = np.nan

        if np.isfinite(l0) and abs(l0) > 1e-15:
            lambda_improvement = (lT - l0) / abs(l0)
        else:
            lambda_improvement = np.nan

        if (
            np.isfinite(k0)
            and abs(k0) > 1e-15
            and T > 0
        ):
            normalized_rate = delta_kappa / (T * abs(k0))
        else:
            normalized_rate = np.nan

        # ----------------------------------------------------
        # Trajectory slope
        # ----------------------------------------------------

        slope, intercept, rvalue = regression_slope(
            sub["step"].to_numpy(dtype=float),
            sub["mean_curvature"].to_numpy(dtype=float)
        )

        # ----------------------------------------------------
        # Variance slope
        # ----------------------------------------------------

        variance_slope, _, variance_r = regression_slope(
            sub["step"].to_numpy(dtype=float),
            sub["curvature_variance"].to_numpy(dtype=float)
        )

        # ----------------------------------------------------
        # Lambda slope
        # ----------------------------------------------------

        lambda_slope, _, lambda_r = regression_slope(
            sub["step"].to_numpy(dtype=float),
            sub["lambda1"].to_numpy(dtype=float)
        )

        # ----------------------------------------------------
        # Monotonicity
        # ----------------------------------------------------

        k_values = sub["mean_curvature"].to_numpy(dtype=float)

        v_values = sub[
            "curvature_variance"
        ].to_numpy(dtype=float)

        lambda_values = sub[
            "lambda1"
        ].to_numpy(dtype=float)

        k_diff = np.diff(k_values)

        v_diff = np.diff(v_values)

        lambda_diff = np.diff(lambda_values)

        if len(k_diff) > 0:
            k_monotonicity = (
                np.sum(
                    np.sign(k_diff)
                    ==
                    np.sign(k_values[-1] - k_values[0])
                )
                / len(k_diff)
            )
        else:
            k_monotonicity = np.nan

        if len(v_diff) > 0:
            v_monotonicity = (
                np.sum(v_diff <= 0)
                / len(v_diff)
            )
        else:
            v_monotonicity = np.nan

        if len(lambda_diff) > 0:
            lambda_monotonicity = (
                np.sum(lambda_diff >= 0)
                / len(lambda_diff)
            )
        else:
            lambda_monotonicity = np.nan

        # ----------------------------------------------------
        # Convergence score
        # ----------------------------------------------------

        # A simple composite score:
        #   1. curvature trajectory monotonicity
        #   2. variance monotonicity
        #   3. lambda monotonicity
        #
        # This is deliberately descriptive rather than
        # claiming a formal mathematical convergence proof.

        score_parts = [
            k_monotonicity,
            v_monotonicity,
            lambda_monotonicity,
        ]

        score_parts = [
            x for x in score_parts
            if np.isfinite(x)
        ]

        if score_parts:
            convergence_score = float(
                np.mean(score_parts)
            )
        else:
            convergence_score = np.nan

        # ----------------------------------------------------
        # Metadata from summary
        # ----------------------------------------------------

        sm = summary[
            summary["dataset"] == dataset
        ]

        if len(sm) > 0:
            sm = sm.iloc[0]

            domain = sm["domain"]

            summary_nodes = sm["nodes"]
            summary_edges = sm["edges"]
            components = sm["components"]
            weight_mode = sm["weight_mode"]
            summary_slope = sm["slope"]
            summary_intercept = sm["intercept"]

        else:
            domain = "Unknown"
            summary_nodes = np.nan
            summary_edges = np.nan
            components = np.nan
            weight_mode = "Unknown"
            summary_slope = np.nan
            summary_intercept = np.nan

        # Prefer raw node/edge values when available
        nodes = first["nodes"]
        edges = first["edges"]

        if not np.isfinite(nodes):
            nodes = summary_nodes

        if not np.isfinite(edges):
            edges = summary_edges

        # ----------------------------------------------------
        # Network density
        # ----------------------------------------------------

        if (
            np.isfinite(nodes)
            and np.isfinite(edges)
            and nodes > 1
        ):
            density = (
                2.0 * edges
                /
                (nodes * (nodes - 1.0))
            )
        else:
            density = np.nan

        rows.append({
            "dataset": dataset,
            "domain": domain,

            "nodes": nodes,
            "edges": edges,
            "components": components,
            "density": density,
            "weight_mode": weight_mode,

            "n_steps": len(sub),
            "initial_step": step0,
            "final_step": stepT,

            "initial_kappa": k0,
            "final_kappa": kT,
            "delta_kappa": delta_kappa,

            "initial_variance": v0,
            "final_variance": vT,
            "variance_reduction": variance_reduction,

            "lambda_initial": l0,
            "lambda_final": lT,
            "lambda_improvement": lambda_improvement,

            "trajectory_slope": slope,
            "trajectory_intercept": intercept,
            "trajectory_r": rvalue,

            "summary_slope": summary_slope,
            "summary_intercept": summary_intercept,

            "variance_slope": variance_slope,
            "variance_r": variance_r,

            "lambda_slope": lambda_slope,
            "lambda_r": lambda_r,

            "normalized_rate": normalized_rate,

            "kappa_monotonicity": k_monotonicity,
            "variance_monotonicity": v_monotonicity,
            "lambda_monotonicity": lambda_monotonicity,

            "convergence_score": convergence_score,
        })

        valid_datasets.append(dataset)

    metrics = pd.DataFrame(rows)

    metrics.to_csv(
        TABLE_DIR / "dataset_level_metrics.csv",
        index=False,
        encoding="utf-8-sig"
    )

    print(
        f"\nValid dataset-level records: {len(metrics)}"
    )

    print(
        f"Domains represented: "
        f"{metrics['domain'].nunique()}"
    )

    return metrics


# ============================================================
# 6. ANALYSIS A
# ============================================================

def analysis_A_universal_evolution(metrics):

    print("\n" + "=" * 72)
    print("ANALYSIS A — UNIVERSAL RICCI FLOW EVOLUTION")
    print("=" * 72)

    cols = [
        "dataset",
        "domain",
        "nodes",
        "edges",
        "density",
        "initial_kappa",
        "final_kappa",
        "delta_kappa",
        "initial_variance",
        "final_variance",
        "variance_reduction",
        "lambda_initial",
        "lambda_final",
        "lambda_improvement",
        "trajectory_slope",
        "normalized_rate",
        "convergence_score",
    ]

    table = metrics[
        [c for c in cols if c in metrics.columns]
    ].copy()

    table.to_csv(
        TABLE_DIR /
        "Table_2_Universal_Ricci_Flow_Evolution.csv",
        index=False,
        encoding="utf-8-sig"
    )

    # --------------------------------------------------------
    # Aggregate summary
    # --------------------------------------------------------

    indicators = [
        "delta_kappa",
        "variance_reduction",
        "lambda_improvement",
        "normalized_rate",
        "convergence_score",
    ]

    rows = []

    for variable in indicators:

        values = pd.to_numeric(
            metrics[variable],
            errors="coerce"
        ).dropna().to_numpy()

        if len(values) == 0:
            continue

        ci_low, ci_high = bootstrap_mean_ci(values)

        median_low, median_high = bootstrap_median_ci(values)

        rows.append({
            "indicator": variable,
            "n": len(values),
            "mean": np.mean(values),
            "median": np.median(values),
            "std": np.std(values, ddof=1)
            if len(values) > 1 else np.nan,
            "min": np.min(values),
            "max": np.max(values),
            "mean_ci95_low": ci_low,
            "mean_ci95_high": ci_high,
            "median_ci95_low": median_low,
            "median_ci95_high": median_high,
        })

    aggregate = pd.DataFrame(rows)

    aggregate.to_csv(
        STAT_DIR /
        "A_universal_evolution_summary.csv",
        index=False,
        encoding="utf-8-sig"
    )

    # --------------------------------------------------------
    # Figure A1 — delta kappa distribution
    # --------------------------------------------------------

    values = metrics["delta_kappa"].dropna()

    if len(values) > 0:

        plt.figure(figsize=(9, 6))

        plt.hist(values, bins=30)

        plt.axvline(
            0,
            linestyle="--",
            linewidth=1.2
        )

        plt.xlabel("Δκ")
        plt.ylabel("Number of datasets")
        plt.title(
            "Distribution of Ricci Curvature Evolution"
        )

        plt.tight_layout()

        plt.savefig(
            FIG_DIR /
            "Fig_A1_delta_kappa_distribution.png",
            dpi=300
        )

        plt.close()


# ============================================================
# 7. ANALYSIS B
# ============================================================

def analysis_B_statistical_significance(metrics):

    print("\n" + "=" * 72)
    print("ANALYSIS B — STATISTICAL SIGNIFICANCE")
    print("=" * 72)

    analyses = []

    test_variables = {
        "delta_kappa": "Δκ",
        "variance_reduction": "Variance reduction",
        "lambda_improvement": "λ1 improvement",
        "normalized_rate": "Normalized evolution rate",
    }

    for variable, label in test_variables.items():

        values = pd.to_numeric(
            metrics[variable],
            errors="coerce"
        ).dropna().to_numpy()

        n = len(values)

        if n == 0:
            continue

        mean_value = np.mean(values)

        median_value = np.median(values)

        std_value = (
            np.std(values, ddof=1)
            if n > 1
            else np.nan
        )

        sem = (
            std_value / math.sqrt(n)
            if n > 1
            else np.nan
        )

        # -----------------------------------------------
        # One-sample t-test
        # -----------------------------------------------

        if n >= 2 and np.isfinite(std_value):

            t_stat, t_p = stats.ttest_1samp(
                values,
                0
            )

        else:
            t_stat = np.nan
            t_p = np.nan

        # -----------------------------------------------
        # Wilcoxon signed-rank test
        # -----------------------------------------------

        nonzero = values[values != 0]

        if len(nonzero) >= 2:

            try:

                w_result = stats.wilcoxon(
                    nonzero,
                    alternative="two-sided",
                    method="auto"
                )

                w_stat = float(w_result.statistic)
                w_p = float(w_result.pvalue)

            except Exception:

                w_stat = np.nan
                w_p = np.nan

        else:

            w_stat = np.nan
            w_p = np.nan

        # -----------------------------------------------
        # Effect sizes
        # -----------------------------------------------

        cohens_d = cohens_d_one_sample(values)

        rank_biserial = (
            rank_biserial_from_wilcoxon(values)
        )

        # -----------------------------------------------
        # Bootstrap CI
        # -----------------------------------------------

        ci_low, ci_high = bootstrap_mean_ci(values)

        median_low, median_high = bootstrap_median_ci(
            values
        )

        # -----------------------------------------------
        # Direction
        # -----------------------------------------------

        positive_fraction = (
            np.mean(values > 0)
        )

        negative_fraction = (
            np.mean(values < 0)
        )

        analyses.append({
            "indicator": variable,
            "label": label,
            "n": n,

            "mean": mean_value,
            "median": median_value,
            "std": std_value,
            "sem": sem,

            "mean_ci95_low": ci_low,
            "mean_ci95_high": ci_high,

            "median_ci95_low": median_low,
            "median_ci95_high": median_high,

            "t_statistic": t_stat,
            "t_pvalue": t_p,

            "wilcoxon_statistic": w_stat,
            "wilcoxon_pvalue": w_p,

            "cohens_d": cohens_d,
            "rank_biserial": rank_biserial,

            "positive_fraction": positive_fraction,
            "negative_fraction": negative_fraction,
        })

    results = pd.DataFrame(analyses)

    # FDR for t-test and Wilcoxon separately

    if len(results) > 0:

        results["t_pvalue_fdr"] = benjamini_hochberg(
            results["t_pvalue"].to_numpy()
        )

        results["wilcoxon_pvalue_fdr"] = benjamini_hochberg(
            results["wilcoxon_pvalue"].to_numpy()
        )

    results.to_csv(
        STAT_DIR /
        "B_statistical_significance.csv",
        index=False,
        encoding="utf-8-sig"
    )

    # --------------------------------------------------------
    # Dataset-level statistical table
    # --------------------------------------------------------

    pair_table = metrics[
        [
            "dataset",
            "domain",
            "delta_kappa",
            "variance_reduction",
            "lambda_improvement",
            "normalized_rate",
        ]
    ].copy()

    pair_table.to_csv(
        TABLE_DIR /
        "Table_3_Statistical_Effects_By_Dataset.csv",
        index=False,
        encoding="utf-8-sig"
    )


# ============================================================
# 8. ANALYSIS C
# ============================================================

def analysis_C_correlation(metrics):

    print("\n" + "=" * 72)
    print("ANALYSIS C — CORRELATION / MECHANISM")
    print("=" * 72)

    outcomes = [
        "delta_kappa",
        "variance_reduction",
        "lambda_improvement",
        "normalized_rate",
    ]

    predictors = [
        "nodes",
        "edges",
        "density",
        "components",
        "initial_kappa",
        "initial_variance",
        "lambda_initial",
    ]

    rows = []

    for predictor in predictors:

        if predictor not in metrics.columns:
            continue

        for outcome in outcomes:

            if outcome not in metrics.columns:
                continue

            x = pd.to_numeric(
                metrics[predictor],
                errors="coerce"
            ).to_numpy(dtype=float)

            y = pd.to_numeric(
                metrics[outcome],
                errors="coerce"
            ).to_numpy(dtype=float)

            pr, pp, pn = safe_pearson(x, y)

            sr, sp, sn = safe_spearman(x, y)

            rows.append({
                "predictor": predictor,
                "outcome": outcome,

                "pearson_r": pr,
                "pearson_pvalue": pp,
                "pearson_n": pn,

                "spearman_rho": sr,
                "spearman_pvalue": sp,
                "spearman_n": sn,
            })

    result = pd.DataFrame(rows)

    if len(result) > 0:

        result["pearson_pvalue_fdr"] = benjamini_hochberg(
            result["pearson_pvalue"].to_numpy()
        )

        result["spearman_pvalue_fdr"] = benjamini_hochberg(
            result["spearman_pvalue"].to_numpy()
        )

    result.to_csv(
        STAT_DIR /
        "C_correlation_analysis.csv",
        index=False,
        encoding="utf-8-sig"
    )

    # --------------------------------------------------------
    # Correlation heatmap
    # --------------------------------------------------------

    matrix_variables = [
        "initial_kappa",
        "initial_variance",
        "delta_kappa",
        "variance_reduction",
        "lambda_improvement",
        "normalized_rate",
        "nodes",
        "edges",
        "density",
    ]

    available = [
        x for x in matrix_variables
        if x in metrics.columns
    ]

    if len(available) >= 2:

        corr = metrics[available].corr(
            method="spearman"
        )

        plt.figure(
            figsize=(11, 9)
        )

        plt.imshow(
            corr,
            aspect="auto"
        )

        plt.colorbar(
            label="Spearman correlation"
        )

        plt.xticks(
            range(len(available)),
            available,
            rotation=70,
            ha="right"
        )

        plt.yticks(
            range(len(available)),
            available
        )

        plt.title(
            "Spearman Correlation Matrix"
        )

        plt.tight_layout()

        plt.savefig(
            FIG_DIR /
            "Fig_C1_correlation_matrix.png",
            dpi=300
        )

        plt.close()


# ============================================================
# 9. ANALYSIS D
# ============================================================

def analysis_D_convergence(raw, metrics):

    print("\n" + "=" * 72)
    print("ANALYSIS D — CONVERGENCE ANALYSIS")
    print("=" * 72)

    # --------------------------------------------------------
    # Aggregate trajectory
    # --------------------------------------------------------

    trajectory_rows = []

    for step, sub in raw.groupby("step", sort=True):

        for variable in [
            "mean_curvature",
            "curvature_variance",
            "lambda1",
        ]:

            values = pd.to_numeric(
                sub[variable],
                errors="coerce"
            ).dropna().to_numpy()

            if len(values) == 0:
                continue

            ci_low, ci_high = bootstrap_mean_ci(
                values,
                iterations=3000
            )

            trajectory_rows.append({
                "step": step,
                "variable": variable,
                "n": len(values),
                "mean": np.mean(values),
                "median": np.median(values),
                "std": (
                    np.std(values, ddof=1)
                    if len(values) > 1
                    else np.nan
                ),
                "ci95_low": ci_low,
                "ci95_high": ci_high,
            })

    trajectory = pd.DataFrame(
        trajectory_rows
    )

    trajectory.to_csv(
        STAT_DIR /
        "D_aggregate_trajectory.csv",
        index=False,
        encoding="utf-8-sig"
    )

    # --------------------------------------------------------
    # Figure D1 — curvature
    # --------------------------------------------------------

    make_trajectory_plot(
        trajectory,
        "mean_curvature",
        "Mean Ricci Curvature",
        "Mean Ricci Curvature Evolution",
        FIG_DIR /
        "Fig_D1_curvature_evolution.png"
    )

    # --------------------------------------------------------
    # Figure D2 — variance
    # --------------------------------------------------------

    make_trajectory_plot(
        trajectory,
        "curvature_variance",
        "Curvature variance",
        "Curvature Variance Evolution",
        FIG_DIR /
        "Fig_D2_variance_evolution.png"
    )

    # --------------------------------------------------------
    # Figure D3 — lambda
    # --------------------------------------------------------

    make_trajectory_plot(
        trajectory,
        "lambda1",
        "λ1",
        "Spectral Evolution",
        FIG_DIR /
        "Fig_D3_lambda_evolution.png"
    )

    # --------------------------------------------------------
    # Convergence classification
    # --------------------------------------------------------

    convergence = metrics[
        [
            "dataset",
            "domain",
            "kappa_monotonicity",
            "variance_monotonicity",
            "lambda_monotonicity",
            "convergence_score",
            "trajectory_r",
        ]
    ].copy()

    convergence["high_convergence"] = (
        convergence["convergence_score"] >= 0.75
    )

    convergence.to_csv(
        TABLE_DIR /
        "Table_4_Convergence_Analysis.csv",
        index=False,
        encoding="utf-8-sig"
    )

    # --------------------------------------------------------
    # Overall convergence statistics
    # --------------------------------------------------------

    score = convergence[
        "convergence_score"
    ].dropna()

    convergence_summary = pd.DataFrame([{
        "n": len(score),
        "mean_convergence_score": (
            score.mean()
            if len(score) else np.nan
        ),
        "median_convergence_score": (
            score.median()
            if len(score) else np.nan
        ),
        "fraction_score_ge_0_75": (
            np.mean(score >= 0.75)
            if len(score) else np.nan
        ),
        "fraction_score_ge_0_90": (
            np.mean(score >= 0.90)
            if len(score) else np.nan
        ),
    }])

    convergence_summary.to_csv(
        STAT_DIR /
        "D_convergence_summary.csv",
        index=False,
        encoding="utf-8-sig"
    )


def make_trajectory_plot(
    trajectory,
    variable,
    ylabel,
    title,
    output_path
):

    sub = trajectory[
        trajectory["variable"] == variable
    ].sort_values("step")

    if len(sub) == 0:
        return

    x = sub["step"].to_numpy(dtype=float)

    mean = sub["mean"].to_numpy(dtype=float)

    low = sub["ci95_low"].to_numpy(dtype=float)

    high = sub["ci95_high"].to_numpy(dtype=float)

    plt.figure(figsize=(9, 6))

    plt.plot(
        x,
        mean,
        marker="o",
        linewidth=2
    )

    plt.fill_between(
        x,
        low,
        high,
        alpha=0.2
    )

    plt.xlabel("Ricci flow step")
    plt.ylabel(ylabel)
    plt.title(title)

    plt.grid(
        alpha=0.25
    )

    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=300
    )

    plt.close()


# ============================================================
# 10. ANALYSIS E
# ============================================================

def analysis_E_robustness(metrics):

    print("\n" + "=" * 72)
    print("ANALYSIS E — ROBUSTNESS")
    print("=" * 72)

    outcomes = [
        "delta_kappa",
        "variance_reduction",
        "lambda_improvement",
        "normalized_rate",
    ]

    # --------------------------------------------------------
    # Leave-one-dataset-out
    # --------------------------------------------------------

    loo_rows = []

    for outcome in outcomes:

        values = pd.to_numeric(
            metrics[outcome],
            errors="coerce"
        )

        valid = metrics.loc[
            values.notna()
        ].copy()

        values = valid[outcome].to_numpy(
            dtype=float
        )

        if len(values) < 3:
            continue

        full_mean = np.mean(values)

        for idx in range(len(valid)):

            remaining = np.delete(
                values,
                idx
            )

            loo_rows.append({
                "outcome": outcome,
                "excluded_dataset": valid.iloc[idx]["dataset"],
                "full_mean": full_mean,
                "loo_mean": np.mean(remaining),
                "difference": (
                    np.mean(remaining)
                    - full_mean
                ),
            })

    loo = pd.DataFrame(loo_rows)

    loo.to_csv(
        STAT_DIR /
        "E_leave_one_dataset_out.csv",
        index=False,
        encoding="utf-8-sig"
    )

    # --------------------------------------------------------
    # Robustness summary
    # --------------------------------------------------------

    summary_rows = []

    for outcome in outcomes:

        sub = loo[
            loo["outcome"] == outcome
        ]

        if len(sub) == 0:
            continue

        full_mean = sub["full_mean"].iloc[0]

        loo_values = sub[
            "loo_mean"
        ].to_numpy(dtype=float)

        deviation = np.abs(
            loo_values - full_mean
        )

        summary_rows.append({
            "outcome": outcome,
            "n_leave_one_out": len(loo_values),
            "full_mean": full_mean,
            "loo_min": np.min(loo_values),
            "loo_max": np.max(loo_values),
            "loo_mean": np.mean(loo_values),
            "max_absolute_deviation": np.max(
                deviation
            ),
            "median_absolute_deviation": np.median(
                deviation
            ),
        })

    robustness_summary = pd.DataFrame(
        summary_rows
    )

    robustness_summary.to_csv(
        TABLE_DIR /
        "Table_5_Robustness_Analysis.csv",
        index=False,
        encoding="utf-8-sig"
    )

    # --------------------------------------------------------
    # Leave-one-domain-out
    # --------------------------------------------------------

    domain_rows = []

    for outcome in outcomes:

        values = pd.to_numeric(
            metrics[outcome],
            errors="coerce"
        )

        valid = metrics.loc[
            values.notna()
        ].copy()

        if len(valid) < 3:
            continue

        full_mean = valid[outcome].mean()

        for domain in sorted(
            valid["domain"].dropna().unique()
        ):

            remaining = valid[
                valid["domain"] != domain
            ]

            if len(remaining) < 2:
                continue

            domain_rows.append({
                "outcome": outcome,
                "excluded_domain": domain,
                "full_mean": full_mean,
                "loo_domain_mean": (
                    remaining[outcome].mean()
                ),
                "difference": (
                    remaining[outcome].mean()
                    - full_mean
                ),
                "remaining_n": len(remaining),
            })

    domain_loo = pd.DataFrame(
        domain_rows
    )

    domain_loo.to_csv(
        STAT_DIR /
        "E_leave_one_domain_out.csv",
        index=False,
        encoding="utf-8-sig"
    )


# ============================================================
# 11. ANALYSIS F
# ============================================================

def analysis_F_cross_domain(metrics):

    print("\n" + "=" * 72)
    print("ANALYSIS F — CROSS-DOMAIN VALIDATION")
    print("=" * 72)

    outcomes = [
        "delta_kappa",
        "variance_reduction",
        "lambda_improvement",
        "normalized_rate",
    ]

    # --------------------------------------------------------
    # Domain-level descriptive statistics
    # --------------------------------------------------------

    domain_rows = []

    for domain, sub in metrics.groupby(
        "domain",
        dropna=False
    ):

        domain_name = str(domain)

        row = {
            "domain": domain_name,
            "n_datasets": len(sub),
        }

        if "nodes" in sub.columns:
            row["median_nodes"] = (
                sub["nodes"].median()
            )

        if "edges" in sub.columns:
            row["median_edges"] = (
                sub["edges"].median()
            )

        if "density" in sub.columns:
            row["median_density"] = (
                sub["density"].median()
            )

        for outcome in outcomes:

            values = pd.to_numeric(
                sub[outcome],
                errors="coerce"
            ).dropna().to_numpy()

            if len(values) == 0:

                row[f"{outcome}_mean"] = np.nan
                row[f"{outcome}_median"] = np.nan
                row[f"{outcome}_ci95_low"] = np.nan
                row[f"{outcome}_ci95_high"] = np.nan

            else:

                ci_low, ci_high = bootstrap_mean_ci(
                    values,
                    iterations=3000
                )

                row[f"{outcome}_mean"] = np.mean(
                    values
                )

                row[f"{outcome}_median"] = np.median(
                    values
                )

                row[f"{outcome}_ci95_low"] = ci_low
                row[f"{outcome}_ci95_high"] = ci_high

        domain_rows.append(row)

    domain_table = pd.DataFrame(
        domain_rows
    )

    domain_table.to_csv(
        TABLE_DIR /
        "Table_6_Cross_Domain_Validation.csv",
        index=False,
        encoding="utf-8-sig"
    )

    # --------------------------------------------------------
    # Within-domain significance
    # --------------------------------------------------------

    significance_rows = []

    for domain, sub in metrics.groupby(
        "domain",
        dropna=False
    ):

        for outcome in outcomes:

            values = pd.to_numeric(
                sub[outcome],
                errors="coerce"
            ).dropna().to_numpy()

            if len(values) < 2:
                continue

            try:

                t_stat, t_p = stats.ttest_1samp(
                    values,
                    0
                )

            except Exception:

                t_stat = np.nan
                t_p = np.nan

            nonzero = values[
                values != 0
            ]

            if len(nonzero) >= 2:

                try:

                    result = stats.wilcoxon(
                        nonzero,
                        method="auto"
                    )

                    w_stat = float(
                        result.statistic
                    )

                    w_p = float(
                        result.pvalue
                    )

                except Exception:

                    w_stat = np.nan
                    w_p = np.nan

            else:

                w_stat = np.nan
                w_p = np.nan

            ci_low, ci_high = bootstrap_mean_ci(
                values,
                iterations=3000
            )

            significance_rows.append({
                "domain": domain,
                "outcome": outcome,
                "n": len(values),
                "mean": np.mean(values),
                "median": np.median(values),
                "ci95_low": ci_low,
                "ci95_high": ci_high,
                "t_statistic": t_stat,
                "t_pvalue": t_p,
                "wilcoxon_statistic": w_stat,
                "wilcoxon_pvalue": w_p,
            })

    domain_significance = pd.DataFrame(
        significance_rows
    )

    if len(domain_significance) > 0:

        domain_significance[
            "t_pvalue_fdr"
        ] = benjamini_hochberg(
            domain_significance[
                "t_pvalue"
            ].to_numpy()
        )

        domain_significance[
            "wilcoxon_pvalue_fdr"
        ] = benjamini_hochberg(
            domain_significance[
                "wilcoxon_pvalue"
            ].to_numpy()
        )

    domain_significance.to_csv(
        STAT_DIR /
        "F_domain_significance.csv",
        index=False,
        encoding="utf-8-sig"
    )

    # --------------------------------------------------------
    # Heterogeneity statistics
    # --------------------------------------------------------

    heterogeneity_rows = []

    for outcome in outcomes:

        domain_effects = []

        for domain, sub in metrics.groupby(
            "domain",
            dropna=False
        ):

            values = pd.to_numeric(
                sub[outcome],
                errors="coerce"
            ).dropna().to_numpy()

            if len(values) >= 2:

                mean_value = np.mean(values)
                variance = np.var(
                    values,
                    ddof=1
                )

                domain_effects.append(
                    (
                        str(domain),
                        len(values),
                        mean_value,
                        variance,
                    )
                )

        if len(domain_effects) < 2:
            continue

        total_n = sum(
            x[1]
            for x in domain_effects
        )

        weighted_mean = (
            sum(
                x[1] * x[2]
                for x in domain_effects
            )
            / total_n
        )

        Q = 0.0

        for (
            domain,
            n,
            mean_value,
            variance
        ) in domain_effects:

            Q += n * (
                mean_value
                - weighted_mean
            ) ** 2

        df_q = len(domain_effects) - 1

        if Q > 0:

            I2 = max(
                0.0,
                (Q - df_q) / Q
            ) * 100.0

        else:

            I2 = 0.0

        try:

            Q_p = stats.chi2.sf(
                Q,
                df_q
            )

        except Exception:

            Q_p = np.nan

        heterogeneity_rows.append({
            "outcome": outcome,
            "domains": len(domain_effects),
            "total_datasets": total_n,
            "weighted_mean": weighted_mean,
            "Q_statistic": Q,
            "Q_df": df_q,
            "Q_pvalue": Q_p,
            "I2_percent": I2,
        })

    heterogeneity = pd.DataFrame(
        heterogeneity_rows
    )

    heterogeneity.to_csv(
        STAT_DIR /
        "F_cross_domain_heterogeneity.csv",
        index=False,
        encoding="utf-8-sig"
    )

    # --------------------------------------------------------
    # Figure F1 — domain effect
    # --------------------------------------------------------

    make_domain_plot(
        domain_table,
        "delta_kappa_mean",
        "Δκ",
        "Cross-Domain Curvature Evolution",
        FIG_DIR /
        "Fig_F1_cross_domain_delta_kappa.png"
    )

    make_domain_plot(
        domain_table,
        "variance_reduction_mean",
        "Variance reduction",
        "Cross-Domain Variance Reduction",
        FIG_DIR /
        "Fig_F2_cross_domain_variance.png"
    )

    make_domain_plot(
        domain_table,
        "lambda_improvement_mean",
        "λ1 improvement",
        "Cross-Domain Spectral Improvement",
        FIG_DIR /
        "Fig_F3_cross_domain_lambda.png"
    )

    make_domain_plot(
        domain_table,
        "normalized_rate_mean",
        "Normalized evolution rate",
        "Cross-Domain Normalized Evolution Rate",
        FIG_DIR /
        "Fig_F4_cross_domain_rate.png"
    )


def make_domain_plot(
    domain_table,
    value_column,
    xlabel,
    title,
    output_path
):

    if value_column not in domain_table.columns:
        return

    plot_df = domain_table[
        [
            "domain",
            value_column,
        ]
    ].dropna()

    if len(plot_df) == 0:
        return

    plot_df = plot_df.sort_values(
        value_column
    )

    plt.figure(
        figsize=(12, max(6, 0.4 * len(plot_df)))
    )

    y = np.arange(
        len(plot_df)
    )

    plt.barh(
        y,
        plot_df[value_column].to_numpy()
    )

    plt.yticks(
        y,
        plot_df["domain"]
    )

    plt.xlabel(xlabel)
    plt.title(title)

    plt.axvline(
        0,
        linestyle="--",
        linewidth=1
    )

    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=300
    )

    plt.close()


# ============================================================
# 12. TABLE 1
# ============================================================

def build_table_1(metrics):

    table = metrics[
        [
            "dataset",
            "domain",
            "nodes",
            "edges",
            "components",
            "density",
            "weight_mode",
        ]
    ].copy()

    table.to_csv(
        TABLE_DIR /
        "Table_1_Dataset_Statistics.csv",
        index=False,
        encoding="utf-8-sig"
    )

    domain_summary = (
        metrics
        .groupby("domain")
        .agg(
            datasets=("dataset", "count"),
            median_nodes=("nodes", "median"),
            median_edges=("edges", "median"),
            median_density=("density", "median"),
        )
        .reset_index()
    )

    domain_summary.to_csv(
        TABLE_DIR /
        "Table_1B_Domain_Statistics.csv",
        index=False,
        encoding="utf-8-sig"
    )


# ============================================================
# 13. DATA QUALITY CHECK
# ============================================================

def data_quality_report(raw, summary, metrics):

    rows = []

    raw_counts = (
        raw
        .groupby("dataset")
        .size()
    )

    for dataset in sorted(
        raw["dataset"].unique()
    ):

        n_raw = int(
            raw_counts.get(
                dataset,
                0
            )
        )

        summary_match = int(
            (
                summary["dataset"]
                == dataset
            ).sum()
        )

        metric_match = int(
            (
                metrics["dataset"]
                == dataset
            ).sum()
        )

        rows.append({
            "dataset": dataset,
            "raw_records": n_raw,
            "summary_records": summary_match,
            "metric_records": metric_match,
            "expected_steps": STEPS_EXPECTED,
            "complete_10_steps": (
                n_raw == STEPS_EXPECTED
            ),
            "has_summary": (
                summary_match > 0
            ),
            "has_metrics": (
                metric_match > 0
            ),
        })

    report = pd.DataFrame(rows)

    report.to_csv(
        STAT_DIR /
        "data_quality_report.csv",
        index=False,
        encoding="utf-8-sig"
    )

    return report


# ============================================================
# 14. FINAL REPORT
# ============================================================

def generate_final_report(
    raw,
    summary,
    metrics,
    quality
):

    report_path = (
        REPORT_DIR /
        "analysis_report.txt"
    )

    with open(
        report_path,
        "w",
        encoding="utf-8"
    ) as f:

        f.write(
            "RICCI FLOW SIX-ANALYSIS REPORT\n"
        )

        f.write(
            "=" * 72 + "\n\n"
        )

        f.write(
            "INPUT DATA\n"
        )

        f.write(
            f"Raw results: {RAW_RESULTS_PATH}\n"
        )

        f.write(
            f"Dataset summary: {DATASET_SUMMARY_PATH}\n\n"
        )

        f.write(
            "DATASET COUNTS\n"
        )

        f.write(
            f"Raw records: {len(raw)}\n"
        )

        f.write(
            f"Raw datasets: "
            f"{raw['dataset'].nunique()}\n"
        )

        f.write(
            f"Summary datasets: "
            f"{summary['dataset'].nunique()}\n"
        )

        f.write(
            f"Valid metric datasets: "
            f"{len(metrics)}\n"
        )

        f.write(
            f"Domains: "
            f"{metrics['domain'].nunique()}\n\n"
        )

        # ----------------------------------------------------
        # Universal results
        # ----------------------------------------------------

        f.write(
            "UNIVERSAL EVOLUTION\n"
        )

        for variable in [
            "delta_kappa",
            "variance_reduction",
            "lambda_improvement",
            "normalized_rate",
            "convergence_score",
        ]:

            values = pd.to_numeric(
                metrics[variable],
                errors="coerce"
            ).dropna()

            if len(values) == 0:
                continue

            ci_low, ci_high = bootstrap_mean_ci(
                values.to_numpy()
            )

            f.write(
                f"{variable}: "
                f"mean={fmt(values.mean())}, "
                f"median={fmt(values.median())}, "
                f"95%CI=["
                f"{fmt(ci_low)}, "
                f"{fmt(ci_high)}], "
                f"n={len(values)}\n"
            )

        f.write("\n")

        # ----------------------------------------------------
        # Domain counts
        # ----------------------------------------------------

        f.write(
            "DOMAIN COUNTS\n"
        )

        domain_counts = (
            metrics["domain"]
            .value_counts()
        )

        for domain, count in domain_counts.items():

            f.write(
                f"{domain}: {count}\n"
            )

        f.write("\n")

        # ----------------------------------------------------
        # Data completeness
        # ----------------------------------------------------

        complete_fraction = (
            quality["complete_10_steps"]
            .mean()
            if len(quality)
            else np.nan
        )

        f.write(
            "DATA COMPLETENESS\n"
        )

        f.write(
            "Fraction with exactly 10 steps: "
            f"{fmt(complete_fraction)}\n"
        )

        f.write("\n")

        f.write(
            "IMPORTANT INTERPRETATION NOTE\n"
        )

        f.write(
            "The dataset is treated as the independent "
            "experimental unit. The ten Ricci-flow steps "
            "within each dataset are treated as repeated "
            "trajectory observations, not as ten independent "
            "replicates.\n"
        )

        f.write(
            "Statistical significance therefore operates "
            "primarily across datasets.\n"
        )

        f.write(
            "Cross-domain analysis uses the domain field "
            "provided in dataset_summary.csv.\n"
        )

    return report_path


# ============================================================
# 15. MAIN
# ============================================================

def main():

    print("\nLoading input data...")

    raw, summary = load_input_data()

    print("\nCalculating dataset-level metrics...")

    metrics = calculate_dataset_metrics(
        raw,
        summary
    )

    if len(metrics) == 0:
        raise RuntimeError(
            "No valid dataset-level records were generated."
        )

    # --------------------------------------------------------
    # Table 1
    # --------------------------------------------------------

    print("\nBuilding Table 1...")

    build_table_1(metrics)

    # --------------------------------------------------------
    # Data quality
    # --------------------------------------------------------

    print("\nRunning data quality checks...")

    quality = data_quality_report(
        raw,
        summary,
        metrics
    )

    # --------------------------------------------------------
    # Analysis A
    # --------------------------------------------------------

    analysis_A_universal_evolution(
        metrics
    )

    # --------------------------------------------------------
    # Analysis B
    # --------------------------------------------------------

    analysis_B_statistical_significance(
        metrics
    )

    # --------------------------------------------------------
    # Analysis C
    # --------------------------------------------------------

    analysis_C_correlation(
        metrics
    )

    # --------------------------------------------------------
    # Analysis D
    # --------------------------------------------------------

    analysis_D_convergence(
        raw,
        metrics
    )

    # --------------------------------------------------------
    # Analysis E
    # --------------------------------------------------------

    analysis_E_robustness(
        metrics
    )

    # --------------------------------------------------------
    # Analysis F
    # --------------------------------------------------------

    analysis_F_cross_domain(
        metrics
    )

    # --------------------------------------------------------
    # Final report
    # --------------------------------------------------------

    report_path = generate_final_report(
        raw,
        summary,
        metrics,
        quality
    )

    print("\n" + "=" * 72)
    print("ALL SIX ANALYSES FINISHED")
    print("=" * 72)

    print(
        f"\nOutput directory:\n{OUTPUT_DIR}"
    )

    print(
        f"\nFinal report:\n{report_path}"
    )

    print("\nTables:")

    for file in sorted(TABLE_DIR.glob("*.csv")):
        print(
            "  ",
            file.name
        )

    print("\nStatistics:")

    for file in sorted(STAT_DIR.glob("*.csv")):
        print(
            "  ",
            file.name
        )

    print("\nFigures:")

    for file in sorted(FIG_DIR.glob("*.png")):
        print(
            "  ",
            file.name
        )

    print("\nDone.")


if __name__ == "__main__":
    main()