"""Analyse the long-horizon Ricci flow run: does curvature variance converge to
zero, plateau at a positive residual, or blow up?

Reads ``ricci_flow_longrun_results`` and writes tables, figures and a text
report into an ``analysis`` sub-directory.
"""

from __future__ import print_function

import os

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


BASE = os.path.dirname(os.path.abspath(__file__))
RESULT_DIR = os.path.join(BASE, "ricci_flow_longrun_results")
OUT_DIR = os.path.join(RESULT_DIR, "analysis")
LEGACY_SUMMARY = os.path.join(BASE, "results_correct_fast", "dataset_summary.csv")


def convergence_diagnostics(raw):
    """Is a run really on a plateau, or was it cut while still falling?

    For every trajectory we fit the local exponential decay rate over the
    last 30% of the recorded steps:  var(t) ~ A exp(-c t)  =>  c > 0 means the
    variance was still shrinking when the run stopped.  Runs stopped by the
    wall-clock budget therefore only give an UPPER bound on the plateau.
    """
    rows = []
    for name, group in raw.groupby("dataset"):
        group = group.dropna(subset=["var_ratio"]).sort_values("step")
        if len(group) < 4:
            continue
        values = group["var_ratio"].values.astype(float)
        times = group["t"].values.astype(float)
        window = max(3, int(0.3 * len(values)))
        tail_v = values[-window:]
        tail_t = times[-window:]
        positive = tail_v > 0
        if positive.sum() >= 2:
            slope = -np.polyfit(tail_t[positive], np.log(tail_v[positive]), 1)[0]
        else:
            slope = float("nan")
        rows.append({
            "dataset": name,
            "steps": int(len(values)),
            "var_ratio_final": float(values[-1]),
            "decay_rate_c": float(slope) if np.isfinite(slope) else np.nan,
            "halving_time": float(np.log(2) / slope) if np.isfinite(slope) and slope > 0 else np.nan,
            "still_falling": bool(np.isfinite(slope) and slope > 0.1)
        })
    return pd.DataFrame(rows)


BUDGET_SECONDS = 600.0
MAX_STEPS = 400
TOL = 1e-4
PATIENCE = 5


def reconstruct_status(group):
    """Recover the stopping reason from a saved trajectory.

    The summary table can lose labels (a restarted worker rewrites finished
    datasets as ``cached``), but every trajectory carries enough information to
    rebuild the decision: step count, per-step wall time and the variance
    series.  This is derived from the data, not remembered.
    """
    group = group.sort_values("step")
    steps = len(group)
    wall = float(pd.to_numeric(group["step_time"], errors="coerce").fillna(0).sum())
    values = pd.to_numeric(group["var_ratio"], errors="coerce").values.astype(float)

    if steps <= 1:
        return "budget_insufficient"
    if np.isfinite(values[-1]) and values[-1] > 10.0:
        return "diverged"
    if steps >= MAX_STEPS:
        return "horizon_reached"
    if np.isfinite(values[-1]) and values[-1] <= 1e-6:
        return "constant_curvature"
    if steps > PATIENCE:
        moves = np.abs(np.diff(values[-(PATIENCE + 1):]))
        if np.all(np.isfinite(moves)) and np.all(moves < TOL):
            return "saturated"
    if wall >= 0.9 * BUDGET_SECONDS:
        return "budget_exceeded"
    return "stopped_early"


def status_map(raw):
    rows = []
    for name, group in raw.groupby("dataset"):
        rows.append({"dataset": name, "status": reconstruct_status(group)})
    return pd.DataFrame(rows)


def covariates(raw):
    """Initial level, initial spread and size, taken from step 0 of each run."""
    first = raw[raw["step"] == 0][["dataset", "mean_curvature",
                                   "curvature_variance", "nodes"]].copy()
    first = first.rename(columns={"mean_curvature": "initial_kappa",
                                  "curvature_variance": "initial_variance"})
    first["log10_nodes"] = np.log10(first["nodes"].astype(float))
    return first[["dataset", "initial_kappa", "initial_variance", "log10_nodes"]]


def stratified_correlations(merged, covars, label):
    out = []
    subset = merged.dropna(subset=["var_ratio_final"])
    for column in covars:
        if column not in subset.columns:
            continue
        pair = subset[["var_ratio_final", column]].dropna()
        if len(pair) < 8:
            continue
        rho = pair.corr(method="spearman").iloc[0, 1]
        out.append({"stratum": label, "n": len(pair), "covariate": column,
                    "spearman": float(rho)})
    return out


def quantile_summary(values, name):
    values = np.asarray([v for v in values if np.isfinite(v)], dtype=float)
    if values.size == 0:
        return {"quantity": name, "n": 0}
    return {
        "quantity": name,
        "n": int(values.size),
        "min": float(np.min(values)),
        "p10": float(np.percentile(values, 10)),
        "p25": float(np.percentile(values, 25)),
        "median": float(np.median(values)),
        "p75": float(np.percentile(values, 75)),
        "p90": float(np.percentile(values, 90)),
        "max": float(np.max(values))
    }


def load_summary():
    """Merge every worker summary file; prefer the richest record per dataset."""
    import glob
    paths = sorted(glob.glob(os.path.join(RESULT_DIR, "longrun_summary*.csv")))
    frames = []
    for path in paths:
        try:
            frame = pd.read_csv(path)
        except Exception:
            continue
        if "dataset" in frame.columns and len(frame):
            frames.append(frame)
    if not frames:
        raise SystemExit("no summary files found in {}".format(RESULT_DIR))
    merged = pd.concat(frames, ignore_index=True)
    # a worker that actually computed a dataset knows more than one that
    # merely saw the cached trajectory, so rank cached last
    merged["_rank"] = (merged["status"].astype(str) == "cached").astype(int)
    merged = merged.sort_values("_rank").drop_duplicates("dataset", keep="first")
    return merged.drop(columns=["_rank"]).reset_index(drop=True)


def main():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)

    summary = load_summary()
    per_dir = os.path.join(RESULT_DIR, "per_dataset")
    frames = []
    for name in sorted(os.listdir(per_dir)):
        if name.endswith(".csv"):
            try:
                frames.append(pd.read_csv(os.path.join(per_dir, name)))
            except Exception:
                pass
    raw = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    partial = not os.path.exists(os.path.join(RESULT_DIR, "longrun_raw.csv"))

    if os.path.exists(LEGACY_SUMMARY):
        legacy = pd.read_csv(LEGACY_SUMMARY)[["dataset", "domain"]]
        summary = summary.merge(legacy, on="dataset", how="left")
        raw = raw.merge(legacy, on="dataset", how="left")
    else:
        summary["domain"] = ""
        raw["domain"] = ""

    # ---- rebuild the stopping reason from the trajectories ------------------
    # A restarted worker rewrites finished datasets as "cached", which destroys
    # the labels.  The trajectories are the ground truth, so derive status from
    # them and keep the summary label only where no trajectory exists.
    recon = status_map(raw)
    original = summary[["dataset", "status"]].copy()
    original = original[original["status"].astype(str) != "cached"]
    original = original.rename(columns={"status": "status_orig"})
    summary = summary.merge(recon, on="dataset", how="left")
    summary["status"] = summary["status_y"].fillna(summary["status_x"]) \
        if "status_y" in summary.columns else summary["status"]
    for column in ("status_x", "status_y"):
        if column in summary.columns:
            summary = summary.drop(columns=[column])

    check = summary.merge(original, on="dataset", how="inner")
    if len(check):
        agree = float((check["status"] == check["status_orig"]).mean())
        disagreements = check[check["status"] != check["status_orig"]]
    else:
        agree = float("nan")
        disagreements = check

    lines = []
    lines.append("LONG-HORIZON RICCI FLOW - ANALYSIS REPORT")
    lines.append("=" * 60)
    if partial:
        lines.append("PARTIAL RESULTS - the batch is still running")
        lines.append("")
    lines.append("datasets in summary: {}".format(len(summary)))
    lines.append("trajectory rows: {}".format(len(raw)))
    lines.append("")

    counts = summary["status"].value_counts()
    lines.append("--- outcome counts ---")
    for status, count in counts.items():
        lines.append("  {:<20s} {}".format(str(status), int(count)))
    lines.append("")

    lines.append("--- status reconstruction check ---")
    lines.append("  stopping reasons rebuilt from saved trajectories: {}".format(
        int(summary["status"].notna().sum())))
    if len(check):
        lines.append("  agreement with labels still present in summary files: {:.1%} (n={})".format(
            agree, len(check)))
        if len(disagreements):
            lines.append("  disagreements (summary label -> rebuilt label):")
            for _, row in disagreements.head(10).iterrows():
                lines.append("    {}: {} -> {}".format(
                    row["dataset"], row["status_orig"], row["status"]))
    lines.append("")

    completed = summary[summary["status"].isin(
        ["saturated", "constant_curvature", "horizon_reached", "budget_exceeded",
         "cached"])].copy()
    completed = completed[completed["var_ratio_final"].notna()]

    lines.append("--- final variance ratio var(T)/var(0) ---")
    rows = []
    for status, group in completed.groupby("status"):
        stats = quantile_summary(group["var_ratio_final"].values, str(status))
        stats["status"] = str(status)
        rows.append(stats)
    ratio_table = pd.DataFrame(rows)
    ratio_table.to_csv(os.path.join(OUT_DIR, "01_plateau_ratio_by_outcome.csv"), index=False)
    lines.append(ratio_table.to_string(index=False))
    lines.append("")

    # how much of the drop had already happened at the original horizon (t=0.5)
    at10 = raw[raw["step"] == 10][["dataset", "var_ratio"]].rename(
        columns={"var_ratio": "var_ratio_at_t0.5"})
    merged = completed.merge(at10, on="dataset", how="left")
    merged["residual_drop_after_t0.5"] = (
        merged["var_ratio_at_t0.5"] - merged["var_ratio_final"])
    merged["fraction_of_drop_in_first_10_steps"] = np.where(
        (1.0 - merged["var_ratio_final"]).abs() > 1e-12,
        (1.0 - merged["var_ratio_at_t0.5"]) / (1.0 - merged["var_ratio_final"]),
        np.nan
    )
    merged.to_csv(os.path.join(OUT_DIR, "02_transient_vs_plateau.csv"), index=False)

    lines.append("--- was the 10-step result only a transient? ---")
    stats = quantile_summary(merged["var_ratio_at_t0.5"].values, "var_ratio at t=0.5")
    lines.append("  var_ratio at t=0.5   median={:.4f} n={}".format(
        stats.get("median", float("nan")), stats.get("n", 0)))
    stats = quantile_summary(merged["var_ratio_final"].values, "var_ratio at stop")
    lines.append("  var_ratio at stop    median={:.4f} n={}".format(
        stats.get("median", float("nan")), stats.get("n", 0)))
    stats = quantile_summary(merged["fraction_of_drop_in_first_10_steps"].values,
                             "share of total drop realised in first 10 steps")
    lines.append("  share of drop in first 10 steps: median={:.3f} p10={:.3f} p90={:.3f}".format(
        stats.get("median", float("nan")), stats.get("p10", float("nan")),
        stats.get("p90", float("nan"))))
    lines.append("")

    # time to saturate
    steps_table = quantile_summary(completed["steps"].values, "steps executed")
    lines.append("--- steps executed (among completed runs) ---")
    lines.append("  median={} p10={} p90={} max={}".format(
        steps_table.get("median"), steps_table.get("p10"),
        steps_table.get("p90"), steps_table.get("max")))
    lines.append("")

    # plateau vs size and initial curvature
    first = raw[raw["step"] == 0][["dataset", "mean_curvature", "curvature_variance"]].rename(
        columns={"mean_curvature": "initial_kappa",
                 "curvature_variance": "initial_variance"})
    # drop the same-named columns coming from the summary table so the merge
    # cannot produce _x / _y duplicates
    merged = merged.drop(columns=[c for c in ("initial_variance", "initial_kappa")
                                  if c in merged.columns])
    merged = merged.merge(first, on="dataset", how="left")
    merged["log10_nodes"] = np.log10(merged["nodes"].astype(float))

    corr_cols = ["var_ratio_final", "log10_nodes", "initial_kappa", "initial_variance", "steps"]
    available = [c for c in corr_cols if c in merged.columns]
    corr = merged[available].corr(method="spearman")
    corr.to_csv(os.path.join(OUT_DIR, "03_spearman_plateau_predictors.csv"))
    lines.append("--- spearman correlations with the plateau ratio ---")
    for column in available:
        if column == "var_ratio_final":
            continue
        value = corr.loc["var_ratio_final", column]
        lines.append("  rho(var_ratio_final, {:<18s}) = {:+.3f}".format(column, value))
    lines.append("")

    # domain table
    if merged["domain"].notna().any():
        domain_rows = []
        for domain, group in merged.groupby("domain"):
            values = group["var_ratio_final"].dropna().values
            if values.size == 0:
                continue
            domain_rows.append({
                "domain": domain,
                "n": int(values.size),
                "median_plateau_ratio": float(np.median(values)),
                "median_steps": float(np.median(group["steps"].dropna().values))
                    if group["steps"].notna().any() else None,
                "n_to_zero": int(np.sum(values < 0.01)),
                "n_plateau_above_0.5": int(np.sum(values > 0.5))
            })
        domain_table = pd.DataFrame(domain_rows).sort_values("median_plateau_ratio")
        domain_table.to_csv(os.path.join(OUT_DIR, "04_domain_plateau.csv"), index=False)
        lines.append("--- plateau ratio by domain (ascending) ---")
        lines.append(domain_table.to_string(index=False))
        lines.append("")

    # divergence diagnostics
    diverged = summary[summary["status"] == "diverged"]
    if len(diverged):
        lines.append("--- diverged runs ---")
        for _, row in diverged.iterrows():
            lines.append("  {}: after {} steps (t={:.2f}) variance ratio reached {:.3g}".format(
                row["dataset"], row["steps"], row["steps"] * 0.05,
                float(row["var_ratio_final"]) if pd.notna(row["var_ratio_final"]) else float("nan")))
        lines.append("  (the original exception text is only in the worker log; "
                     "cause is weight collapse -> d -> 0 -> kappa blow-up)")
        lines.append("")

    # convergence diagnostics: separate real plateaus from truncated runs
    diag = convergence_diagnostics(raw)
    diag = diag.merge(summary[["dataset", "status", "domain"]], on="dataset", how="left")
    CONVERGED = ("saturated", "constant_curvature", "horizon_reached")
    diag["converged"] = diag["status"].isin(CONVERGED)
    # save only after the flag exists, otherwise downstream tools cannot split
    diag.to_csv(os.path.join(OUT_DIR, "05_convergence_diagnostics.csv"), index=False)
    conv = diag[diag["converged"]]
    trunc = diag[~diag["converged"] & diag["var_ratio_final"].notna()]
    trunc_falling = trunc[trunc["still_falling"]]

    lines.append("--- convergence diagnostics ---")
    lines.append("  runs that stopped on a criterion (saturated / constant / horizon): {}".format(len(conv)))
    lines.append("  runs cut by the wall-clock budget:                                {}".format(len(trunc)))
    lines.append("  of those, still decaying fast at the cut (c > 0.1 per unit t):     {}".format(len(trunc_falling)))
    if len(trunc_falling):
        lines.append("    their reported plateau is an UPPER BOUND; median = {:.4f}".format(
            float(np.median(trunc_falling["var_ratio_final"].values))))
        lines.append("    median decay rate among them: c = {:.3f} per unit t "
                     "(halving time {:.2f})".format(
                         float(np.nanmedian(trunc_falling["decay_rate_c"].values)),
                         float(np.nanmedian(trunc_falling["halving_time"].values))))
    lines.append("")

    # the size effect is confounded with truncation: bigger graphs are slower
    # per step, so they get cut more often.  Test the correlation inside the
    # stratum of runs that genuinely converged.
    covars = ["log10_nodes", "initial_kappa", "initial_variance", "steps"]
    strat_rows = stratified_correlations(merged, covars, "all completed")
    diag_merged = diag.merge(covariates(raw), on="dataset", how="left")
    strat_rows += stratified_correlations(
        diag_merged[diag_merged["converged"]], covars, "converged only")
    strat_rows += stratified_correlations(
        diag_merged[~diag_merged["converged"]], covars, "truncated only")
    strat = pd.DataFrame(strat_rows)
    strat.to_csv(os.path.join(OUT_DIR, "06_stratified_correlations.csv"), index=False)
    lines.append("--- plateau predictors, stratified by whether the run converged ---")
    lines.append("  (truncation is confounded with size: slow = big = cut early)")
    lines.append(strat.to_string(index=False))
    lines.append("")

    # ---- one authoritative per-dataset table -------------------------------
    # Downstream tools must read this, not the worker summary files: the summary
    # files lose their stopping reasons whenever a worker is restarted.
    master = diag_merged.copy()
    if "nodes" not in master.columns:
        master = master.merge(summary[["dataset", "nodes"]], on="dataset", how="left")
    missing = summary[~summary["dataset"].isin(master["dataset"])]
    if len(missing):
        extra = pd.DataFrame({
            "dataset": missing["dataset"].values,
            "status": missing["status"].values,
            "steps": 0,
            "var_ratio_final": np.nan,
            "converged": False,
        })
        if "domain" in missing.columns:
            extra["domain"] = missing["domain"].values
        master = pd.concat([master, extra], ignore_index=True)
    master.to_csv(os.path.join(OUT_DIR, "00_dataset_table.csv"), index=False)

    # verdict
    lines.append("--- verdict ---")
    plateau_values = completed["var_ratio_final"].dropna().values
    if plateau_values.size:
        near_zero = float(np.mean(plateau_values < 0.05))
        positive = float(np.mean(plateau_values > 0.3))
        lines.append("  among all completed long runs:")
        lines.append("    fraction with plateau ratio < 0.05 (variance -> 0): {:.2f}".format(near_zero))
        lines.append("    fraction with plateau ratio > 0.30 (partial only): {:.2f}".format(positive))
        lines.append("    median plateau ratio: {:.4f}".format(float(np.median(plateau_values))))
        lines.append("")
        if len(conv):
            degenerate = diag[diag["var_ratio_final"].isna()
                              | (diag["var_ratio_final"].abs() > 1e3)]
            zero_var = raw[raw["step"] == 0]
            zero_var = zero_var[pd.to_numeric(zero_var["curvature_variance"],
                                              errors="coerce").fillna(0) <= 0]
            lines.append("")
            lines.append("  degenerate cases that cannot carry a variance-ratio claim:")
            lines.append("    datasets with zero/undefined initial curvature variance: {} ({})".format(
                len(zero_var), ", ".join(sorted(zero_var["dataset"].tolist())[:8])))
            lines.append("    datasets whose ratio is NaN or exploded:               {} ({})".format(
                len(degenerate), ", ".join(sorted(degenerate["dataset"].tolist())[:8])))
            lines.append("    these must be excluded from any 'variance shrinks' statement, "
                         "not silently dropped")
            frozen = diag[(diag["converged"]) & (diag["var_ratio_final"] >= 0.99)]
            increased = diag[diag["var_ratio_final"] > 1.01]
            lines.append("")
            lines.append("  counterexamples to a naive 'everything homogenises' claim:")
            lines.append("    runs ending with variance >= 99% of the initial value: {} ({})".format(
                len(frozen), ", ".join(sorted(frozen["dataset"].tolist())[:8])))
            lines.append("    runs ending with variance ABOVE the initial value:    {} ({})".format(
                len(increased), ", ".join(sorted(increased["dataset"].tolist())[:8])))
            conv_values = conv["var_ratio_final"].dropna().values
            lines.append("  restricted to the {} runs that genuinely converged:".format(len(conv_values)))
            lines.append("    fraction with plateau ratio < 0.05: {:.2f}".format(
                float(np.mean(conv_values < 0.05))))
            lines.append("    fraction with plateau ratio > 0.30: {:.2f}".format(
                float(np.mean(conv_values > 0.30))))
            lines.append("    median plateau ratio: {:.4f}".format(float(np.median(conv_values))))
            lines.append("    quartiles: p25={:.4f} p75={:.4f}".format(
                float(np.percentile(conv_values, 25)),
                float(np.percentile(conv_values, 75))))
    else:
        lines.append("  no run stopped on a convergence criterion yet.")
    lines.append("")

    # figures
    sample = raw[raw["dataset"].isin(completed["dataset"].head(25))]
    if len(sample):
        plt.figure(figsize=(7.5, 5))
        for name, group in sample.groupby("dataset"):
            plt.plot(group["t"], group["var_ratio"], linewidth=0.9, alpha=0.6)
        plt.xlabel("flow time t")
        plt.ylabel("var(kappa, t) / var(kappa, 0)")
        plt.title("Curvature variance trajectories (first 25 completed datasets)")
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, "fig_trajectories.png"), dpi=130)
        plt.close()

    if "var_ratio" in raw.columns and raw["var_ratio"].notna().any():
        grid = np.linspace(0, raw["t"].max(), 60)
        curves = []
        for name, group in raw.groupby("dataset"):
            group = group.dropna(subset=["var_ratio"])
            if len(group) < 3:
                continue
            curves.append(np.interp(grid, group["t"], group["var_ratio"]))
        if curves:
            matrix = np.asarray(curves)
            plt.figure(figsize=(7.5, 5))
            plt.plot(grid, np.median(matrix, axis=0), color="#185FA5",
                     linewidth=2.2, label="median")
            plt.fill_between(grid,
                             np.percentile(matrix, 25, axis=0),
                             np.percentile(matrix, 75, axis=0),
                             color="#185FA5", alpha=0.18, label="IQR")
            plt.axhline(0.0, color="#888780", linewidth=0.8, linestyle="--")
            plt.xlabel("flow time t")
            plt.ylabel("var(kappa, t) / var(kappa, 0)")
            plt.title("Median variance trajectory over all datasets")
            plt.legend()
            plt.grid(alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(OUT_DIR, "fig_median_envelope.png"), dpi=130)
            plt.close()

    if merged["var_ratio_final"].notna().any():
        plt.figure(figsize=(7.5, 5))
        plt.scatter(merged["log10_nodes"], merged["var_ratio_final"],
                    s=18, alpha=0.7, color="#185FA5")
        plt.xlabel("log10(nodes)")
        plt.ylabel("final variance ratio")
        plt.title("Plateau level versus network size")
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, "fig_plateau_vs_size.png"), dpi=130)
        plt.close()

    report_path = os.path.join(OUT_DIR, "LONG_REPORT.txt")
    with open(report_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))

    print("\n".join(lines))
    print("\nreport written to {}".format(report_path))


if __name__ == "__main__":
    main()
