# -*- coding: utf-8 -*-
"""
Convergence-rate probe: does the Ricci-flow homogenization rate track the
spectral gap?

Model to test
-------------
    V(t) = var(t)/var(0)  ~=  V_inf + (1 - V_inf) * exp(-r * t)

Hypothesis H1 (needed before any convergence-rate theorem):
    r is an increasing function of lambda1(0), the ALGEBRAIC CONNECTIVITY of
    the initial (unweighted) graph Laplacian; ideally r ~ c * lambda1^a with
    a ~ 1.

Definitions used here (must match the flow code):
  * lambda1 = 2nd smallest eigenvalue of the weighted graph Laplacian
    (weight="weight"); returns 0.0 for disconnected graphs, so those are
    excluded by construction (they carry no usable gap).
  * one flow step = dt = 0.05 in the update w <- w*exp(-dt*kappa).
  * V_inf is estimated from the tail of the trajectory (see below).

Falsification rules fixed in advance
------------------------------------
  R1. Primary statistic: Spearman(r, lambda1_0) on the CONVERGED subset
      (status == saturated). horizon_reached runs are reported separately
      and must never be pooled into the headline.
  R2. r must survive partial correlation controlling log10(nodes) and
      log10(mean degree); lambda1 is itself size-dependent.
  R3. Within-family check on ENZYMES (largest family). If the pooled
      correlation is driven entirely by between-family contrast, report it
      as ecological, not predictive.
  R4. Exponent a is reported with its CI; a ~ 0 (no dependence) is a
      negative result and kills the theorem route.
"""

import os
import json
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
RESULT_DIR = os.path.join(HERE, "ricci_flow_longrun_results")
PER_DIR = os.path.join(RESULT_DIR, "per_dataset")
OUT_DIR = os.path.join(HERE, "convergence_results")
os.makedirs(OUT_DIR, exist_ok=True)

MIN_STEPS = 40          # need enough trajectory to fit a rate
TAIL_FRAC = 0.15        # last 15% of steps define the plateau estimate
FIT_LOW = 0.03          # fit window: (V - Vinf) / (1 - Vinf) in [FIT_LOW, FIT_HIGH]
FIT_HIGH = 0.60
EARLY_FRAC = 0.20       # early-time window for the plateau-free rate estimate


def early_rate(step, vr):
    """Plateau-free rate: slope of log V over the leading EARLY_FRAC of the
    run. No V_inf needed, so truncation cannot bias it."""
    n = len(step)
    k = max(5, int(np.floor(n * EARLY_FRAC)))
    x = step[:k]
    y = np.log(np.clip(vr[:k], 1e-12, None))
    slope, _, rr, _, se = stats.linregress(x, y)
    return {"rate_early": float(-slope), "rate_early_r2": float(rr ** 2),
            "rate_early_se": float(se)}


def joint_fit(step, vr):
    """Three-parameter fit V = Vinf + (1-Vinf) exp(-r t); Vinf is free, so a
    truncated run is not forced to use its last point as the plateau."""
    from scipy.optimize import curve_fit
    def model(t, vinf, r):
        return vinf + (1.0 - vinf) * np.exp(-r * t)
    try:
        popt, _ = curve_fit(model, step, vr, p0=[max(vr[-1], 1e-6), 0.05],
                            bounds=([0.0, 1e-6], [1.0, 5.0]), maxfev=20000)
        resid = vr - model(step, *popt)
        ss = 1.0 - np.sum(resid ** 2) / np.sum((vr - vr.mean()) ** 2)
        return {"v_inf_joint": float(popt[0]), "rate_joint": float(popt[1]),
                "rate_joint_r2": float(ss)}
    except Exception:
        return {"v_inf_joint": float("nan"), "rate_joint": float("nan"),
                "rate_joint_r2": float("nan")}


def load_summary():
    """Prefer the worker-2 summary (real stop reasons); fall back to w1."""
    p2 = os.path.join(RESULT_DIR, "longrun_summary_w2.csv")
    p1 = os.path.join(RESULT_DIR, "longrun_summary.csv")
    p = p2 if os.path.exists(p2) else p1
    return pd.read_csv(p), os.path.basename(p)


def trajectory_converged(step, vr, frac=0.25, tol=0.02):
    """Self-declared plateau test: relative drift of V over the last `frac`
    of the trajectory is below `tol`. Avoids relying on the 'cached' status
    whose stop reason was not recorded."""
    n = len(step)
    k = max(3, int(np.floor(n * frac)))
    tail = vr[-k:]
    start = tail[0]
    if not np.isfinite(start) or abs(start) < 1e-12:
        return False, float("nan")
    drift = abs(tail[-1] - start) / abs(start)
    return bool(drift < tol), float(drift)


def fit_rate(step, var_ratio):
    """Return dict with plateau estimate, per-step decay rate, t_half."""
    n = len(step)
    tail_start = max(1, int(np.floor(n * (1.0 - TAIL_FRAC))))
    v_inf = float(np.median(var_ratio[tail_start:]))
    span = var_ratio[0] - v_inf
    if not np.isfinite(span) or span <= 0:
        return None
    excess = (var_ratio - v_inf) / span
    # monotone-ish decay region
    mask = (excess > FIT_LOW) & (excess < FIT_HIGH)
    # keep only the leading contiguous block (avoid re-crossings)
    idx = np.where(mask)[0]
    if idx.size >= 5:
        # cut at first gap
        breaks = np.where(np.diff(idx) > 1)[0]
        if breaks.size:
            idx = idx[: breaks[0] + 1]
        if idx.size >= 5:
            x = step[idx].astype(float)
            y = np.log(excess[idx])
            slope, intercept, rr, pv, se = stats.linregress(x, y)
            n_fit = idx.size
        else:
            return None
    else:
        return None
    # t_half: first step where excess drops below 0.5
    below = np.where(excess <= 0.5)[0]
    t_half = float(step[below[0]]) if below.size else float("nan")
    return {
        "v_inf": v_inf,
        "rate_per_step": float(-slope),
        "rate_se": float(se),
        "rate_r2": float(rr ** 2),
        "t_half": t_half,
        "n_fit": int(n_fit),
    }


def main():
    summary, summary_name = load_summary()
    rows = []
    for _, s in summary.iterrows():
        name = s["dataset"]
        fpath = os.path.join(PER_DIR, name + ".csv")
        if not os.path.exists(fpath):
            continue
        try:
            tr = pd.read_csv(fpath)
        except Exception:
            continue
        if len(tr) < MIN_STEPS:
            continue
        # initial variance taken from the trajectory's own step 0, because the
        # summary column is mostly missing
        init_var = float(tr["curvature_variance"].astype(float).values[0])
        if not np.isfinite(init_var) or init_var <= 1e-30:
            continue
        lam = tr["lambda1"].astype(float).values
        lam0 = lam[0]
        if not np.isfinite(lam0) or lam0 <= 0:
            continue
        step = tr["step"].astype(float).values
        vr = tr["var_ratio"].astype(float).values
        vr = vr / vr[0]  # renormalise defensively (should already start at 1)
        fit = fit_rate(step, vr)
        if fit is None:
            continue
        conv_flag, drift = trajectory_converged(step, vr)
        # nodes/edges taken from the trajectory itself (the summary columns
        # are NaN for part of the sample)
        nodes = float(tr["nodes"].astype(float).iloc[0])
        edges = float(tr["edges"].astype(float).iloc[0])
        er = early_rate(step, vr)
        jf = joint_fit(step, vr)
        rows.append({
            "dataset": name,
            "status": s.get("status"),
            "stop_reason": s.get("stop_reason"),
            "traj_converged": conv_flag,
            "tail_drift": drift,
            "steps": int(s.get("steps", 0)),
            "nodes": nodes,
            "edges": edges,
            "mean_degree": (2.0 * edges / nodes) if nodes else np.nan,
            "initial_variance": init_var,
            "lambda1_0": float(lam0),
            "lambda1_final": float(lam[-1]) if np.isfinite(lam[-1]) else np.nan,
            "var_ratio_final": float(s.get("var_ratio_final", np.nan)),
            "v_inf_fit": fit["v_inf"],
            "rate_per_step": fit["rate_per_step"],
            "rate_se": fit["rate_se"],
            "rate_r2": fit["rate_r2"],
            "t_half": fit["t_half"],
            "n_fit": fit["n_fit"],
            "rate_early": er["rate_early"],
            "rate_early_r2": er["rate_early_r2"],
            "rate_joint": jf["rate_joint"],
            "v_inf_joint": jf["v_inf_joint"],
            "rate_joint_r2": jf["rate_joint_r2"],
        })

    df = pd.DataFrame(rows)
    df["family"] = df["dataset"].str.split("_").str[0]
    df["log10_nodes"] = np.log10(df["nodes"])
    df["log10_md"] = np.log10(df["mean_degree"])
    df["log_rate"] = np.log10(df["rate_per_step"].clip(lower=1e-12))
    df["log_lambda"] = np.log10(df["lambda1_0"])
    df["log_thalf"] = np.log10(df["t_half"])
    df.to_csv(os.path.join(OUT_DIR, "convergence_rates.csv"), index=False)

    lines = []
    def w(x=""):
        lines.append(str(x))
        print(x)

    w("CONVERGENCE-RATE PROBE:  rate vs spectral gap")
    w("=" * 66)
    w("model: V(t) = Vinf + (1-Vinf) exp(-r t);  lambda1 = 2nd smallest")
    w("       eigenvalue of the weighted Laplacian (0 for disconnected).")
    w("")
    w("summary used: %s" % summary_name)
    w("units with usable fit: %d   (of %d datasets in the summary)" % (len(df), len(summary)))
    w("status breakdown:")
    for k, v in df["stop_reason"].value_counts().items():
        w("    %-18s %d" % (k, v))
    w("trajectory plateau test (last 25%% drift < 2%%): %d / %d converged"
      % (int(df["traj_converged"].sum()), len(df)))
    w("")

    conv = df[df["traj_converged"]].copy()
    w("--- PRIMARY (trajectory-converged subset, n=%d) ---" % len(conv))
    if len(conv) >= 8:
        rho, pv = stats.spearmanr(conv["rate_per_step"], conv["lambda1_0"])
        w("  Spearman(rate, lambda1_0)      rho = %+.3f  p = %.3g" % (rho, pv))
        rho_t, pv_t = stats.spearmanr(conv["t_half"], conv["lambda1_0"])
        w("  Spearman(t_half, lambda1_0)    rho = %+.3f  p = %.3g  (expect negative)" % (rho_t, pv_t))
        sl, ic, rr, pp, se = stats.linregress(conv["log_lambda"], conv["log_rate"])
        w("  log-log fit  log10 r = %+.3f + %.3f * log10 lambda1" % (ic, sl))
        w("     exponent a = %.3f  (95%% CI %.3f .. %.3f), R2 = %.3f, p = %.3g"
          % (sl, sl - 1.96 * se, sl + 1.96 * se, rr ** 2, pp))
        w("     => r = %.4g * lambda1^%.3f" % (10 ** ic, sl))
        # R2: partial correlation controlling size + mean degree
        for ctrl in (["log10_nodes"], ["log10_md"], ["log10_nodes", "log10_md"]):
            pr, ppv = partial_corr(conv, "rate_per_step", "lambda1_0", ctrl)
            w("  partial (rate, lambda1_0 | %-28s) rho = %+.3f p = %.3g"
              % (",".join(ctrl), pr, ppv))
    else:
        w("  too few converged units (n=%d) for the primary test" % len(conv))

    w("")
    w("--- SECONDARY (all statuses pooled; headline must NOT use this) ---")
    rho2, pv2 = stats.spearmanr(df["rate_per_step"], df["lambda1_0"])
    w("  Spearman(rate, lambda1_0)      rho = %+.3f  p = %.3g  (n=%d)" % (rho2, pv2, len(df)))

    w("")
    w("--- WITHIN FAMILY (R3: is the effect ecological?) ---")
    fam_counts = df["family"].value_counts()
    for fam in fam_counts.index:
        sub = df[df["family"] == fam]
        if len(sub) < 8:
            continue
        r_, p_ = stats.spearmanr(sub["rate_per_step"], sub["lambda1_0"])
        w("  %-12s n=%-3d rho = %+.3f  p = %.3g" % (fam, len(sub), r_, p_))

    w("")
    w("--- RATE DEFINITIONS x PREDICTORS (full usable sample, n=%d) ---" % len(df))
    df["lam_over_md"] = df["lambda1_0"] / df["mean_degree"]
    df["lam_times_n"] = df["lambda1_0"] * df["nodes"]
    rates = ["rate_early", "rate_joint", "rate_per_step"]
    preds = ["lambda1_0", "lam_over_md", "lam_times_n", "lambda1_final"]
    w("  %-14s %s" % ("", "".join("%14s" % p for p in preds)))
    for rt in rates:
        cells = []
        for p in preds:
            sub = df[[rt, p]].replace([np.inf, -np.inf], np.nan).dropna()
            sub = sub[sub[rt] > 0]
            if len(sub) >= 8:
                r_, p_ = stats.spearmanr(sub[rt], sub[p])
                cells.append("%+7.3f(%.2f)" % (r_, p_))
            else:
                cells.append("%14s" % "n/a")
        w("  %-14s %s" % (rt, "".join("%14s" % c for c in cells)))
    w("  (cells are Spearman rho with p in parentheses)")
    w("")
    w("  lambda1_0 range: %.4g .. %.4g   (median %.4g)"
      % (df["lambda1_0"].min(), df["lambda1_0"].max(), df["lambda1_0"].median()))
    for rt in rates:
        v = df[rt].dropna()
        v = v[v > 0]
        w("  %-14s range %.4g .. %.4g (median %.4g)" % (rt, v.min(), v.max(), v.median()))
    w("")
    w("  partial correlations (rate_early, lambda1_0 | controls):")
    for ctrl in (["log10_nodes"], ["log10_md"], ["log10_nodes", "log10_md"]):
        pr, ppv = partial_corr(df, "rate_early", "lambda1_0", ctrl)
        w("     | %-28s rho = %+.3f p = %.3g" % (",".join(ctrl), pr, ppv))
    w("")
    w("  is the lambda1*n signal just size? (must be checked before quoting it)")
    r_, p_ = stats.spearmanr(df["lam_times_n"], df["log10_nodes"])
    w("     Spearman(lambda1*n, log10_nodes) rho = %+.3f p = %.3g" % (r_, p_))
    pr, ppv = partial_corr(df, "rate_early", "lam_times_n", ["log10_nodes"])
    w("     partial (rate_early, lambda1*n | log10_nodes) rho = %+.3f p = %.3g" % (pr, ppv))
    w("")
    w("  secondary: faster decay vs lower plateau (both plateau-dependent)")
    pr, ppv = partial_corr(df, "rate_early", "var_ratio_final", ["log10_nodes"])
    r_, p_ = stats.spearmanr(df["rate_early"], df["var_ratio_final"])
    w("     raw    rho = %+.3f p = %.3g" % (r_, p_))
    w("     | log10_nodes  rho = %+.3f p = %.3g" % (pr, ppv))
    pr, ppv = partial_corr(df, "rate_early", "var_ratio_final", ["log10_nodes", "log10_md"])
    w("     | nodes+degree rho = %+.3f p = %.3g" % (pr, ppv))

    w("")
    w("  within-family (ENZYMES only, guards against ecological correlation):")
    for rt in rates:
        sub = df[df["family"] == "ENZYMES"][[rt, "lambda1_0"]].dropna()
        if len(sub) >= 8 and sub[rt].gt(0).all():
            r_, p_ = stats.spearmanr(sub[rt], sub["lambda1_0"])
            w("     %-14s n=%d rho = %+.3f p = %.3g" % (rt, len(sub), r_, p_))
    w("")
    w("  POWER NOTE: with n=%d, |rho| >= %.2f is detectable at alpha=0.05;"
      % (len(df), 2.0 / np.sqrt(len(df))))
    w("  a null here is informative, not merely underpowered.")

    w("")
    w("--- SCALING: log10(rate_early) ~ log10(lambda1_0) + log10(nodes) + log10(mean_degree) ---")
    import itertools
    sub = df[["log_rate", "log_lambda", "log10_nodes", "log10_md"]].replace(
        [np.inf, -np.inf], np.nan).dropna()
    sub = sub[pd.notna(df.loc[sub.index, "rate_early"])]
    sub = sub[df.loc[sub.index, "rate_early"] > 0]
    if len(sub) >= 15:
        X = np.column_stack([np.ones(len(sub)), sub["log_lambda"], sub["log10_nodes"],
                             sub["log10_md"]])
        y = sub["log_rate"].values
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta
        dof = len(sub) - X.shape[1]
        s2 = resid @ resid / dof
        cov = s2 * np.linalg.pinv(X.T @ X)
        se = np.sqrt(np.diag(cov))
        names = ["const", "log10 lambda1", "log10 nodes", "log10 mean_degree"]
        for nm, b, s in zip(names, beta, se):
            t = b / s if s > 0 else float("nan")
            pv = 2 * (1 - stats.t.cdf(abs(t), dof)) if np.isfinite(t) else float("nan")
            w("     %-20s beta = %+.3f  se = %.3f  p = %.3g" % (nm, b, s, pv))
        w("     R2 = %.3f  (n=%d)" % (1 - resid @ resid / np.sum((y - y.mean()) ** 2), len(sub)))
        w("     => rate ~ lambda1^%.2f * nodes^%.2f * mean_degree^%.2f"
          % (beta[1], beta[2], beta[3]))

    w("")
    w("--- SANITY: what DOES the rate track? ---")
    san = conv if len(conv) >= 8 else df
    w("  (using %s, n=%d)" % ("converged subset" if len(conv) >= 8 else "full sample", len(san)))
    for col in ("log10_nodes", "log10_md", "mean_degree", "initial_variance",
                "lambda1_0", "var_ratio_final"):
        if col not in san.columns:
            continue
        sub = san[["rate_early", col]].dropna()
        if len(sub) < 8:
            continue
        r_, p_ = stats.spearmanr(sub["rate_early"], sub[col])
        w("  Spearman(rate_early, %-16s) rho = %+.3f p = %.3g" % (col, r_, p_))

    w("")
    w("--- FIT QUALITY ---")
    w("  median R2 of the decay fits: per-step %.3f / early %.3f / joint %.3f"
      % (df["rate_r2"].median(), df["rate_early_r2"].median(), df["rate_joint_r2"].median()))
    w("  joint-fit Vinf vs trajectory tail V: median |diff| = %.4f"
      % float((df["v_inf_joint"] - df["v_inf_fit"]).abs().median()))
    w("")
    w("--- INTERPRETATION RULE ---")
    w("  a ~ 1 with a surviving partial correlation => convergence-rate")
    w("  theorem route is viable (predict r from a static spectral quantity).")
    w("  a ~ 0 or partial collapses to 0 => rate is set by something other")
    w("  than the Laplacian gap; do NOT build the theorem on lambda1.")

    with open(os.path.join(OUT_DIR, "CONVERGENCE_RATE_REPORT.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def partial_corr(df, x, y, controls):
    """Spearman-based partial correlation via rank regression residuals."""
    sub = df[[x, y] + controls].replace([np.inf, -np.inf], np.nan).dropna()
    sub = sub[sub[x] > 0]
    if len(sub) < len(controls) + 6:
        return float("nan"), float("nan")
    R = sub.rank()
    def resid(target, others):
        X = np.column_stack([np.ones(len(R))] + [R[o].values for o in others])
        beta, *_ = np.linalg.lstsq(X, R[target].values, rcond=None)
        return R[target].values - X @ beta
    rx = resid(x, controls)
    ry = resid(y, controls)
    return stats.spearmanr(rx, ry)


if __name__ == "__main__":
    main()
