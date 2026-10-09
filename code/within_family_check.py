# -*- coding: utf-8 -*-
"""Is the assortativity -> plateau effect a WITHIN-family effect?

The leave-one-family-out test came back essentially null (pooled oos R2 ~ 0.04,
p = 0.11) and even negative inside ENZYMES. That is the signature of an
ECOLOGICAL correlation: a relation that exists between group means but not
between individuals inside a group. This script separates the two, because
they license completely different claims.

  * within-family  : does a network with higher assortativity than its
                     siblings also retain more curvature variance?
  * between-family : do families with higher mean assortativity retain more?

If the effect is only between-family and rests on ~7 points, the paper may
say "a between-family tendency" and MUST NOT say "assortativity predicts the
residual of a network". The synthetic sweep (beta = +0.34, p = 3e-10) stays
valid either way: there the degree sequence is literally held fixed, so the
causal statement is about the synthetic ensemble, not about real networks.
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


def rho_p(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 5 or np.std(x[m]) == 0 or np.std(y[m]) == 0:
        return np.nan, np.nan, int(m.sum())
    r, p = spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


def main():
    frame = pd.read_csv(SRC)
    lines = []
    lines.append("WITHIN-FAMILY vs BETWEEN-FAMILY DECOMPOSITION")
    lines.append("=" * 62)

    if "initial_variance" in frame.columns:
        frame = frame[frame["initial_variance"].astype(float) > MIN_VARIANCE]
    if "nodes" in frame.columns:
        frame = frame[frame["nodes"].astype(float) >= MIN_NODES]
    frame = frame[np.isfinite(frame["var_ratio_final"].astype(float))]
    frame = frame.copy()
    frame["log_ratio"] = np.log10(frame["var_ratio_final"].astype(float))
    lines.append("sample after exclusion chain: {} networks / {} families".format(
        len(frame), frame["family"].nunique()))
    lines.append("")

    # ---- 1. within-family
    lines.append("--- 1. WITHIN-FAMILY Spearman(assortativity, plateau) ---")
    lines.append("{:12s} {:>4s} {:>9s} {:>9s} {:>9s} {:>9s}".format(
        "family", "n", "rho", "p", "assort_sd", "log10sd"))
    rows = []
    for fam, sub in frame.groupby("family"):
        sub = sub.dropna(subset=["assortativity", "log_ratio"])
        if len(sub) < 5:
            lines.append("{:12s} {:>4d}  (n<5, not estimable)".format(
                fam, len(sub)))
            continue
        r, p, n = rho_p(sub["assortativity"], sub["log_ratio"])
        rows.append({"family": fam, "n": n, "rho": r, "p": p})
        lines.append("{:12s} {:>4d} {:>9.3f} {:>9.3f} {:>9.3f} {:>9.3f}".format(
            fam, n, r if r == r else float("nan"),
            p if p == p else float("nan"),
            float(np.std(sub["assortativity"])),
            float(np.std(sub["log_ratio"]))))
    wf = pd.DataFrame(rows)
    if len(wf):
        lines.append("  median within-family rho = {:+.3f} ({} families)".format(
            float(np.nanmedian(wf["rho"])), len(wf)))
        pos = int((wf["rho"] > 0).sum())
        lines.append("  positive within-family: {}/{}".format(pos, len(wf)))
    lines.append("")

    # same for modularity, for completeness
    lines.append("--- 1b. WITHIN-FAMILY Spearman(modularity, plateau) ---")
    rows = []
    for fam, sub in frame.groupby("family"):
        sub = sub.dropna(subset=["modularity", "log_ratio"])
        if len(sub) < 5:
            continue
        r, p, n = rho_p(sub["modularity"], sub["log_ratio"])
        rows.append({"family": fam, "n": n, "rho": r, "p": p})
        lines.append("{:12s} {:>4d} {:>9.3f} {:>9.3f}".format(
            fam, n, r if r == r else float("nan"),
            p if p == p else float("nan")))
    if rows:
        wm = pd.DataFrame(rows)
        lines.append("  median within-family rho = {:+.3f}".format(
            float(np.nanmedian(wm["rho"]))))
    lines.append("")

    # ---- 2. between-family
    lines.append("--- 2. BETWEEN-FAMILY (family means) ---")
    fam_means = frame.groupby("family").agg(
        n=("log_ratio", "size"),
        mean_log_ratio=("log_ratio", "mean"),
        mean_assort=("assortativity", "mean"),
        mean_modul=("modularity", "mean"),
        mean_trans=("transitivity", "mean"),
    ).reset_index()
    fam_means = fam_means[fam_means["n"] >= 2]
    lines.append(fam_means.to_string(index=False,
                                     float_format=lambda v: "{:+.3f}".format(v)))
    for cov in ("mean_assort", "mean_modul", "mean_trans"):
        r, p, n = rho_p(fam_means[cov], fam_means["mean_log_ratio"])
        lines.append("  Spearman({:12s}, mean plateau) rho={:+.3f} p={:.3f} "
                     "n={}".format(cov, r, p, n))
    lines.append("")

    # ---- 3. jackknife: how fragile is the between-family rho?
    lines.append("--- 3. JACKKNIFE on the between-family rho (drop one family) ---")
    base_r, base_p, _ = rho_p(fam_means["mean_assort"],
                              fam_means["mean_log_ratio"])
    lines.append("  full-sample family-level rho = {:+.3f} (p={:.3f}, n={})".format(
        base_r, base_p, len(fam_means)))
    jr = []
    for i in range(len(fam_means)):
        sub = fam_means.drop(fam_means.index[i])
        r, p, n = rho_p(sub["mean_assort"], sub["mean_log_ratio"])
        jr.append(r)
        lines.append("  drop {:12s} -> rho={:+.3f}".format(
            str(fam_means.iloc[i]["family"]), r))
    jr = np.asarray(jr, dtype=float)
    lines.append("  jackknife range: {:+.3f} .. {:+.3f}".format(
        float(np.nanmin(jr)), float(np.nanmax(jr))))
    lines.append("")

    # ---- 4. variance decomposition: how much of log-plateau is between-family?
    grand = float(frame["log_ratio"].mean())
    ss_total = float(((frame["log_ratio"] - grand) ** 2).sum())
    ss_between = float(sum(len(sub) * (sub["log_ratio"].mean() - grand) ** 2
                           for _, sub in frame.groupby("family")))
    icc = ss_between / ss_total if ss_total > 0 else np.nan
    lines.append("--- 4. where does the variance live? ---")
    lines.append("  ICC (fraction of log-plateau variance BETWEEN families) "
                 "= {:.3f}".format(icc))
    lines.append("  => {:.0f}% between-family / {:.0f}% within-family".format(
        100 * icc, 100 * (1 - icc)))

    lines.append("")
    lines.append("--- VERDICT RULE ---")
    lines.append("  Effect is usable at network level only if within-family rho")
    lines.append("  is consistently positive AND the LOFO oos R2 > 0. Otherwise")
    lines.append("  the correlation is ecological: describe it as a between-")
    lines.append("  family tendency over n families, and never as a predictor.")

    with open(os.path.join(OUT_DIR, "WITHIN_FAMILY_REPORT.txt"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    wf.to_csv(os.path.join(OUT_DIR, "17_within_family.csv"), index=False)
    fam_means.to_csv(os.path.join(OUT_DIR, "17_family_means.csv"), index=False)
    print("ok")


if __name__ == "__main__":
    main()
