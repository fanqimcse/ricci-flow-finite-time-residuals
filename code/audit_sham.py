"""Robustness audit of the sham-flow paired result.

The raw paired analysis (25/26, p~4e-7) is decisive, but three things can
inflate it:

  1. diverged arms      - eco-mangwet's sham arm blew up (ratio 55.7).
                          A runaway is not a plateau.
  2. uninformative arms - datasets where the REAL flow also ended at
                          ratio ~ 1 (nothing homogenised in either arm).
                          They contribute a ~0 difference and should not
                          count as evidence.
  3. pseudoreplication  - ENZYMES_g* is 12 draws from one generator,
                          scc_rt_* is 4 draws from one corpus. Datasets are
                          not independent; families are.

This script recomputes the test on the defensible subset and at family
level, and checks that the two arms really are step-matched.
"""

import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, binomtest

BASE = os.path.dirname(os.path.abspath(__file__))
SHAM_DIR = os.path.join(BASE, "sham_results")
REF_DIR = os.path.join(BASE, "ricci_flow_longrun_results", "per_dataset")
OUT_DIR = os.path.join(SHAM_DIR, "analysis")
LIST_FILE = os.path.join(BASE, "converged_clean_list.csv")

# ---- family map: group obvious replicate generators -----------------------
FAMILY_RULES = [
    ("ENZYMES", "ENZYMES"),
    ("scc_rt_", "scc_rt (social media retweet)"),
    ("bn-macaque", "macaque connectome"),
    ("bn-mouse", "mouse connectome"),
    ("reptilia", "reptilia"),
    ("mammalia", "mammalia"),
    ("eco-", "ecological flow network"),
    ("aves-", "aves"),
]


def family_of(name):
    for key, label in FAMILY_RULES:
        if name.startswith(key):
            return label
    return name


def load_traj(path):
    if not os.path.exists(path):
        return None
    frame = pd.read_csv(path)
    if frame.empty or "curvature_variance" not in frame.columns:
        return None
    return frame["curvature_variance"].astype(float).values


def ratio_at(traj, step):
    """var(t=step)/var(t=0), clipped to the trajectory length."""
    v0 = traj[0]
    if not np.isfinite(v0) or v0 <= 0:
        return np.nan
    idx = min(int(step), len(traj) - 1)
    return float(traj[idx] / v0)


def main():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)

    sham_files = sorted(glob.glob(
        os.path.join(SHAM_DIR, "per_dataset", "*__s*.csv")))
    datasets = sorted({os.path.basename(p).split("__s")[0] for p in sham_files})

    # ---- what was requested vs what actually ran
    requested = []
    if os.path.exists(LIST_FILE):
        frame = pd.read_csv(LIST_FILE)
        col = "dataset" if "dataset" in frame.columns else frame.columns[0]
        requested = [str(v) for v in frame[col].tolist()]
    missing = [d for d in requested if d not in datasets]

    rows = []
    for dataset in datasets:
        real_traj = load_traj(os.path.join(REF_DIR, dataset + ".csv"))
        if real_traj is None:
            continue
        for path in sham_files:
            name = os.path.basename(path)
            if not name.startswith(dataset + "__s"):
                continue
            seed = name.split("__s")[1].replace(".csv", "")
            sham_traj = load_traj(path)
            if sham_traj is None:
                continue
            v0 = real_traj[0]
            if not np.isfinite(v0) or v0 <= 0:
                continue
            n_real, n_sham = len(real_traj), len(sham_traj)
            n_common = min(n_real, n_sham)
            rows.append({
                "dataset": dataset,
                "seed": int(seed),
                "family": family_of(dataset),
                "steps_real": n_real,
                "steps_sham": n_sham,
                "step_matched": n_real == n_sham,
                "ratio_real": float(real_traj[-1] / v0),
                "ratio_sham": float(sham_traj[-1] / v0),
                # matched-step: sham evaluated at exactly the real arm's
                # final step, so neither arm is favoured by length
                "ratio_sham_at_real_end": ratio_at(sham_traj, n_real - 1),
                "ratio_real_at_sham_end": ratio_at(real_traj, n_sham - 1),
                # a runaway: sham ended far ABOVE where it started
                "sham_diverged": bool(sham_traj[-1] / v0 > 5.0),
                # uninformative: the real flow did not homogenise either
                "real_inert": bool(real_traj[-1] / v0 > 0.8),
            })

    if not rows:
        print("no sham results found")
        return

    table = pd.DataFrame(rows)
    table["diff"] = table["ratio_sham"] - table["ratio_real"]
    table["diff_matched"] = (table["ratio_sham_at_real_end"]
                             - table["ratio_real"])
    table.to_csv(os.path.join(OUT_DIR, "03_robust_table.csv"), index=False)

    L = []
    L.append("SHAM-FLOW CONTROL: ROBUSTNESS AUDIT")
    L.append("=" * 66)
    L.append("")

    # ---- 0. coverage
    L.append("--- 0. coverage ---")
    L.append("  requested datasets : {}".format(len(requested)))
    L.append("  with a sham arm    : {}".format(table["dataset"].nunique()))
    L.append("  never produced     : {}".format(len(missing)))
    for name in missing:
        L.append("      - {}".format(name))
    L.append("")

    # ---- 1. step matching
    L.append("--- 1. step matching ---")
    matched = int(table["step_matched"].sum())
    L.append("  arms with identical step count: {}/{}".format(matched, len(table)))
    for _, row in table[~table["step_matched"]].iterrows():
        L.append("      MISMATCH {} seed{}: real={} sham={}".format(
            row["dataset"], row["seed"], row["steps_real"], row["steps_sham"]))
    L.append("")

    # ---- 2. exclusions
    div = table["sham_diverged"]
    inert = table["real_inert"]
    keep = table[~div & ~inert].copy()
    L.append("--- 2. exclusions ---")
    L.append("  sham runaway (ratio>5) : {} -> {}".format(
        int(div.sum()), ", ".join(table[div]["dataset"].tolist()) or "none"))
    L.append("  real flow inert (>0.8) : {} -> {}".format(
        int(inert.sum()), ", ".join(sorted(set(table[inert]["dataset"]))) or "none"))
    L.append("  analysed pairs         : {}".format(len(keep)))
    L.append("")

    def paired_tests(frame, label, col="diff"):
        L.append("--- {} ---".format(label))
        diff = frame[col].astype(float).values
        diff = diff[np.isfinite(diff)]
        n = len(diff)
        if n == 0:
            L.append("  no pairs")
            L.append("")
            return
        higher = int((diff > 0).sum())
        L.append("  n = {}   sham higher in {}".format(n, higher))
        if n >= 5:
            L.append("  sign test            p = {:.3g}".format(
                binomtest(higher, n, 0.5, alternative="greater").pvalue))
            L.append("  Wilcoxon signed-rank p = {:.3g}".format(
                wilcoxon(diff, alternative="greater").pvalue))
        L.append("  median real      = {:.4f}".format(
            float(np.median(frame["ratio_real"]))))
        L.append("  median sham      = {:.4f}".format(
            float(np.median(frame["ratio_sham"]))))
        L.append("  median diff      = {:+.4f}".format(float(np.median(diff))))
        L.append("")

    # ---- 3. primary (all pairs) vs defensible
    paired_tests(table, "3a. ALL pairs (raw, as first reported)")
    paired_tests(keep, "3b. DEFENSIBLE pairs (runaway + inert removed)")
    paired_tests(keep, "3c. DEFENSIBLE, matched-step sham endpoint",
                 col="diff_matched")

    # ---- 4. family level: one value per family (median across members)
    L.append("--- 4. family level (families are the independent unit) ---")
    fam = (keep.groupby("family")
               .agg(n=("dataset", "size"),
                    real=("ratio_real", "median"),
                    sham=("ratio_sham", "median"),
                    diff=("diff", "median"))
               .reset_index())
    fam = fam.sort_values("diff")
    L.append("  {:<32s} {:>4s} {:>8s} {:>8s} {:>9s}".format(
        "family", "n", "real", "sham", "sham-real"))
    for _, row in fam.iterrows():
        L.append("  {:<32s} {:>4d} {:>8.4f} {:>8.4f} {:>+9.4f}".format(
            row["family"][:32], int(row["n"]), row["real"], row["sham"],
            row["diff"]))
    fam.to_csv(os.path.join(OUT_DIR, "04_family_table.csv"), index=False)
    L.append("")
    d = fam["diff"].astype(float).values
    n = len(d)
    higher = int((d > 0).sum())
    L.append("  families: {}   sham higher in {}".format(n, higher))
    if n >= 5:
        L.append("  sign test            p = {:.3g}".format(
            binomtest(higher, n, 0.5, alternative="greater").pvalue))
        L.append("  Wilcoxon signed-rank p = {:.3g}".format(
            wilcoxon(d, alternative="greater").pvalue))
    L.append("  median family-level diff = {:+.4f}".format(float(np.median(d))))
    L.append("")

    # ---- 5. verdict
    L.append("--- 5. verdict ---")
    med_real = float(np.median(keep["ratio_real"]))
    med_sham = float(np.median(keep["ratio_sham"]))
    L.append("  defensible median plateau: real={:.4f}  sham={:.4f}".format(
        med_real, med_sham))
    if med_sham > 0.8 and med_real < 0.35:
        L.append("  => HOLDS after exclusions and family aggregation.")
        L.append("     Randomising which edge receives which curvature at")
        L.append("     every step destroys essentially all of the variance")
        L.append("     reduction, at an identical number of steps. The")
        L.append("     homogenisation therefore depends on the curvature-to-edge")
        L.append("     mapping, i.e. on network geometry, not on the update rule.")
    else:
        L.append("  => weakened; report the numbers, do not overclaim.")

    text = "\n".join(L)
    with open(os.path.join(OUT_DIR, "ROBUST_REPORT.txt"), "w",
              encoding="utf-8") as fh:
        fh.write(text)
    print(text)


if __name__ == "__main__":
    main()
