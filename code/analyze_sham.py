"""Paired analysis of the sham-flow control.

Compares, for every dataset, the plateau reached by the real Ricci flow
against the plateau reached by the sham flow (curvature values randomly
permuted across edges at every step, everything else identical).

Both arms run the same number of steps, so the comparison is exactly paired.

H1 (structural)  predicts sham >> real  - destroying the curvature-to-edge
                 mapping destroys the homogenisation
H2 (mechanical)  predicts sham ~= real  - the update rule compresses the
                 weight spread regardless of which edge gets which curvature

Outputs: sham_results/analysis/{01_paired_table.csv,02_summary.txt,SHAM_REPORT.txt}
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


def load_arm(dataset, suffix):
    path = os.path.join(SHAM_DIR, "per_dataset", dataset + suffix)
    if not os.path.exists(path):
        return None
    frame = pd.read_csv(path)
    if frame.empty or "curvature_variance" not in frame.columns:
        return None
    v0 = float(frame["curvature_variance"].iloc[0])
    if not np.isfinite(v0) or v0 <= 0:
        return None
    return {
        "steps": int(len(frame)),
        "initial_variance": v0,
        "final_variance": float(frame["curvature_variance"].iloc[-1]),
        "var_ratio": float(frame["curvature_variance"].iloc[-1]) / v0,
        "ratio_at_half": float(
            frame["curvature_variance"].iloc[len(frame) // 2]) / v0,
        "min_ratio": float(frame["curvature_variance"].min()) / v0,
    }


def main():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)

    sham_files = sorted(glob.glob(
        os.path.join(SHAM_DIR, "per_dataset", "*__s*.csv")))
    datasets = sorted({os.path.basename(p).split("__s")[0] for p in sham_files})

    rows = []
    for dataset in datasets:
        real = load_arm(dataset, "__real.csv")
        if real is None:
            ref = os.path.join(REF_DIR, dataset + ".csv")
            if not os.path.exists(ref):
                continue
            frame = pd.read_csv(ref)
            v0 = float(frame["curvature_variance"].iloc[0])
            if not np.isfinite(v0) or v0 <= 0:
                continue
            real = {
                "steps": int(len(frame)),
                "initial_variance": v0,
                "final_variance": float(frame["curvature_variance"].iloc[-1]),
                "var_ratio": float(frame["curvature_variance"].iloc[-1]) / v0,
                "ratio_at_half": float(
                    frame["curvature_variance"].iloc[len(frame) // 2]) / v0,
                "min_ratio": float(frame["curvature_variance"].min()) / v0,
            }
        for path in sham_files:
            name = os.path.basename(path)
            if not name.startswith(dataset + "__s"):
                continue
            seed = name.split("__s")[1].replace(".csv", "")
            sham = load_arm(dataset, "__s{}.csv".format(seed))
            if sham is None:
                continue
            rows.append({
                "dataset": dataset, "seed": seed,
                "steps_real": real["steps"], "steps_sham": sham["steps"],
                "ratio_real": real["var_ratio"], "ratio_sham": sham["var_ratio"],
                "log_ratio_real": np.log10(max(real["var_ratio"], 1e-12)),
                "log_ratio_sham": np.log10(max(sham["var_ratio"], 1e-12)),
                "half_real": real["ratio_at_half"],
                "half_sham": sham["ratio_at_half"],
                "min_real": real["min_ratio"], "min_sham": sham["min_ratio"],
                "initial_variance": real["initial_variance"],
            })

    if not rows:
        print("no sham results found yet")
        return

    table = pd.DataFrame(rows)
    table.to_csv(os.path.join(OUT_DIR, "01_paired_table.csv"), index=False)

    lines = []
    lines.append("SHAM-FLOW CONTROL: PAIRED ANALYSIS")
    lines.append("=" * 62)
    lines.append("")
    lines.append("paired comparisons: {} (datasets {})".format(
        len(table), table["dataset"].nunique()))
    lines.append("")

    sortable = table.sort_values("ratio_real")
    lines.append("--- per dataset (sorted by real-flow plateau) ---")
    lines.append("  {:<44s} {:>8s} {:>8s} {:>9s}".format(
        "dataset", "real", "sham", "sham-real"))
    for _, row in sortable.iterrows():
        lines.append("  {:<44s} {:>8.4f} {:>8.4f} {:>+9.4f}".format(
            row["dataset"][:44], row["ratio_real"], row["ratio_sham"],
            row["ratio_sham"] - row["ratio_real"]))
    lines.append("")

    # ---- distribution summary
    lines.append("--- distribution ---")
    for label in ["ratio_real", "ratio_sham", "log_ratio_real", "log_ratio_sham"]:
        v = table[label].astype(float).values
        lines.append("  {:<16s} min={:>7.4f} p25={:>7.4f} med={:>7.4f} "
                     "p75={:>7.4f} max={:>7.4f}".format(
                         label, *np.nanpercentile(v, [0, 25, 50, 75, 100])))
    lines.append("")

    # ---- paired tests
    lines.append("--- paired tests (H1: sham plateau is HIGHER than real) ---")
    diff = (table["ratio_sham"] - table["ratio_real"]).values
    n = len(diff)
    n_higher = int((diff > 0).sum())
    lines.append("  sham higher than real in {}/{} datasets".format(n_higher, n))
    if n >= 5:
        try:
            bt = binomtest(n_higher, n, 0.5, alternative="greater")
            lines.append("  sign test                 p = {:.3g}".format(bt.pvalue))
        except Exception as exc:
            lines.append("  sign test failed: {}".format(exc))
        try:
            w = wilcoxon(diff, alternative="greater")
            lines.append("  Wilcoxon signed-rank      p = {:.3g}".format(w.pvalue))
        except Exception as exc:
            lines.append("  Wilcoxon failed: {}".format(exc))
        try:
            wl = wilcoxon(table["log_ratio_sham"] - table["log_ratio_real"],
                          alternative="greater")
            lines.append("  Wilcoxon on log10 ratio   p = {:.3g}".format(wl.pvalue))
        except Exception as exc:
            lines.append("  Wilcoxon (log) failed: {}".format(exc))
    lines.append("  median difference (sham-real) = {:+.4f}".format(
        float(np.median(diff))))
    lines.append("  median log10 fold change      = {:+.3f}".format(
        float(np.median(table["log_ratio_sham"] - table["log_ratio_real"]))))
    lines.append("")

    # ---- how far does sham get at the halfway point
    lines.append("--- trajectory midpoint (same step count for both arms) ---")
    lines.append("  median ratio at half of the steps: real={:.4f}  sham={:.4f}".format(
        float(np.median(table["half_real"])), float(np.median(table["half_sham"]))))
    lines.append("  median lowest ratio reached:       real={:.4f}  sham={:.4f}".format(
        float(np.median(table["min_real"])), float(np.median(table["min_sham"]))))
    lines.append("")

    # ---- verdict
    lines.append("--- verdict ---")
    med_real = float(np.median(table["ratio_real"]))
    med_sham = float(np.median(table["ratio_sham"]))
    lines.append("  median plateau: real={:.4f}   sham={:.4f}".format(
        med_real, med_sham))
    if med_sham > 0.8 and med_real < 0.35:
        lines.append("  => H1 supported: the sham flow barely homogenises at all.")
        lines.append("     Homogenisation REQUIRES the true curvature-to-edge mapping,")
        lines.append("     so the residual is a property of the network's geometry,")
        lines.append("     not an artefact of the update rule or of stopping too early.")
    elif abs(med_sham - med_real) < 0.1:
        lines.append("  => H2 supported: the sham flow reaches essentially the same")
        lines.append("     plateau. The variance drop is mechanical; the 'residual is")
        lines.append("     structural' claim must be withdrawn.")
    else:
        lines.append("  => intermediate: sham is worse than real but still drops.")
        lines.append("     Report the partial contribution explicitly; the structural")
        lines.append("     reading holds only for the sham-real gap.")

    text = "\n".join(lines)
    with open(os.path.join(OUT_DIR, "SHAM_REPORT.txt"), "w", encoding="utf-8") as fh:
        fh.write(text)
    with open(os.path.join(OUT_DIR, "02_summary.txt"), "w", encoding="utf-8") as fh:
        fh.write(text)
    print(text)


if __name__ == "__main__":
    main()
