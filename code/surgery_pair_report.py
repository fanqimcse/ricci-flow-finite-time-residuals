# -*- coding: utf-8 -*-
"""Pair the no-surgery plateau against the surgery plateau.

Rules that keep this honest
---------------------------
1. A run with 0 removals is bit-identical to the no-surgery arm. Including it
   dilutes the comparison towards "no difference" for a purely mechanical
   reason, so those runs are reported separately and excluded from the test.
2. A run that hit its wall-clock budget is still falling. Its endpoint is a
   truncation artefact, not a plateau, so completed arms are analysed on their
   own.
3. Removing 1-2 edges is a perturbation, not surgery: the graph barely changed
   but the variance trajectory is restarted. Those are excluded from the
   headline and shown separately.
"""

import os
import glob

import numpy as np
import pandas as pd
from scipy.stats import binomtest

BASE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE, "ricci_flow_longrun_results", "analysis")
LONGRUN = os.path.join(OUT_DIR, "07_merged_with_structure.csv")
SURG_GLOB = os.path.join(BASE, "surgery_results", "*", "per_dataset", "*.csv")
DONE = {"saturated", "constant_curvature", "horizon_reached"}
MIN_REMOVED = 3


def main():
    base = pd.read_csv(LONGRUN)[
        ["dataset", "var_ratio_final", "status", "nodes", "edges", "family"]]
    base = base.rename(columns={"var_ratio_final": "no_surgery_ratio",
                                "status": "no_surgery_status"})
    base = base.drop_duplicates(subset=["dataset"], keep="first")

    rows = []
    for path in sorted(glob.glob(SURG_GLOB)):
        head = os.path.relpath(path, os.path.join(BASE, "surgery_results"))
        if head.split(os.sep)[0] in ("test", "test2"):
            continue
        name = os.path.basename(path)[:-4]
        try:
            frame = pd.read_csv(path)
        except Exception:
            continue
        if frame.empty:
            continue
        tol = (float(frame["surgery_tol"].iloc[0])
               if "surgery_tol" in frame.columns else np.nan)
        ratio = frame["ratio"].dropna()
        rows.append({
            "dataset": name, "tol": tol,
            "surgery_ratio": float(ratio.iloc[-1]) if len(ratio) else np.nan,
            "surgery_steps": int(frame["step"].max()) + 1,
            "surgery_status": str(frame["status"].iloc[-1]),
            "edges_removed": int(frame["removed"].max()),
            "edges_before": int(frame["edges"].iloc[0]),
        })
    surg = pd.DataFrame(rows)
    if surg.empty:
        raise SystemExit("no surgery results found")
    surg = (surg.sort_values("surgery_steps", ascending=False)
                .drop_duplicates(subset=["dataset", "tol"], keep="first"))

    m = surg.merge(base, on="dataset", how="left")
    m["fired"] = m["edges_removed"] > 0
    m["completed"] = m["surgery_status"].isin(DONE)
    m["substantial"] = m["edges_removed"] >= MIN_REMOVED
    m["drop_factor"] = m["no_surgery_ratio"] / m["surgery_ratio"]

    L = []
    A = L.append
    A("SURGERY vs NO-SURGERY   (paired, identical numerics)")
    A("=" * 66)

    for tol, grp in m.groupby("tol"):
        A("")
        A("### tol = {}".format(tol))
        A("{:32s} {:>9s} {:>9s} {:>6s} {:>8s} {:>10s}".format(
            "dataset", "no_surg", "surgery", "rm", "drop x", "status"))
        for _, r in grp.sort_values(
                ["substantial", "completed"], ascending=False).iterrows():
            A("{:32s} {:9.4f} {:9.4f} {:6d} {:8.1f} {:>10s}".format(
                str(r["dataset"])[:32], float(r["no_surgery_ratio"]),
                float(r["surgery_ratio"]), int(r["edges_removed"]),
                float(r["drop_factor"]) if np.isfinite(r["drop_factor"])
                else float("nan"),
                str(r["surgery_status"])[:10]))

    A("")
    A("--- HEADLINE SET: removed >= {} AND surgery arm completed ---".format(
        MIN_REMOVED))
    head = m[m["substantial"] & m["completed"] & np.isfinite(m["drop_factor"])]
    A("{:32s} {:>9s} {:>9s} {:>6s} {:>8s}".format(
        "dataset", "no_surg", "surgery", "rm", "drop x"))
    for _, r in head.sort_values("drop_factor", ascending=False).iterrows():
        A("{:32s} {:9.4f} {:9.4f} {:6d} {:8.1f}".format(
            str(r["dataset"])[:32], float(r["no_surgery_ratio"]),
            float(r["surgery_ratio"]), int(r["edges_removed"]),
            float(r["drop_factor"])))
    if len(head):
        better = int((head["surgery_ratio"] < head["no_surgery_ratio"]).sum())
        p = binomtest(better, len(head), 0.5, alternative="greater").pvalue
        A("")
        A("  surgery arm lower in {}/{} datasets".format(better, len(head)))
        A("  one-sided exact sign test p = {:.4f} (n={})".format(p, len(head)))
        A("  drop factor: median {:.1f}x, range {:.1f}x .. {:.1f}x".format(
            float(head["drop_factor"].median()),
            float(head["drop_factor"].min()),
            float(head["drop_factor"].max())))
        A("  families covered: {}".format(
            sorted(head["family"].dropna().unique().tolist())))
        A("  edges removed: {}..{} of {}..{} total".format(
            int(head["edges_removed"].min()), int(head["edges_removed"].max()),
            int(head["edges_before"].min()), int(head["edges_before"].max())))

    A("")
    A("--- excluded, and why ---")
    A("  surgery never fired (arm identical to no-surgery): {}".format(
        ", ".join(sorted(m.loc[~m["fired"], "dataset"].tolist())) or "none"))
    A("  fired but < {} edges (perturbation, not surgery): {}".format(
        MIN_REMOVED,
        ", ".join(sorted(m.loc[m["fired"] & ~m["substantial"],
                               "dataset"].tolist())) or "none"))
    A("  surgery arm truncated by the wall-clock budget: {}".format(
        ", ".join(sorted(m.loc[m["substantial"] & ~m["completed"],
                               "dataset"].tolist())) or "none"))

    A("")
    A("--- reading ---")
    if len(head) and float(head["drop_factor"].median()) > 2:
        A("  ANSWER B. Where surgery actually removes a meaningful number of")
        A("  edges and the arm runs to completion, the residual collapses.")
        A("  The plateau is largely the cost of NOT doing surgery.")
        A("")
        A("  Consequence for the paper: the claim cannot be 'a structural")
        A("  obstruction to Ricci flow'. It must be 'the uniformisation limit")
        A("  of the TOPOLOGY-PRESERVING flow' - how far a graph is from")
        A("  constant curvature if its edge set may not change. Surgery-based")
        A("  flows (Ma-Yang 2025) destroy exactly this quantity, so the result")
        A("  is complementary to, not in conflict with, that literature.")
    else:
        A("  ANSWER A (or inconclusive). The residual survives surgery.")

    out = os.path.join(BASE, "surgery_results", "SURGERY_REPORT_v2.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))
    m.to_csv(os.path.join(BASE, "surgery_results", "surgery_paired_v2.csv"),
             index=False)
    print("\n".join(L))


if __name__ == "__main__":
    main()
