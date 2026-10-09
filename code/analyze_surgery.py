# -*- coding: utf-8 -*-
"""Pair the no-surgery plateau against the surgery plateau.

Only datasets where surgery actually fired (edges_removed > 0) carry
information: a run with 0 removals is bit-identical to the no-surgery arm and
must be excluded from the paired comparison, otherwise it dilutes the test
towards "no difference" for a purely mechanical reason.
"""

import os
import glob

import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE, "ricci_flow_longrun_results", "analysis")
LONGRUN = os.path.join(OUT_DIR, "07_merged_with_structure.csv")
SURG_GLOB = os.path.join(BASE, "surgery_results", "*", "per_dataset", "*.csv")


def main():
    base = pd.read_csv(LONGRUN)
    base = base[["dataset", "var_ratio_final", "status", "nodes", "edges"]]
    base = base.rename(columns={"var_ratio_final": "no_surgery_ratio",
                                "status": "no_surgery_status"})

    rows = []
    for path in sorted(glob.glob(SURG_GLOB)):
        name = os.path.basename(path)[:-4]
        # scratch dirs from smoke tests must not enter the analysis
        head = os.path.relpath(path, os.path.join(BASE, "surgery_results"))
        if head.split(os.sep)[0] in ("test", "test2"):
            continue
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
            "surgery_min_ratio": float(ratio.min()) if len(ratio) else np.nan,
            "surgery_steps": int(frame["step"].max()) + 1,
            "surgery_status": str(frame["status"].iloc[-1]),
            "edges_removed": int(frame["removed"].max()),
            "edges_before": int(frame["edges"].iloc[0]),
            "edges_after": int(frame["edges"].iloc[-1]),
        })
    surg = pd.DataFrame(rows)
    if surg.empty:
        raise SystemExit("no surgery results found")
    # keep the longest run per (dataset, tol) so truncated scratch runs
    # never overwrite a completed one
    surg = (surg.sort_values("surgery_steps", ascending=False)
                .drop_duplicates(subset=["dataset", "tol"], keep="first")
                .sort_values(["tol", "dataset"]))

    base = base.drop_duplicates(subset=["dataset"], keep="first")
    merged = surg.merge(base, on="dataset", how="left")
    merged["fired"] = merged["edges_removed"] > 0
    merged["drop_factor"] = merged["no_surgery_ratio"] / merged["surgery_ratio"]

    lines = []
    lines.append("SURGERY vs NO-SURGERY  (paired, same numerics)")
    lines.append("=" * 62)
    for tol, grp in merged.groupby("tol"):
        lines.append("")
        lines.append("### surgery tol = {}".format(tol))
        lines.append("{:30s} {:>8s} {:>8s} {:>7s} {:>9s} {:>9s}".format(
            "dataset", "no_surg", "surgery", "removed", "drop x", "status"))
        for _, r in grp.sort_values("fired", ascending=False).iterrows():
            lines.append("{:30s} {:8.4f} {:8.4f} {:7d} {:9.1f} {:>9s}".format(
                str(r["dataset"])[:30],
                float(r["no_surgery_ratio"]), float(r["surgery_ratio"]),
                int(r["edges_removed"]),
                float(r["drop_factor"]) if np.isfinite(r["drop_factor"])
                else np.nan,
                str(r["surgery_status"])[:9]))
        fired = grp[grp["fired"] & np.isfinite(grp["drop_factor"])]
        lines.append("  fired: {}/{}   drop factor median {:.1f}x "
                     "(range {:.1f}..{:.1f})   surgery lower in {}/{}".format(
                         len(fired), len(grp),
                         float(fired["drop_factor"].median()) if len(fired)
                         else np.nan,
                         float(fired["drop_factor"].min()) if len(fired)
                         else np.nan,
                         float(fired["drop_factor"].max()) if len(fired)
                         else np.nan,
                         int((fired["surgery_ratio"]
                              < fired["no_surgery_ratio"]).sum())
                         if len(fired) else 0, len(fired)))

    fired = merged[merged["fired"] & np.isfinite(merged["drop_factor"])]
    lines.append("")
    lines.append("datasets where surgery FIRED (all tolerances): {}".format(
        len(fired)))
    if len(fired):
        lines.append("  drop factor (no-surgery / surgery): "
                     "median {:.1f}x, range {:.1f}x .. {:.1f}x".format(
                         float(fired["drop_factor"].median()),
                         float(fired["drop_factor"].min()),
                         float(fired["drop_factor"].max())))
        lines.append("  surgery arm lower in {}/{} datasets".format(
            int((fired["surgery_ratio"] < fired["no_surgery_ratio"]).sum()),
            len(fired)))
        lines.append("  edges removed: {} .. {} of {} .. {}".format(
            int(fired["edges_removed"].min()), int(fired["edges_removed"].max()),
            int(fired["edges_before"].min()), int(fired["edges_before"].max())))
    lines.append("")
    lines.append("datasets where surgery did NOT fire: {} "
                 "(excluded: arm is identical to no-surgery)".format(
                     int((~merged["fired"]).sum())))
    done = fired[fired["surgery_status"].isin(
        {"saturated", "constant_curvature", "horizon_reached"})]
    lines.append("")
    lines.append("--- completed surgery arms only (truncated runs excluded) ---")
    lines.append("{:30s} {:>8s} {:>8s} {:>7s} {:>9s}".format(
        "dataset", "no_surg", "surgery", "removed", "drop x"))
    for _, r in done.iterrows():
        lines.append("{:30s} {:8.4f} {:8.4f} {:7d} {:9.1f}".format(
            str(r["dataset"])[:30], float(r["no_surgery_ratio"]),
            float(r["surgery_ratio"]), int(r["edges_removed"]),
            float(r["drop_factor"])))
    if len(done):
        better = int((done["surgery_ratio"] < done["no_surgery_ratio"]).sum())
        lines.append("  {}/{} lower with surgery; median drop {:.1f}x "
                     "(range {:.1f}..{:.1f})".format(
                         better, len(done),
                         float(done["drop_factor"].median()),
                         float(done["drop_factor"].min()),
                         float(done["drop_factor"].max())))
        lines.append("  one-sided sign test p = {:.3f} (n={}) -- with n this "
                     "small the sign test cannot carry the claim alone;".format(
                         0.5 ** len(done), len(done)))
        lines.append("  the claim rests on the effect SIZE, not on the p value.")
    lines.append("")
    lines.append("--- reading ---")
    if len(fired) and float(fired["drop_factor"].median()) > 2:
        lines.append("  Answer B: the plateau is largely the cost of NOT doing")
        lines.append("  surgery. Rewrite the claim as: 'for the topology-")
        lines.append("  preserving flow, the residual measures how far the graph")
        lines.append("  is from being uniformisable WITHOUT changing its edge")
        lines.append("  set'. Cite Ma-Yang 2025 as the reason it vanishes with")
        lines.append("  surgery, and position the residual as the quantity that")
        lines.append("  surgery-based flows destroy.")
    else:
        lines.append("  Answer A: the residual survives surgery; it is a genuine")
        lines.append("  obstruction and sits against the convergence literature.")

    with open(os.path.join(BASE, "surgery_results", "SURGERY_REPORT.txt"),
              "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    merged.to_csv(os.path.join(BASE, "surgery_results", "surgery_paired.csv"),
                  index=False)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
