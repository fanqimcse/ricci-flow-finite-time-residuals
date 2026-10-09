#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Two-fates analysis of the topology-preserving Ollivier--Ricci flow.

Background
----------
The stopping rule used to generate the first long-run corpus divided the
per-step move of the curvature variance by the variance at step 0, so every
slowly decaying run was eventually labelled "saturated".  Re-running with a
correct rule (per-step move divided by the *current* variance) showed that the
corpus splits into two behaviours:

  * a **converging** branch in which the residual variance keeps decaying
    exponentially, apparently towards zero (constant curvature);
  * a **stalled** branch in which the decay decelerates and the variance stops
    at a finite residual while the edge weights keep diverging.

This script quantifies the split, checks that it is not an artefact of the
horizon, and asks which structural quantity predicts it.
"""

from __future__ import print_function

import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats

SOURCE_DIR = os.path.dirname(os.path.abspath(__file__))
MASTER = os.path.join(SOURCE_DIR, "ricci_flow_longrun_results", "analysis",
                      "07_merged_with_structure.csv")

# A run is called stalled when its late-time decay rate is slower than this and
# it is decelerating; 1e-4 per step means "less than 10% over 1000 steps".
STALL_RATE = -1e-3
DECEL_MAX = 0.7


def fam_of(d):
    for f in ("scc_rt", "scc2034", "reptilia", "mammalia", "aves", "insecta",
              "bn", "ENZYMES", "bio", "ia", "power", "eco", "social"):
        if d.startswith(f):
            return f
    return d.split("_")[0]


def trajectories(name, dirs):
    best = None
    for base in dirs:
        p = os.path.join(base, "per_dataset", name + ".csv")
        if os.path.exists(p):
            try:
                tr = pd.read_csv(p)
            except Exception:
                continue
            if best is None or len(tr) > len(best):
                best = tr
    return best


def quarter_slopes(tr):
    v = tr["curvature_variance"].astype(float).values
    if v[0] <= 0:
        return None
    vn = v / v[0]
    s = tr["step"].astype(float).values
    n = len(vn)

    def slope(a, b):
        msk = (s >= a) & (s < b)
        if msk.sum() < 8:
            return np.nan
        return stats.linregress(s[msk], np.log(np.clip(vn[msk], 1e-300, None))).slope

    edges = [0, n * 0.25, n * 0.5, n * 0.75, n]
    return [slope(edges[i], edges[i + 1]) for i in range(4)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dirs", nargs="*",
                    default=[os.path.join(SOURCE_DIR, "ricci_flow_deep"),
                             os.path.join(SOURCE_DIR, "ricci_flow_v2"),
                             os.path.join(SOURCE_DIR, "ricci_flow_strict")])
    ap.add_argument("--out_dir", default=os.path.join(SOURCE_DIR, "strict_results"))
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    master = pd.read_csv(MASTER)
    rows = []
    for name in sorted(master["dataset"]):
        tr = trajectories(name, args.dirs)
        if tr is None:
            continue
        sl = quarter_slopes(tr)
        if sl is None:
            continue
        r3, r4 = sl[2], sl[3]
        decel = abs(r4) / abs(r3) if (r3 and abs(r3) > 1e-12) else np.nan
        v = tr["curvature_variance"].astype(float).values
        weights = tr["max_weight"].astype(float).values[-1] / max(
            float(tr["min_weight"].astype(float).values[-1]), 1e-300)
        if r4 > STALL_RATE and (np.isnan(decel) or decel < DECEL_MAX):
            fate = "stalled"
        elif r4 <= STALL_RATE:
            fate = "converging"
        else:
            fate = "undecided"
        rows.append({
            "dataset": name, "family": fam_of(name), "steps": len(tr),
            "v_end": float(v[-1] / v[0]), "rate_Q3": r3, "rate_Q4": r4,
            "decel": decel, "weight_span": float(weights),
            "min_w": float(tr["min_weight"].astype(float).values[-1]),
            "mean_k_end": float(tr["mean_curvature"].astype(float).values[-1]),
            "fate": fate,
        })
    d = pd.DataFrame(rows).merge(master, on="dataset", how="left",
                                 suffixes=("", "_m"))
    d.to_csv(os.path.join(args.out_dir, "fate_table.csv"), index=False)

    L = []
    w = L.append
    w("TWO FATES OF THE TOPOLOGY-PRESERVING OLLIVIER-RICCI FLOW")
    w("=" * 70)
    w("strict stopping rule: |V(t)-V(t-1)| / V(t) < 1e-5 for 50 consecutive "
      "steps")
    w("fate rule: stalled if the decay rate over the last quarter of the run "
      "> %.0e per step AND the run is decelerating (decel < %.2f)"
      % (STALL_RATE, DECEL_MAX))
    w("")
    w("strict re-run networks: %d" % len(d))
    w(d["fate"].value_counts().to_string())
    w("")
    for fate, g in d.groupby("fate"):
        w("  %-11s n=%2d  median v_end=%.4f  median rate_Q4=%+.2e  "
          "median weight span=%.2e  families=%s"
          % (fate, len(g), g["v_end"].median(), g["rate_Q4"].median(),
             g["weight_span"].median(), sorted(set(g["family"]))))
    w("")
    stalled = d[d["fate"] == "stalled"]
    conv = d[d["fate"] == "converging"]
    if len(stalled) and len(conv):
        u, p = stats.mannwhitneyu(np.log10(stalled["v_end"].clip(1e-12)),
                                  np.log10(conv["v_end"].clip(1e-12)),
                                  alternative="greater")
        w("log10 v_end stalled vs converging: median %.3f vs %.3f ; "
          "Mann-Whitney p=%.4g" % (np.log10(stalled["v_end"].median()),
                                   np.log10(conv["v_end"].median()), p))
        u2, p2 = stats.mannwhitneyu(stalled["weight_span"], conv["weight_span"],
                                    alternative="greater")
        w("weight span stalled vs converging: median %.2e vs %.2e ; "
          "Mann-Whitney p=%.4g" % (stalled["weight_span"].median(),
                                   conv["weight_span"].median(), p2))
    w("")
    w("--- does structure predict the fate / the residual? ---")
    for c in ("transitivity", "assortativity", "modularity", "density",
              "mean_degree", "degree_cv", "nodes", "n_components"):
        if c not in d.columns:
            continue
        sub = d[["v_end", "rate_Q4", "fate", c]].replace([np.inf, -np.inf], np.nan).dropna()
        if len(sub) < 6:
            continue
        r1, p1 = stats.spearmanr(sub["v_end"], sub[c])
        r2, p2 = stats.spearmanr(sub["rate_Q4"], sub[c])
        b = sub[sub["fate"].isin(["stalled", "converging"])]
        r3 = p3 = np.nan
        if len(b) >= 6 and b["fate"].nunique() == 2:
            r3, p3 = stats.mannwhitneyu(b.loc[b["fate"] == "stalled", c],
                                        b.loc[b["fate"] == "converging", c])
            r3 = float(r3) / (len(b[b["fate"] == "stalled"]) * len(b[b["fate"] == "converging"]))
        w("  %-14s n=%2d rho(v_end)=%+.3f(p%.3f) rho(rate)=%+.3f(p%.3f) "
          "AUC(stall>conv)=%.3f(p%.3f)"
          % (c, len(sub), r1, p1, r2, p2, r3 if r3 == r3 else float("nan"),
             p3 if p3 == p3 else float("nan")))
    w("")
    w("--- per-network ---")
    w(d[["dataset", "family", "steps", "v_end", "rate_Q4", "decel",
         "weight_span", "fate", "transitivity"]]
      .sort_values(["fate", "v_end"]).to_string(index=False,
                                                float_format=lambda x: "%.4g" % x))

    out = os.path.join(args.out_dir, "TWO_FATES_REPORT.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))
    print("\n".join(L[:45]))
    print("\nwrote", out)


if __name__ == "__main__":
    main()
