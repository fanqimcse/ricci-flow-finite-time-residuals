# -*- coding: utf-8 -*-
"""Round 7 (reviewer point #1).

Recompute, from the per-step trajectories of the 16 strict re-runs
(ricci_flow_v2/per_dataset, tol=1e-5, tol_mode=current, patience=50,
max_steps=400, budget=700 s), whether the operational stopping criterion
could ever have fired:

    rel_move_t = |V_t - V_{t-1}| / V_t  <  1e-5   for 50 consecutive steps

and, if it did not, how large the tolerance would have had to be.

Produces:
  Supplementary_Table_S15_stopping_rule_diagnostic.csv
  Supplementary_Table_S16_descriptive_plateau_three_networks.csv  (exact transport)
and prints the numbers used in the manuscript.
"""
import glob
import json
import os

import numpy as np
import pandas as pd

ROOT = r"D:\2026codes\2026821100data"
V2 = os.path.join(ROOT, "ricci_flow_v2")
OUT = os.path.join(ROOT, "paper_a", "manuscript")
TOL = 1e-5
PATIENCE = 50


def streaks(mask):
    """longest run of True, and the run of True that ends the series."""
    best = cur = 0
    for m in mask:
        cur = cur + 1 if m else 0
        best = max(best, cur)
    tail = 0
    for m in reversed(mask):
        if m:
            tail += 1
        else:
            break
    return best, tail


def diagnose(v, name, planned_horizon=None, wall=None):
    v = np.asarray(v, dtype=float)
    n = len(v) - 1                       # number of flow steps recorded
    dv = np.diff(v)
    rel = np.abs(dv) / np.maximum(v[1:], 1e-300)
    r = -np.log(v[1:] / v[:-1])
    mask = rel < TOL
    best, tail = streaks(mask)
    last50 = rel[-50:] if n >= 50 else rel
    r50 = r[-50:] if n >= 50 else r
    # smallest tolerance at which the rule would fire at least once in the
    # terminal 50 steps == the largest per-step relative move there
    tol_needed_last50 = float(np.max(last50))
    tol_needed_anywhere = float(np.min(rel))
    fired = bool(tail >= PATIENCE)
    return dict(
        dataset=name,
        n_steps=n,
        planned_horizon=planned_horizon if planned_horizon is not None else n,
        wall_seconds=wall,
        V_end=float(v[-1]),
        V_end_over_V0=float(v[-1] / v[0]) if v[0] else np.nan,
        n_steps_rel_lt_tol=int(mask.sum()),
        frac_steps_rel_lt_tol=float(mask.mean()),
        longest_streak_rel_lt_tol=int(best),
        trailing_streak_rel_lt_tol=int(tail),
        rule_fired_patience50=fired,
        rel_move_min=float(np.min(rel)),
        rel_move_median=float(np.median(rel)),
        rel_move_last50_median=float(np.median(last50)),
        rel_move_last50_max=tol_needed_last50,
        tol_that_would_fire_last50=tol_needed_last50,
        r_median_last50=float(np.median(r50)),
        dV_rel_last50=float((v[-1] - v[-51]) / v[-51]) if n >= 51 else np.nan,
    )


def main():
    prog = {}
    pj = os.path.join(V2, "batch_progress.json")
    if os.path.exists(pj):
        prog = json.load(open(pj))

    rows = []
    for fp in sorted(glob.glob(os.path.join(V2, "per_dataset", "*.csv"))):
        ds = os.path.basename(fp)[:-4]
        tr = pd.read_csv(fp)
        v = tr["curvature_variance"].astype(float).values
        info = prog.get(ds, {})
        rows.append(diagnose(
            v, ds,
            planned_horizon=int(info.get("steps", len(v) - 1)),
            wall=float(info.get("wall", np.nan)) if info else np.nan,
        ))
    df = pd.DataFrame(rows).sort_values("dataset").reset_index(drop=True)
    df.to_csv(os.path.join(
        OUT, "Supplementary_Table_S15_stopping_rule_diagnostic.csv"), index=False)

    fates = pd.read_csv(os.path.join(ROOT, "strict_results", "fate_table.csv"))
    fm = dict(zip(fates["dataset"], fates["fate"]))
    df.insert(1, "pattern", [fm.get(d, "?") for d in df.dataset])

    print("=" * 78)
    print("C4 strict re-runs (n=%d): did the operational rule fire?" % len(df))
    print("=" * 78)
    print(df[["dataset", "pattern", "n_steps", "rule_fired_patience50",
              "longest_streak_rel_lt_tol", "frac_steps_rel_lt_tol",
              "rel_move_last50_median", "rel_move_last50_max",
              "r_median_last50"]].to_string(index=False))

    n_fired = int(df.rule_fired_patience50.sum())
    print()
    print("runs in which the rule fired (>=50 consecutive steps): %d / %d"
          % (n_fired, len(df)))
    print("runs with ANY step below the tolerance: %d / %d"
          % (int((df.n_steps_rel_lt_tol > 0).sum()), len(df)))
    print("median over runs of (max rel move in last 50 steps): %.3g"
          % float(np.median(df.rel_move_last50_max)))
    print("range of that quantity: %.3g .. %.3g"
          % (float(df.rel_move_last50_max.min()), float(df.rel_move_last50_max.max())))
    print("median over runs of (median rel move in last 50 steps): %.3g"
          % float(np.median(df.rel_move_last50_median)))
    print("median r_t last50: %.3g ; range %.3g .. %.3g"
          % (float(np.median(df.r_median_last50)),
             float(df.r_median_last50.min()), float(df.r_median_last50.max())))
    print("n runs terminal median rate > 1e-4: %d"
          % int((df.r_median_last50 > 1e-4).sum()))
    print("n runs terminal median rate < 0: %d"
          % int((df.r_median_last50 < 0).sum()))
    print("median |dV| over last 50 steps: %.3g ; max %.3g"
          % (float(np.median(np.abs(df.dV_rel_last50))),
             float(np.max(np.abs(df.dV_rel_last50)))))
    print("total wall clock recorded in batch_progress: %.1f s = %.2f h"
          % (float(np.nansum(df.wall_seconds)),
             float(np.nansum(df.wall_seconds)) / 3600.0))

    # ---- the three 5,000-step exact-transport runs, for contrast ----------
    print()
    print("=" * 78)
    print("C5 three 5,000-step exact-transport runs (the only runs that DO meet it)")
    print("=" * 78)
    rows2 = []
    for fp in sorted(glob.glob(os.path.join(
            ROOT, "exp_round4", "traj_long", "*__exact5000.csv"))):
        ds = os.path.basename(fp).split("__")[0]
        tr = pd.read_csv(fp)
        col = "V_within" if "V_within" in tr.columns else "curvature_variance"
        v = tr[col].astype(float).values
        rows2.append(diagnose(v, ds))
    df2 = pd.DataFrame(rows2).sort_values("dataset").reset_index(drop=True)
    df2.to_csv(os.path.join(
        OUT, "Supplementary_Table_S16_descriptive_plateau_three_networks.csv"),
        index=False)
    print(df2[["dataset", "n_steps", "rule_fired_patience50",
               "longest_streak_rel_lt_tol", "frac_steps_rel_lt_tol",
               "rel_move_last50_max", "r_median_last50"]].to_string(index=False))

    # ---- ratio of tolerances -------------------------------------------
    print()
    med_c4 = float(np.median(df.rel_move_last50_max))
    med_c5 = float(np.median(df2.rel_move_last50_max))
    print("median tol needed by C4 (400-step Sinkhorn)      : %.3g" % med_c4)
    print("median max rel move C5 (5,000-step exact)        : %.3g" % med_c5)
    print("ratio                                            : %.3g x" % (med_c4 / med_c5))


if __name__ == "__main__":
    main()
