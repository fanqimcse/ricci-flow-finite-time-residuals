#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Summary statistics for the round-6 weight-semantics control (S12)."""
from __future__ import print_function

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
S12 = pd.read_csv(os.path.join(HERE, "..", "paper_a", "manuscript",
                               "Supplementary_Table_S12_weight_semantics.csv"))

raw = pd.read_csv(os.path.join(HERE, "exp_round6",
                               "raw_weight_column_stats.csv"))
if "raw_cv" not in S12.columns:
    S12 = S12.merge(raw[["dataset", "raw_cv", "year_like"]], on="dataset",
                    how="left")

done = S12.dropna(subset=["V0_asis", "V0_recip", "V0_unit"]).copy()
print("datasets with all three conventions: %d of %d"
      % (len(done), len(S12)))

v0 = done[["V0_asis", "V0_recip", "V0_unit"]].to_numpy(dtype=float)
done["v0_spread"] = v0.max(axis=1) / v0.min(axis=1)
done["v0_rel_rng"] = (v0.max(axis=1) - v0.min(axis=1)) / v0.min(axis=1)

degen = done[done["raw_cv"] < 1e-3]
var = done[done["raw_cv"] >= 1e-3]
print("\nnear-constant third column (cv < 1e-3): n=%d" % len(degen))
print("  max V0 spread across conventions: %.4f (%.3f%%)"
      % (degen["v0_spread"].max(), 100 * degen["v0_rel_rng"].max()))
print("  " + ", ".join(sorted(degen["dataset"])))
print("\nvariable third column (cv >= 1e-3): n=%d" % len(var))
print("  V0 spread across conventions: median %.2f x, max %.2f x"
      % (var["v0_spread"].median(), var["v0_spread"].max()))
print("  " + ", ".join(sorted(var["dataset"])))

print("\nper-dataset V0 spread (variable-column datasets):")
for _, r in var.sort_values("v0_spread", ascending=False).iterrows():
    print("  %-38s cv=%6.3f  V0 spread=%.2f x" % (r["dataset"], r["raw_cv"],
                                                  r["v0_spread"]))

# residual ratio
for c in ("asis", "recip", "unit"):
    col = "ratio_within_%s" % c
    if col in done.columns:
        v = done[col].dropna()
        print("\n%-6s n=%d  ratio<1 in %d/%d  median=%.3f  min=%.4f  max=%.3f"
              % (c, len(v), int((v < 1).sum()), len(v), v.median(), v.min(),
                 v.max()))

allr = []
for c in ("asis", "recip", "unit"):
    col = "ratio_within_%s" % c
    if col in done.columns:
        allr.extend(done[col].dropna().tolist())
allr = np.asarray(allr, dtype=float)
print("\npooled over convention x dataset: n=%d, ratio<1 in %d, max=%.4f"
      % (allr.size, int((allr < 1).sum()), allr.max()))

# how many of the 28 are documented as a distance
print("\ndocumented as distance: %d of %d"
      % (int((S12["documented_as_distance"] == "yes").sum()), len(S12)))
print(S12["documented_semantics"].value_counts().to_string())
print("\nyear-like third column: %d datasets" % int(S12["year_like"].sum()))
print("cv < 1e-3: %d datasets" % int((S12["raw_cv"] < 1e-3).sum()))
print("in analysis set: %d" % int(S12["in_analysis_set"].sum()))
