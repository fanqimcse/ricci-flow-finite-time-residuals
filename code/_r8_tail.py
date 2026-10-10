# -*- coding: utf-8 -*-
"""
Round 8, reviewer point #2.

The manuscript says all 28 runs labelled `saturated` by the OLD rule were
still falling; Fig. 2d prints "%d/%d still falling" from a stricter
threshold (net change over the final quarter < -2%), and its caption says
"all but one".  Recompute both statistics from the raw per-step
trajectories so the paper can quote one number per definition.

Definitions:
  S1  last-step decrease       : V[-1] < V[-2]
  S2  final-quarter net change : V[-1]/V[k] - 1 < 0 ,  k = floor(0.75*n)
  S2' same but stricter (< -2%)
  S3  monotone non-increase over the final quarter: every step <= previous
"""
import os
import glob
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))          # 2026821100data
MASTER = os.path.join(ROOT, "ricci_flow_longrun_results", "analysis",
                      "07_merged_with_structure.csv")
PER_DIR = os.path.join(ROOT, "ricci_flow_longrun_results", "per_dataset")

master = pd.read_csv(MASTER)
sat = master[master["status"].astype(str) == "saturated"]
print("master rows:", len(master), "| status == saturated:", len(sat))

rows = []
skip = []
for fp in sorted(glob.glob(os.path.join(PER_DIR, "*.csv"))):
    name = os.path.basename(fp)[:-4]
    r = master.loc[master["dataset"] == name]
    if len(r) == 0 or str(r["status"].iloc[0]) != "saturated":
        continue
    try:
        tr = pd.read_csv(fp)
    except Exception as e:
        skip.append((name, "read"))
        continue
    v = tr["curvature_variance"].astype(float).values
    if len(v) < 40:
        skip.append((name, "len%d" % len(v)))
        continue
    if v[0] <= 0:
        skip.append((name, "v0<=0"))
        continue
    v = v / v[0]
    k = int(len(v) * 0.75)
    q = v[k:]
    rows.append(dict(
        dataset=name,
        family=r["family"].iloc[0] if "family" in r.columns else "",
        n_steps=len(v),
        V_end=v[-1],
        last_step_down=bool(v[-1] < v[-2]),
        last_step_rel=(v[-1] / v[-2] - 1.0),
        q4_net=(v[-1] / v[k] - 1.0),
        q4_monotone=bool(np.all(np.diff(q) <= 0)),
        q4_any_up=bool(np.any(np.diff(q) > 0)),
        min_in_q4=v[k:].min(),
        end_over_q4min=(v[-1] / v[k:].min()),
    ))

d = pd.DataFrame(rows)
print("usable trajectories:", len(d))
print("skipped:", skip)
print()
d["q4_net_lt0"] = d["q4_net"] < 0
d["q4_net_lt_m2"] = d["q4_net"] < -0.02

print("S1 last step decreased          : %d / %d" % (int(d.last_step_down.sum()), len(d)))
print("S2 final-quarter net change < 0 : %d / %d" % (int(d.q4_net_lt0.sum()), len(d)))
print("S2' same, stricter than -2%%     : %d / %d" % (int(d.q4_net_lt_m2.sum()), len(d)))
print("S3 monotone over final quarter  : %d / %d" % (int(d.q4_monotone.sum()), len(d)))
print()
print("median final-quarter net change : %+.2f%%" % (100 * d.q4_net.median()))
print("worst (most negative)           : %+.2f%%" % (100 * d.q4_net.min()))
print("most positive                   : %+.2f%%" % (100 * d.q4_net.max()))
print("median last-step relative move  : %+.3e" % d.last_step_rel.median())
print()
print(d[["dataset", "n_steps", "V_end", "last_step_rel", "q4_net",
         "q4_monotone"]].to_string(float_format=lambda x: "%.4g" % x))

out = os.path.join(HERE, "Supplementary_Table_S17_old_rule_trajectory_audit.csv")
d.to_csv(out, index=False)
print("\nwrote", out)
