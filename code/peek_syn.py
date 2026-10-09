# -*- coding: utf-8 -*-
import os
import pandas as pd
import numpy as np

BASE = r"D:\2026codes\2026821100data"
frames = []
for i in range(4):
    p = os.path.join(BASE, "synthetic_results", "assort_base{}.csv".format(i))
    if os.path.exists(p):
        frames.append(pd.read_csv(p))
d = pd.concat(frames, ignore_index=True)
lines = []
lines.append("rows={}".format(len(d)))
lines.append("steps: median={:.0f} mean={:.1f} min={} max={}".format(
    d["steps"].median(), d["steps"].mean(), d["steps"].min(), d["steps"].max()))
lines.append(d["steps"].describe().to_string())
lines.append("")
lines.append("var_ratio: median={:.5f}".format(d["var_ratio"].median()))
lines.append("nodes={} edges median={:.0f}".format(
    d["nodes"].iloc[0], d["edges"].median()))
with open(os.path.join(BASE, "_synsteps.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines))
print("ok")
