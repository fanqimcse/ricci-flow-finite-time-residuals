# -*- coding: utf-8 -*-
import os
import pandas as pd

BASE = r"D:\2026codes\2026821100data"
OUT_DIR = os.path.join(BASE, "ricci_flow_longrun_results", "analysis")
lines = []
for f in sorted(os.listdir(OUT_DIR)):
    if not f.endswith(".csv"):
        continue
    p = os.path.join(OUT_DIR, f)
    try:
        df = pd.read_csv(p)
    except Exception as exc:
        lines.append("{}: unreadable ({})".format(f, exc))
        continue
    lines.append("{}: rows={} cols={}".format(f, len(df), list(df.columns)))
    lines.append("")

with open(os.path.join(BASE, "_cols.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines))
print("ok")
