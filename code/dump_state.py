# -*- coding: utf-8 -*-
"""Dump surgery results + synthetic progress into one readable text file."""
import os
import glob

D = r"D:\2026codes\2026821100data"
OUT = []


def read_tail(path, n=25):
    if not os.path.exists(path):
        return ["(missing)"]
    raw = open(path, "rb").read()
    for enc in ("utf-8", "utf-16", "utf-8-sig", "gbk", "latin-1"):
        try:
            t = raw.decode(enc)
            break
        except Exception:
            continue
    else:
        t = raw.decode("utf-8", errors="replace")
    lines = [ln.rstrip() for ln in t.splitlines() if ln.strip()]
    return lines[-n:]


OUT.append("=== surgery_results tree ===")
sr = os.path.join(D, "surgery_results")
for root, dirs, files in os.walk(sr):
    for f in sorted(files):
        p = os.path.join(root, f)
        OUT.append("  {}  ({} bytes)".format(os.path.relpath(p, sr), os.path.getsize(p)))

OUT.append("")
OUT.append("=== surgery per_dataset contents ===")
pd = os.path.join(sr, "per_dataset")
if os.path.isdir(pd):
    for f in sorted(os.listdir(pd)):
        p = os.path.join(pd, f)
        OUT.append("--- {} ---".format(f))
        OUT.extend(read_tail(p, 8))

OUT.append("")
OUT.append("=== surgery logs ===")
for name in ("_surg_a.txt", "_surg_b.txt", "_surg_a2.txt", "_surg_debug.txt", "_dbg_run.txt"):
    p = os.path.join(D, name)
    if os.path.exists(p):
        OUT.append("--- {} ---".format(name))
        OUT.extend(read_tail(p, 12))

OUT.append("")
OUT.append("=== synthetic csv row counts ===")
syn = os.path.join(D, "synthetic_results")
for i in range(4):
    c = os.path.join(syn, "assort_base{}.csv".format(i))
    if os.path.exists(c):
        with open(c, encoding="utf-8", errors="replace") as fh:
            n = sum(1 for _ in fh) - 1
        OUT.append("  base{} rows={}".format(i, n))
    else:
        OUT.append("  base{} missing".format(i))

with open(os.path.join(D, "_state.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT))
print("ok")
