# -*- coding: utf-8 -*-
import os

D = r"D:\2026codes\2026821100data"
lines = []
for w in ("w0", "w1"):
    p = os.path.join(D, "surgery_results", w, "progress.log")
    lines.append("=== {} ===".format(w))
    if os.path.exists(p):
        with open(p, encoding="utf-8", errors="replace") as fh:
            lines.extend(ln.rstrip() for ln in fh if ln.strip())
    else:
        lines.append("(no progress.log)")
    pd_dir = os.path.join(D, "surgery_results", w, "per_dataset")
    if os.path.isdir(pd_dir):
        lines.append("files: " + ", ".join(sorted(os.listdir(pd_dir))))
    else:
        lines.append("files: none")
    lines.append("")

with open(os.path.join(D, "_surgprog.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines))
print("ok")
