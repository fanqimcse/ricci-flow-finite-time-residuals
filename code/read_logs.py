# -*- coding: utf-8 -*-
import os

D = r"D:\2026codes\2026821100data"
OUT = []
for i in range(4):
    p = os.path.join(D, "_surg_w{}.txt".format(i))
    OUT.append("=== _surg_w{}.txt ({} bytes) ===".format(i, os.path.getsize(p)))
    raw = open(p, "rb").read()
    t = None
    for enc in ("utf-8", "utf-16", "utf-16-le", "gbk", "latin-1"):
        try:
            t = raw.decode(enc)
            break
        except Exception:
            continue
    if t is None:
        t = raw.decode("utf-8", errors="replace")
    lines = [ln.rstrip() for ln in t.splitlines() if ln.strip()]
    OUT.extend(lines[-30:])
    OUT.append("")

with open(os.path.join(D, "_surg_logs.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT))
print("ok")
