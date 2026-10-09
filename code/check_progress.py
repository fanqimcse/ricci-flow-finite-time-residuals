"""Print progress of the synthetic sweep workers + rough ETA."""

import glob
import os
import time

BASE = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.join(BASE, "_syn_b{}.txt")
CSV = os.path.join(BASE, "synthetic_results", "assort_base{}.csv")
OUT = os.path.join(BASE, "_progress.txt")

lines = []
total_done = 0
total_expected = 60
for i in range(4):
    last = ""
    if os.path.exists(LOG.format(i)):
        # PowerShell's '*>' redirect writes UTF-16LE; fall back cleanly
        raw = open(LOG.format(i), "rb").read()
        for enc in ("utf-16", "utf-8"):
            try:
                text = raw.decode(enc)
                break
            except Exception:
                text = raw.decode("utf-8", errors="replace")
        for row in text.splitlines():
            if row.strip():
                last = row.strip()
    last = "".join(ch for ch in last if ch.isprintable())
    rows = 0
    if os.path.exists(CSV.format(i)):
        with open(CSV.format(i), encoding="utf-8", errors="replace") as fh:
            rows = max(sum(1 for _ in fh) - 1, 0)
    total_done += rows
    lines.append("worker {}: csv_rows={}  last={}".format(i, rows, last[:110]))

lines.append("")
lines.append("total completed graphs: {} / {}".format(total_done, total_expected))
lines.append("checked at {}".format(
    time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

with open(OUT, "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines))
print("\n".join(lines))
