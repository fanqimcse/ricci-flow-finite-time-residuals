#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Round-6f: subheading lengths + reference audit."""

import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PATH = os.path.join(HERE, "manuscript.tex")
raw = io.open(PATH, "r", encoding="utf-8", newline="").read()
crlf = "\r\n" in raw
text = raw.replace("\r\n", "\n")

for old, new in [
    (r"\subsection*{Relaxing the edge-set constraint largely removes the residual}",
     r"\subsection*{Relaxing the edge-set constraint}"),
    (r"\subsection*{No evidence that the initial spectral gap predicts early decay}",
     r"\subsection*{No evidence that the spectral gap predicts decay}"),
]:
    assert text.count(old) == 1, old
    text = text.replace(old, new)
print("subheadings fixed")

out = text.replace("\n", "\r\n") if crlf else text
io.open(PATH, "w", encoding="utf-8", newline="").write(out)

# ------------------------------------------------------------- reference audit
tex = io.open(PATH, encoding="utf-8").read().replace("\r\n", "\n")
tex = re.sub(r"(?<!\\)%.*", "", tex)
cited = set()
for m in re.finditer(r"\\cite[a-z]*\{([^}]*)\}", tex):
    for k in m.group(1).split(","):
        cited.add(k.strip())

bib = io.open(os.path.join(HERE, "references.bib"), encoding="utf-8").read()
keys = re.findall(r"@\w+\{([^,]+),", bib)
dup = [k for k in set(keys) if keys.count(k) > 1]
missing = sorted(cited - set(keys))
orphan = sorted(set(keys) - cited)
print("bib entries      :", len(keys))
print("duplicate keys   :", dup)
print("cited-but-missing:", missing)
print("uncited entries  :", orphan)
titles = re.findall(r"\btitle\s*=\s*\{(.*?)\}", bib, re.S)
print("titles parsed    :", len(titles))
if missing or dup:
    sys.exit(1)
