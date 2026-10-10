# -*- coding: utf-8 -*-
"""Insert the Zenodo DOI into Data/Code availability and cite the archive."""
import io, os, re

MS = r"D:\2026codes\2026821100data\paper_a\manuscript"
URL = "https://github.com/fanqimcse/ricci-flow-finite-time-residuals"
DOI_V = "10.5281/zenodo.23277162"
DOI_C = "10.5281/zenodo.23277161"


def sub(path, old, new, tag):
    p = os.path.join(MS, path)
    s = io.open(p, encoding="utf-8").read()
    r = re.compile(r"\s+".join(map(re.escape, old.split())))
    if not r.search(s):
        print("SKIP", tag); return
    s = r.sub(lambda m: new, s, count=1)
    io.open(p, "w", encoding="utf-8").write(s)
    print("ok  ", tag)


sub("manuscript.tex",
    """The frozen release is available at \\url{%s}
(tag \\texttt{v1.0.0}). An archival request has been submitted to {Software
Heritage} (save request \\texttt{2549855}, accepted); its permanent intrinsic
identifier is issued once the archival visit completes, and a {Zenodo} {DOI}
will be added when the {GitHub}--{Zenodo} integration is enabled.""" % URL,
    """The frozen release is available at \\url{%s}
(tag \\texttt{v1.0.0}) and is archived with a permanent {DOI},
\\url{https://doi.org/%s} (version) and
\\url{https://doi.org/%s} (concept, resolving to the latest version)
\\cite{Fan2026Code}. An archival request has also been submitted to {Software
Heritage} (save request \\texttt{2549855}, accepted), whose intrinsic
identifier is issued once the archival visit completes.""" % (URL, DOI_V, DOI_C),
    "code availability + DOI")

# Data availability: point at the same archive for the derived data
sub("manuscript.tex",
    """are provided with this submission
as Supplementary Data, together with Supplementary Tables~S1--S10.""",
    """are provided with this submission
as Supplementary Data, together with Supplementary Tables~S1--S10, and are
archived alongside the code at \\url{https://doi.org/%s}.""" % DOI_V,
    "data availability + DOI")

s = io.open(os.path.join(MS, "cover_letter.md"), encoding="utf-8").read()
s = s.replace("every one of its 44 references is\ncited in the text and resolves to a primary record.",
              "every one of its 45 references is\ncited in the text and resolves to a primary record.")
s = s.replace("""under the MIT licence at
https://github.com/fanqimcse/ricci-flow-finite-time-residuals and frozen at
`v1.0.0`; no number in the paper is hard-coded, and a single script
regenerates every figure and statistic from the raw result files, emitting a
log of the number behind each plotted point. An archival request has been
submitted to Software Heritage (save request 2549855, accepted).""",
              """under the MIT licence at
https://github.com/fanqimcse/ricci-flow-finite-time-residuals and frozen at
`v1.0.0`; it has a permanent Zenodo DOI, 10.5281/zenodo.23277162 (version) /
10.5281/zenodo.23277161 (concept). No number in the paper is hard-coded, and a
single script regenerates every figure and statistic from the raw result files,
emitting a log of the number behind each plotted point. An archival request has
also been submitted to Software Heritage (save request 2549855, accepted).""")
io.open(os.path.join(MS, "cover_letter.md"), "w", encoding="utf-8").write(s)
print("ok   cover letter")

p = os.path.join(MS, "SUBMISSION_CHECKLIST.md")
s = io.open(p, encoding="utf-8").read()
s = s.replace(
    "（551 文件 / 25 MB，tag `v1.0.0`，Release 已发布；Software Heritage 归档请求 id 2549855 已 accepted）。",
    "（551 文件 / 25 MB，tag `v1.0.0`，Release 已发布）。**Zenodo DOI 已签发**："
    "`10.5281/zenodo.23277162`（v1.0.0）/ `10.5281/zenodo.23277161`（concept，已写入正文）；"
    "Software Heritage 归档请求 id 2549855 已 accepted。")
io.open(p, "w", encoding="utf-8").write(s)
print("ok   checklist")
