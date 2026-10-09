# -*- coding: utf-8 -*-
"""Remove the fabricated GitHub URL (user 'fanqimai' does not exist) from ms + cover letter."""
import io, re

MS = r"D:\2026codes\2026821100data\paper_a\manuscript\manuscript.tex"
CL = r"D:\2026codes\2026821100data\paper_a\manuscript\cover_letter.md"


def wsrep(path, old, new, tag):
    s = io.open(path, encoding="utf-8").read()
    r = re.compile(r"\s+".join(map(re.escape, old.split())))
    if not r.search(s):
        print("SKIP", tag)
        return
    s = r.sub(lambda m: new, s, count=1)
    io.open(path, "w", encoding="utf-8").write(s)
    print("ok  ", tag)


NEW_SEC = r"""\section*{Code availability}
All analysis code is released under the {MIT} licence and accompanies this
submission as a single archive (Supplementary Data), frozen at version
\texttt{v1.0.0}. The archive contains the flow implementation, every
experiment and analysis script, the machine-readable configuration
(\texttt{config.yaml}, which records the unsafe legacy stopping rule as well,
so that the artefact itself can be reproduced), the network manifest, the
fixed random seeds, the raw result files, the supplementary tables and the
figure-generation script. No number in the paper is hard-coded: every figure
and statistic is regenerated from the raw result files by
\texttt{make\_paper\_figures.py}, which also emits a log giving the number
behind each plotted point. The archive is available from the author on
request and will be deposited in a public version-controlled repository with
a permanent identifier on publication."""

wsrep(MS, r"""\section*{Code availability}
All analysis code is released under the {MIT} licence at
\url{https://github.com/fanqimai/ricci-flow-finite-time-residuals} and frozen
at version \texttt{v1.0.0} for this submission. The release contains the flow
implementation, every experiment and analysis script, the machine-readable
configuration (\texttt{config.yaml}, which records the unsafe legacy stopping
rule as well, so that the artefact itself can be reproduced), the network
manifest, the fixed random seeds, the raw result files, the supplementary
tables and the figure-generation script. No number in the paper is hard-coded:
every figure and statistic is regenerated from the raw result files by
\texttt{make\_paper\_figures.py}, which also emits a log giving the number
behind each plotted point. The release is additionally archived by {Software
Heritage}, which assigns it a permanent intrinsic identifier; a {Zenodo} {DOI}
will be added if the author enables the {GitHub}--{Zenodo} integration.""",
      NEW_SEC, "code availability")

wsrep(MS, """provided with this submission
as Supplementary Data and in the repository named below, together with
Supplementary Tables~S1--S10.""",
      """provided with this submission
as Supplementary Data, together with Supplementary Tables~S1--S10.""",
      "data availability cross-ref")

# cover letter
s = io.open(CL, encoding="utf-8").read()
s = s.replace("https://github.com/fanqimai/ricci-flow-finite-time-residuals",
              "https://github.com/<account>/ricci-flow-finite-time-residuals")
io.open(CL, "w", encoding="utf-8").write(s)
print("ok   cover letter url -> placeholder (cover letter is not published)")
