# -*- coding: utf-8 -*-
"""Finalise the repository URL: manuscript wording + cover letter + checklist."""
import io, os, re

MS = r"D:\2026codes\2026821100data\paper_a\manuscript"
URL = "https://github.com/fanqimcse/ricci-flow-finite-time-residuals"


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
    """The frozen release is available at
\\url{https://github.com/fanqimcse/ricci-flow-finite-time-residuals}; it is additionally archived by {Software
Heritage}, which assigns it a permanent intrinsic identifier, and a {Zenodo} {DOI} will be added once the
{GitHub}--{Zenodo} integration is enabled.""",
    """The frozen release is available at \\url{%s}
(tag \\texttt{v1.0.0}). An archival request has been submitted to {Software
Heritage} (save request \\texttt{2549855}, accepted); its permanent intrinsic
identifier is issued once the archival visit completes, and a {Zenodo} {DOI}
will be added when the {GitHub}--{Zenodo} integration is enabled.""" % URL,
    "code availability wording")

sub("cover_letter.md",
    """**Data and code.** All networks are from the Network Data Repository
(http://networkrepository.com), distributed openly for research use; we follow
its acknowledgement policy and ship the derived data --- manifest, trajectories
and intermediate results --- with the submission. All analysis code is released
under the MIT licence and accompanies the submission as a single archive,
frozen at `v1.0.0`; no number in the paper is hard-coded, and a single script
regenerates every figure and statistic from the raw result files, emitting a
log of the number behind each plotted point. The public repository address and
its permanent identifier will be supplied on request and at proof stage.""",
    """**Data and code.** All networks are from the Network Data Repository
(http://networkrepository.com), distributed openly for research use; we follow
its acknowledgement policy and ship the derived data --- manifest, trajectories
and intermediate results --- with the submission. All analysis code is released
under the MIT licence at %s, frozen at `v1.0.0`; no number
in the paper is hard-coded, and a single script regenerates every figure and
statistic from the raw result files, emitting a log of the number behind each
plotted point. An archival request has been submitted to Software Heritage
(save request 2549855, accepted).""" % URL,
    "cover letter")

s = io.open(os.path.join(MS, "SUBMISSION_CHECKLIST.md"), encoding="utf-8").read()
s = s.replace("https://github.com/<account>/ricci-flow-finite-time-residuals", URL)
s = s.replace("https://github.com/fanqimai/ricci-flow-finite-time-residuals", URL)
io.open(os.path.join(MS, "SUBMISSION_CHECKLIST.md"), "w", encoding="utf-8").write(s)
print("ok   checklist")
