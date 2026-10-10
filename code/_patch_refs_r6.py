"""Round-6b: fix the four bibliographic errors and the two over-stated claims
found by the three-way reference audit.  Idempotent: every edit asserts its
anchor is found exactly once; re-running reports 0 replacements."""
import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BIB = os.path.join(HERE, "references.bib")
TEX = os.path.join(HERE, "manuscript.tex")

edits = 0


def sub1(text, old, new, tag):
    global edits
    n = text.count(old)
    if n != 1:
        raise SystemExit("ANCHOR FAIL (%s): found %d, expected 1\n%r"
                         % (tag, n, old[:120]))
    edits += 1
    return text.replace(old, new)


def load(p):
    t = io.open(p, encoding="utf-8", newline="").read()
    return t.replace("\r\n", "\n")


def save(p, t):
    io.open(p, "w", encoding="utf-8", newline="").write(t.replace("\n", "\r\n"))


# ---------------------------------------------------------------- bib ----
bib = load(BIB)

bib = sub1(bib,
           "  volume  = {27}, number = {2}, pages = {343--348}, year = {1967}",
           "  volume  = {21}, number = {2}, pages = {343--348}, year = {1967}",
           "SinkhornKnopp1967 volume 27->21")

bib = sub1(bib,
           "  volume  = {17}, number = {2}, pages = {345--358}, year = {2010}",
           "  volume  = {17}, number = {2}, pages = {343--356}, year = {2010}",
           "LinYau2010 pages 345-358 -> 343-356")

bib = sub1(bib,
           "  doi     = {10.1090/proc/13404}",
           "  doi     = {10.1090/proc/13145}",
           "BauerChungLinLiu2017 DOI")

# Feydy2019: AISTATS 2019 = PMLR 89:2681-2690, not ICML 2019 = PMLR 97:1958-1967
bib = sub1(
    bib,
    "  booktitle = {Proceedings of the 36th International Conference on Machine\n"
    "               Learning (ICML)},\n"
    "  volume    = {97}, pages = {1958--1967}, year = {2019}",
    "  booktitle = {Proceedings of the 22nd International Conference on Artificial\n"
    "               Intelligence and Statistics (AISTATS)},\n"
    "  volume    = {89}, pages = {2681--2690}, year = {2019}",
    "Feydy2019 venue/volume/pages")

save(BIB, bib)

# ---------------------------------------------------------------- tex ----
tex = load(TEX)

# (1) Samal et al. (2018) report that the Ollivier and Forman discretisations
#     are *highly correlated* across many networks; "disagree" misreads it.
#     LinYau2010 is added as the Bakry-Emery source (the four original
#     citations do not treat Bakry-Emery at all).
tex = sub1(
    tex,
    "the Bakry--\\'Emery, Ollivier and Forman discretisations disagree on the same\n"
    "graph \\cite{BauerChungLinLiu2017,Sreejith2016,Samal2018,Azarhooshang2020},",
    "the {Bakry--\\'Emery}, Ollivier and Forman discretisations give different\n"
    "values on the same graph\n"
    "\\cite{LinYau2010,BauerChungLinLiu2017,Sreejith2016,Samal2018,Azarhooshang2020},",
    "claim: discretisations disagree")

# (2) Ma-Yang's surgery deletes edges; nothing in that paper adds them.
tex = sub1(tex,
           "flows with \\emph{surgery}---deleting or adding edges---converge",
           "flows with \\emph{surgery}---deleting edges---converge",
           "claim: surgery adds edges")

save(TEX, tex)

print("replacements applied: %d" % edits)
print("bib ->", BIB)
print("tex ->", TEX)
sys.exit(0)
