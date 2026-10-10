#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Round-6d: final word-budget trim (main text <= 5,000, abstract <= 150)."""

import io
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PATH = os.path.join(HERE, "manuscript.tex")
raw = io.open(PATH, "r", encoding="utf-8", newline="").read()
crlf = "\r\n" in raw
text = raw.replace("\r\n", "\n")


def repl(tag, old, new, count=1):
    n = text.count(old)
    if n != count:
        print("FAIL %-28s found %d expected %d" % (tag, n, count))
        sys.exit(1)
    print("ok   %-28s" % tag)
    return text.replace(old, new, count)


# abstract: 155 -> <=150
text = repl("abstract a", r"""curvature; a finite computation sees something else. Here we separate four
conflated quantities: the finite-horizon residual, its numerical counterpart,
the trend, and the geometric infimum $V_G^*$, which we do not compute. First,""",
            r"""curvature; finite computations differ. Here we separate four
conflated quantities: the finite-horizon residual, its numerical counterpart,
the trend, and the geometric infimum $V_G^*$ (not computed here). First,""")

text = repl("abstract b", r"""residual is therefore a measurement, not a
convergence certificate: a sham control destroying the curvature--edge map
yields larger residuals ($8/8$ families), and structure does not predict the
residual across families.""",
            r"""residual is a measurement, not a
convergence certificate: a sham control destroying the curvature--edge map
yields larger residuals ($8/8$ families); structure does not predict the
residual across families.""")

# stopping-criterion subsection
text = repl("stop 1", r"""The reason is worth stating precisely. Suppose $V_t=V_\infty+Ae^{-\lambda t}$.
Then""", r"""Suppose $V_t=V_\infty+Ae^{-\lambda t}$. Then""")

text = repl("stop 2", r"""the numerators are the same and both decay, but the denominators are not. With
$V_0$ fixed the ratio decays like $e^{-\lambda t}$ whatever $V_\infty$ is, so
the rule reports elapsed time. With $V_t$ in the denominator the ratio tends to
$1-e^{-\lambda}$ when $V_\infty=0$---it does not tend to zero---and to zero
like $e^{-\lambda t}/V_\infty$ when $V_\infty>0$. The corrected rule is
triggered by runs that are settling onto a positive value and is not triggered
by runs still decaying towards zero. That is the diagnostic we want, and it is
still not a convergence theorem: it says the variance has stopped moving
relative to itself over the observed window, nothing about the limit.""",
            r"""the numerators are the same and both decay, but the denominators are not.
With $V_0$ fixed the ratio decays like $e^{-\lambda t}$ whatever $V_\infty$ is,
so the rule reports elapsed time; with $V_t$ in the denominator it tends to
$1-e^{-\lambda}$ when $V_\infty=0$---not to zero---and to zero like
$e^{-\lambda t}/V_\infty$ when $V_\infty>0$. The corrected rule fires on runs
settling onto a positive value and not on runs still decaying towards zero.
That is the diagnostic we want, and it is still not a convergence theorem: it
says the variance has stopped moving relative to itself over the observed
window, nothing about the limit.""")

# analysis populations
text = repl("pop 1", r"""The second distinction is the one that organises the paper. A number produced
by a finite computation can be any of four things, and empirical work reports
them interchangeably:""",
            r"""The second distinction organises the paper. A number produced by a finite
computation can be any of four things, and empirical work reports them
interchangeably:""")

out = text.replace("\n", "\r\n") if crlf else text
io.open(PATH, "w", encoding="utf-8", newline="").write(out)
print("written", len(out))
