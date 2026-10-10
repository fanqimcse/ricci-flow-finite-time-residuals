#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Round-6e: last 30 words."""

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


text = repl("table1 caption", r"""\caption{\textbf{Pooled versus within-component curvature variance.} For every
disconnected network of $C_4$ the within-component variance falls by one to
three orders of magnitude while the pooled variance does not: the flow has
uniformised each component and the estimator reports their disagreement. The
two connected networks at the bottom are shown for contrast---for them the
decomposition is trivial. Horizons differ because each run had a wall-clock
budget, so the connected residuals are finite-horizon values, not bounds.}""",
            r"""\caption{\textbf{Pooled versus within-component curvature variance.} For every
disconnected network of $C_4$ the within-component variance falls by one to
three orders of magnitude while the pooled variance does not. The two connected
networks are shown for contrast. Horizons differ because each run had a
wall-clock budget, so the connected residuals are finite-horizon values.}""")

text = repl("g10 a", r"""They are not interchangeable, and the differences are the reason two of them
carry the case study.""",
            r"""They are not interchangeable, which is why two of them carry the case
study.""")

text = repl("g10 b", r"""A
finite-horizon residual is a value on a non-monotone trajectory, and none of
these numbers is a statement about the infimum.""",
            r"""A
finite-horizon residual is a value on a non-monotone trajectory, not a
statement about the infimum.""")

text = repl("practitioners", r"""For practitioners: stop when the variance stops moving relative to
\emph{itself}; report the within-component variance whenever the graph is
disconnected; with a fixed-$\varepsilon$ solver report $\varepsilon/s_t$ or use
a scale-equivariant schedule; treat neither the residual as predictable across
families nor networks from one source as independent replicates. The design
limitations that bound every statement above are listed in Methods.""",
            r"""For practitioners: stop when the variance stops moving relative to
\emph{itself}; report the within-component variance whenever the graph is
disconnected; with a fixed-$\varepsilon$ solver report $\varepsilon/s_t$ or use
a scale-equivariant schedule; treat neither the residual as predictable across
families nor networks from one source as independent replicates. The
limitations bounding every statement above are listed in Methods.""")

out = text.replace("\n", "\r\n") if crlf else text
io.open(PATH, "w", encoding="utf-8", newline="").write(out)
print("written", len(out))
