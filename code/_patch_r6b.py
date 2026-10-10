#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Round-6b: compliance pass.

  * main text (Introduction+Results+Discussion) to <= 5,000 words
    -> Discussion rewritten to remove overlap with Results (the journal
       checklist asks for no overlap) and the detailed limitations list moved
       to Methods (Methods does not count towards the limit);
  * abstract to <= 150 words;
  * Results subheadings shortened to < 60 characters each;
  * duplicate-free Methods gains the limitations list and the weight-semantics
    result.
"""

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
        print("FAIL %-30s found %d expected %d" % (tag, n, count))
        sys.exit(1)
    print("ok   %-30s" % tag)
    return text.replace(old, new, count)


# ------------------------------------------------------------------ abstract
OLD_ABS = r"""Discrete {Ricci} flow on a connected weighted graph converges to constant
curvature; what a finite computation sees differs. We separate four quantities
empirical work conflates: the finite-horizon residual $V(T)$, its numerical
counterpart, the trend of $V(t)$, and the geometric infimum $V_G^*$, which we
do not compute. First, the stopping rule used in empirical work normalises the
per-step change of variance by its initial value, which any decaying run
passes: all $28$ runs it labelled saturated were still falling. Second, on a
disconnected graph the pooled variance contains a between-component term that
need not vanish even when the within-component variance approaches zero;
within components it falls on $34/34$ networks, median factor
$5.7$. Third, under exact transport the residual is largely stable to
step-size variation at matched physical time for the tested non-converged
trajectories and shows negligible decay over the tested $5{,}000$-step
horizon; fixed-$\varepsilon$ regularisation breaks scale invariance. Destroying the
curvature--edge correspondence in a sham control is associated with larger
residuals ($8/8$ families); the structure-to-residual map does not transfer."""

NEW_ABS = r"""Discrete {Ricci} flow on a connected weighted graph converges to constant
curvature; a finite computation sees something else. Here we separate four
quantities that empirical work conflates: the finite-horizon residual, its
numerical counterpart, the trend, and the geometric infimum $V_G^*$, which we
do not compute. First, the stopping rule used in empirical work normalises the
per-step variance change by its initial value, which any decaying run passes:
all $28$ runs it labelled saturated were still falling. Second, on a
disconnected graph the pooled variance contains a between-component term that
need not vanish even when the within-component variance approaches zero; within
components it falls on $34/34$ networks, median factor $5.7$. Third, under
exact transport the residual is largely stable to step-size variation at
matched physical time and barely decays over $5{,}000$ steps;
fixed-$\varepsilon$ regularisation breaks scale invariance. We conclude that a
finite-time residual is a measurement, not a convergence certificate: destroying
the curvature--edge correspondence in a sham control is associated with larger
residuals ($8/8$ families), and structure does not predict the residual across
families."""

text = repl("abstract", OLD_ABS, NEW_ABS)

# ------------------------------------------------------- Results subheadings
for old, new in [
    (r"""\subsection*{Analysis populations, evidence tiers, and four quantities that
are routinely conflated}""",
     r"""\subsection*{Analysis populations and four conflated quantities}"""),
    (r"""\subsection*{The stopping criterion used in previous work is an operational
diagnostic, not a convergence test}""",
     r"""\subsection*{The stopping criterion is not a convergence test}"""),
    (r"""\subsection*{The regularised curvature is not scale invariant---and the
residual survives that being fixed}""",
     r"""\subsection*{Regularisation breaks scale invariance}"""),
    (r"""\subsection*{The plateau is associated with the curvature--edge
correspondence}""",
     r"""\subsection*{The plateau tracks the curvature--edge map}"""),
    (r"""\subsection*{In controlled rewiring, transitivity changes the relaxation
time, not the endpoint}""",
     r"""\subsection*{Transitivity changes the relaxation time}"""),
    (r"""\subsection*{The structure-to-residual mapping does not generalise across
families}""",
     r"""\subsection*{The structure-to-residual map does not transfer}"""),
]:
    text = repl("subheading", old, new)

# ------------------------------------------------------- Discussion rewrite
start = text.index("Four statements survive the controls above.")
end = text.index(r"\section*{Methods}")
old_len = len(text[start:end])

NEW_DISC = r"""Four statements survive the controls above, in order of evidential
strength, each stated at the level the evidence supports.

\emph{1. A finite-time residual is not a convergence certificate.} The stopping
rule used in empirical work divides the per-step change of the variance by its
initial value; Eq.~(\ref{eq:stoprate}) shows that this ratio decays like
$e^{-\lambda t}$ whatever the limit is, so the rule reports elapsed time. All
$28$ runs it labelled saturated were still falling. The \emph{operational
stopping criterion} used instead says only that the variance has stopped moving
relative to itself over the observed window; every residual here is reported
with its horizon.

\emph{2. On a disconnected graph the pooled curvature variance can
substantially overstate within-component non-uniformity.} Eq.~(\ref{eq:union})
exhibits the case $V_{\mathrm{within}}=0$ with $V_{\mathrm{pooled}}>0$: a graph
can be perfectly uniformised within every component and still report a nonzero
residual. The between-component term is a function of the component mean
curvatures, which themselves evolve, so we do not claim that it cannot
decay---only that it need not, and that on this corpus it does not. Across all
$34$ disconnected networks of the analysis set the within-component variance
falls, by a median factor of $5.7$; our earlier claim of three to four orders of
magnitude on every disconnected network is withdrawn. A third of the corpus is
disconnected and the pooled variance is the standard reporting choice, so the
correction applies to previous empirical work as much as to ours.

\emph{3. Entropic regularisation changes the numerical curvature dynamics
unless its scale is controlled.} With $\varepsilon$ held fixed the regularised
curvature is not invariant under $w\mapsto cw$, so the reading of the per-step
renormalisation as a gauge fixing is unavailable and we withdraw it. Exact
transport and a scale-equivariant schedule restore invariance, and the residual
survives both: at a matched horizon it moves by $\le9.7\%$ over a $20$-fold
$\varepsilon$ reduction and by $\le4.9\%$ over a ten-fold $\Delta t$ range at
$T=20$, the larger figure occurring on residuals already at the uniformisation
threshold. We claim that the residual is not explained by the tested controls;
we do not claim immunity to every regularisation, conditioning or precision
effect.

\emph{4. After these controls some connected networks still relax very slowly,
but the evidence does not establish a positive geometric infimum.} Three
connected {ENZYMES} networks---three graphs from one family of
seventeen---changed their within-component residual by $+0.009\%$, $+0.0002\%$
and $+0.17\%$ over $3{,}000$ further steps of exact transport, to a
$5{,}000$-step horizon, while their weight span grew to $10^{65}$--$10^{72}$;
at those weights the residual survives every numerical diagnostic we applied
(Methods). Two of the three carry the full control battery and are the case
study of Fig.~\ref{fig:8}; the third flattened later and from a higher residual
and is reported on equal terms. We treat all three as one case study, have not
computed $V_G^*$, and note that a metastable transient longer than any horizon
we could afford fits the data equally well.

Homogenisation itself is not in doubt and is tied to the curvature--edge
correspondence: destroying it in the sham control is associated with suppressed
homogenisation in $8/8$ families---under one intervention, not a proof of
necessity. What remains unpredicted is which real networks relax quickly: the
structure-to-residual map does not transfer across families, transitivity moves
the relaxation time and not the endpoint, and a rate measured in one window
carries no information about the next, which is why Fig.~\ref{fig:8} reports
the rate trajectory rather than a classification.

Table~\ref{tab:assumptions} sets out, hypothesis by hypothesis, why the
convergence theorem most relevant to us does not transfer to the scheme studied
here---a statement about the hypotheses, not a criticism of the theorem. Where
the hypotheses hold we observe what the theorem predicts, and the remaining
differences---exponential versus Euler update, regularised versus exact
transport---we tested directly rather than leaving them open (Methods); neither
produces the residual.

\begin{table}[t]
\centering
\begin{tabular}{p{2.6cm}p{3.0cm}p{3.0cm}}
\toprule
Hypothesis / ingredient & Li~\&~M\"unch \cite{LiMunch2026} & Scheme used here \\
\midrule
Graph & finite, weighted, connected (irreducible) & finite, weighted; $34$ of $93$ disconnected \\
Curvature & exact {Ollivier} $W_1$ & entropy-regularised {Sinkhorn}, $\varepsilon=0.05$ \\
Update & $d_{n+1}=(1-\alpha\kappa)d_n$ & $w_e\leftarrow w_e e^{-\Delta t\,\kappa_e}$, $\Delta t=0.05$ \\
Renormalisation & not part of the statement & mean weight renormalised every step \\
Step size & admissible range & fixed; $\alpha=0.5$ and $0.05$ both tested \\
Edge set & fixed, no surgery & fixed; surgery arm adds a control only \\
\bottomrule
\end{tabular}
\caption{\textbf{Why the convergence theorem does not transfer.} Each row is an
ingredient on which the scheme studied differs from the one proved to
converge.}
\label{tab:assumptions}
\end{table}

For practitioners: stop when the variance stops moving relative to
\emph{itself}; report the within-component variance whenever the graph is
disconnected; with a fixed-$\varepsilon$ solver report $\varepsilon/s_t$ or use
a scale-equivariant schedule; treat neither the residual as predictable across
families nor networks from one source as independent replicates. The design
limitations that bound every statement above are listed in Methods.
"""

text = text[:start] + NEW_DISC + text[end:]
print("ok   %-30s (%d -> %d chars)" % ("discussion", old_len, len(NEW_DISC)))

# --------------------------------------------------- Methods: limitations
anchor = r"\paragraph{Computational budget.} All runs were executed on an 8-core"
LIMITS = r"""\paragraph{Limitations of the design.} Five limitations bound every
statement above. \emph{(i) Estimator.} $V_{\mathrm{pooled}}$ and
$V_{\mathrm{within}}$ are not interchangeable on a disconnected graph; both are
reported for every network, and the within-component re-run reached a median of
only $86$ steps, so those residuals are themselves finite-horizon quantities.
\emph{(ii) Transport.} The production curvature is a {Sinkhorn} surrogate at
fixed $\varepsilon=0.05$ and idleness $\alpha=0.5$; the extreme weight ratios
reached in Pattern~II make the transport subproblem ill conditioned and put the
surrogate out of range altogether (the kernel underflows), so the long-horizon
residuals rest on exact transport alone, and we have tested but not eliminated
a regularisation or precision origin for the residual. \emph{(iii) Initial
weights.} A positive third column read as a length is a modelling choice, not a
verified data dictionary; the sensitivity control above bounds, but does not
remove, its influence. \emph{(iv) Independence.} The unit of independence is
the source family, and within the synthetic sweeps the base graph---two of
which is not a basis for a general claim. \emph{(v) Attribution.} The
manipulations hold the degree sequence fixed but not motif structure, shortest
paths or spectra, and the assortativity and community sweeps inherit the unsafe
stopping rule, so their magnitudes carry horizon bias. Population $C_4$ covers
$16$ networks because a step costs $0.03$--$12$~s and the corpus would need
$\sim80$~h to take all $93$ to $400$ steps: the sample is cost-stratified, not
random, and biased towards small graphs. A lower bound on $V_G^*$ in terms of a
conflict between curvature-feasible intervals is conjectured, not proven.

"""
text = repl("methods limitations", anchor, LIMITS + anchor)

out = text.replace("\n", "\r\n") if crlf else text
io.open(PATH, "w", encoding="utf-8", newline="").write(out)
print("written", len(out))
