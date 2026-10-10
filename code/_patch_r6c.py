#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Round-6c: word-budget pass (main text <= 5,000 words, abstract <= 150)."""

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


# --------------------------------------------------------------- abstract ---
OLD = r"""Discrete {Ricci} flow on a connected weighted graph converges to constant
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

NEW = r"""Discrete {Ricci} flow on a connected weighted graph converges to constant
curvature; a finite computation sees something else. Here we separate four
conflated quantities: the finite-horizon residual, its numerical counterpart,
the trend, and the geometric infimum $V_G^*$, which we do not compute. First,
the stopping rule used in empirical work normalises the per-step variance
change by its initial value, which any decaying run passes: all $28$ runs it
labelled saturated were still falling. Second, on a disconnected graph the
pooled variance contains a between-component term that need not vanish; within
components it falls on $34/34$ networks, median factor $5.7$. Third, under
exact transport the residual is largely stable to step-size variation and
barely decays over $5{,}000$ steps; fixed-$\varepsilon$ regularisation breaks
scale invariance. A finite-time residual is therefore a measurement, not a
convergence certificate: a sham control destroying the curvature--edge map
yields larger residuals ($8/8$ families), and structure does not predict the
residual across families."""
text = repl("abstract", OLD, NEW)

# ----------------------------------------------- analysis populations -------
OLD = r"""Two further distinctions organise the evidence. The first is an \emph{evidence tier}: a result is \emph{primary} when it
rests on at least eight independent source families and \emph{exploratory}
when it rests on fewer than five, in which case its $p$-value is an upper bound
on the evidence rather than a measurement of it. The sham control, the
corpus-level component decomposition and the main variance contraction are
primary; the surgery control ($2$ families), the transitivity manipulation
($2$ base graphs) and the assortativity and community sweeps ($4$ base graphs)
are exploratory and are labelled as such."""
NEW = r"""Two further distinctions organise the evidence. The first is an
\emph{evidence tier}: a result is \emph{primary} when it rests on at least
eight independent source families and \emph{exploratory} when it rests on fewer
than five, in which case its $p$-value is an upper bound on the evidence rather
than a measurement of it. The sham control, the corpus-level component
decomposition and the main variance contraction are primary; the surgery
control ($2$ families), the transitivity manipulation ($2$ base graphs) and the
assortativity and community sweeps ($4$ base graphs) are exploratory."""
text = repl("populations", OLD, NEW)

# --------------------------------------------------- pooled estimator -------
OLD = r"""Over the corpus rather than the strict sample the picture changes
quantitatively. Across all $34$ disconnected networks the within-component
variance falls in $34/34$, but by a median factor of only $5.7$
(interquartile range $2.2$--$10$), exceeding $10^{2}$ on two networks, with
final values spanning $1.9\times10^{-32}$ to $8.1\times10^{-2}$
(Supplementary Table~S3); the between-component term accounts for a median
$65\%$ of the pooled variance and more than half of it on $22$ of $34$. The
direction is universal; the magnitude on which we earlier based the
claim---three to four orders of magnitude---is the behaviour of the two most
extreme networks, not the typical one, and we withdraw it."""
NEW = r"""Across all $34$ disconnected networks the within-component variance falls in
$34/34$, but by a median factor of only $5.7$ (interquartile range $2.2$--$10$,
exceeding $10^{2}$ on two networks), with final values spanning
$1.9\times10^{-32}$ to $8.1\times10^{-2}$ (Supplementary Table~S3); the
between-component term accounts for a median $65\%$ of the pooled variance and
more than half of it on $22$ of $34$. The direction is universal; the magnitude
on which we earlier based the claim---three to four orders of magnitude---is the
behaviour of the two most extreme networks, not the typical one, and we
withdraw it."""
text = repl("pooled corpus", OLD, NEW)

OLD = r"""A graph can therefore be perfectly uniformised within every component and still
report a nonzero pooled residual. 
"""
NEW = r"""A graph can therefore be perfectly uniformised within every component and still
report a nonzero pooled residual.
"""
text = repl("whitespace", OLD, NEW)

# ---------------------------------------------------- scale invariance ------
OLD = r"""\emph{Scale invariance.} For the exact curvature
$\kappa(G,cw)=\kappa(G,w)$, which is what would license the per-step
renormalisation as removal of a null direction. For the entropy-regularised
cost this fails: with $\varepsilon$ fixed, $C\mapsto cC$ changes the kernel
$\exp(-C/\varepsilon)$ and hence the plan. We measured the violation on nine
networks and six rescaling factors (Methods): under exact $W_1$ and under a
scale-\emph{equivariant} schedule $\varepsilon_t\propto s_t$ the curvature is
invariant to machine precision, while under the fixed-$\varepsilon=0.05$
solver used elsewhere in the paper $\max_e|\kappa_e(cw)-\kappa_e(w)|$ reaches
$1.40$ at the initial weights and $3.86$ after $150$ steps. The claim in our
previous version---that the renormalisation removes a null direction---is
therefore false as stated, and we withdraw it; $\sum_e\log w_e=0$ in the
independent optimisation below is a normalisation of the search space, not a
gauge fixing."""
NEW = r"""\emph{Scale invariance.} For the exact curvature
$\kappa(G,cw)=\kappa(G,w)$, which is what would license the per-step
renormalisation as removal of a null direction; for the entropy-regularised
cost this fails, because with $\varepsilon$ fixed $C\mapsto cC$ changes the
kernel $\exp(-C/\varepsilon)$. Measured on nine networks and six rescaling
factors (Methods), the violation reaches
$\max_e|\kappa_e(cw)-\kappa_e(w)|=1.40$ at the initial weights and $3.86$
after $150$ steps, while exact $W_1$ and a scale-\emph{equivariant} schedule
$\varepsilon_t\propto s_t$ are invariant to machine precision. The claim in our
previous version---that the renormalisation removes a null direction---is false
as stated and we withdraw it; $\sum_e\log w_e=0$ below is a normalisation of
the search space, not a gauge fixing."""
text = repl("scale invariance", OLD, NEW)

OLD = r"""Two things bound the damage: $s_t$ drifts by only $0.98$--$1.29\times$ over
$200$ steps because the mean weight is renormalised every step, so the
violation does not compound; and repairing it barely changes the numbers. Under
the scale-equivariant schedule at a matched horizon of $200$ steps the pooled
residual moves by a median $-4.6\%$ (at most $-14\%$) with no change of
ordering, and replacing regularised transport by exact transport
\emph{raises} the residual on the two connected networks ($+5\%$ on
\texttt{g105}, $+193\%$ on \texttt{g106})---slower relaxation, not a lower
floor. The conclusion survives the repair; only the interpretation had to go."""
NEW = r"""Two things bound the damage: $s_t$ drifts by only $0.98$--$1.29\times$ over
$200$ steps, so the violation does not compound; and repairing it barely
changes the numbers. Under the scale-equivariant schedule at a matched horizon
of $200$ steps the pooled residual moves by a median $-4.6\%$ (at most $-14\%$)
with no change of ordering, and replacing regularised transport by exact
transport \emph{raises} the residual on the two connected networks ($+5\%$ on
\texttt{g105}, $+193\%$ on \texttt{g106})---slower relaxation, not a lower
floor."""
text = repl("bound damage", OLD, NEW)

# ----------------------------------------------- regularisation strength -----
OLD = r"""\emph{Regularisation strength.} Entropic smoothing shrinks curvature
magnitudes, and could hide or invent a spread. Re-run at a \emph{matched}
horizon---the stopping step itself moves with $\varepsilon$---over a
$20$-fold reduction of $\varepsilon$ with periodic exact-transport
checkpoints, the exact-transport residual spans $5.1\%$ on \texttt{g105}
and $9.7\%$ on \texttt{g106}. Where the estimators disagree, {Sinkhorn} is
\emph{below} exact, by $3.6\times$ on \texttt{g106}, so the artefact hides
the residual; at $\varepsilon=0.002$ the solve fails and we exclude those
runs (Supplementary Table~S9)."""
NEW = r"""\emph{Regularisation strength.} Entropic smoothing shrinks curvature
magnitudes, and could hide or invent a spread. Re-run at a \emph{matched}
horizon over a $20$-fold reduction of $\varepsilon$ with periodic
exact-transport checkpoints, the exact-transport residual spans $5.1\%$ on
\texttt{g105} and $9.7\%$ on \texttt{g106}; where the estimators disagree,
{Sinkhorn} is \emph{below} exact ($3.6\times$ on \texttt{g106}$\!\!$), so
regularisation hides rather than creates the residual. At
$\varepsilon=0.002$ the solve fails and those runs are excluded
(Supplementary Table~S9)."""
text = repl("reg strength", OLD, NEW)

# ------------------------------------------------- relaxation patterns ------
OLD = r"""Because a rate estimated in one window says nothing about the next, we also
report the rate trajectory. Figure~\ref{fig:8} plots the local relaxation rate
$r_t=-\log(V_t/V_{t-1})$ for all $16$ runs of $C_4$: $r_t>0$ is decay,
$r_t\approx0$ a plateau, $r_t<0$ a rebound. The patterns appear as two
bands; the $10^{-3}$ per step line between them is a descriptive threshold,
not a theoretical boundary: runs cross it, so it carries no information about
the limit."""
NEW = r"""Because a rate estimated in one window says nothing about the next,
Fig.~\ref{fig:8} plots the local rate $r_t=-\log(V_t/V_{t-1})$ for all $16$
runs of $C_4$: $r_t>0$ is decay, $r_t\approx0$ a plateau, $r_t<0$ a rebound.
The patterns appear as two bands; the $10^{-3}$ per step line between them is a
descriptive threshold, not a theoretical boundary---runs cross it."""
text = repl("patterns", OLD, NEW)

out = text.replace("\n", "\r\n") if crlf else text
io.open(PATH, "w", encoding="utf-8", newline="").write(out)
print("written", len(out))
