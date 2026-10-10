#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Round-6 manuscript patch (second re-review of Communications Physics).

Applies, in order:
  1  abstract:  between-component term "need not vanish"; step-size claim softened
  2  intro:     same two downgrades
  3  results:   g10 reported alongside g105/g106
  4  results:   extended-horizon + step-size paragraphs rewritten
  5  discussion: between-component term evolves
  6  discussion: "not explained by the tested controls"; three networks
  7  methods:    initial-weight semantics flagged as a modelling choice
  8  methods:    delete duplicated "Sample sizes / Further numerical controls"
  9  methods:    delete duplicated "In controlled rewiring" subsection
 10  methods:    new weight-semantics sensitivity control
 11  methods:    numerical reliability at extreme weight span
 12  preamble:   numbered citations, single References heading
 13  back matter: Funding declaration, data-availability table range
"""

import io
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PATH = os.path.join(HERE, "manuscript.tex")

raw = io.open(PATH, "r", encoding="utf-8", newline="").read()
crlf = "\r\n" in raw
text = raw.replace("\r\n", "\n")

edits = []


def add(tag, old, new, count=1):
    edits.append((tag, old, new, count))


# ---------------------------------------------------------------- 1 --------
add("1 abstract", r"""disconnected graph the pooled variance contains a between-component term that
cannot decay; within components it falls on $34/34$ networks, median factor
$5.7$. Third, under exact transport the residual is step-size converged and
shows negligible decay over the tested $5{,}000$-step horizon;
fixed-$\varepsilon$ regularisation breaks scale invariance.""",
    r"""disconnected graph the pooled variance contains a between-component term that
need not vanish even when the within-component variance approaches zero;
within components it falls on $34/34$ networks, median factor
$5.7$. Third, under exact transport the residual is largely stable to
step-size variation at matched physical time for the tested non-converged
trajectories and shows negligible decay over the tested $5{,}000$-step
horizon; fixed-$\varepsilon$ regularisation breaks scale invariance.""")

# ---------------------------------------------------------------- 2 --------
add("2 intro", r"""pooled variance on a disconnected graph is dominated by a between-component
term that cannot decay, so a third of the corpus looks like a failure when it
is not. Fixed-$\varepsilon$ entropic regularisation breaks the scale
invariance that would license the standard renormalisation as a gauge fixing.
After all of these are controlled, two connected networks from a single family
still relax unusually slowly---a case study, not a demonstration that $V_G^*$
is positive.""",
    r"""pooled variance on a disconnected graph is dominated by a between-component
term that need not vanish, so a third of the corpus looks like a failure when
it is not. Fixed-$\varepsilon$ entropic regularisation breaks the scale
invariance that would license the standard renormalisation as a gauge fixing.
After all of these are controlled, connected networks from one single family
still relax unusually slowly---a case study, not a demonstration that $V_G^*$
is positive.""")

# ---------------------------------------------------------------- 3 --------
add("3 g10", r"""We treat these two as \emph{case studies}, not as evidence about connected
networks in general: two graphs from one family out of seventeen support no
prevalence claim.""",
    r"""We treat these two as \emph{case studies}, not as evidence about connected
networks in general: two graphs from one family out of seventeen support no
prevalence claim. A third connected graph of the same family,
\texttt{ENZYMES\_g10}, was classified Pattern~I at the $400$-step horizon of
the {Sinkhorn} re-run and is therefore not one of the two; under exact
transport it also stops moving, but later and from a higher residual
(Section~\emph{Extended horizon}, Supplementary Table~S11). We report all
three on equal terms there and count all three as one family.""")

# ---------------------------------------------------------------- 4 --------
add("4 horizon+stepsize", r"""\emph{Extended horizon.} Extending the three connected ENZYMES runs from
$2{,}000$ to $5{,}000$ steps under exact transport changes their residual by
$+0.009\%$ (\texttt{g105}), $+0.0002\%$ (\texttt{g106}) and $+0.17\%$
(\texttt{g10}), while the weight span grows to $2\times10^{71}$
(Supplementary Table~S7). Over three thousand further steps of an unregularised
solver these networks show no appreciable decay---evidence about the
trajectory, not about the infimum.

\emph{Step size at fixed physical time.} Holding physical time fixed at $T=20$
and varying $\Delta t$ ten-fold, from $0.01$ ($2{,}000$ steps) to $0.1$ ($200$
steps) under exact transport on $C_5$, the residual at $T$ changes by at most
$0.05\%$ on the seven networks that have not converged by then and by at most
$4.9\%$ on the two that have converged by then
(Supplementary Table~S7): not a time-discretisation artefact.""",
    r"""\emph{Extended horizon.} Extending the three connected {ENZYMES} runs from
$2{,}000$ to $5{,}000$ steps under exact transport changes their
within-component residual by $+0.009\%$ (\texttt{g105}), $+0.0002\%$
(\texttt{g106}) and $+0.17\%$ (\texttt{g10}); in absolute terms
$1.1\times10^{-6}$, $6.8\times10^{-8}$ and $8.7\times10^{-5}$, while the
weight span grows to $10^{65}$--$10^{72}$
(Supplementary Table~S11). All three therefore meet the same descriptive
near-plateau criterion used elsewhere in the paper---a relative per-step move
of the variance below $10^{-5}$ sustained for at least $50$ consecutive
steps---and all three sustain it for $3{,}900$--$4{,}200$ of their $5{,}000$
steps. We report the three together and count them as one family, not as three
independent observations.

They are not interchangeable, and the differences are the reason two of them
carry the case study. \texttt{g10} is still relaxing at step $1{,}000$: its
residual falls $14.3\%$ between steps $1{,}000$ and $1{,}500$ before
flattening, so its plateau begins later and at a value
($5.2\times10^{-2}$) that is not comparable with a $2{,}000$-step reading of
the other two. Its terminal drift is correspondingly larger---a median local
rate $r_t = -1.4\times10^{-7}$ per step over the final $500$ steps, against
$-1.3\times10^{-10}$ (\texttt{g105}) and $-2.4\times10^{-13}$
(\texttt{g106})---but still three orders of magnitude below the descriptive
$10^{-3}$ per step line. In all three the sign of $r_t$ at the end of the run
is negative: the variance is rising, not falling, and the residual at the
horizon lies above the minimum reached earlier in the same run by $0.5\%$
(\texttt{g10}), $2.8\%$ (\texttt{g105}) and $11\%$ (\texttt{g106}). A
finite-horizon residual is a value on a non-monotone trajectory, and none of
these numbers is a statement about the infimum.

\emph{Step size at fixed physical time.} Holding physical time fixed at $T=20$
and varying $\Delta t$ ten-fold, from $0.01$ ($2{,}000$ steps) to $0.1$ ($200$
steps) under exact transport on $C_5$, the pooled residual at $T$ changes by
at most $0.033\%$ on the seven networks that have not reached the operational
uniformisation threshold by then, and by up to $4.9\%$ on the two that have
(\texttt{aves-weaver-social-14}, $8.0\times10^{-6}\to7.8\times10^{-6}$, and
\texttt{scc\_rt\_http}, $1.06\times10^{-5}\to1.01\times10^{-5}$;
Supplementary Table~S7). Both of the larger changes are on residuals already
at the uniformisation threshold, so $4.9\%$ of $10^{-5}$ is $5\times10^{-7}$
in absolute units and affects no conclusion drawn from them. The same caveat
applies to the disconnected networks, whose within-component residual is
already $\sim10^{-7}$ and moves by up to $18.6\%$ in relative terms on
\texttt{scc\_rt\_oman} ($7.1\times10^{-8}\to5.8\times10^{-8}$), although by
only $1.3\times10^{-8}$ in absolute terms. What we can claim is therefore that
the residual is \emph{largely stable to step-size variation at matched physical
time for the tested non-converged trajectories}. That is weaker than
step-size convergence---one physical time, one idleness, one initial weighting
and nine networks were tested---and we do not claim the latter.""")

# ---------------------------------------------------------------- 5 --------
add("5 discussion 2", r"""uniformised within every component and still report a nonzero residual. Across
all $34$ disconnected networks""",
    r"""uniformised within every component and still report a nonzero residual. The
between-component term is a function of the component mean curvatures, which
themselves evolve, so we do not claim that it cannot decay---only that it need
not, and that on this corpus it does not. Across
all $34$ disconnected networks""")

# ---------------------------------------------------------------- 6 --------
add("6 discussion 3-4", r"""invariance, and the residual survives both. Nor is it a regularisation or discretisation artefact: at a matched
horizon it moves by $\le9.7\%$ over a $20$-fold $\varepsilon$ reduction and
by $\le0.05\%$ over a ten-fold $\Delta t$ range at $T=20$.

\emph{4. After these controls some connected networks still relax very slowly,
but the present evidence does not establish a positive geometric infimum.} Two
connected ENZYMES networks---two graphs from one family of seventeen---showed
no appreciable decay over $3{,}000$ further steps of exact transport, to a
$5{,}000$-step horizon, while their weight span grew to $10^{71}$. We treat it as a case study.""",
    r"""invariance, and the residual survives both. Nor is it explained by the
numerical controls we applied: at a matched horizon it moves by $\le9.7\%$
over a $20$-fold $\varepsilon$ reduction and by $\le4.9\%$ over a ten-fold
$\Delta t$ range at $T=20$, the larger figure occurring on the two networks
that had already reached the uniformisation threshold, where it amounts to
$5\times10^{-7}$ in absolute units. ``Not explained by the tested controls''
is what we claim; ``not a regularisation or discretisation artefact'' would
require exhausting every regularisation, conditioning and precision effect,
which we have not done.

\emph{4. After these controls some connected networks still relax very slowly,
but the present evidence does not establish a positive geometric infimum.}
Three connected {ENZYMES} networks---three graphs, all from one family of
seventeen---showed no appreciable decay over $3{,}000$ further steps of exact
transport, to a $5{,}000$-step horizon, changing their within-component
residual by $+0.009\%$, $+0.0002\%$ and $+0.17\%$
(Supplementary Table~S11), while their weight span grew to
$10^{65}$--$10^{72}$. At those weights the residual is insensitive to the
linear-programming backend, to rescaling the cost matrices out of that range
and to a $1\%$ perturbation of the initial weights (Methods). Two of the three
carry the full control battery and are the case study of Fig.~\ref{fig:8};
the third flattened later and from a higher residual and is reported on equal
terms. We treat all three as one case study.""")

# ---------------------------------------------------------------- 7 --------
add("7 initial weights", r"""from the production scheme is the update rule. Initial weights follow dataset metadata: where a dataset
supplies positive edge weights they are interpreted as lengths
(\texttt{weight\_mode = auto} with the \texttt{distance} convention), and
where the dataset is unweighted or carries non-positive weights---the majority
of the corpus---every edge is initialised to unit length. ``Weighted
network'' in this paper therefore refers to the \emph{evolving} weight
function, not to an assumption that the source data carry physical weights.
Runs use a step horizon and a wall-clock budget, whichever is reached first;
per-step trajectories are stored.""",
    r"""from the production scheme is the update rule. \emph{Initial weights.} Where
a dataset supplies finite, strictly positive edge values those values are used
as the initial edge lengths (\texttt{weight\_mode = auto} with the
\texttt{distance} convention); where a dataset is unweighted or carries
non-positive values---the majority of the corpus---every edge is initialised to
unit length. ``Weighted network'' in this paper therefore refers to the
\emph{evolving} weight function, not to an assumption that the source data
carry physical weights. Reading a positive third column as a length is a
\emph{modelling choice}: in public repositories such a value usually encodes an
affinity, a count, a flow or a timestamp, and we hold no verified data
dictionary for any of the $28$ weighted datasets. We therefore do not assert a
physical meaning for any of them; Supplementary Table~S12 lists every weighted
dataset with the observed range and type of its third column, and the control
below re-runs them under two alternative readings.
Runs use a step horizon and a wall-clock budget, whichever is reached first;
per-step trajectories are stored.""")

# ---------------------------------------------------------------- 8 --------
add("8 delete duplicate block", r"""\paragraph{Sample sizes.} \emph{Further numerical controls.} Re-running the flow with""",
    r"""@@DELETE@@\paragraph{Sample sizes.} \emph{Further numerical controls.} Re-running the flow with""")

# ---------------------------------------------------------------- 10 -------
add("10 weight semantics", r"""\paragraph{Sham flow.} Identical numerics, but at each step the vector of""",
    r"""\paragraph{Weight-semantics sensitivity.} Twenty-eight of the $117$
assembled datasets carry a finite, strictly positive third column and were
therefore initialised with it as the edge length. To test how much of the
result depends on that reading, each of them is re-run under three initial
conventions with identical numerics and a matched step horizon:
\emph{as-is} ($l_e = w_e$, the production choice), \emph{reciprocal}
($l_e = 1/w_e$, the reading appropriate if the value is an affinity, a count
or a flow) and \emph{unit} ($l_e = 1$, topology only). The three conventions
give the same topology and the same flow, and differ only in the initial
metric. For each network we compare the initial curvature variance $V_0$, the
residual at the matched horizon and the late-window decay rate. Two of the
twenty-eight (\texttt{bio-DM-LC}, whose third column is constant, and every
network whose third column is a year) have a coefficient of variation below
$10^{-3}$, so for them the three conventions are numerically the same input and
the control is vacuous; they are reported as such rather than counted. The
outcome is summarised in Results and given per dataset in Supplementary
Table~S12. No convention is designated correct: the control bounds how much of
the corpus result is a property of the data rather than of the reading we gave
them.

\paragraph{Sham flow.} Identical numerics, but at each step the vector of""")

# ---------------------------------------------------------------- 11 -------
add("11 extreme span", r"""\paragraph{What the extended-horizon runs do and do not establish.} Three
propositions are kept separate. \emph{Observed:} for a fixed exact-transport
algorithm, step size, idleness and initial weighting, the residual on two
connected networks changes by less than $0.2\%$ over $3{,}000$ further steps of
a $5{,}000$-step run.""",
    r"""\paragraph{Numerical reliability at extreme weight span.} The
$5{,}000$-step exact-transport runs end with edge weights spanning
$10^{65}$--$10^{72}$, so we report what the residual at that span does and does
not rest on (Supplementary Table~S13). \emph{Solvers.} Every transport problem
is solved by network simplex (\texttt{ot.emd} from POT $0.9.7$,
\texttt{numItermax}$=10^{6}$); the \texttt{scipy.linprog} {HiGHS} fallback
exists only for the case in which {POT} raises, and it was never invoked---$0$
of $9{,}120$ transport problems across the $60$ checkpoints of the three
extended runs fell back. Where we forced it, it agreed with the network
simplex to a maximum relative cost difference of $1.2\times10^{-15}$ while the
weight span was below $10^{5}$, degraded to $3\times10^{-9}$--$6\times10^{-8}$
at spans of $10^{6}$--$10^{7}$, and to $O(1)$--$O(10)$ on a minority of edges
beyond $10^{8}$ (the \emph{median} difference stays below $10^{-9}$
throughout). The fallback is not trustworthy at extreme span, and no reported
number depends on it. \emph{Conditioning.} Dividing each cost matrix by its own
mean before solving---which multiplies the exact $W_1$ cost by the same
factor and therefore leaves $\kappa_e$ unchanged in exact arithmetic---reproduces
every curvature to $\le2.2\times10^{-16}$ and the variance to a relative
$3.6\times10^{-16}$ at a span of $10^{26}$--$10^{29}$. The magnitude of the
costs, as opposed to their ratio structure, therefore costs nothing in double
precision, and the residual is not a floating-point artefact of the span.
\emph{Consistency checks.} At every checkpoint the curvature of two reference
graphs whose symmetry forces constant curvature ($C_6$ and $K_6$, unit
weights) is recomputed and required to have zero variance: $60$ of $60$
checkpoints passed. (This check is applied to the reference graphs, not to the
network under study, whose curvature is not constant at unit weights.)
\emph{The regularised solver.} At these spans the fixed-$\varepsilon=0.05$
{Sinkhorn} kernel underflows and the surrogate is not usable:
$\max_e|\kappa_e^{\rm sinkhorn}-\kappa_e^{\rm exact}|$ is of order
$10^{16}$--$10^{25}$. No {Sinkhorn} number is therefore reported beyond
$2{,}000$ steps, and the extended-horizon residual rests on exact transport
alone.

\paragraph{What the extended-horizon runs do and do not establish.} Three
propositions are kept separate. \emph{Observed:} for a fixed exact-transport
algorithm, step size, idleness and initial weighting, the residual on three
connected networks changes by less than $0.2\%$ over $3{,}000$ further steps of
a $5{,}000$-step run.""")

# ---------------------------------------------------------------- 12 -------
add("12 numbered citations", r"""\usepackage[round]{natbib}""",
    r"""\usepackage[numbers]{natbib}""")

add("12 single References heading", r"""\section*{References}
\bibliographystyle{unsrtnat}""",
    r"""\bibliographystyle{unsrtnat}""")

# ---------------------------------------------------------------- 13 -------
add("13 funding", r"""\section*{Author contributions}""",
    r"""\section*{Funding}
This work received no external funding.

\section*{Author contributions}""")

add("13 data availability tables", r"""as Supplementary Data, together with Supplementary Tables~S1--S10, and are
archived alongside the code at \url{https://doi.org/10.5281/zenodo.23277162}.""",
    r"""as Supplementary Data, together with Supplementary Tables~S1--S13, and are
archived alongside the code at \url{https://doi.org/10.5281/zenodo.23277162}.""")


# ---------------------------------------------------------------- apply -----
for tag, old, new, count in edits:
    n = text.count(old)
    if n != count:
        print("FAIL  %-28s found %d, expected %d" % (tag, n, count))
        sys.exit(1)
    if new.startswith("@@DELETE@@"):
        # deletion: remove from this anchor up to the next paragraph marker
        marker = new[len("@@DELETE@@"):]
        i = text.index(marker)
        # find the following \paragraph{Two experimental stages.}
        j = text.index(r"\paragraph{Two experimental stages.}", i)
        print("  delete block: %d chars" % (j - i))
        text = text[:i] + text[j:]
    else:
        text = text.replace(old, new, count)
    print("ok    %-28s" % tag)

# second duplication: the Methods copy of the transitivity subsection
marker = r"\subsection*{In controlled rewiring, transitivity changes the relaxation" + "\n" + r"time, not the endpoint}"
first = text.index(marker)
second = text.index(marker, first + 1)
endmark = r"\paragraph{Monotonicity of the variance (control).}"
j = text.index(endmark, second)
print("  delete duplicated Methods subsection: %d chars" % (j - second))
text = text[:second] + text[j:]

out = text.replace("\n", "\r\n") if crlf else text
io.open(PATH, "w", encoding="utf-8", newline="").write(out)
print("\nwritten:", PATH, len(out), "chars")
