# Cover Letter

**Manuscript type:** Article
**Title:** *Finite-time residuals in discrete Ollivier–Ricci flow: a measurement
framework, its controls, and its limits*
**Author:** Qi Fan (Huaibei Normal University, Huaibei 235000, China;
e-mail 8073592@qq.com; ORCID 0000-0003-3591-3363)

Dear Editors,

We submit the above manuscript for consideration as an Article in
*Communications Physics*. This is a substantially restructured revision: the
paper no longer asks "which networks retain a residual?", because that
question cannot be answered by a finite computation. It asks instead what a
finite computation can and cannot establish, and it answers by
distinguishing four quantities that empirical work on this flow has
routinely conflated, and by measuring the artefact associated with each.

**The four quantities.** (L1) $V(T)$, the curvature variance at the horizon
$T$ that a run actually reaches; (L2) $V_{\rm num}(T;\theta)$, the same
quantity as produced by the solver, which differs from L1 by discretisation
and regularisation error; (L3) the trend of $V(t)$, the only quantity for
which "still falling" and "plateaued" are well-defined; and (L4)
$V_G^*=\inf_w \mathrm{Var}[\kappa(G,w)]$, the geometric infimum over all
positive weightings. We compute L1–L3 and explicitly do **not** compute L4.
Nothing in the manuscript is claimed about L4, and every statement about a
"residual" is indexed to a named level.

**Finding 1 — the stopping rule is an operational diagnostic, not a
convergence test.** The rule used in empirical work divides the per-step
change of the variance by its *initial* value. Writing $V_t=V_\infty+Ae^{-\lambda t}$,
the ratio tends to $1-e^{-\lambda}$ if $V_\infty=0$ but to $0$ like
$e^{-\lambda t}/V_\infty$ if $V_\infty>0$: the rule is satisfied by *any*
exponentially decaying run, whatever its limit. All 28 runs it labelled
saturated in our corpus were still falling, by a median of 24.7% over the
final quarter. We relabel the two observed behaviours Pattern I (continued
relaxation) and Pattern II (near-plateau over the observed horizon), and we
report them as trajectory properties, not as asymptotic fates.

**Finding 2 — the estimator.** On a disconnected graph the pooled
edge-curvature variance contains a between-component term that cannot decay
under any weight update. Decomposing it exactly, the within-component
variance falls on 34/34 networks (median factor 5.7, up to $10^{30}$), while
the pooled value does not. We prove a disjoint-union identity that isolates
this term exactly and verify it numerically to machine precision. This also
resolves an apparent conflict with the convergence theorem of Li & Münch
(*J. Funct. Anal.* 290, 111367, 2026), whose connectedness hypothesis our
disconnected networks violate: their conclusion holds precisely where their
hypotheses hold, component by component.

**Finding 3 — the solver.** The entropy-regularised Sinkhorn curvature used
in practice is not invariant under $w\mapsto cw$, because a fixed
$\varepsilon$ does not scale with the cost matrix. We measure the violation
(curvatures move by up to 3.86, against $4\times10^{-16}$ for exact $W_1$)
and we isolate a scale-equivariant schedule $\varepsilon_t\propto s_t$ that
restores it. Two new experiments close the remaining gap. A full
$\varepsilon$-trajectory sweep with exact-transport checkpoints, read at a
*matched* horizon rather than each run's own stopping step, shows that the
exact-transport residual moves by only 5.1% and 9.7% on the two test networks
across a 20-fold reduction of $\varepsilon$. Where the regularised and exact
estimators disagree they do so in the direction that *hides* the residual: the
Sinkhorn value on one network is 3.6× below the exact value on the same
weights, so removing the regularisation makes the residual larger, not
smaller — it cannot be the source. Only at the smallest $\varepsilon$ does the
regularised solve itself break down, and those runs are excluded rather than
interpreted. A step-size sweep at *fixed physical time* shows the residual
changes by at most 0.05% over a ten-fold change in $\Delta t$, so it is not a
discretisation artefact either. Extending three networks to 5,000 steps
under exact transport — a 25-fold horizon, with weights spanning $10^{65}$
to $10^{72}$ — changes the residual by +0.009%, +0.0002% and +0.17%. We
report this as trajectory evidence and explicitly not as evidence that
$V_G^*>0$.

**Finding 4 — what survives, and what does not.** Homogenisation requires
the curvature–edge correspondence: an information-destroying sham flow
reverses it on 24/25 matched pairs and on 8/8 source families, so the effect
is carried by the pairing of curvature to edge and not by the weight
dynamics alone. A constraint-relaxation control (surgery) removes most of
the residual, median 6.6-fold, but rests on two families and is reported
with the effect size as primary. The structure-to-residual mapping does
**not** generalise: leave-one-family-out out-of-sample $R^2$ is at or below
the intercept baseline for every predictor we tried. One correlation we
previously reported is withdrawn, one is demoted, and the residual itself is
reported for two connected graphs from a single family of seventeen — which
supports no prevalence claim, and is written as a case study.

**Methodological advance.** Every headline inference is clustered at the
level of the source family rather than the network, with exact family-level
permutation tests, cluster bootstrap and leave-one-family-out jackknife. Six
analysis populations $C_0$–$C_5$ are defined once and used consistently, and
each claim carries an evidence tier (primary: $\ge 8$ families; exploratory:
$<5$). A reproducibility package ships with the submission — a 93-row
network manifest, the population definitions, a machine-readable
configuration including the unsafe legacy stopping rule so that the artefact
itself can be reproduced, and a script index — rather than being promised
after acceptance.

**Why *Communications Physics*.** The paper takes a tool with a strong
geometric promise and tests it where the promise has been asserted more
often than checked, answering with controls rather than with argument: an
information-destroying sham, a constraint-relaxation control,
degree-preserving synthetic manipulation, exact transport, a
scale-equivariant solver, a second update rule, an $\varepsilon$-trajectory
sweep, a fixed-physical-time step-size sweep, and a 5,000-step extended
horizon. It delivers three corrections — to the stopping rule, the
estimator, and the reading of scale — that affect any prior or future
empirical use of this flow, and it states its negative results plainly.

**Compliance.** Main text 4,997 words (limit $\sim$5,000); abstract 149
words (limit 150); 10 display items (8 figures, 2 tables; limit 10); Methods
is not word-limited and carries the numerical detail. Eleven supplementary
tables S1–S10 (S7 has a companion S7b), including a per-network provenance
table listing the landing page and verified download URL of all 93 networks.
The manuscript compiles without error and every one of its 45 references is
cited in the text and resolves to a primary record.

**Data and code.** All networks are from the Network Data Repository
(http://networkrepository.com), distributed openly for research use; we follow
its acknowledgement policy and ship the derived data — manifest, trajectories
and intermediate results — with the submission. All analysis code is released
under the MIT licence at
https://github.com/fanqimcse/ricci-flow-finite-time-residuals and frozen at
`v1.0.1`; it has a permanent Zenodo DOI, 10.5281/zenodo.23330125 (version) /
10.5281/zenodo.23277161 (concept). No number in the paper is hard-coded, and a
single script regenerates every figure and statistic from the raw result files,
emitting a log of the number behind each plotted point. An archival request has
also been submitted to Software Heritage (save request 2549855, accepted).

Suggested reviewers:
- [Reviewer 1 — discrete curvature on graphs; name/affiliation/e-mail to be completed]
- [Reviewer 2 — network geometry / complex networks; to be completed]
- [Reviewer 3 — optimal transport and numerical analysis; to be completed]

We thank you for your consideration.

Sincerely,

Qi Fan
Huaibei Normal University, Huaibei 235000, China
8073592@qq.com — ORCID 0000-0003-3591-3363
