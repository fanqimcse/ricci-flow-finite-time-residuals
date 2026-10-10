# Response to the fourth-round review

Round 8 of the internal review cycle. The fourth-round document contains **4 substantive items, 4 theoretical/wording items, 2 reproducibility/format items**, plus a prioritised table. All are addressed below. Two of them were verified by recomputation rather than by rewording, and one turned out to expose a real defect in the released archive.

The reviewer's overall verdict is unchanged and we agree with it: **do not submit yet.** We treat items 1–3 as the submission blockers; items 4–8 as required before submission; the remainder as format work. Item 8 was not answered by assertion — we downloaded the archive and tried to reproduce from it.

---

## Priority 1 — the limit in Eq. (2) (mathematical error)

**Correct. Fixed.** The reviewer's algebra is right, and so is the diagnosis: we had mixed two denominators.

With $V_t=V_\infty+Ae^{-\lambda t}$ and the *current-value* relative change

$$R_t=\frac{|V_t-V_{t-1}|}{V_t},$$

the numerator is $A(1-e^{-\lambda})e^{-\lambda(t-1)}$ in both cases, so

$$R_t=\frac{A(1-e^{-\lambda})e^{-\lambda(t-1)}}{V_\infty+Ae^{-\lambda t}}
\;\xrightarrow[V_\infty=0]{}\;e^{\lambda}-1 .$$

The $1-e^{-\lambda}$ limit belongs to the *previous-value* denominator $|V_t-V_{t-1}|/V_{t-1}$, not to $R_t$.

Eq. (2) now displays both ratios side by side over the same numerator, states that the numerators agree and decay while the denominators do not, and gives the limits separately: $e^{\lambda}-1$ when $V_\infty=0$ (**not** zero), and decay like $e^{-\lambda t}/V_\infty$ when $V_\infty>0$. The downstream condition is written as the reviewer proposes, $\lambda<\ln(1+10^{-5})$.

We retained the reviewer's own caveat, which we agree with and have made explicit in the text: for small $\lambda$ the two expressions agree to first order, so this does not overturn the slow-decay argument — but a formal derivation error cannot stay in the paper.

**One place still uses the old form and is now deliberately framed as such.** The old rule divides by $V_0$, whose limit genuinely *is* $1-e^{-\lambda}$; that identity is retained where the old rule is being analysed, and the two are no longer presented as the same quantity.

## Priority 2 — "28/28" against "27/28"

**Correct, and the reviewer's instruction not to simply pick one number was the right one.** We recomputed both statistics from the 28 raw trajectories rather than choosing a side. They are genuinely different questions and give different answers:

| statistic | count | median | extreme |
|---|---|---|---|
| variance fell over the **final quarter** | **27 / 28** | $-24.7\%$ | $-84.9\%$ |
| variance fell at the **terminal step** | **26 / 28** | $<0.02\%$ | one rose $+1.1\%$ over the quarter |

The one network rising at the last step did so by less than $0.02\%$ while its quarter-net change was strongly negative — so neither count is a rounding artefact of the other. Every layer now reports both, with the quarter statistic first, because that is the one the claim rests on:

- **abstract**: "$27$ of the $28$ runs it labelled saturated were still falling over the final quarter";
- **Results**: both counts, with the $<0.02\%$ and $+1.1\%$ extremes;
- **Discussion** statement 1: $27/28$ over the quarter;
- **Methods**: both counts;
- **Fig. 2d in-plot annotation**: "n=28: 27/28 falling over the final quarter, 26/28 at the last step" (was "27/28 still falling");
- **new Supplementary Table S17**: per-network `last_step_down`, `q4_net`, `q4_monotone`, `q4_any_up`, `end_over_q4min` for all 28.

The previous "28/28" in the abstract was wrong, and the previous figure annotation was right for the wrong reason — it did not distinguish the two statistics. The caption no longer says "all but one were still falling".

## Priority 3 — $V_\infty$ and "plateau" labels in figures

**Fixed systematically.** The paper's central point is that a finite-time residual is not an asymptotic claim, so labels asserting a limit were self-undermining. Adopted the reviewer's terminology table:

| was | now |
|---|---|
| plateau $V_\infty$ (axis/in-plot) | **finite-horizon residual $V(T)$** |
| converged (reaching a numeric threshold) | **reached the operational threshold** |
| plateau (small terminal change) | **near-plateau over the observed window** |
| convergence criterion (our relative-change rule) | **operational stopping criterion** |

Applied to Figs. 2b/2c, 4, 6 axis labels and to every in-plot annotation; $V_\infty$ now survives only inside the analytic model of Eq. (2), where it is a genuine asymptotic quantity. Verified by extracting the text layer of all eight PDFs: no $V_\infty$, "plateau" or "converge" axis label remains in any figure.

Two wording points the reviewer's table implies, applied in the text: "near-plateau" is now reserved for the three 5,000-step runs that actually satisfy the criterion (S16), and no run is described as "converged".

## Priority 4 — abstract conflated step-size stability with longer runs

**Correct. Split.** The single clause "stable to step size, moving less than 0.2% from 2,000 to 5,000 steps" answered two different questions. The abstract now reports them as two separate results, each with its own scope:

> at matched physical time the residual is stable to step size; extending three runs from $2{,}000$ to $5{,}000$ steps moves it by $<0.2\%$

Step-size sensitivity is a discretisation check at fixed $T=20$; the horizon extension is a statement about the observed window on three trajectories. Neither is presented as evidence about the other, and both keep the "tested trajectories" limitation.

## Priority 5 — "any decaying trajectory passes the old rule"

**Fixed as the reviewer phrases it.** The claim was too wide. It now reads: with $V_0$ fixed in the denominator while the per-step change decays, **any run whose increment tends to zero may eventually satisfy it — exponential decay is the textbook case**. "Any" is gone, and the mechanism (constant denominator, decaying numerator) is stated rather than assumed.

## Priority 6 — "converged" / "convergence criterion" in Methods and labels

**Fixed.** All instances replaced:

- "convergence criterion" → **operational stopping criterion**, stated explicitly as "never a convergence criterion";
- "41 runs satisfy the convergence criterion" → "41 runs satisfy the operational stopping criterion — $34$ on the old rule and $7$ on their step horizon — and a stricter tail-drift test leaves 31; **neither count is a convergence claim**";
- "$V\le 10^{-6}V_0$" → "**reached the operational threshold** $V=10^{-6}$" (three places);
- figures: no convergence label remains (see item 3).

## Priority 7 — extreme weight spans

**Kept at the reviewer's ceiling, no further.** The manuscript now claims only that "the tested numerical controls do not" account for the residual, and immediately states what this does **not** establish: no claim that the residual is artefact-free, no claim of global optimality (two optimisers failed to beat the flow), and no claim that untested controls would also fail. We did not strengthen this in any direction.

## Priority 8 — archive verification (performed, not asserted)

The reviewer asked us to confirm the DOIs resolve, check the package matches v1.0.0, and rebuild at least the key tables and Figs. 2, 3, 8 in a clean environment. We did this rather than restating the claim.

**DOIs resolve.** `10.5281/zenodo.23277162` (version) and `10.5281/zenodo.23277161` (concept) both resolve through `doi.org`; the version record is titled *"Finite-time residuals in discrete Ollivier–Ricci flow: reproducibility package"*, version `v1.0.0`, one file, 10.79 MB, 603 entries.

**The package was downloaded and run.** Extracted into an empty directory outside the working tree; `make_paper_figures.py` was run there against the archived result files only.

- **All eight figures regenerate.** The data log emitted by the clean-room run is **byte-identical to the working copy across all 84 lines** — one line per plotted quantity. Every number behind Figs. 1–8 reproduces from the archive alone.
- **The check found a real defect.** Four paths in the released script still resolved against the package root instead of `results/`, so the run **stopped at the third figure**. This was invisible from the working tree, where results sit at the top level. Fixed; the corrected script is what the release now carries, and the clean-room verification was **repeated with the corrected script and passed**.
- **The archive was one round behind the manuscript.** It carried Tables S1–S10 while the paper cites S1–S17, and it did not contain the scripts that generate S15 and S17. Both are now in the release; all 17 tables and both audit scripts are present.

We report the defect rather than quietly shipping the fix, because a reproducibility statement that has not been executed is exactly what the reviewer was right to refuse to accept. The Data/Code availability statements now say what was run, what matched, and what the check found.

## Format items

- **Abstract** compressed to **145 words** (target 140–145).
- **Title** changed to *"What finite-time computations reveal about graph curvature flows"* — 9 words, no colon, no technical compound noun, and it states the paper's actual contribution rather than its method.
- **Display items**: 8 figures + 2 tables = 10, within the journal's range.
- **Source files**: editable LaTeX plus separate figure files (PDF + PNG) are in the release.
- **References**: all 45 entries were independently verified in round 6 against Crossref/OpenAlex and, where applicable, the PMLR/arXiv record of record — title, authors, volume, pages, DOI — with 4 bibliographic corrections and 2 over-claimed statements fixed. This round we re-checked that no citation is orphaned, uncited or duplicated: 0/0/0.

---

## What we did not do, and why

- **We did not make the residual look asymptotic.** Items 3 and 6 remove labels rather than add caveats. If anything, the paper now *states less* about limits than before.
- **We did not add experiments.** The reviewer's own priority is arithmetic, alignment and labelling, not more data; we agree, and the compute went into verification (item 8) rather than new runs.
- **We did not change the sample, the estimator or the stopping tolerance.** Nothing in this round touched the 93-network analysis set or the flow implementation.

## Current status

0 LaTeX errors, 0 undefined citations/references, 27 pages, main text **4,998 words**, abstract **145 words**, 10 display items, 45 references, Supplementary Tables S1–S17.

**Outstanding before submission** (needs author action, not analysis): push the updated release and tag **v1.0.1**; request the new Zenodo version DOI (the manuscript currently cites v1.0.0, whose archive predates S11–S17 — the corrected archive must be re-archived before submission, and the DOI updated); add the {Software Heritage} intrinsic identifier once the archival visit completes; revoke the personal access token exposed earlier in this project.
