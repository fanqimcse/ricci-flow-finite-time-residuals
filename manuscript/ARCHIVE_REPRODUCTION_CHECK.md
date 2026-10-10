# Clean-room reproduction of the v1.0.0 reproducibility archive

Verification performed 2026-10-11 in response to item 8 of the fourth-round
internal review: *"Please actually confirm in a normal browser that both DOIs
resolve, check that the download matches the version described in the
manuscript, and rebuild at least the key tables and Figs. 2, 3 and 8 in a
clean environment."*

This file records what was run and what came out. It is deliberately not a
claim; it is a log.

---

## 1. DOI resolution

| identifier | resolves | record |
|---|---|---|
| `10.5281/zenodo.23277162` (version) | yes | title *"Finite-time residuals in discrete Ollivier–Ricci flow: reproducibility package"*, version `v1.0.0`, 1 file, 10.79 MB, 603 entries |
| `10.5281/zenodo.23277161` (concept) | yes | concept DOI, resolves to the latest version |

Both resolve through `doi.org`. Package metadata: title, version and licence
match the manuscript's Code availability statement.

## 2. Package contents (as downloaded)

```
_tools/       3 entries
code/        65 scripts
data/         6 files   (manifest, config.yaml, sources, download script)
figures/     19 files   (8 x PDF + 8 x PNG + make_paper_figures.py + data log)
manuscript/ 14 files   (manuscript.tex, references.bib, Tables S1-S10)
results/    487 files   (per-dataset trajectories, experiment outputs)
environment.yml, requirements.txt, CITATION.cff, LICENSE, README.md
```

## 3. Clean-room rebuild

Procedure: downloaded the archive through the DOI, extracted it into an
empty directory **outside the working tree**, and ran
`figures/make_paper_figures.py` there against the archived result files only.
No file from the working tree was present except the figure script itself.

### Result

**All eight figures regenerate, and the emitted data log is identical to the
working copy line for line.**

```
archive-rebuilt lines: 84
local lines:           84
IDENTICAL: every number in the data log reproduces from the Zenodo package
```

The log carries one line per plotted quantity (medians, ratios, Spearman
$\rho$ and $p$, LOFO scores, partial correlations, the 27/28 annotation).
Every one reproduces exactly.

### Defect found by this check

The released script contained **four paths that resolved against the package
root instead of `results/`**:

| line | path |
|---|---|
| 671 | `os.path.join(ROOT, "exp_round3")` |
| 747 | `os.path.join(ROOT, "ricci_flow_v2", "per_dataset")` |
| 748 | `os.path.join(ROOT, "strict_results", "fate_table.csv")` |
| 750 | `os.path.join(ROOT, "exp_round4", "traj_long")` |

Consequence: the clean-room run **raised `FileNotFoundError` at the third
figure** (Fig. 8). This is invisible in the working tree, where raw results
sit at the top level; it only appears in the released layout.

Fix: all four now resolve against `RES`, the result root that the script
already computes for its other paths. The corrected script regenerated all
eight figures in the same clean directory and its data log was again
identical to the working copy. A residual assertion in the patch script
fails if any data path still references `ROOT`.

### Version gap

The downloaded archive carried **Supplementary Tables S1–S10** while the
manuscript cites **S1–S17**, and did not include the scripts generating
S15 and S17. The release has been updated to carry all 17 tables plus
`code/_r7_stop.py` and `code/_r8_tail.py`.

## 4. Not verified here

- **The flow integration itself was not re-run.** Rebuilding the trajectories
  from raw datasets takes ~34 h of the wall-clock budget and is a separate
  exercise; this check covers the analysis and figure layer, which is where
  the paper's numbers are produced.
- **The Supplementary Tables were not regenerated in the clean room.** They
  are shipped as CSVs; the scripts that produce S15 and S17 are now in the
  archive, and regenerating them requires the per-step trajectories, which are
  present under `results/`.
- **Software Heritage identifier** not yet issued (save request 2549855,
  accepted; the archival visit has not completed).

## 5. Files

- working copy of the corrected script: `paper_a/make_paper_figures.py`
- release copy: `github_release/figures/make_paper_figures.py`
- clean-room copy and logs: `_verify_zenodo/` (not part of the archive)
