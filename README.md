# Finite-time residuals in discrete Ollivier–Ricci flow

**Reproducibility package** (version `v1.0.1`) for the manuscript

> **Finite-time residuals in discrete Ollivier–Ricci flow: a measurement
> framework, its controls, and its limits**
> Qi Fan (Huaibei Normal University, 8073592@qq.com,
> ORCID 0000-0003-3591-3363) — submitted to *Communications Physics*.

This repository contains every script, configuration file, network manifest and
raw result file needed to regenerate **all** numbers, tables and figures in the
paper. No value in the manuscript is hard-coded: `figures/make_paper_figures.py`
reads the files under `results/` and writes both the figure files and
`FIGURE_DATA_LOG.txt`, which lists the number behind every panel.

---

## 1. What the paper measures

Discrete Ollivier–Ricci flow on a weighted graph updates every edge weight by
$w_e \leftarrow w_e \exp(-\Delta t\,\kappa_e)$, with the mean edge weight
renormalised after each step, idleness $\alpha = 0.5$, step size
$\Delta t = 0.05$, and curvature $\kappa_e = 1 - W_{1,e}/d_e$ computed either
by exact $W_1$ optimal transport (`ot.emd`, with a HiGHS `linprog` fallback) or
by entropy-regularised Sinkhorn (`reg = 0.05`).

The paper separates four quantities that empirical work on this flow routinely
conflates:

| level | symbol | meaning | computed here? |
|---|---|---|---|
| L1 | $V(T)$ | curvature variance at the horizon actually reached | yes |
| L2 | $V_{\rm num}(T;\theta)$ | the same quantity as produced by the solver | yes |
| L3 | trend of $V(t)$ | the only level at which "still falling" is defined | yes |
| L4 | $V_G^*=\inf_w \mathrm{Var}[\kappa(G,w)]$ | geometric infimum over all positive weightings | **no** |

Nothing in the paper is claimed about L4.

---

## 2. Repository layout

```
code/        analysis and experiment scripts (Python)
data/        network_manifest.csv       93-network analysis set + all settings
             analysis_populations.csv   definitions of C0-C5 and control samples
             dataset_sources.csv        per-network landing page + download URL
             download_networks.py       fetches the corpus from networkrepository.com
config.yaml  machine-readable run configuration (including the unsafe legacy
             stopping rule, so the artefact can be reproduced deliberately)
results/     raw result files, one directory per experiment
manuscript/  manuscript.tex, references.bib, Supplementary_Table_S1-S10
figures/     make_paper_figures.py, Fig1-Fig8 (pdf + png), FIGURE_DATA_LOG.txt
```

## 3. Data

All 93 networks come from the **Network Data Repository**,
<http://networkrepository.com> (Rossi & Ahmed, AAAI 2015). They are
distributed by that repository for research use; please follow the terms stated
on its website and cite it when reusing the data.

`data/dataset_sources.csv` gives, for every network in the analysis set, its
family, size, component count, landing page
(`https://networkrepository.com/<id>.php`) and verified direct download URL.
To fetch the corpus:

```bash
python data/download_networks.py --out data/networks
```

## 4. Reproducing the figures

```bash
python -m venv .venv && .venv/Scripts/pip install -r requirements.txt
python figures/make_paper_figures.py
```

This reads `results/` and rewrites `figures/Fig1`–`Fig8` and
`figures/FIGURE_DATA_LOG.txt` from scratch (a few minutes; no network re-run
is required).

## 5. Reproducing the experiments

| step | script | output |
|---|---|---|
| first corpus pass | `code/ricci_flow_program_100data.py` | `results/ricci_flow_longrun_results/` |
| within-component decomposition | `code/exp_components.py`, `code/exp_corpus_within.py` | `results/exp_comp/`, Table S3 |
| strict re-run ($C_4$) | `code/strict_batch.py`, `code/strict_plateau.py`, `code/fate_analysis.py` | `results/strict_results/`, `results/ricci_flow_v2/`, Table S2 |
| family-level inference | `code/exp_family_inference.py` | `results/exp_family/` |
| sham control | `code/sham_flow.py`, `code/analyze_sham.py` | `results/sham_results/` |
| constraint-relaxation control | `code/analyze_surgery.py`, `code/surgery_pair_report.py` | `results/surgery_results/` |
| synthetic rewiring | `code/synthetic_assortativity.py`, `code/synthetic_clustering.py`, `code/synthetic_community2.py`, `code/clust_continue.py` | `results/synthetic_results/` |
| regularisation probes | `code/exp_sinkhorn_eps.py`, `code/exp_indep_opt.py` | `results/exp_eps/`, `results/exp_indep_opt/` |
| exact transport, scale equivariance | `code/exp_round3.py` | `results/exp_round3/` |
| extended horizon, step size, $\varepsilon$ and $V_G^*$ checks | `code/exp_round4.py` | `results/exp_round4/`, Tables S7–S9 |
| predictability | `code/lofo_cv.py`, `code/final_predictive_audit.py` | `results/convergence_results/` |

Random seeds are fixed (seed 0) and recorded in `data/network_manifest.csv`.

## 6. Computational budget

Runs were executed on an 8-core x86-64 CPU with 9.9 GB RAM, one thread per
process (`OMP/OPENBLAS/MKL_NUM_THREADS=1`) and at most two concurrent
processes. Measured wall clock (summed per-step timers, single process, BLAS pinned to one
thread): first pass over the corpus $34$~h; strict re-runs $2.0$~h
($1.5$~h for the $16$-network set reported in the paper, $0.5$~h for an auxiliary
set with extended horizons); exact-transport runs $0.11$~h; step-size sweep
$0.08$~h; $\varepsilon$-trajectory sweep $1.1$~h. The exact-transport figure is
not comparable with the corpus median per-step cost: those graphs have 5--33
nodes and 6--69 edges, so a step costs 5--21~ms rather than seconds. A further
$34$ networks were re-run disconnected-only and were not timed per run, so no
total is reported for them.

## 7. Version

This is the release accompanying the submission (`v1.0.1`). The version number in
the manuscript's Code Availability statement matches this tag.

`v1.0.1` exists because the archive was verified rather than assumed. The
`v1.0.0` package was downloaded from its DOI into an empty directory and
`make_paper_figures.py` was run there against the archived result files alone.
All eight figures regenerated and the emitted data log matched this tree
line for line, but the run stopped at the third figure: four paths still
resolved against the package root instead of `results/`, which is invisible
in a working tree. Those paths are fixed here, and the clean-room run was
repeated with the corrected script and passed. The release also now carries
Supplementary Tables S11--S17 and the scripts that generate S15 and S17
(`code/_r7_stop.py`, `code/_r8_tail.py`), which the manuscript cites and which
`v1.0.0` did not contain.

See `manuscript/ARCHIVE_REPRODUCTION_CHECK.md` for the full log, including what
was *not* verified (the flow integration itself was not re-run).

## 8. Licence

Code: MIT (see `LICENSE`). Manuscript text and figures: CC BY 4.0. The network
data are governed by the terms of the Network Data Repository.

## Publishing this repository

The tree is versioned locally (annotated tag `v1.0.1`). To create the public
repository and upload everything in one step:

```bash
python push_github.py <GITHUB_TOKEN>          # classic PAT with the `repo` scope
```

The script resolves the account login from the token (so the repository URL is
never guessed), creates the repository if needed, pushes `main` and the
`v1.0.1` tag (falling back to the GitHub Contents REST API if `git push` cannot
reach `github.com`), creates the GitHub Release, requests a Software Heritage
archival, and writes the real URL back into `manuscript.tex`,
`cover_letter.md` and `SUBMISSION_CHECKLIST.md`.
