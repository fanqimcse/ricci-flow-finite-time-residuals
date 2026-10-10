#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Build Supplementary Tables S11, S12, S13 for round 6 of the review.

  S11  long-horizon diagnostics for the three connected ENZYMES networks
       (why g10 is reported next to g105/g106)
  S12  weight-semantics table: what the third column of every weighted
       dataset is, and the three-convention sensitivity sweep
  S13  numerical reliability at extreme weight span: per-checkpoint
       consistency checks, final-weight solver/cross-check probes, and the
       initial-weight perturbation test

Outputs are written next to the other supplementary tables.
"""
from __future__ import print_function

import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "paper_a", "manuscript")
EXP = os.path.join(HERE, "exp_round6")
FOCUS = os.path.join(HERE, "exp_round6_focus")
REST = os.path.join(HERE, "exp_round6_rest")

os.makedirs(OUT, exist_ok=True)


def wide(df, index, column, values, suffix=""):
    """Pivot and flatten to <value>_<column> columns."""
    p = df.pivot_table(index=index, columns=column, values=values, aggfunc="first")
    p.columns = ["%s_%s" % (v, c) for v, c in p.columns]
    return p.reset_index()


# ------------------------------------------------------------------ S11 -----
s11 = pd.read_csv(os.path.join(EXP, "longrun_three_networks.csv"))
s11 = s11.rename(columns={"dataset": "network"})
s11.to_csv(os.path.join(OUT, "Supplementary_Table_S11_long_horizon_three_networks.csv"),
           index=False)
print("S11 rows:", len(s11))

# ------------------------------------------------------------------ S12 -----
raw = pd.concat([pd.read_csv(os.path.join(d, "raw_weight_column_stats.csv"))
                 for d in (EXP, FOCUS, REST)
                 if os.path.exists(os.path.join(d, "raw_weight_column_stats.csv"))],
                ignore_index=True).drop_duplicates(subset=["dataset"], keep="first")
sweeps = []
for d in (EXP, FOCUS, REST):
    f = os.path.join(d, "weight_semantics_sweep.csv")
    if os.path.exists(f):
        s = pd.read_csv(f)
        s["source"] = os.path.basename(d)
        sweeps.append(s)
sw = pd.concat(sweeps, ignore_index=True)
# keep the run with most steps per (dataset, convention)
sw = sw.sort_values("steps", ascending=False)
sw = sw.drop_duplicates(subset=["dataset", "convention"], keep="first")

V0 = wide(sw, "dataset", "convention", ["V0"])
RATIO = wide(sw, "dataset", "convention", ["ratio_within"])
RATE = wide(sw, "dataset", "convention", ["rate_Q4"])
STEP = wide(sw, "dataset", "convention", ["steps"])

s12 = raw.merge(V0, on="dataset", how="outer").merge(RATIO, on="dataset", how="outer")
s12 = s12.merge(RATE, on="dataset", how="outer").merge(STEP, on="dataset", how="outer")

manifest = pd.read_csv(os.path.join(HERE, "github_release", "data",
                                    "network_manifest.csv"))
s12["in_analysis_set"] = s12["dataset"].isin(set(manifest["dataset_id"]))
fam = dict(zip(manifest["dataset_id"], manifest["family"]))
s12["family"] = s12["dataset"].map(fam)

# ---- documented semantics (retrieved from the repository metadata pages) ---
VERIFIED = {
    "reptilia-tortoise-network-hw": ("unweighted; third column = timestamp (year)",
                                     "repository metadata: Edge weights = Unweighted"),
    "reptilia-tortoise-network-mc": ("unweighted; third column = timestamp (year)",
                                     "repository metadata: Edge weights = Unweighted"),
    "reptilia-tortoise-network-sg": ("unweighted; third column = timestamp (year)",
                                     "repository metadata (same study, same collection)"),
    "reptilia-tortoise-network-bsv": ("unweighted; third column = timestamp (year)",
                                      "inferred: same study and collection as -hw/-mc"),
    "reptilia-tortoise-network-cs": ("unweighted; third column = timestamp (year)",
                                     "inferred: same study and collection as -hw/-mc"),
    "reptilia-tortoise-network-fi": ("unweighted; third column = timestamp (year)",
                                     "inferred: same study and collection as -hw/-mc"),
    "reptilia-tortoise-network-pv": ("unweighted; third column = timestamp (year)",
                                     "inferred: same study and collection as -hw/-mc"),
    "aves-barn-swallow-non-physical": ("frequency (number of interactions)",
                                       "repository metadata: Edge weight type = frequency"),
    "aves-sparrow-social-2009": ("association index / frequency",
                                 "inferred: same Animal Networks collection"),
    "aves-sparrow-social": ("association index / frequency",
                            "inferred: same Animal Networks collection"),
    "aves-sparrowlyon-flock-season2": ("association index / frequency",
                                       "inferred: same Animal Networks collection"),
    "aves-sparrowlyon-flock-season3": ("association index / frequency",
                                       "inferred: same Animal Networks collection"),
    "aves-thornbill-farine": ("association index / frequency",
                              "inferred: same Animal Networks collection"),
    "mammalia-raccoon-proximity": ("duration of proximity",
                                   "repository metadata: Edge weight type = duration"),
    "mammalia-primate-association": ("association index",
                                     "inferred: same Animal Networks collection"),
    "mammalia-voles-bhp-trapping": ("count of co-captures",
                                    "inferred: same Animal Networks collection"),
    "mammalia-voles-kcs-trapping": ("count of co-captures",
                                    "inferred: same Animal Networks collection"),
    "mammalia-voles-plj-trapping": ("count of co-captures",
                                    "inferred: same Animal Networks collection"),
    "insecta-ant-trophallaxis-colony1": ("count of trophallaxis events",
                                         "inferred: same Animal Networks collection"),
    "insecta-ant-trophallaxis-colony2": ("count of trophallaxis events",
                                         "inferred: same Animal Networks collection"),
    "bio-CE-HT": ("gene functional-association score (WormNet)",
                  "repository metadata: Weighted, WormNet Bayesian score"),
    "bio-CE-LC": ("gene functional-association score (WormNet)",
                  "inferred: same collection as bio-CE-HT"),
    "bio-DM-LC": ("gene functional-association score (WormNet)",
                  "inferred: same collection as bio-CE-HT"),
    "bio-SC-TS": ("gene functional-association score (WormNet)",
                  "inferred: same collection as bio-CE-HT"),
    "eco-mangwet": ("undocumented; values span 5.6e-9 to 720 (flow-like)",
                    "no weight metadata on the repository page; cited source is a trophic-flow food web"),
    "eco-stmarks": ("undocumented; values span 1e-3 to 265 (flow-like)",
                    "no weight metadata on the repository page; cited source is a trophic-flow food web"),
    "ia-workplace-contacts": ("undocumented; values 2.9e4 to 1.0e6 (duration-like)",
                              "not retrieved"),
    "ia-dnc-corecipient": ("count of shared emails",
                           "in-file header comment (dataset preamble)"),
}

s12["documented_semantics"] = s12["dataset"].map(
    lambda d: VERIFIED.get(d, ("unknown", "not retrieved"))[0])
s12["basis"] = s12["dataset"].map(
    lambda d: VERIFIED.get(d, ("unknown", "not retrieved"))[1])
s12["documented_as_distance"] = "no"
s12["conventions_tested"] = "as-is (l=w), reciprocal (l=1/w), unit (l=1)"

STOP = wide(sw, "dataset", "convention", ["stop"])
s12 = s12.merge(STOP, on="dataset", how="outer")

cols = ["dataset", "family", "in_analysis_set", "n_values", "raw_min", "raw_max",
        "raw_median", "raw_cv", "all_integer", "year_like",
        "documented_semantics", "basis", "documented_as_distance",
        "conventions_tested",
        "steps_asis", "steps_recip", "steps_unit",
        "stop_asis", "stop_recip", "stop_unit",
        "V0_asis", "V0_recip", "V0_unit",
        "ratio_within_asis", "ratio_within_recip", "ratio_within_unit",
        "rate_Q4_asis", "rate_Q4_recip", "rate_Q4_unit"]
s12 = s12[[c for c in cols if c in s12.columns]]
s12 = s12.sort_values("dataset")
s12.to_csv(os.path.join(OUT, "Supplementary_Table_S12_weight_semantics.csv"),
           index=False)
print("S12 rows:", len(s12))

# ------------------------------------------------------------------ S13 -----
ck = pd.read_csv(os.path.join(EXP, "num_checkpoints.csv"))
fp = pd.read_csv(os.path.join(EXP, "num_final_probes.csv"))
init_path = os.path.join(EXP, "init_weight_perturbation.csv")
init = pd.read_csv(init_path) if os.path.exists(init_path) else None

ck_out = ck.copy()
ck_out.to_csv(os.path.join(OUT, "Supplementary_Table_S13a_checkpoints.csv"),
              index=False)
fp.to_csv(os.path.join(OUT, "Supplementary_Table_S13b_final_probes.csv"),
          index=False)
if init is not None:
    init.to_csv(os.path.join(OUT, "Supplementary_Table_S13c_initial_weight_perturbation.csv"),
                index=False)
    print("S13c rows:", len(init))
print("S13a rows:", len(ck_out), " S13b rows:", len(fp))

# ------------------------------------------------------- console summary ----
print("\n--- S12: change in V0 and in the residual across conventions ---")
sub = s12.dropna(subset=["V0_asis", "V0_recip", "V0_unit"])
if len(sub):
    for _, r in sub.iterrows():
        print("  %-38s V0: %.3e / %.3e / %.3e   ratio: %s / %s / %s" % (
            r["dataset"], r["V0_asis"], r["V0_recip"], r["V0_unit"],
            ("%.3f" % r["ratio_within_asis"]) if r["ratio_within_asis"] == r["ratio_within_asis"] else "nan",
            ("%.3f" % r["ratio_within_recip"]) if r["ratio_within_recip"] == r["ratio_within_recip"] else "nan",
            ("%.3f" % r["ratio_within_unit"]) if r["ratio_within_unit"] == r["ratio_within_unit"] else "nan"))
print("\nsign of the residual ratio (contraction) per convention:")
for c in ("asis", "recip", "unit"):
    col = "ratio_within_%s" % c
    if col in sub.columns:
        v = sub[col].dropna()
        print("  %-6s n=%2d  median=%.3f  max=%.3f  fraction<1: %.2f" % (
            c, len(v), v.median(), v.max(), float((v < 1).mean())))

if init is not None:
    print("\n--- initial-weight perturbation ---")
    for ds, g in init.groupby("dataset"):
        base = g[g.rep == 0]["V_within"].iloc[0]
        pert = g[g.rep > 0]["V_within"].values
        print("  %-16s baseline=%.6e perturbed range=[%.6e, %.6e] max rel dev=%.3e" % (
            ds, base, pert.min(), pert.max(),
            float(np.max(np.abs(pert / base - 1.0)))))
