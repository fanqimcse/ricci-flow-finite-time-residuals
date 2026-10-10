#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Extend Supplementary Table S1 with the weight provenance columns.

Review point 5 asks whether S1 records, for every dataset, where the weights
come from and how they were transformed.  It did not.  This adds

  weight_mode        auto | unweighted | strength  (production run)
  weight_reason      why that mode was chosen by the loader
  has_third_column   whether the file carries a numeric third column
  raw_min/max/median/cv/distinct   statistics of the raw third column
  raw_all_integer / raw_year_like
  documented_semantics             what the third column is, per source
  semantics_basis                  how that was established
  documented_as_distance           yes / no / n-a
  length_transform                 the map raw value -> edge length
  conventions_tested               which alternative readings were run
"""

from __future__ import print_function

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
MS = os.path.join(HERE, "paper_a", "manuscript")

S1 = pd.read_csv(os.path.join(MS, "Supplementary_Table_S1_dataset_sources.csv"))
SUM = pd.read_csv(os.path.join(HERE, "results_correct_fast", "dataset_summary.csv"))
S12 = pd.read_csv(os.path.join(MS, "Supplementary_Table_S12_weight_semantics.csv"))

mode = dict(zip(SUM["dataset"], SUM["weight_mode"]))
reason = dict(zip(SUM["dataset"], SUM["weight_reason"]))

_RAW = pd.read_csv(os.path.join(HERE, "exp_round6",
                                "raw_weight_column_stats.csv"))
raw_cols = ["n_values", "raw_min", "raw_max", "raw_median", "raw_cv",
            "raw_distinct", "all_integer", "year_like"]
rawmap = {c: dict(zip(_RAW["dataset"], _RAW[c])) for c in raw_cols}
sem = dict(zip(S12["dataset"], S12["documented_semantics"]))
basis = dict(zip(S12["dataset"], S12["basis"]))


def _mode(d):
    return mode.get(d, "")


S1["weight_mode"] = S1["dataset_id"].map(_mode)
S1["weight_reason"] = S1["dataset_id"].map(lambda d: reason.get(d, ""))
RENAME = {"n_values": "n_values", "raw_min": "raw_min", "raw_max": "raw_max",
          "raw_median": "raw_median", "raw_cv": "raw_cv",
          "raw_distinct": "raw_distinct", "all_integer": "raw_all_integer",
          "year_like": "raw_year_like"}
for c in raw_cols:
    S1[RENAME[c]] = S1["dataset_id"].map(lambda d: rawmap[c].get(d, np.nan))

S1["has_third_column"] = np.where(S1["n_values"].notna(), "yes", "no")
S1["documented_semantics"] = S1["dataset_id"].map(
    lambda d: sem.get(d, "no numeric third column"))
S1["semantics_basis"] = S1["dataset_id"].map(lambda d: basis.get(d, ""))
S1["documented_as_distance"] = np.where(
    S1["has_third_column"] == "no", "n/a",
    np.where(S1["dataset_id"].isin(set(S12["dataset"])), "no", "unchecked"))

TRANSFORM = {
    "unweighted": "l_e = 1 for every edge (no third column used)",
    "auto": "l_e = w_e (raw third column, distance reading)",
    "strength": "l_e = 1 / w_e",
}
S1["length_transform"] = S1["weight_mode"].map(
    lambda m: TRANSFORM.get(m, "unknown"))
S1["conventions_tested"] = np.where(
    S1["dataset_id"].isin(set(S12["dataset"])),
    "as-is (l=w); reciprocal (l=1/w); unit (l=1)", "n/a")

ordered = ["dataset_id", "family", "domain", "n_nodes", "n_edges",
           "n_components", "connected", "has_third_column", "weight_mode",
           "length_transform", "weight_reason", "n_values", "raw_min",
           "raw_max", "raw_median", "raw_cv", "raw_distinct", "raw_all_integer",
           "raw_year_like", "documented_semantics", "semantics_basis",
           "documented_as_distance", "conventions_tested", "landing_page",
           "download_url_zip", "page_status"]
S1 = S1[[c for c in ordered if c in S1.columns]]
S1.to_csv(os.path.join(MS, "Supplementary_Table_S1_dataset_sources.csv"),
          index=False)

print("S1 rows: %d" % len(S1))
print("with numeric third column: %d" % int((S1["has_third_column"] == "yes").sum()))
print("documented as distance: %d"
      % int((S1["documented_as_distance"] == "yes").sum()))
print("weight modes:")
print(S1["weight_mode"].value_counts().to_string())
print("year-like third column: %d" % int((S1["raw_year_like"] == True).sum()))  # noqa: E712
