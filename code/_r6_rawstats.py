#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Fast extraction of the raw third-column statistics (no graph loading).

Writes ``exp_round6/raw_weight_column_stats.csv``, which
``make_supplementary_r6.py`` merges into Supplementary Table S12.
"""
from __future__ import print_function

import os
import sys

import numpy as np
import pandas as pd

SOURCE_DIR = os.path.dirname(os.path.abspath(__file__))
if SOURCE_DIR not in sys.path:
    sys.path.insert(0, SOURCE_DIR)

import exp_round6 as E

OUT = os.path.join(SOURCE_DIR, "exp_round6")


def main():
    rows = []
    for ds in E.WEIGHTED:
        st = E.raw_third_column_stats(ds)
        print("[raw] %-42s %s" % (ds, "MISSING" if st is None else
                                  "min=%.4g max=%.4g cv=%.3f year_like=%s"
                                  % (st["raw_min"], st["raw_max"],
                                     st["raw_cv"], st["year_like"])),
              flush=True)
        if st:
            rows.append(st)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "raw_weight_column_stats.csv"), index=False)
    print("\nwrote %d rows" % len(df))
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
