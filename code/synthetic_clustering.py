# -*- coding: utf-8 -*-
"""Experiment D: manipulate transitivity with the degree sequence FIXED.

Why
---
The pooled audit (real + synthetic, n=97) ranked transitivity above every
other candidate: Spearman(log10 plateau, transitivity) = +0.821 pooled,
+0.621 within real networks, +0.736 within synthetic graphs, and real
networks carry 5.8x more transitivity than the configuration-model graphs
(0.463 vs 0.079). That pattern is exactly what a driver of the ~90x
magnitude gap should look like.

But transitivity has never been *manipulated*. In experiment B it entered
only as a covariate (and dropped out under base fixed effects), and in
experiment C it was collinear with modularity by construction (rho = 0.965).
Correlation across heterogeneous graphs cannot separate it from size, family,
or degree sequence. This experiment can.

Design
------
  * same configuration-model base graphs as experiments B and C, so the three
    experiments share a common reference point
  * degree-preserving double-edge swaps, accepted only when they move
    transitivity towards the target -- the degree sequence is invariant
  * targets: 0.10 (baseline-ish) / 0.30 / 0.50
  * every graph run to saturation (not a matched horizon): the run-to-
    saturation check in experiment C showed that matched-horizon comparisons
    overstate effects that are really differences in convergence speed

The number that matters is whether the SATURATED plateau rises towards the
real-network median (0.197) as transitivity rises towards real-network values
(0.46). A rise from ~0.002 to ~0.1 would close most of the gap.
"""

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
import networkx as nx

SOURCE_DIR = os.path.dirname(os.path.abspath(__file__))
if SOURCE_DIR not in sys.path:
    sys.path.insert(0, SOURCE_DIR)

from synthetic_assortativity import (base_graph, try_swap, run_flow)  # noqa


def rewire_to_transitivity(graph, target, rng, max_attempts=6000, tol=0.025):
    g = graph.copy()
    try:
        current = float(nx.transitivity(g))
    except Exception:
        current = 0.0
    best, best_val = g.copy(), current
    for _ in range(max_attempts):
        if abs(current - target) <= tol:
            break
        candidate = g.copy()
        if not try_swap(candidate, rng):
            continue
        try:
            new_val = float(nx.transitivity(candidate))
        except Exception:
            continue
        if abs(new_val - target) < abs(current - target):
            g, current = candidate, new_val
            if abs(current - target) < abs(best_val - target):
                best, best_val = g.copy(), current
    return best, best_val


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    parser.add_argument("--n", type=int, default=50)
    parser.add_argument("--mean_degree", type=int, default=6)
    parser.add_argument("--levels", default="0.10,0.30,0.50")
    parser.add_argument("--base_list", default="0")
    parser.add_argument("--reps", type=int, default=1)
    parser.add_argument("--max_steps", type=int, default=200)
    parser.add_argument("--tol", type=float, default=1e-5)
    parser.add_argument("--patience", type=int, default=50)
    parser.add_argument("--tol_mode", default="current",
                        choices=("current", "initial"))
    parser.add_argument("--seed0", type=int, default=3000)
    args = parser.parse_args()

    base_indices = [int(x) for x in args.base_list.split(",") if x.strip()]
    levels = [float(x) for x in args.levels.split(",")]
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    rows = []
    if os.path.exists(args.out):
        try:
            rows = pd.read_csv(args.out).to_dict("records")
        except Exception:
            rows = []
    seen = {(int(r["base"]), float(r["target_T"]), int(r["rep"]))
            for r in rows}

    t_start = time.time()
    for base_index in base_indices:
        graph = base_graph(args.n, args.mean_degree, args.seed0 + base_index)
        degree_seq = sorted(d for _, d in graph.degree())
        for level in levels:
            for rep in range(args.reps):
                if (base_index, level, rep) in seen:
                    continue
                rng = np.random.default_rng(
                    args.seed0 + 97 * base_index + 13 * rep
                    + int(level * 100) + 700)
                rewired, achieved = rewire_to_transitivity(graph, level, rng)
                try:
                    trans = float(nx.transitivity(rewired))
                    r_final = float(
                        nx.degree_assortativity_coefficient(rewired))
                except Exception:
                    trans, r_final = np.nan, np.nan
                res = run_flow(rewired.copy(), max_steps=args.max_steps,
                               tol=args.tol, patience=args.patience,
                               tol_mode=args.tol_mode)
                rows.append({
                    "base": base_index, "target_T": level, "rep": rep,
                    "achieved_T": achieved, "transitivity": trans,
                    "final_r": r_final,
                    "nodes": rewired.number_of_nodes(),
                    "edges": rewired.number_of_edges(),
                    "degree_seq_hash": hash(tuple(degree_seq)) & 0xFFFF,
                    "var_ratio": res["var_ratio"], "steps": res["steps"],
                    "status": res["status"], "message": res["message"],
                    "mean_kappa_final": res["mean_kappa_final"],
                })
                pd.DataFrame(rows).to_csv(args.out, index=False)
                print("base={} T={:.2f}->{:.3f} rep={} r={:+.3f} "
                      "ratio={:.5f} ({} steps, {:.0f}s)".format(
                          base_index, level, achieved, rep,
                          r_final if r_final == r_final else float("nan"),
                          res["var_ratio"] if res["var_ratio"] ==
                          res["var_ratio"] else float("nan"),
                          res["steps"], time.time() - t_start), flush=True)

    pd.DataFrame(rows).to_csv(args.out, index=False)
    print("DONE {} rows -> {}".format(len(rows), args.out), flush=True)


if __name__ == "__main__":
    main()
