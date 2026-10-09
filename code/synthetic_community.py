# -*- coding: utf-8 -*-
"""Experiment C: does community structure explain the 92x magnitude gap?

The puzzle
----------
The assortativity sweep (synthetic_assortativity.py) reproduces the DIRECTION
of the real-network effect but not its SIZE: configuration-model graphs plateau
at var_ratio ~ 0.002 while real networks plateau at ~ 0.20, a 92x gap.
Configuration-model graphs have no community structure. Real networks do.

Hypothesis
----------
Community structure -- not degree mixing -- sets the MAGNITUDE of the residual;
assortativity sets the ORDERING within a comparably-structured family.
If so, adding communities to a configuration-model graph at fixed degree
sequence should push the plateau from ~0.002 towards ~0.2.

Design
------
  * configuration model base, degree sequence FIXED (same generator as the
    assortativity sweep, so the two experiments are directly comparable)
  * nodes assigned to k groups at random
  * degree-preserving double-edge swaps driven towards a target cross-group
    edge fraction mu  (mu = 1 - 1/k means "no structure" for k equal groups)
  * mu targets: 0.75 (none, k=4) / 0.45 / 0.15 (strong)
  * 3 bases x 3 levels x 2 replicates = 18 graphs
  * every graph: Ricci flow to saturation, same numerics as everywhere else

Caveat to report: pushing towards strong community structure with
heterogeneous degrees also moves assortativity (hubs cluster with hubs), so mu
and r are not independent by construction. Both are measured and both enter
the joint model, exactly as in the assortativity sweep.
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


def cross_fraction(graph, groups):
    n_cross = sum(1 for u, v in graph.edges() if groups[u] != groups[v])
    return n_cross / max(graph.number_of_edges(), 1)


def modularity_planted(graph, groups):
    """Modularity of the PLANTED partition (not a detected one)."""
    try:
        communities = [[] for _ in range(max(groups.values()) + 1)]
        for node, g in groups.items():
            communities[g].append(node)
        return float(nx.community.modularity(graph, communities))
    except Exception:
        return np.nan


def rewire_to_mu(graph, groups, target_mu, rng, max_attempts=8000, tol=0.03):
    g = graph.copy()
    current = cross_fraction(g, groups)
    best, best_mu = g.copy(), current
    for _ in range(max_attempts):
        if abs(current - target_mu) <= tol:
            break
        candidate = g.copy()
        if not try_swap(candidate, rng):
            continue
        new_mu = cross_fraction(candidate, groups)
        if abs(new_mu - target_mu) < abs(current - target_mu):
            g, current = candidate, new_mu
            if abs(current - target_mu) < abs(best_mu - target_mu):
                best, best_mu = g.copy(), current
    return best, best_mu


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=os.path.join(
        SOURCE_DIR, "synthetic_results", "community_sweep.csv"))
    parser.add_argument("--n", type=int, default=50)
    parser.add_argument("--mean_degree", type=int, default=6)
    parser.add_argument("--groups", type=int, default=4)
    parser.add_argument("--levels", default="0.75,0.45,0.15")
    parser.add_argument("--bases", type=int, default=3)
    parser.add_argument("--reps", type=int, default=2)
    parser.add_argument("--max_steps", type=int, default=200)
    parser.add_argument("--seed0", type=int, default=2000)
    parser.add_argument("--base_list", default="")
    args = parser.parse_args()

    if args.base_list.strip():
        base_indices = [int(x) for x in args.base_list.split(",") if x.strip()]
    else:
        base_indices = list(range(args.bases))

    levels = [float(x) for x in args.levels.split(",")]
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    rows = []
    total = len(base_indices) * len(levels) * args.reps
    done = 0
    t_start = time.time()

    for base_index in base_indices:
        graph = base_graph(args.n, args.mean_degree, args.seed0 + base_index)
        rng0 = np.random.default_rng(args.seed0 + 31 * base_index)
        order = list(graph.nodes())
        rng0.shuffle(order)
        groups = {}
        for position, node in enumerate(order):
            groups[node] = position % args.groups
        degree_seq = sorted(d for _, d in graph.degree())

        for level in levels:
            for rep in range(args.reps):
                rng = np.random.default_rng(
                    args.seed0 + 97 * base_index + 13 * rep
                    + int(level * 100) + 500)
                rewired, achieved = rewire_to_mu(graph, groups, level, rng)
                try:
                    trans = nx.transitivity(rewired)
                    r_final = nx.degree_assortativity_coefficient(rewired)
                except Exception:
                    trans, r_final = np.nan, np.nan
                q = modularity_planted(rewired, groups)
                res = run_flow(rewired.copy(), max_steps=args.max_steps)
                rows.append({
                    "base": base_index, "target_mu": level, "rep": rep,
                    "achieved_mu": achieved,
                    "modularity": q,
                    "final_r": r_final, "transitivity": trans,
                    "nodes": rewired.number_of_nodes(),
                    "edges": rewired.number_of_edges(),
                    "degree_seq_hash": hash(tuple(degree_seq)) & 0xFFFF,
                    "var_ratio": res["var_ratio"], "steps": res["steps"],
                    "status": res["status"], "message": res["message"],
                    "mean_kappa_final": res["mean_kappa_final"],
                })
                done += 1
                pd.DataFrame(rows).to_csv(args.out, index=False)
                print("[{}/{}] base={} mu={:.2f}->{:.2f} rep={} Q={:.3f} "
                      "r={:+.3f} ratio={:.5f} ({} steps, {:.0f}s)".format(
                          done, total, base_index, level, achieved, rep,
                          q if q == q else float("nan"),
                          r_final if r_final == r_final else float("nan"),
                          res["var_ratio"] if res["var_ratio"] == res["var_ratio"]
                          else float("nan"), res["steps"],
                          time.time() - t_start), flush=True)

    pd.DataFrame(rows).to_csv(args.out, index=False)
    print("written: {}".format(args.out), flush=True)


if __name__ == "__main__":
    main()
