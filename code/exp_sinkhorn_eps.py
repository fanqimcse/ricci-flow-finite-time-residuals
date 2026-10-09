"""Experiment A + B: is the stalled residual a numerical artefact?

Two questions the reviewer says must be answered before any claim about a
"geometric limit":

  A. Sinkhorn regularisation sensitivity.  We do not compute the exact
     Wasserstein-1 distance; we compute an entropy-regularised surrogate with
     a *fixed* epsilon = 0.05.  In the stalled branch the edge weights span
     10^5 - 10^11, so the transport cost matrix is badly scaled and the
     regularisation bias can be large and weight-dependent.  If the plateau
     value V_infty moves with epsilon, the stalled branch is a regularisation
     artefact.

  B. Solver-precision sensitivity.  At each horizon we re-score the final
     weights with a much tighter Sinkhorn solver.  If V changes appreciably
     when the solver is tightened, the value is set by the solver, not by the
     graph.

Design: the flow update is byte-identical to the production run
(w <- w exp(-dt kappa), renormalised by the mean weight, dt = 0.05,
alpha = 0.5); only the transport sub-solver settings change.

Usage
-----
    python exp_sinkhorn_eps.py --out exp_eps --steps 200 \
        --regs 0.1,0.05,0.02,0.01 --only scc_rt_oman,scc2034_kilifi...
"""

from __future__ import print_function

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

SOURCE_DIR = os.path.dirname(os.path.abspath(__file__))
if SOURCE_DIR not in sys.path:
    sys.path.insert(0, SOURCE_DIR)

import networkx as nx
import ricci_flow_program_100data as R
import ricci_flow_longrun as L
from fast_sinkhorn import batched_sinkhorn_cost

DATA_DIR = r"D:\100data"
ALPHA = 0.5
DT = 0.05

TARGETS = [
    ("scc_rt_oman", "stalled"),
    ("scc2034_kilifi_all_contacts_across_households", "stalled"),
    ("ENZYMES_g106", "stalled"),
    ("scc_rt_ksa", "converging"),
]


def find_file(dataset_id):
    for path in R.discover_graph_files(DATA_DIR):
        if R.safe_dataset_id(path, DATA_DIR) == dataset_id:
            return path
    return None


def make_workspace(graph):
    return R.CurvatureWorkspace(graph, ALPHA, "auto", 256.0)


def curvature_variance(ws, graph, reg, max_iter, stop_thr):
    ws.prepare_step()
    mass_u, mass_v, costs, dists = [], [], [], []
    for i in range(len(ws.edges)):
        _, _, mu, mv, cost, d = ws.edge_payload(i)
        mass_u.append(mu)
        mass_v.append(mv)
        costs.append(cost)
        dists.append(d)
    values = batched_sinkhorn_cost(mass_u, mass_v, costs, reg, max_iter,
                                   stop_thr)
    k = np.array([1.0 - float(v) / float(d)
                  for v, d in zip(values, dists)], dtype=float)
    return k


def run_one(dataset_id, fate, reg, steps, out_dir, budget_s,
            tight_thr=1e-13, tight_iter=50000):
    path = find_file(dataset_id)
    if path is None:
        print("MISSING %s" % dataset_id)
        return None
    base = L.build_base_args(["--data_dir", DATA_DIR])
    graph = R.load_graph(path, base)
    ws = make_workspace(graph)
    edges = list(ws.edges)

    out_csv = os.path.join(out_dir, "%s_reg%g.csv" % (dataset_id, reg))
    if os.path.exists(out_csv):
        try:
            prior = pd.read_csv(out_csv)
            if len(prior) >= steps:
                print("   [skip] %s reg=%g already has %d steps"
                      % (dataset_id, reg, len(prior)))
                return prior
        except Exception:
            pass

    rows = []
    t_start = time.time()
    step = 0
    try:
        k = curvature_variance(ws, graph, reg, 5000, 1e-9)
        v0 = float(k.var())
        while step < steps:
            rows.append(dict(step=step, var_ratio=float(k.var() / v0),
                             variance=float(k.var()),
                             mean_k=float(k.mean()),
                             max_abs_k=float(np.abs(k).max()),
                             w_min=float(min(d["weight"]
                                             for _, _, d in graph.edges(data=True))),
                             w_max=float(max(d["weight"]
                                             for _, _, d in graph.edges(data=True)))))
            # ---- production update, identical to the long-run driver ----
            curv = dict(zip(edges, k))
            R.ricci_flow_step(graph, curv, DT)
            R.normalize_weights(graph)
            k = curvature_variance(ws, graph, reg, 5000, 1e-9)
            step += 1
            if time.time() - t_start > budget_s:
                break
    except Exception as exc:
        print("   [stop] %s reg=%g at step %d: %r" % (dataset_id, reg, step, exc))

    df = pd.DataFrame(rows)
    if len(df):
        # ---- Experiment B: re-score the final weights with a tighter solver --
        try:
            k_tight = curvature_variance(ws, graph, reg, tight_iter, tight_thr)
            v_std = float(df["variance"].iloc[-1])
            v_tight = float(k_tight.var())
            df["v_tight"] = np.nan
            df.loc[df.index[-1], "v_tight"] = v_tight
            df["rel_solver_shift"] = np.nan
            df.loc[df.index[-1], "rel_solver_shift"] = (abs(v_tight - v_std)
                                                        / max(v_std, 1e-300))
        except Exception as exc:
            print("   [tighten failed] %r" % exc)
        df["dataset"] = dataset_id
        df["fate"] = fate
        df["reg"] = reg
        df.to_csv(out_csv, index=False)
        shift = df["rel_solver_shift"].iloc[-1]
        print("   %-52s reg=%-5g steps=%d  V/V0=%.4e  span=%.1e  "
              "solver_shift=%.2e"
              % (dataset_id[:52], reg, len(df),
                 df["var_ratio"].iloc[-1],
                 df["w_max"].iloc[-1] / max(df["w_min"].iloc[-1], 1e-300),
                 shift if shift == shift else float("nan")))
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="exp_eps")
    ap.add_argument("--steps", type=int, default=200)
    ap.add_argument("--regs", default="0.1,0.05,0.02,0.01")
    ap.add_argument("--budget", type=float, default=150.0)
    ap.add_argument("--only", default="")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    regs = [float(x) for x in args.regs.split(",") if x]
    wanted = [t for t in args.only.split(",") if t]
    for dataset_id, fate in TARGETS:
        if wanted and dataset_id not in wanted:
            continue
        print("[%s] %s" % (fate, dataset_id))
        for reg in regs:
            run_one(dataset_id, fate, reg, args.steps, args.out, args.budget)


if __name__ == "__main__":
    main()
