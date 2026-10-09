# -*- coding: utf-8 -*-
"""Chunked, resumable continuation of experiment D.

Why
---
Long single-process runs are unreliable in this environment (they are killed
at unpredictable times), so the flow is advanced in short chunks and every
chunk checkpoints the edge weights, the step counter and the full history.
Any chunk can be re-run; the trajectory is bit-for-bit the one
``synthetic_assortativity.run_flow`` would produce, because the update, the
config and the seeding are copied from it exactly.

Usage
-----
    python clust_continue.py --base 0 --level 0.50 --total 700 \
        --max_seconds 240 --out synthetic_results/clust_strict_b0_t050_700.csv

Run it repeatedly until the CSV reports a terminal status or the target step
count is reached.
"""

import argparse
import json
import os
import sys
import time

import numpy as np
import networkx as nx

SOURCE_DIR = os.path.dirname(os.path.abspath(__file__))
if SOURCE_DIR not in sys.path:
    sys.path.insert(0, SOURCE_DIR)

from synthetic_assortativity import base_graph, CFG                    # noqa
from synthetic_clustering import rewire_to_transitivity                # noqa
from ricci_flow_longrun import compute_curvature, build_base_args      # noqa
import ricci_flow_program_100data as R                                 # noqa


def build_initial_graph(base_index, level, n=50, mean_degree=6, seed0=3000,
                        rep=0):
    """The exact graph synthetic_clustering.py would build for this cell."""
    graph = base_graph(n, mean_degree, seed0 + base_index)
    rng = np.random.default_rng(seed0 + 97 * base_index + 13 * rep
                               + int(level * 100) + 700)
    rewired, achieved = rewire_to_transitivity(graph, level, rng)
    return rewired, achieved


def state_dir(base_index, level, root):
    return os.path.join(root, "ckpt_base%d_t%03d" % (base_index,
                                                     int(round(level * 100))))


def load_state(g, sd):
    meta_p = os.path.join(sd, "meta.json")
    w_p = os.path.join(sd, "weights.csv")
    if os.path.exists(meta_p) and os.path.exists(w_p):
        meta = json.load(open(meta_p, encoding="utf-8"))
        import csv as _csv
        with open(w_p, newline="", encoding="utf-8") as fh:
            for row in _csv.DictReader(fh):
                g[int(row["u"])][int(row["v"])]["weight"] = float(row["w"])
        return meta
    return {"step": 0, "var0": None, "prev": None, "stable": 0,
            "status": "running", "mean_k": np.nan}


def save_state(g, sd, meta):
    os.makedirs(sd, exist_ok=True)
    with open(os.path.join(sd, "weights.csv"), "w", encoding="utf-8") as fh:
        fh.write("u,v,w\n")
        for u, v, d in g.edges(data=True):
            fh.write("%d,%d,%.17g\n" % (u, v, float(d.get("weight", 1.0))))
    json.dump({k: (None if v is None else v) for k, v in meta.items()},
              open(os.path.join(sd, "meta.json"), "w", encoding="utf-8"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", type=int, required=True)
    ap.add_argument("--level", type=float, required=True)
    ap.add_argument("--total", type=int, default=700)
    ap.add_argument("--max_seconds", type=float, default=240.0)
    ap.add_argument("--tol", type=float, default=1e-5)
    ap.add_argument("--patience", type=int, default=50)
    ap.add_argument("--var_floor", type=float, default=1e-6)
    ap.add_argument("--tol_mode", default="current",
                    choices=("current", "initial"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--root", default=os.path.join(SOURCE_DIR,
                                                   "synthetic_results"))
    args = ap.parse_args()

    sd = state_dir(args.base, args.level, args.root)
    g, achieved = build_initial_graph(args.base, args.level)
    meta = load_state(g, sd)
    if meta["status"] not in ("running",):
        print("cell already terminal: %s at step %d" % (meta["status"],
                                                        meta["step"]))
        return

    os.makedirs(sd, exist_ok=True)
    hist_p = os.path.join(sd, "hist.csv")
    args_ = build_base_args([])
    workspace = R.CurvatureWorkspace(
        g, alpha=CFG["alpha"], backend=args_.distance_backend,
        max_distance_mb=args_.max_distance_mb)

    t0 = time.time()
    step = int(meta["step"])
    var0 = meta["var0"]
    prev = meta["prev"]
    stable = int(meta["stable"])
    status = "running"
    ratio = np.nan
    mean_k = np.nan
    done = 0
    fh = open(hist_p, "a", encoding="utf-8")
    if step == 0:
        fh.write("step,variance,var_ratio,min_w,max_w,mean_kappa\n")

    while step < args.total:
        mean_k, var_k, edge_k, _, _ = compute_curvature(
            workspace, CFG, CFG["solver"])
        var_k = float(var_k)
        if step == 0:
            var0 = var_k
        ratio = var_k / var0 if var0 and var0 > 0 else np.nan
        ws = np.array([d.get("weight", 1.0) for _, _, d in g.edges(data=True)])
        fh.write("%d,%.17g,%.17g,%.6g,%.6g,%.10g\n"
                 % (step, var_k, ratio, ws.min(), ws.max(), float(mean_k)))
        fh.flush()
        if step > 0 and var0 and var0 > 0:
            denom = var_k if args.tol_mode == "current" else var0
            move = abs(var_k - prev) / max(denom, 1e-300)
            stable = stable + 1 if move < args.tol else 0
            if stable >= args.patience:
                status = ("saturated_current" if args.tol_mode == "current"
                          else "saturated")
                step += 1
                break
        if var0 and var0 > 0 and var_k <= args.var_floor * var0:
            status = "constant_curvature"
            step += 1
            break
        prev = var_k
        R.ricci_flow_step(g, edge_k, dt=CFG["dt"])
        step += 1
        done += 1
        if time.time() - t0 > args.max_seconds:
            break
    fh.close()

    if status == "running" and step >= args.total:
        status = "horizon_reached"

    meta = {"step": step, "var0": var0, "prev": prev, "stable": stable,
            "status": status, "mean_k": float(mean_k) if mean_k == mean_k else None}
    save_state(g, sd, meta)

    try:
        trans = float(nx.transitivity(g))
        r_final = float(nx.degree_assortativity_coefficient(g))
    except Exception:
        trans, r_final = np.nan, np.nan
    row = {
        "base": args.base, "target_T": args.level, "rep": 0,
        "achieved_T": achieved, "transitivity": trans, "final_r": r_final,
        "nodes": g.number_of_nodes(), "edges": g.number_of_edges(),
        "degree_seq_hash": hash(tuple(sorted(d for _, d in g.degree()))) & 0xFFFF,
        "var_ratio": ratio, "steps": step, "status": status, "message": "",
        "mean_kappa_final": float(mean_k) if mean_k == mean_k else np.nan,
    }
    import pandas as pd
    pd.DataFrame([row]).to_csv(args.out, index=False)
    print("base=%d T=%.2f achieved_T=%.6f step=%d V=%.6g status=%s "
          "(+%d steps this chunk, %.0fs)"
          % (args.base, args.level, achieved, step, ratio, status, done,
             time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
