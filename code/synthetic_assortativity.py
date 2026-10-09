"""Synthetic test: does degree assortativity *cause* the residual variance?

Why this is needed
------------------
Across 25 converged real networks the plateau ratio correlates with
assortativity (rho = +0.45, and +0.886 at family level). But real networks
vary in everything at once: size, density, clustering, modularity and
assortativity are entangled, so the correlation cannot separate them.

Here the degree sequence is held FIXED and only the wiring is changed, by
degree-preserving double-edge swaps driven towards a target assortativity.
Anything that then moves the plateau is attributable to the rewiring.

Caveat that must be reported: double-edge swaps also change transitivity
(clustering). So the experiment does not isolate assortativity from
clustering by construction; instead every generated graph has both r and
transitivity measured, and the analysis asks which one tracks the plateau
once both are in the model.

Design
------
  * 4 base degree sequences (configuration model, different seeds)
  * x 5 target assortativity levels (-0.4, -0.2, 0.0, +0.2, +0.4)
  * x 3 replicate rewirings  =>  60 graphs
  * each graph: Ricci flow to saturation (or 200 steps), record plateau ratio
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

import ricci_flow_program_100data as R                    # noqa: E402
from ricci_flow_longrun import compute_curvature, build_base_args  # noqa: E402

CFG = {
    "alpha": 0.5,
    "dt": 0.05,
    "sinkhorn_reg": 0.05,
    "sinkhorn_max_iter": 5000,
    "stop_thr": 1e-9,
    "solver": "batched",
}


def base_graph(n, mean_degree, seed):
    """Configuration model with a mild degree spread, largest component kept."""
    rng = np.random.default_rng(seed)
    # shifted Poisson-ish degrees, minimum 2 so the LCC stays dense enough
    degrees = np.maximum(rng.poisson(mean_degree, n), 2)
    if degrees.sum() % 2:
        degrees[0] += 1
    graph = nx.configuration_model(degrees.tolist(), seed=seed)
    graph = nx.Graph(graph)
    graph.remove_edges_from(list(nx.selfloop_edges(graph)))
    comps = list(nx.connected_components(graph))
    graph = graph.subgraph(max(comps, key=len)).copy()
    return nx.convert_node_labels_to_integers(graph)


def try_swap(graph, rng):
    edges = list(graph.edges())
    if len(edges) < 4:
        return False
    i, j = rng.choice(len(edges), 2, replace=False)
    u, v = edges[i]
    x, y = edges[j]
    if len({u, v, x, y}) < 4:
        return False
    if graph.has_edge(u, y) or graph.has_edge(x, v):
        return False
    graph.remove_edges_from([(u, v), (x, y)])
    graph.add_edges_from([(u, y), (x, v)])
    return True


def rewire_to_target(graph, target_r, rng, max_attempts=6000, tol=0.03):
    g = graph.copy()
    try:
        current = nx.degree_assortativity_coefficient(g)
    except Exception:
        return g, np.nan
    if not np.isfinite(current):
        return g, np.nan
    best = g.copy()
    best_r = current
    for _ in range(max_attempts):
        if abs(current - target_r) <= tol:
            break
        candidate = g.copy()
        if not try_swap(candidate, rng):
            continue
        try:
            new_r = nx.degree_assortativity_coefficient(candidate)
        except Exception:
            continue
        if not np.isfinite(new_r):
            continue
        if abs(new_r - target_r) < abs(current - target_r):
            g = candidate
            current = new_r
            if abs(current - target_r) < abs(best_r - target_r):
                best = g.copy()
                best_r = current
    return best, best_r


def run_flow(graph, max_steps=200, tol=1e-4, patience=5, var_floor=1e-6,
             tol_mode="current"):
    """Run the flow until it stops changing.

    ``tol_mode="initial"`` reproduces the historical (and unsafe) rule, which
    divides the per-step move by the variance at step 0; because that number
    is fixed while the variance decays, ANY slowly decaying run is eventually
    declared saturated.  ``tol_mode="current"`` divides by the current
    variance, so ``tol`` becomes a bound on the per-step exponential decay
    rate (tol=1e-5: the variance moves by <1% over the next 1000 steps).
    """
    args = build_base_args([])
    workspace = R.CurvatureWorkspace(
        graph, alpha=CFG["alpha"], backend=args.distance_backend,
        max_distance_mb=args.max_distance_mb)
    var0 = None
    prev = None
    stable = 0
    status = "horizon_reached"
    for step in range(max_steps):
        try:
            mean_k, var_k, edge_k, _, _ = compute_curvature(
                workspace, CFG, CFG["solver"])
        except Exception as exc:
            return {"steps": step, "var_ratio": np.nan, "status": "failed",
                    "message": repr(exc)[:100], "mean_kappa_final": np.nan}
        if step == 0:
            var0 = float(var_k)
        ratio = float(var_k) / var0 if var0 and var0 > 0 else np.nan
        if step > 0 and var0 and var0 > 0:
            denom = float(var_k) if tol_mode == "current" else var0
            move = abs(float(var_k) - prev) / max(denom, 1e-300)
            stable = stable + 1 if move < tol else 0
            if stable >= patience:
                status = ("saturated_current" if tol_mode == "current"
                          else "saturated")
                return {"steps": step + 1, "var_ratio": ratio, "status": status,
                        "message": "", "mean_kappa_final": float(mean_k)}
        if var0 and var0 > 0 and float(var_k) <= var_floor * var0:
            status = "constant_curvature"
            return {"steps": step + 1, "var_ratio": ratio, "status": status,
                    "message": "", "mean_kappa_final": float(mean_k)}
        prev = float(var_k)
        try:
            R.ricci_flow_step(graph, edge_k, dt=CFG["dt"])
        except Exception as exc:
            return {"steps": step + 1, "var_ratio": ratio, "status": "diverged",
                    "message": repr(exc)[:100], "mean_kappa_final": float(mean_k)}
    return {"steps": max_steps, "var_ratio": ratio, "status": status,
            "message": "", "mean_kappa_final": float(mean_k)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=os.path.join(
        SOURCE_DIR, "synthetic_results", "assortativity_sweep.csv"))
    parser.add_argument("--n", type=int, default=80)
    parser.add_argument("--mean_degree", type=int, default=6)
    parser.add_argument("--levels", default="-0.4,-0.2,0.0,0.2,0.4")
    parser.add_argument("--bases", type=int, default=4)
    parser.add_argument("--reps", type=int, default=3)
    parser.add_argument("--max_steps", type=int, default=200)
    parser.add_argument("--seed0", type=int, default=1000)
    parser.add_argument("--base_list", default="",
                        help="comma separated base indices for this worker")
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
        degree_seq = sorted(d for _, d in graph.degree())
        for level in levels:
            for rep in range(args.reps):
                rng = np.random.default_rng(
                    args.seed0 + 97 * base_index + 13 * rep + int(level * 10) + 1000)
                rewired, achieved_r = rewire_to_target(graph, level, rng)
                try:
                    trans = nx.transitivity(rewired)
                    r_final = nx.degree_assortativity_coefficient(rewired)
                except Exception:
                    trans, r_final = np.nan, np.nan
                res = run_flow(rewired.copy(), max_steps=args.max_steps)
                rows.append({
                    "base": base_index, "target_r": level, "rep": rep,
                    "achieved_r": achieved_r if achieved_r == achieved_r else np.nan,
                    "final_r": r_final, "transitivity": trans,
                    "nodes": rewired.number_of_nodes(),
                    "edges": rewired.number_of_edges(),
                    "degree_seq_hash": hash(tuple(degree_seq)) & 0xFFFF,
                    "var_ratio": res["var_ratio"], "steps": res["steps"],
                    "status": res["status"], "message": res["message"],
                    "mean_kappa_final": res["mean_kappa_final"],
                })
                done += 1
                # incremental write: a long sweep must survive an interrupt
                if done % 1 == 0:
                    pd.DataFrame(rows).to_csv(args.out, index=False)
                elapsed = time.time() - t_start
                print("[{}/{}] base={} target_r={:+.1f} rep={} -> r={:+.3f} "
                      "Q={:.3f} ratio={:.4f} ({} steps, {:.1f}s elapsed)".format(
                          done, total, base_index, level, rep,
                          r_final if r_final == r_final else float("nan"),
                          trans if trans == trans else float("nan"),
                          res["var_ratio"] if res["var_ratio"] == res["var_ratio"]
                          else float("nan"), res["steps"], elapsed), flush=True)

    frame = pd.DataFrame(rows)
    frame.to_csv(args.out, index=False)
    print("\nwritten: {}".format(args.out), flush=True)

    # quick summary
    ok = frame[np.isfinite(frame["var_ratio"])]
    print("\n--- grouped by achieved assortativity quartile ---", flush=True)
    if len(ok) > 8:
        ok = ok.copy()
        ok["r_bin"] = pd.qcut(ok["final_r"], 4, duplicates="drop")
        print(ok.groupby("r_bin", observed=True).agg(
            n=("var_ratio", "count"),
            median_ratio=("var_ratio", "median"),
            median_transitivity=("transitivity", "median")).to_string(), flush=True)


if __name__ == "__main__":
    main()
