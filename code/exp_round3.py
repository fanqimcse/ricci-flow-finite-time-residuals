"""Round-3 experiments (reviewer round 2), four decisive numerics.

R2-1  fixed-eps Sinkhorn breaks the scale-invariance argument used in the
      manuscript (kappa invariant under w -> c w).  We *measure* the violation
      for three estimators:  fixed-eps Sinkhorn, scale-equivariant Sinkhorn
      (eps proportional to the cost scale), and exact W_1.
R2-2  exact optimal transport (network simplex) instead of Sinkhorn: does the
      residual survive when the regularisation is removed entirely?
R2-3  scale-equivariant Sinkhorn schedule eps_t = eps_0 * s_t / s_0.
R2-4  Li-Muench-compatible update d <- (1 - alpha kappa) d  (first order Euler
      on the metric) versus the exponential update used throughout the paper.

Modes:  scale | exact | eqreg | linear
"""

from __future__ import print_function

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import networkx as nx

SOURCE_DIR = os.path.dirname(os.path.abspath(__file__))
if SOURCE_DIR not in sys.path:
    sys.path.insert(0, SOURCE_DIR)

import ricci_flow_program_100data as R
import ricci_flow_longrun as L
from fast_sinkhorn import batched_sinkhorn_cost

DATA_DIR = r"D:\100data"
ALPHA = 0.5
DT = 0.05
REG = 0.05

TARGETS = [
    ("ENZYMES_g105", "stalled"),
    ("ENZYMES_g106", "stalled"),
    ("scc2034_kilifi_all_contacts_across_households", "stalled"),
    ("scc_rt_oman", "stalled"),
    ("scc_rt_ksa", "converging"),
    ("reptilia-tortoise-network-hw", "converging"),
    ("aves-weaver-social-14", "converging"),
    ("scc_rt_http", "converging"),
    ("ENZYMES_g10", "converging"),
]


# --------------------------------------------------------------------------
# transport solvers
# --------------------------------------------------------------------------
def _exact_costs(mu_list, mv_list, cost_list):
    """Exact W_1 by network simplex (POT ot.emd); falls back to HiGHS LP."""
    try:
        import ot
    except Exception:
        ot = None
    out = []
    for a, b, M in zip(mu_list, mv_list, cost_list):
        a = np.asarray(a, dtype=np.float64)
        b = np.asarray(b, dtype=np.float64)
        M = np.asarray(M, dtype=np.float64)
        a = a / a.sum()
        b = b / b.sum()
        if ot is not None:
            try:
                G = ot.emd(a, b, M, numItermax=1000000)
                out.append(float((G * M).sum()))
                continue
            except Exception:
                pass
        n, m = M.shape
        rows = []
        for i in range(n):
            row = np.zeros((n, m))
            row[i, :] = 1.0
            rows.append(row.ravel())
        for j in range(m):
            row = np.zeros((n, m))
            row[:, j] = 1.0
            rows.append(row.ravel())
        from scipy.optimize import linprog
        res = linprog(M.ravel(), A_eq=np.asarray(rows),
                      b_eq=np.concatenate([a, b]), bounds=(0, None),
                      method="highs")
        out.append(float(res.fun))
    return np.asarray(out, dtype=np.float64)


def _scale_of(cost_list):
    """Mean entry of the cost matrices: the natural scale of the problems."""
    vals = [float(np.mean(np.asarray(c, dtype=float))) for c in cost_list
            if np.size(c) > 0]
    return float(np.mean(vals)) if vals else 1.0


class Solver(object):
    """Curvature evaluator with a selectable transport backend."""

    def __init__(self, ws, mode, reg=REG, s0=None):
        self.ws = ws
        self.mode = mode
        self.reg = reg
        self.s0 = s0 if s0 else 1.0

    def payload(self):
        self.ws.prepare_step()
        mu, mv, costs, dists = [], [], [], []
        for i in range(len(self.ws.edges)):
            _, _, a, b, c, d = self.ws.edge_payload(i)
            mu.append(a)
            mv.append(b)
            costs.append(c)
            dists.append(d)
        return mu, mv, costs, np.asarray(dists, dtype=float)

    def curvature(self, return_scale=False):
        mu, mv, costs, dists = self.payload()
        scale = _scale_of(costs)
        if self.mode == "sinkhorn":
            reg = self.reg
        elif self.mode == "exact":
            reg = None
        elif self.mode == "eqreg":
            reg = self.reg * scale / self.s0
        else:
            raise ValueError(self.mode)
        if reg is None:
            vals = _exact_costs(mu, mv, costs)
        else:
            vals = batched_sinkhorn_cost(mu, mv, costs, reg, 5000, 1e-9)
        k = np.array([1.0 - float(v) / float(d)
                      for v, d in zip(vals, dists)], dtype=float)
        if return_scale:
            return k, scale
        return k


# --------------------------------------------------------------------------
# graph helpers
# --------------------------------------------------------------------------
def find_file(dataset_id):
    for path in R.discover_graph_files(DATA_DIR):
        if R.safe_dataset_id(path, DATA_DIR) == dataset_id:
            return path
    return None


def setup(dataset_id, base):
    path = find_file(dataset_id)
    if path is None:
        return None
    graph = R.load_graph(path, base)
    ws = R.CurvatureWorkspace(graph, ALPHA, "auto", 256.0)
    edges = list(ws.edges)
    comp = {}
    for i, c in enumerate(nx.connected_components(graph)):
        for n in c:
            comp[n] = i
    comp_of_edge = np.array([comp[u] for u, v in edges], dtype=int)
    return graph, ws, edges, comp_of_edge


def decompose(k, comp_of_edge):
    k = np.asarray(k, dtype=float)
    pooled = float(k.var())
    labels = sorted(set(comp_of_edge.tolist()))
    counts, means, wsum, n = [], [], 0.0, 0
    for c in labels:
        kc = k[comp_of_edge == c]
        if kc.size == 0:
            continue
        counts.append(int(kc.size))
        means.append(float(kc.mean()))
        wsum += float(kc.var()) * kc.size
        n += kc.size
    within = wsum / max(n, 1)
    counts = np.asarray(counts, dtype=float)
    means = np.asarray(means, dtype=float)
    if counts.sum():
        gm = float(np.average(means, weights=counts))
        between = float(np.average((means - gm) ** 2, weights=counts))
    else:
        between = 0.0
    return pooled, within, between


def weight_span(graph):
    w = np.array([float(d.get("weight", 1.0))
                  for _, _, d in graph.edges(data=True)], dtype=float)
    if w.size == 0 or np.min(w) <= 0:
        return float("nan")
    return float(np.max(w) / np.min(w))


def rescale(graph, factor):
    for _, _, d in graph.edges(data=True):
        d["weight"] = float(d.get("weight", 1.0)) * factor


# --------------------------------------------------------------------------
# mode: scale-invariance probe
# --------------------------------------------------------------------------
def run_scale(args):
    base = L.build_base_args(["--data_dir", DATA_DIR])
    rows = []
    factors = [0.25, 0.5, 2.0, 4.0, 10.0, 100.0]
    for dataset_id, fate in TARGETS:
        if args.only and dataset_id not in args.only.split(","):
            continue
        st = setup(dataset_id, base)
        if st is None:
            print("MISSING", dataset_id)
            continue
        graph, ws, edges, comp_of_edge = st
        sol_s = Solver(ws, "sinkhorn")
        sol_e = Solver(ws, "exact")
        sol_q = Solver(ws, "eqreg")
        for phase, warmup in (("initial", 0), ("late", args.warmup)):
            st2 = setup(dataset_id, base)
            graph, ws, edges, comp_of_edge = st2
            sol_s = Solver(ws, "sinkhorn")
            sol_e = Solver(ws, "exact")
            sol_q = Solver(ws, "eqreg")
            if warmup:
                k = sol_s.curvature()
                for _ in range(warmup):
                    try:
                        R.ricci_flow_step(graph, dict(zip(edges, k)), DT)
                    except Exception as exc:
                        print("   [warmup stop]", repr(exc))
                        break
                    k = sol_s.curvature()
            _, s0 = sol_s.curvature(return_scale=True)
            sol_q = Solver(ws, "eqreg", s0=s0)
            k_ref = {"sinkhorn": sol_s.curvature(),
                     "exact": sol_e.curvature(),
                     "eqreg": sol_q.curvature()}
            for name in ("sinkhorn", "exact", "eqreg"):
                base_k = k_ref[name]
                for c in factors:
                    rescale(graph, c)
                    sol = {"sinkhorn": Solver(ws, "sinkhorn"),
                           "exact": Solver(ws, "exact"),
                           "eqreg": Solver(ws, "eqreg", s0=s0)}[name]
                    kk = sol.curvature()
                    dk = float(np.max(np.abs(kk - base_k)))
                    dv = float(abs(kk.var() - base_k.var()))
                    rows.append(dict(dataset=dataset_id, fate=fate,
                                     phase=phase, warmup=warmup,
                                     estimator=name, factor=c,
                                     max_abs_dkappa=dk,
                                     d_variance=dv,
                                     var_ref=float(base_k.var())))
                    rescale(graph, 1.0 / c)
            print("  %-52s %-8s span=%.2e" %
                  (dataset_id[:52], phase, weight_span(graph)))
    df = pd.DataFrame(rows)
    os.makedirs(args.out, exist_ok=True)
    df.to_csv(os.path.join(args.out, "scale_probe.csv"), index=False)
    print(df.groupby(["estimator", "phase"])[["max_abs_dkappa", "d_variance"]]
          .max().to_string())
    print("wrote", os.path.join(args.out, "scale_probe.csv"))


# --------------------------------------------------------------------------
# generic integrator used by exact / eqreg / linear
# --------------------------------------------------------------------------
def integrate(dataset_id, fate, base, mode, steps, budget, alpha_lin=0.5,
              renormalise=True, out_dir=None):
    st = setup(dataset_id, base)
    if st is None:
        return None
    graph, ws, edges, comp_of_edge = st
    sol = Solver(ws, "sinkhorn")
    _, s0 = sol.curvature(return_scale=True)
    if mode == "sinkhorn":
        solver = Solver(ws, "sinkhorn")
    elif mode == "exact":
        solver = Solver(ws, "exact")
    elif mode == "eqreg":
        solver = Solver(ws, "eqreg", s0=s0)
    else:
        solver = Solver(ws, "sinkhorn")

    k, scale = solver.curvature(return_scale=True)
    p0, w0, b0 = decompose(k, comp_of_edge)
    hist = [dict(step=0, V=k.var(), V_within=w0, span=weight_span(graph),
                 scale=scale, mean_k=float(k.mean()))]
    t0 = time.time()
    step = 0
    stop = "horizon"
    try:
        while step < steps and time.time() - t0 < budget:
            if mode == "linear":
                for (u, v), kv in zip(edges, k):
                    old = float(graph[u][v].get("weight", 1.0))
                    new = old * (1.0 - alpha_lin * float(kv))
                    if not np.isfinite(new) or new <= 0.0:
                        raise RuntimeError("non-positive weight in linear update")
                    graph[u][v]["weight"] = new
                if renormalise:
                    R.normalize_weights(graph)
            else:
                R.ricci_flow_step(graph, dict(zip(edges, k)), DT)
            k, scale = solver.curvature(return_scale=True)
            p, w, b = decompose(k, comp_of_edge)
            hist.append(dict(step=step + 1, V=p, V_within=w,
                             span=weight_span(graph), scale=scale,
                             mean_k=float(k.mean())))
            step += 1
    except Exception as exc:
        stop = "error: " + repr(exc)[:120]
    p1, w1, b1 = decompose(k, comp_of_edge)
    rec = dict(dataset=dataset_id, fate=fate, mode=mode,
               alpha_linear=(alpha_lin if mode == "linear" else ""),
               renormalise=renormalise, steps=step, stop=stop,
               V0=p0, V0_within=w0,
               V_end=p1, V_within=w1, V_between=b1,
               ratio=p1 / max(p0, 1e-300),
               ratio_within=w1 / max(w0, 1e-300),
               between_share=b1 / max(p1, 1e-300),
               span=weight_span(graph), mean_k_end=float(k.mean()),
               seconds=time.time() - t0)
    if out_dir:
        os.makedirs(os.path.join(out_dir, "traj"), exist_ok=True)
        pd.DataFrame(hist).to_csv(
            os.path.join(out_dir, "traj", "%s__%s.csv" % (dataset_id, mode)),
            index=False)
    return rec


def run_modes(args):
    base = L.build_base_args(["--data_dir", DATA_DIR])
    modes = args.modes.split(",")
    rows = []
    for mode in modes:
        for dataset_id, fate in TARGETS:
            if args.only and dataset_id not in args.only.split(","):
                continue
            if mode == "linear":
                for a_l in [float(x) for x in args.alpha_linear.split(",")]:
                    for ren in (True, False):
                        rec = integrate(dataset_id, fate, base, mode,
                                        args.steps, args.budget,
                                        alpha_lin=a_l, renormalise=ren,
                                        out_dir=args.out)
                        if rec:
                            rows.append(rec)
                            print("  %-46s %-8s a=%.2f ren=%d steps=%3d "
                                  "V=%.3e Vw=%.3e ratio=%.2e" %
                                  (dataset_id[:46], mode, a_l, ren,
                                   rec["steps"], rec["V_end"], rec["V_within"],
                                   rec["ratio"]))
            else:
                rec = integrate(dataset_id, fate, base, mode, args.steps,
                                args.budget, out_dir=args.out)
                if rec:
                    rows.append(rec)
                    print("  %-46s %-8s steps=%3d V=%.3e Vw=%.3e "
                          "ratio=%.2e span=%.1e" %
                          (dataset_id[:46], mode, rec["steps"], rec["V_end"],
                           rec["V_within"], rec["ratio"], rec["span"]))
    df = pd.DataFrame(rows)
    os.makedirs(args.out, exist_ok=True)
    tag = "_".join(modes)
    df.to_csv(os.path.join(args.out, "flow_%s.csv" % tag), index=False)
    print(df.to_string(index=False))
    print("wrote", os.path.join(args.out, "flow_%s.csv" % tag))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="scale")
    ap.add_argument("--modes", default="exact")
    ap.add_argument("--out", default="exp_round3")
    ap.add_argument("--only", default="")
    ap.add_argument("--steps", type=int, default=250)
    ap.add_argument("--budget", type=float, default=240.0)
    ap.add_argument("--warmup", type=int, default=200)
    ap.add_argument("--alpha_linear", default="0.5,0.05")
    args = ap.parse_args()
    if args.mode == "scale":
        run_scale(args)
    else:
        run_modes(args)


if __name__ == "__main__":
    main()
