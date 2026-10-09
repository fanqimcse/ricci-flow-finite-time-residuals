"""Independent-optimisation control for the stalled residual.

Question posed by the reviewer (their central point):

    Is the observed stalled curvature variance a *dynamical failure* of the
    Ricci flow, or the *geometric minimum* imposed by the fixed topology?

The Ricci flow is one particular (and possibly badly conditioned) dynamical
route towards uniform curvature.  We therefore ask, independently of the flow,

    V*_G = inf_{w : E -> R_+}  Var_E[ kappa(G, w) ]

by direct numerical minimisation over log-weights theta = log w, using exactly
the same curvature oracle (same alpha, same Sinkhorn solver, same regularisation)
as the flow itself.  The oracle is identical, so the comparison isolates the
*search procedure*, not the discretisation:

  * V_opt  ~=  V_flow  > 0  ->  the residual is a property of the edge set
                                (geometric infeasibility of uniformisation)
  * V_opt  <<  V_flow       ->  the flow itself failed to reach the achievable
                                minimum (dynamical / numerical failure)

kappa is invariant under a global rescaling of w (distances and transport costs
scale together), so the objective has one flat direction; we remove it by
projecting theta onto sum(theta) = 0.

Usage
-----
    python exp_indep_opt.py --probe            # timing only
    python exp_indep_opt.py --out exp_indep/   # full run
"""

from __future__ import print_function

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

SOURCE_DIR = os.path.dirname(os.path.abspath(__file__))
if SOURCE_DIR not in sys.path:
    sys.path.insert(0, SOURCE_DIR)

import ricci_flow_program_100data as R
import ricci_flow_longrun as L
from fast_sinkhorn import batched_sinkhorn_cost

DATA_DIR = r"D:\100data"

# Networks covering both branches.  All are small enough that a full |E|-dim
# optimisation is affordable.  fate / v_end come from strict_results/fate_table.csv.
TARGETS = [
    # stalled branch
    ("ENZYMES_g106", "stalled"),
    ("ENZYMES_g105", "stalled"),
    ("scc_rt_oman", "stalled"),
    ("scc2034_kilifi_all_contacts_across_households", "stalled"),
    # converging branch (controls: here the flow already succeeded)
    ("ENZYMES_g10", "converging"),
    ("ENZYMES_g111", "converging"),
    ("scc_rt_ksa", "converging"),
    ("aves-weaver-social-14", "converging"),
]

REG = 0.05
MAX_ITER = 5000
STOP_THR = 1e-9
ALPHA = 0.5


def run_root_find(oracle, n_dim, starts, max_evals, clip=30.0):
    """Independent test of a *different algorithmic class*.

    The flow descends Var(kappa) by gradient-like steps.  Here we instead solve
    the uniformisation equation directly,

        F_e(theta) = kappa_e(exp theta) - mean_f kappa_f(exp theta) = 0,

    with a damped Newton / trust-region least-squares solver.  This asks the
    sharper question "does a constant-curvature metric exist on this edge set",
    which is exactly what the convergence literature asserts.  A solver that
    drives ||F|| to ~0 implies uniformisation is attainable and the stalled
    residual is dynamical; a solver that stalls at a finite ||F|| is evidence
    (though not proof) that the edge set obstructs it.
    """
    used = [0]
    limit = [0]

    def resid(x):
        used[0] += 1
        if used[0] > limit[0]:
            raise _Budget
        k = oracle.curvature(project(x))
        if k is None:
            return np.full(n_dim, 1e3)
        r = k - k.mean()
        if oracle.best[0] > float(k.var()):
            oracle.best = (float(k.var()), np.asarray(x, dtype=float).copy())
        return r

    out = []
    per_start = max(60, int(max_evals // max(1, len(starts))))
    for label, x0 in starts:
        limit[0] = min(max_evals, used[0] + per_start)
        try:
            from scipy.optimize import least_squares
            res = least_squares(resid, project(x0),
                                bounds=(-clip, clip),
                                xtol=1e-12, ftol=1e-12, gtol=1e-12,
                                max_nfev=per_start)
            out.append((label, float(np.abs(res.fun).max()),
                        float(res.fun.var()), project(res.x).copy()))
        except _Budget:
            break
        except Exception as exc:
            print("      root start %s failed: %r" % (label, exc))
    if not out:
        if np.isfinite(oracle.best[0]) and oracle.best[1] is not None:
            return ("partial", np.nan, float(oracle.best[0]),
                    project(oracle.best[1]).copy()), used[0]
        return None
    out.sort(key=lambda t: t[2])
    return out[0], used[0]


class _Budget(Exception):
    pass


def find_file(dataset_id):
    for path in R.discover_graph_files(DATA_DIR):
        if R.safe_dataset_id(path, DATA_DIR) == dataset_id:
            return path
    return None


class Oracle(object):
    """Var_E[kappa(G, exp(theta))] with the flow's own curvature oracle.

    ``search`` mode relaxes the Sinkhorn tolerance so a |E|-dimensional
    optimisation is affordable; every candidate optimum is re-scored with the
    *exact* solver settings used by the flow before it is reported.
    """

    def __init__(self, graph, clip=30.0, stop_thr=STOP_THR, max_iter=MAX_ITER,
                 search=False):
        self.graph = graph
        self.clip = float(clip)
        self.stop_thr = float(stop_thr)
        self.max_iter = int(max_iter)
        self.search = search
        self.ws = R.CurvatureWorkspace(graph, ALPHA, "auto", 256.0)
        self.edges = list(self.ws.edges)
        self.n_evals = 0
        self.best = (np.inf, None)

    def _set_weights(self, theta):
        w = np.exp(np.clip(theta, -self.clip, self.clip))
        for (u, v), value in zip(self.edges, w):
            self.graph[u][v]["weight"] = float(value)

    def curvature(self, theta):
        self._set_weights(theta)
        try:
            self.ws.prepare_step()
            mass_u, mass_v, costs, dists = [], [], [], []
            for i in range(len(self.edges)):
                _, _, mu, mv, cost, d = self.ws.edge_payload(i)
                mass_u.append(mu)
                mass_v.append(mv)
                costs.append(cost)
                dists.append(d)
            values = batched_sinkhorn_cost(mass_u, mass_v, costs,
                                           REG, self.max_iter, self.stop_thr)
            k = np.array([1.0 - float(val) / float(d)
                          for val, d in zip(values, dists)], dtype=float)
        except Exception:
            return None
        if not np.all(np.isfinite(k)):
            return None
        return k

    def objective(self, theta):
        """Returns (variance, mean_curvature) or (inf, nan) on failure."""
        self.n_evals += 1
        k = self.curvature(theta)
        if k is None:
            return np.inf, np.nan
        var = float(k.var())
        if var < self.best[0]:
            self.best = (var, np.asarray(theta, dtype=float).copy())
        return var, float(k.mean())


def project(theta):
    """Remove the scale-invariance null direction."""
    theta = np.asarray(theta, dtype=float)
    return theta - theta.mean()


def run_optimisation(oracle, n_dim, starts, max_evals, seed=0, clip=30.0):
    """Multi-start L-BFGS-B on log-weights, projected to sum(theta)=0.

    The evaluation budget is split evenly across starts; if a start exhausts it
    mid-flight we keep the best point the oracle has seen so far instead of
    discarding the whole cell.
    """
    rng = np.random.default_rng(seed)
    results = []
    used = [0]
    limit = [0]

    def f(x):
        used[0] += 1
        if used[0] > limit[0]:
            raise _Budget
        return oracle.objective(project(x))[0]

    per_start = max(40, int(max_evals // max(1, len(starts))))
    # L-BFGS-B costs (n_dim + 1) evaluations per iteration when the gradient
    # is finite-differenced
    maxiter = max(5, int(per_start // max(1, n_dim + 1)))

    for label, x0 in starts:
        limit[0] = min(max_evals, used[0] + per_start)
        try:
            from scipy.optimize import minimize
            res = minimize(f, project(x0), method="L-BFGS-B",
                           bounds=[(-clip, clip)] * n_dim,
                           options=dict(maxiter=maxiter, maxfun=per_start,
                                        maxcor=15, ftol=1e-15, gtol=1e-12,
                                        eps=1e-5))
            results.append((label, float(res.fun), project(res.x).copy()))
        except _Budget:
            break
        except Exception as exc:
            print("      start %s failed: %r" % (label, exc))

    # derivative-free polish from the best point (kappa is only piecewise
    # smooth in w because shortest paths switch)
    if results:
        results.sort(key=lambda t: t[1])
        best_label, best_v, best_x = results[0]
        limit[0] = min(max_evals, used[0] + max(60, int(max_evals * 0.3)))
        try:
            from scipy.optimize import minimize
            res = minimize(f, best_x, method="Powell",
                           bounds=[(-clip, clip)] * n_dim,
                           options=dict(maxiter=max(2, int(max_evals * 0.3)),
                                        maxfev=max(60, int(max_evals * 0.3)),
                                        xtol=1e-6, ftol=1e-12))
            if res.fun < best_v:
                results.append((best_label + "+powell", float(res.fun),
                                project(res.x).copy()))
                results.sort(key=lambda t: t[1])
        except _Budget:
            pass
        except Exception:
            pass

    if not results:
        # budget exhausted inside the very first start: keep the best point seen
        if np.isfinite(oracle.best[0]) and oracle.best[1] is not None:
            return ("partial", float(oracle.best[0]),
                    project(oracle.best[1]).copy()), used[0]
        return None
    return results[0], used[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="exp_indep_opt")
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--clip", type=float, default=30.0)
    ap.add_argument("--max_evals", type=int, default=6000)
    ap.add_argument("--n_random", type=int, default=6)
    ap.add_argument("--fast_thr", type=float, default=1e-6)
    ap.add_argument("--fast_iter", type=int, default=200)
    ap.add_argument("--only", type=str, default="")
    ap.add_argument("--budget", type=float, default=200.0,
                    help="wall-clock seconds per network")
    ap.add_argument("--mode", choices=["min", "root"], default="min")
    ap.add_argument("--tag", default="")
    args = ap.parse_args()

    base = L.build_base_args(["--data_dir", DATA_DIR])
    os.makedirs(args.out, exist_ok=True)
    rows = []
    wanted = [t for t in args.only.split(",") if t]

    for dataset_id, fate in TARGETS:
        if wanted and dataset_id not in wanted:
            continue
        path = find_file(dataset_id)
        if path is None:
            print("MISSING %s" % dataset_id)
            continue
        try:
            graph = R.load_graph(path, base)
        except Exception as exc:
            print("LOAD FAIL %s: %r" % (dataset_id, exc))
            continue
        if graph is None:
            print("LOAD FAIL %s" % dataset_id)
            continue

        # two oracles over the SAME graph: cheap one for the search, exact
        # one (identical to the flow) for scoring
        fast = Oracle(graph, clip=args.clip, stop_thr=args.fast_thr,
                      max_iter=args.fast_iter, search=True)
        exact = Oracle(graph, clip=args.clip)
        oracle = fast
        n = len(oracle.edges)
        t0 = time.time()
        v0, m0 = exact.objective(np.zeros(n))
        t_eval = time.time() - t0
        t0 = time.time()
        fast.objective(np.zeros(n))
        t_fast = time.time() - t0
        print("[%s] %s  n=%d m=%d  V_exact(w=1)=%.6e  eval_exact=%.3fs eval_fast=%.3fs"
              % (fate[:4], dataset_id, graph.number_of_nodes(), n, v0,
                 t_eval, t_fast))
        if args.probe:
            rows.append(dict(dataset=dataset_id, fate=fate, nodes=graph.number_of_nodes(),
                             edges=n, V_unit=v0, mean_k_unit=m0, eval_s=t_eval))
            continue

        starts = [("unit", np.zeros(n))]
        rng = np.random.default_rng(12345)
        for i in range(args.n_random):
            starts.append(("rand%d" % i, rng.normal(0.0, 3.0, n)))
        # a start that mimics the strongly heterogeneous weights the flow finds
        starts.append(("spread", rng.normal(0.0, 8.0, n)))

        # keep every network inside a wall-clock budget
        if args.budget > 0 and t_fast > 0:
            max_evals = int(max(60, min(args.max_evals,
                                        args.budget / max(t_fast, 1e-3))))
        else:
            max_evals = args.max_evals
        print("      budget: %d evals (%.0fs)" % (max_evals, max_evals * t_fast))

        t0 = time.time()
        if args.mode == "root":
            out = run_root_find(oracle, n, starts, max_evals, clip=args.clip)
        else:
            out = run_optimisation(oracle, n, starts, max_evals,
                                   clip=args.clip)
        elapsed = time.time() - t0
        if out is None:
            print("      optimisation failed")
            continue
        (label, v_search, res_norm, x_opt), used = out
        # re-score the optimum with the EXACT solver settings of the flow
        v_opt, mean_k = exact.objective(project(x_opt))
        v_unit_exact = v0
        span = float(np.exp(np.clip(project(x_opt), -args.clip, args.clip)).max() /
                     np.exp(np.clip(project(x_opt), -args.clip, args.clip)).min())
        print("      best start=%-8s max|F|=%s  V_opt=%.6e  (V_unit=%.6e,"
              " ratio=%.3e)  mean_k=%+.4f  span=%.2e  evals=%d  %.1fs"
              % (label, ("%.3e" % res_norm) if res_norm == res_norm else "n/a",
                 v_opt, v_unit_exact, v_opt / v_unit_exact, mean_k,
                 span, used, elapsed))
        rows.append(dict(dataset=dataset_id, fate=fate,
                         nodes=graph.number_of_nodes(), edges=n,
                         V_unit=v_unit_exact, mean_k_unit=m0, V_opt=v_opt,
                         ratio_opt_over_unit=v_opt / v_unit_exact,
                         mean_k_opt=mean_k,
                         weight_span_opt=span, best_start=label,
                         max_abs_resid=res_norm,
                         n_evals=used, seconds=elapsed, clip=args.clip,
                         mode=args.mode))

        # incremental checkpoint: a killed process must not lose finished cells
        if not args.probe:
            _save(pd.DataFrame(rows), args.out, args.tag)
            print("      [checkpoint written]")

    out_csv = _save(pd.DataFrame(rows), args.out, args.tag)
    print("\nwrote %s" % out_csv)
    print(pd.DataFrame(rows).to_string(index=False))


def _save(df, out_dir, tag=""):
    out_csv = os.path.join(out_dir, "indep_opt%s.csv" % ("_" + tag if tag else ""))
    key = ["dataset", "mode"] if "mode" in df.columns else ["dataset"]
    if os.path.exists(out_csv):
        try:
            prior = pd.read_csv(out_csv)
            cols = set(prior.columns) & set(df.columns)
            if "V_opt" in cols:
                df = pd.concat([prior, df], ignore_index=True)
                df = df.sort_values("V_opt").drop_duplicates(key, keep="first")
        except Exception:
            pass
    df.to_csv(out_csv, index=False)
    return out_csv


if __name__ == "__main__":
    main()
