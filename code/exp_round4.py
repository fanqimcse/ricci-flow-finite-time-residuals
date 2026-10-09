"""Round-4 experiments (reviewer round 4).

R4-P3   extended-horizon exact-transport run:   5,000 steps of exact W_1
        transport on the slowest connected networks (g105, g106, g10).
R4-P4   epsilon trajectory sensitivity:  eps in {0.1 ... 0.002} on g105/g106,
        recording the *whole* trajectory V(t), the local rate r_t, the weight
        span, and periodic exact-OT checkpoints (transport residual).
R4-P9   step-size convergence test at fixed *physical* time
        T = n_steps * dt, dt in {0.01, 0.025, 0.05, 0.1}, so that the
        comparison is not confounded by comparing different physical times.
R4-P15  V_G^* sanity checks on graph classes where symmetry forces
        constant curvature (complete graph, cycle) and on a disjoint union,
        where V_within = 0 but V_pooled = |E1||E2|/(|E1|+|E2|)^2 (k1-k2)^2.
"""

from __future__ import print_function

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

import ricci_flow_program_100data as R
import ricci_flow_longrun as L

from exp_round3 import (Solver, setup, decompose, weight_span,
                        _exact_costs, TARGETS, ALPHA, REG)

DT = 0.05


def _wmin_wmax(graph):
    w = np.array([float(d.get("weight", 1.0))
                  for _, _, d in graph.edges(data=True)], dtype=float)
    if w.size == 0:
        return float("nan"), float("nan")
    return float(w.min()), float(w.max())


def run_flow(dataset_id, base, mode="sinkhorn", reg=REG, steps=200,
             budget=300.0, dt=DT, ckpt=0, traj_path=None, tag=""):
    """Integrate the flow and return (summary dict, trajectory DataFrame)."""
    st = setup(dataset_id, base)
    if st is None:
        return None, None
    graph, ws, edges, comp_of_edge = st
    probe = Solver(ws, "sinkhorn")
    _, s0 = probe.curvature(return_scale=True)
    if mode == "sinkhorn":
        solver = Solver(ws, "sinkhorn", reg=reg)
    elif mode == "eqreg":
        solver = Solver(ws, "eqreg", reg=reg, s0=s0)
    else:
        solver = Solver(ws, "exact")
    exact_solver = Solver(ws, "exact")

    k, scale = solver.curvature(return_scale=True)
    p0, w0, b0 = decompose(k, comp_of_edge)
    wmin, wmax = _wmin_wmax(graph)
    hist = [dict(step=0, t=0.0, V=p0, V_within=w0, span=weight_span(graph),
                 w_min=wmin, w_max=wmax, mean_k=float(k.mean()))]
    t0 = time.time()
    step = 0
    stop = "horizon"
    try:
        while step < steps and (time.time() - t0) < budget:
            R.ricci_flow_step(graph, dict(zip(edges, k)), dt)
            k, scale = solver.curvature(return_scale=True)
            p, w, b = decompose(k, comp_of_edge)
            wmin, wmax = _wmin_wmax(graph)
            row = dict(step=step + 1, t=(step + 1) * dt, V=p, V_within=w,
                       span=weight_span(graph), w_min=wmin, w_max=wmax,
                       mean_k=float(k.mean()))
            if ckpt and (step + 1) % ckpt == 0:
                ke = exact_solver.curvature()
                pe, we, be = decompose(ke, comp_of_edge)
                row["V_exact"] = pe
                row["V_within_exact"] = we
                row["max_dkappa_vs_exact"] = float(np.max(np.abs(ke - k)))
            hist.append(row)
            step += 1
    except Exception as exc:
        stop = "error: " + repr(exc)[:120]
    p1, w1, b1 = decompose(k, comp_of_edge)
    ke = exact_solver.curvature()
    pe, we, be = decompose(ke, comp_of_edge)
    rec = dict(dataset=dataset_id, mode=mode, reg=reg, dt=dt, tag=tag,
               steps=step, T_phys=step * dt, stop=stop,
               V0=p0, V0_within=w0, V_end=p1, V_within=w1, V_between=b1,
               ratio=p1 / max(p0, 1e-300),
               ratio_within=w1 / max(w0, 1e-300),
               between_share=b1 / max(p1, 1e-300),
               span=weight_span(graph), mean_k_end=float(k.mean()),
               V_end_exact_check=pe, V_within_exact_check=we,
               max_dkappa_vs_exact=float(np.max(np.abs(ke - k))),
               seconds=time.time() - t0)
    traj = pd.DataFrame(hist)
    if traj_path:
        os.makedirs(os.path.dirname(traj_path), exist_ok=True)
        traj.to_csv(traj_path, index=False)
    return rec, traj


# ---------------------------------------------------------------- R4-P3 ----
def mode_long(args):
    base = L.build_base_args(["--data_dir", r"D:\100data"])
    out = args.out
    rows = []
    for ds in args.only.split(","):
        print("[long] %s  %d steps, exact transport" % (ds, args.steps),
              flush=True)
        rec, _ = run_flow(ds, base, mode="exact", steps=args.steps,
                          budget=args.budget, dt=DT, ckpt=args.ckpt,
                          traj_path=os.path.join(
                              out, "traj_long", "%s__exact%d.csv" % (ds, args.steps)),
                          tag="long%d" % args.steps)
        if rec:
            rows.append(rec)
            print("    steps=%d V=%.4e Vw=%.4e ratio=%.2e span=%.2e "
                  "(%.0fs)" % (rec["steps"], rec["V_end"], rec["V_within"],
                               rec["ratio"], rec["span"], rec["seconds"]),
                  flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out, "flow_exact_long%d.csv" % args.steps),
              index=False)
    print(df.to_string(index=False))


# ---------------------------------------------------------------- R4-P4 ----
def mode_eps(args):
    base = L.build_base_args(["--data_dir", r"D:\100data"])
    out = args.out
    eps_list = [float(x) for x in args.eps.split(",")]
    rows = []
    for ds in args.only.split(","):
        for eps in eps_list:
            print("[eps] %s eps=%.4g steps<=%d" % (ds, eps, args.steps),
                  flush=True)
            rec, _ = run_flow(
                ds, base, mode="sinkhorn", reg=eps, steps=args.steps,
                budget=args.budget, dt=DT, ckpt=args.ckpt,
                traj_path=os.path.join(
                    out, "traj_eps", "%s__reg%g.csv" % (ds, eps)),
                tag="eps%g" % eps)
            if rec:
                rows.append(rec)
                print("    steps=%3d V=%.4e Vw=%.4e |dk_vs_exact|=%.3e "
                      "span=%.1e (%.0fs)" %
                      (rec["steps"], rec["V_end"], rec["V_within"],
                       rec["max_dkappa_vs_exact"], rec["span"],
                       rec["seconds"]), flush=True)
            pd.DataFrame(rows).to_csv(os.path.join(out, "eps_sweep.csv"),
                                      index=False)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out, "eps_sweep.csv"), index=False)
    print(df.to_string(index=False))


# ---------------------------------------------------------------- R4-P9 ----
def mode_dt(args):
    base = L.build_base_args(["--data_dir", r"D:\100data"])
    out = args.out
    dts = [float(x) for x in args.dts.split(",")]
    T = args.T_phys
    rows = []
    for ds, fate in TARGETS:
        if args.only and ds not in args.only.split(","):
            continue
        for dt in dts:
            steps = int(round(T / dt))
            print("[dt] %s dt=%.3g -> %d steps (T=%.1f, %s)" %
                  (ds, dt, steps, T, args.solver), flush=True)
            rec, _ = run_flow(
                ds, base, mode=args.solver, steps=steps, budget=args.budget,
                dt=dt, ckpt=0,
                traj_path=os.path.join(
                    out, "traj_dt", "%s__%s_dt%g.csv" % (ds, args.solver, dt)),
                tag="dt%g" % dt)
            if rec:
                rows.append(rec)
                print("    steps=%4d V=%.4e ratio=%.2e (%.0fs)" %
                      (rec["steps"], rec["V_end"], rec["ratio"],
                       rec["seconds"]), flush=True)
            pd.DataFrame(rows).to_csv(
                os.path.join(out, "dt_sweep_%s.csv" % args.solver), index=False)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out, "dt_sweep_%s.csv" % args.solver), index=False)
    print(df.to_string(index=False))


# --------------------------------------------------------------- R4-P15 ----
def mode_sym(args):
    """V_G^* sanity checks: symmetry forces constant curvature."""
    import ricci_flow_program_100data as RR
    base = L.build_base_args(["--data_dir", r"D:\100data"])
    out = args.out
    rows = []

    def evaluate(graph, name, note):
        ws = RR.CurvatureWorkspace(graph, ALPHA, "auto", 256.0)
        comp = {}
        for i, c in enumerate(nx.connected_components(graph)):
            for n in c:
                comp[n] = i
        comp_of_edge = np.array([comp[u] for u, v in ws.edges], dtype=int)
        sol = Solver(ws, "exact")
        k = sol.curvature()
        p, w, b = decompose(k, comp_of_edge)
        rows.append(dict(case=name, note=note,
                         n=graph.number_of_nodes(), m=graph.number_of_edges(),
                         n_components=len(set(comp_of_edge.tolist())),
                         V_pooled=p, V_within=w, V_between=b,
                         kappa_min=float(k.min()), kappa_max=float(k.max())))
        print("  %-28s n=%3d m=%4d comp=%d  V_pooled=%.3e V_within=%.3e"
              % (name, graph.number_of_nodes(), graph.number_of_edges(),
                 len(set(comp_of_edge.tolist())), p, w))

    for n in (6, 9, 12):
        evaluate(nx.complete_graph(n), "K_%d" % n, "edge-transitive")
    for n in (6, 12, 20):
        evaluate(nx.cycle_graph(n), "C_%d" % n, "edge-transitive")
    for n in (6, 12):
        evaluate(nx.path_graph(n), "P_%d" % n, "not edge-transitive")
    for n in (4, 6):
        evaluate(nx.star_graph(n), "S_%d" % (n + 1), "not edge-transitive")

    # disjoint union of two cycles with different sizes -> constant curvature
    # inside each component, different constants between components
    for (n1, n2) in ((5, 9), (6, 20)):
        g = nx.disjoint_union(nx.cycle_graph(n1), nx.cycle_graph(n2))
        evaluate(g, "C_%d + C_%d" % (n1, n2),
                 "within=0 expected; pooled>0 if constants differ")
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out, "vstar_sanity.csv"), index=False)
    print(df.to_string(index=False))
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="long")
    ap.add_argument("--out", default="exp_round4")
    ap.add_argument("--only", default="ENZYMES_g105,ENZYMES_g106")
    ap.add_argument("--steps", type=int, default=5000)
    ap.add_argument("--budget", type=float, default=1800.0)
    ap.add_argument("--ckpt", type=int, default=100)
    ap.add_argument("--eps", default="0.1,0.05,0.02,0.01,0.005,0.002")
    ap.add_argument("--dts", default="0.01,0.025,0.05,0.1")
    ap.add_argument("--T_phys", type=float, default=20.0)
    ap.add_argument("--solver", default="exact")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    if args.mode == "long":
        mode_long(args)
    elif args.mode == "eps":
        mode_eps(args)
    elif args.mode == "dt":
        mode_dt(args)
    elif args.mode == "sym":
        mode_sym(args)


if __name__ == "__main__":
    main()
