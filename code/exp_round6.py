#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Round-6 experiments (second re-review of the Communications Physics draft).

R6-W   weight-semantics sensitivity (reviewer point 5).
       Twenty-eight datasets of the corpus carry a positive third column.  In
       the production run those values were interpreted as edge *lengths*
       (``positive_weight_interpretation = distance``).  That is a modelling
       choice, not a verified property of the data: in public repositories a
       positive edge value usually encodes an affinity, a count, a flow or a
       timestamp.  Here every weighted dataset is re-run under three initial
       conventions --

         asis     l_e = w_e                      (production)
         recip    l_e = 1 / w_e                  (affinity -> distance)
         unit     l_e = 1                        (topology only)

       with identical numerics and a matched step horizon, and the initial
       geometry (V0), the residual at the horizon and the late-window decay
       rate are compared across conventions.

R6-N   numerical reliability at extreme weight span (reviewer point 7).
       The 5,000-step exact-transport runs end with weight spans of
       1e65-1e72.  We measure, at those weights, how much of the reported
       residual is attributable to (a) the linear-programming backend,
       (b) floating-point conditioning of cost matrices of that magnitude,
       (c) the Sinkhorn tolerance, and (d) the initial weighting, and we
       record the per-checkpoint consistency checks with their pass/fail
       status.

R6-G   g10 reconciliation (reviewer point 2): per-network long-horizon
       diagnostics for the three connected ENZYMES graphs, so that the
       "two networks" case study can be reported next to the third run.
"""

from __future__ import print_function

import argparse
import copy
import json
import os
import platform
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

from exp_round3 import (_exact_costs, Solver, setup, decompose,
                        weight_span, DATA_DIR, ALPHA, REG)

DT = 0.05

# --------------------------------------------------------------------------
# weighted datasets of the corpus (positive third column, distance reading)
# --------------------------------------------------------------------------
WEIGHTED = [
    "aves-barn-swallow-non-physical",
    "aves-sparrow-social-2009",
    "aves-sparrow-social",
    "aves-sparrowlyon-flock-season2",
    "aves-sparrowlyon-flock-season3",
    "aves-thornbill-farine",
    "bio-CE-HT",
    "bio-CE-LC",
    "bio-DM-LC",
    "bio-SC-TS",
    "eco-mangwet",
    "eco-stmarks",
    "ia-dnc-corecipient",
    "ia-workplace-contacts",
    "insecta-ant-trophallaxis-colony1",
    "insecta-ant-trophallaxis-colony2",
    "mammalia-primate-association",
    "mammalia-raccoon-proximity",
    "mammalia-voles-bhp-trapping",
    "mammalia-voles-kcs-trapping",
    "mammalia-voles-plj-trapping",
    "reptilia-tortoise-network-bsv",
    "reptilia-tortoise-network-cs",
    "reptilia-tortoise-network-fi",
    "reptilia-tortoise-network-hw",
    "reptilia-tortoise-network-mc",
    "reptilia-tortoise-network-pv",
    "reptilia-tortoise-network-sg",
]

CONVENTIONS = [("asis", "auto", "distance"),
               ("recip", "strength", "distance"),
               ("unit", "unweighted", "distance")]


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
_FILES = None


def all_files():
    """Cached dataset listing (scanning D:\\100data costs ~11 s)."""
    global _FILES
    if _FILES is None:
        _FILES = R.discover_graph_files(DATA_DIR)
    return _FILES


def find_path(dataset_id):
    for p in all_files():
        if R.safe_dataset_id(p, DATA_DIR) == dataset_id:
            return p
    return None


def raw_third_column_stats(dataset_id):
    """Statistics of the raw third column *before* any interpretation."""
    path = find_path(dataset_id)
    if path is None:
        return None
    vals = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("%") or s.startswith("#"):
                continue
            parts = s.replace(",", " ").split()
            if len(parts) >= 3:
                try:
                    vals.append(float(parts[2]))
                except ValueError:
                    pass
    if not vals:
        return None
    v = np.asarray(vals, dtype=float)
    out = dict(dataset=dataset_id, n_values=int(v.size),
               raw_min=float(v.min()), raw_max=float(v.max()),
               raw_median=float(np.median(v)),
               raw_cv=float(v.std() / max(v.mean(), 1e-300)),
               raw_distinct=int(np.unique(np.round(v, 9)).size))
    out["all_integer"] = bool(np.all(np.isclose(v, np.round(v))))
    out["year_like"] = bool(out["all_integer"] and v.min() >= 1900
                            and v.max() <= 2100)
    return out


def load_with_convention(dataset_id, convention):
    base = L.build_base_args(["--data_dir", DATA_DIR])
    if convention == "asis":
        pass
    elif convention == "recip":
        base.weight_mode = "strength"
    elif convention == "unit":
        base.weight_mode = "unweighted"
    else:
        raise ValueError(convention)
    path = find_path(dataset_id)
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


def late_slopes(v, frac=0.25):
    """Per-step slope of log V over the final quarter (as in fate_analysis.py)."""
    v = np.asarray(v, dtype=float)
    n = len(v)
    a = int(np.floor(n * (1.0 - frac)))
    if n - a < 8:
        return float("nan")
    s = np.arange(n, dtype=float)
    msk = s >= a
    good = np.isfinite(np.log(np.clip(v[msk], 1e-300, None)))
    if good.sum() < 8:
        return float("nan")
    from scipy import stats
    return float(stats.linregress(s[msk][good],
                                  np.log(np.clip(v[msk], 1e-300, None))[good]
                                  ).slope)


# --------------------------------------------------------------------------
# R6-W  weight-semantics sensitivity
# --------------------------------------------------------------------------
def run_one_convention(dataset_id, convention, steps, budget):
    st = load_with_convention(dataset_id, convention)
    if st is None:
        return None
    graph, ws, edges, comp_of_edge = st
    solver = Solver(ws, "sinkhorn", reg=REG)
    k = solver.curvature()
    p0, w0, b0 = decompose(k, comp_of_edge)
    hist = [dict(step=0, V=p0, V_within=w0, V_between=b0,
                 span=weight_span(graph))]
    t0 = time.time()
    step = 0
    stop = "horizon"
    try:
        while step < steps and (time.time() - t0) < budget:
            R.ricci_flow_step(graph, dict(zip(edges, k)), DT)
            k = solver.curvature()
            p, w, b = decompose(k, comp_of_edge)
            hist.append(dict(step=step + 1, V=p, V_within=w, V_between=b,
                             span=weight_span(graph)))
            step += 1
    except Exception as exc:
        stop = "error: " + repr(exc)[:120]
    if (time.time() - t0) >= budget and step < steps:
        stop = "budget"
    traj = pd.DataFrame(hist)
    v = traj["V_within"].astype(float).values
    return dict(dataset=dataset_id, convention=convention, steps=step,
                stop=stop, V0=p0, V0_within=w0,
                V_end=float(traj["V"].iloc[-1]),
                V_within_end=float(traj["V_within"].iloc[-1]),
                ratio=float(traj["V"].iloc[-1]) / max(p0, 1e-300),
                ratio_within=float(traj["V_within"].iloc[-1]) / max(w0, 1e-300),
                rate_Q4=late_slopes(v),
                span=float(traj["span"].iloc[-1]),
                seconds=round(time.time() - t0, 2))


def mode_sem(args):
    out = args.out
    os.makedirs(os.path.join(out, "traj_sem"), exist_ok=True)
    rows = []
    meta = []
    datasets = args.only.split(",") if args.only else WEIGHTED
    for ds in datasets:
        st = raw_third_column_stats(ds)
        if st:
            meta.append(st)
        for conv, _, _ in CONVENTIONS:
            print("[sem] %-42s %-6s steps<=%d" % (ds, conv, args.steps),
                  flush=True)
            rec = run_one_convention(ds, conv, args.steps, args.budget)
            if rec is None:
                print("    MISSING", flush=True)
                continue
            rows.append(rec)
            print("    steps=%3d V0=%.4e ratio=%.4e rateQ4=%+.2e (%.0fs)"
                  % (rec["steps"], rec["V0"], rec["ratio"],
                     rec["rate_Q4"], rec["seconds"]), flush=True)
            pd.DataFrame(rows).to_csv(
                os.path.join(out, "weight_semantics_sweep.csv"), index=False)
    pd.DataFrame(rows).to_csv(
        os.path.join(out, "weight_semantics_sweep.csv"), index=False)
    if meta:
        pd.DataFrame(meta).to_csv(
            os.path.join(out, "raw_weight_column_stats.csv"), index=False)
    print("\n=== summary ===")
    df = pd.DataFrame(rows)
    piv = df.pivot_table(index="dataset", columns="convention",
                         values=["ratio_within", "rate_Q4"])
    print(piv.to_string())
    return df


# --------------------------------------------------------------------------
# R6-N  numerical reliability at extreme weight span
# --------------------------------------------------------------------------
def _linprog_cost(a, b, M):
    from scipy.optimize import linprog
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    M = np.asarray(M, dtype=np.float64)
    a = a / a.sum()
    b = b / b.sum()
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
    res = linprog(M.ravel(), A_eq=np.asarray(rows),
                  b_eq=np.concatenate([a, b]), bounds=(0, None),
                  method="highs")
    return float(res.fun), res


def symmetry_check():
    """Reference graphs whose symmetry forces constant curvature."""
    out = []
    for name, g in (("C_6", nx.cycle_graph(6)), ("K_6", nx.complete_graph(6))):
        for u, v in g.edges():
            g[u][v]["weight"] = 1.0
        ws = R.CurvatureWorkspace(g, ALPHA, "auto", 256.0)
        sol = Solver(ws, "exact")
        k = sol.curvature()
        out.append(dict(reference=name, var_kappa=float(k.var()),
                        max_dev=float(np.max(np.abs(k - k.mean())))))
    return out


def mode_num(args):
    import ot
    out = args.out
    os.makedirs(out, exist_ok=True)
    rows_ckpt = []
    rows_final = []
    for ds in args.only.split(","):
        st = setup(ds, L.build_base_args(["--data_dir", DATA_DIR]))
        if st is None:
            print("MISSING", ds)
            continue
        graph, ws, edges, comp_of_edge = st
        solver = Solver(ws, "exact")
        k = solver.curvature()
        p0, w0, b0 = decompose(k, comp_of_edge)
        t0 = time.time()
        step = 0
        while step < args.steps and (time.time() - t0) < args.budget:
            R.ricci_flow_step(graph, dict(zip(edges, k)), DT)
            k = solver.curvature()
            step += 1
            if step % args.ckpt == 0 or step == args.steps:
                p, w, b = decompose(k, comp_of_edge)
                # ---- checkpoint consistency battery ------------------------
                mu, mv, costs, dists = solver.payload()
                gaps = []
                lp_fail = 0
                for a, b_, M in zip(mu, mv, costs):
                    a = np.asarray(a, dtype=np.float64)
                    b_ = np.asarray(b_, dtype=np.float64)
                    M = np.asarray(M, dtype=np.float64)
                    a = a / a.sum()
                    b_ = b_ / b_.sum()
                    try:
                        G = ot.emd(a, b_, M, numItermax=1000000)
                        w_emd = float((G * M).sum())
                    except Exception:
                        w_emd = float("nan")
                    try:
                        w_lp, _ = _linprog_cost(a, b_, M)
                    except Exception:
                        w_lp = float("nan")
                        lp_fail += 1
                    gaps.append(abs(w_lp - w_emd) / max(abs(w_emd), 1e-300))
                gaps = np.asarray([g for g in gaps if np.isfinite(g)])
                sym = symmetry_check()
                rows_ckpt.append(dict(
                    dataset=ds, step=step, V_within=w,
                    span=weight_span(graph),
                    n_edges=len(edges),
                    n_transport_problems=len(mu),
                    lp_failures=lp_fail,
                    max_rel_cost_gap_emd_vs_highs=(float(gaps.max())
                                                   if gaps.size else float("nan")),
                    median_rel_cost_gap=(float(np.median(gaps))
                                         if gaps.size else float("nan")),
                    sym_C6_var=float(sym[0]["var_kappa"]),
                    sym_K6_var=float(sym[1]["var_kappa"]),
                    sym_pass=bool(sym[0]["var_kappa"] < 1e-20
                                  and sym[1]["var_kappa"] < 1e-20),
                    checkpoint_pass=bool(
                        (not lp_fail)
                        and (gaps.size == 0 or float(gaps.max()) < 1e-9)
                        and sym[0]["var_kappa"] < 1e-20
                        and sym[1]["var_kappa"] < 1e-20)))
            if step % args.ckpt == 0:
                print("  [%s] step %5d V_within=%.4e span=%.2e" %
                      (ds, step, decompose(k, comp_of_edge)[1],
                       weight_span(graph)), flush=True)
            pd.DataFrame(rows_ckpt).to_csv(
                os.path.join(out, "num_checkpoints.csv"), index=False)
        # ---------------- final-weight probes ----------------
        p, w, b = decompose(k, comp_of_edge)
        span = weight_span(graph)
        mu, mv, costs, dists = solver.payload()
        k_emd = k.copy()

        # (a) HiGHS instead of network simplex
        k_lp = []
        for a, b_, M, d in zip(mu, mv, costs, dists):
            w_lp, _ = _linprog_cost(a, b_, M)
            k_lp.append(1.0 - w_lp / float(d))
        k_lp = np.asarray(k_lp, dtype=float)

        # (b) cost matrices normalised to unit mean (removes the 1e72 scale)
        k_scaled = []
        for a, b_, M, d in zip(mu, mv, costs, dists):
            M = np.asarray(M, dtype=np.float64)
            sc = float(np.mean(M))
            Mn = M / sc if sc > 0 else M
            aa = np.asarray(a, dtype=np.float64)
            bb = np.asarray(b_, dtype=np.float64)
            G = ot.emd(aa / aa.sum(), bb / bb.sum(), Mn, numItermax=1000000)
            k_scaled.append(1.0 - (float((G * Mn).sum()) * sc) / float(d))
        k_scaled = np.asarray(k_scaled, dtype=float)

        # (c) Sinkhorn at two tolerances
        ks = []
        for tol, it in ((1e-9, 5000), (1e-13, 50000)):
            vals = batched_sinkhorn_cost(mu, mv, costs, REG, it, tol)
            ks.append(np.array([1.0 - float(v) / float(d)
                                for v, d in zip(vals, dists)], dtype=float))

        rows_final.append(dict(
            dataset=ds, steps=step, span=span,
            V_end=float(p), V_within=float(w),
            var_k_emd=float(k_emd.var()), var_k_lp=float(k_lp.var()),
            max_abs_dk_emd_vs_highs=float(np.max(np.abs(k_emd - k_lp))),
            rel_dvar_emd_vs_highs=float(abs(k_emd.var() - k_lp.var())
                                        / max(k_emd.var(), 1e-300)),
            max_abs_dk_emd_vs_scaledcosts=float(
                np.max(np.abs(k_emd - k_scaled))),
            rel_dvar_emd_vs_scaledcosts=float(
                abs(k_emd.var() - k_scaled.var()) / max(k_emd.var(), 1e-300)),
            max_abs_dk_sink_tol9=float(np.max(np.abs(k_emd - ks[0]))),
            max_abs_dk_sink_tol13=float(np.max(np.abs(k_emd - ks[1]))),
            rel_dvar_sink_tol9=float(abs(k_emd.var() - ks[0].var())
                                     / max(k_emd.var(), 1e-300)),
            rel_dvar_sink_tol13=float(abs(k_emd.var() - ks[1].var())
                                      / max(k_emd.var(), 1e-300)),
            seconds=round(time.time() - t0, 2)))
        print("[num] %s done: steps=%d span=%.2e V_within=%.4e |dk(HiGHS)|=%.2e"
              " |dk(scaled)|=%.2e" %
              (ds, step, span, w,
               rows_final[-1]["max_abs_dk_emd_vs_highs"],
               rows_final[-1]["max_abs_dk_emd_vs_scaledcosts"]), flush=True)
        pd.DataFrame(rows_final).to_csv(
            os.path.join(out, "num_final_probes.csv"), index=False)
    pd.DataFrame(rows_ckpt).to_csv(
        os.path.join(out, "num_checkpoints.csv"), index=False)
    pd.DataFrame(rows_final).to_csv(
        os.path.join(out, "num_final_probes.csv"), index=False)
    ck = pd.DataFrame(rows_ckpt)
    if len(ck):
        print("\ncheckpoints: %d, passed %d (%.1f%%)" %
              (len(ck), int(ck["checkpoint_pass"].sum()),
               100.0 * ck["checkpoint_pass"].mean()))
        print(ck.groupby("dataset")["checkpoint_pass"].agg(["size", "sum"])
              .to_string())
    print(pd.DataFrame(rows_final).to_string(index=False))


# --------------------------------------------------------------------------
# R6-N(b) initial-weight perturbation
# --------------------------------------------------------------------------
def mode_init(args):
    out = args.out
    os.makedirs(out, exist_ok=True)
    rows = []
    for ds in args.only.split(","):
        for rep in range(args.reps + 1):
            rng = np.random.default_rng(1000 + rep)
            st = setup(ds, L.build_base_args(["--data_dir", DATA_DIR]))
            if st is None:
                continue
            graph, ws, edges, comp_of_edge = st
            if rep > 0:
                for _, _, d in graph.edges(data=True):
                    d["weight"] = float(d.get("weight", 1.0)) * float(
                        1.0 + rng.uniform(-args.jitter, args.jitter))
            solver = Solver(ws, "exact")
            k = solver.curvature()
            p0, w0, b0 = decompose(k, comp_of_edge)
            t0 = time.time()
            step = 0
            while step < args.steps and (time.time() - t0) < args.budget:
                R.ricci_flow_step(graph, dict(zip(edges, k)), DT)
                k = solver.curvature()
                step += 1
            p, w, b = decompose(k, comp_of_edge)
            rows.append(dict(dataset=ds, rep=rep, jitter=(0.0 if rep == 0
                                                          else args.jitter),
                             steps=step, V0=p0, V_end=p, V_within=w,
                             ratio=p / max(p0, 1e-300),
                             span=weight_span(graph),
                             seconds=round(time.time() - t0, 2)))
            print("[init] %-16s rep=%d steps=%d ratio=%.6e (%.0fs)" %
                  (ds, rep, step, rows[-1]["ratio"], rows[-1]["seconds"]),
                  flush=True)
            pd.DataFrame(rows).to_csv(
                os.path.join(out, "init_weight_perturbation.csv"), index=False)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out, "init_weight_perturbation.csv"), index=False)
    print("\n=== spread of the residual under 1%% initial-weight jitter ===")
    for ds, g in df.groupby("dataset"):
        base = g[g.rep == 0]["V_within"].iloc[0]
        pert = g[g.rep > 0]["V_within"].values
        print("  %-16s baseline=%.6e  perturbed range=[%.6e, %.6e]  "
              "max rel dev=%.3e" %
              (ds, base, pert.min(), pert.max(),
               float(np.max(np.abs(pert / base - 1.0)))))


# --------------------------------------------------------------------------
# R6-G  long-horizon diagnostics for the three connected ENZYMES graphs
# --------------------------------------------------------------------------
def mode_g10(args):
    """Per-network diagnostics from the stored 5,000-step trajectories."""
    out = args.out
    rows = []
    for ds in args.only.split(","):
        path = os.path.join(args.traj_dir, "%s__exact%d.csv" % (ds, args.steps))
        if not os.path.exists(path):
            print("missing", path)
            continue
        tr = pd.read_csv(path)
        V = tr["V_within"].astype(float).values
        span = tr["span"].astype(float).values
        n = len(tr) - 1
        rec = dict(dataset=ds, steps=n,
                   V0=float(V[0]),
                   V_2000=float(V[2000]) if n >= 2000 else float("nan"),
                   V_end=float(V[-1]),
                   abs_change_2000_end=float(V[-1] - V[2000]) if n >= 2000
                   else float("nan"),
                   rel_change_2000_end_pct=float(100 * (V[-1] / V[2000] - 1))
                   if n >= 2000 else float("nan"),
                   span_2000=float(span[2000]) if n >= 2000 else float("nan"),
                   span_end=float(span[-1]))
        i = int(np.argmin(V))
        rec["V_min"] = float(V[i])
        rec["step_of_min"] = i
        rec["rise_after_min_pct"] = float(100 * (V[-1] / V[i] - 1))
        for window in (500, 1000, 2000, 3000):
            seg = V[-(window + 1):]
            r = -np.diff(np.log(seg))
            rec["median_rate_last_%d" % window] = float(np.median(r))
            rec["rel_change_last_%d_pct" % window] = float(
                100 * (V[-1] / V[-window - 1] - 1))
            rec["frac_steps_rising_last_%d" % window] = float(np.mean(r <= 0))
        rel = np.abs(np.diff(V)) / V[1:]
        ok = rel < 1e-5
        best = cur = 0
        for b in ok:
            cur = cur + 1 if b else 0
            best = max(best, cur)
        rec["longest_streak_rel_move_lt_1e-5"] = int(best)
        rec["frac_steps_rel_move_lt_1e-5"] = float(np.mean(ok))
        rec["meets_descriptive_plateau"] = bool(best >= 50)
        rows.append(rec)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out, "longrun_three_networks.csv"), index=False)
    print(df.to_string(index=False))


def versions():
    import scipy
    import networkx
    try:
        import ot
        pot = ot.__version__
    except Exception:
        pot = "n/a"
    return dict(python=sys.version.split()[0], numpy=np.__version__,
                scipy=scipy.__version__, networkx=networkx.__version__,
                pandas=pd.__version__, POT=pot,
                platform=platform.platform())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="sem")
    ap.add_argument("--out", default="exp_round6")
    ap.add_argument("--only", default="")
    ap.add_argument("--steps", type=int, default=400)
    ap.add_argument("--budget", type=float, default=120.0)
    ap.add_argument("--ckpt", type=int, default=100)
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--jitter", type=float, default=0.01)
    ap.add_argument("--traj_dir", default=os.path.join("exp_round4", "traj_long"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "versions.json"), "w",
              encoding="utf-8") as fh:
        json.dump(versions(), fh, indent=2)
    print(json.dumps(versions(), indent=2))
    if args.mode == "sem":
        mode_sem(args)
    elif args.mode == "num":
        mode_num(args)
    elif args.mode == "init":
        mode_init(args)
    elif args.mode == "g10":
        mode_g10(args)


if __name__ == "__main__":
    main()
