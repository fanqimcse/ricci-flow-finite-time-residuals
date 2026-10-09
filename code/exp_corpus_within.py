"""Reviewer round 2, points 18-19: recompute the corpus-level variance on the
within-component definition.

The paper shows that on a disconnected graph the *pooled* edge-curvature
variance is the wrong convergence estimator, but then continues to use it for
every corpus-level statistic.  That is inconsistent.  Here we re-run every
disconnected network of the analysis set to the horizon it reached in the
first pass and record both estimators:

    V_pooled  = V_within + V_between
    V_within  = sum_c (m_c/m) Var_c
    V_between = sum_c (m_c/m) (kbar_c - kbar)^2

For connected networks V_within == V_pooled identically, so only the
disconnected members of the corpus need re-running.  Results are written
incrementally per dataset so an interrupted batch can be resumed.
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
from exp_round3 import Solver, decompose, weight_span

DATA_DIR = r"D:\100data"
ALPHA = 0.5
DT = 0.05
REG = 0.05


def find_file(dataset_id):
    for path in R.discover_graph_files(DATA_DIR):
        if R.safe_dataset_id(path, DATA_DIR) == dataset_id:
            return path
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="exp_corpus_within")
    ap.add_argument("--cap", type=int, default=200)
    ap.add_argument("--budget", type=float, default=900.0)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    out_csv = os.path.join(args.out, "corpus_within.csv")

    struct = pd.read_csv(
        os.path.join("ricci_flow_longrun_results", "analysis",
                     "07_merged_with_structure.csv"))
    sizes = pd.read_csv("_perdataset_sizes.csv").rename(
        columns={"steps": "nsteps"})
    plan = struct.merge(sizes, left_on="dataset", right_on="ds", how="inner")
    plan = plan[plan["n_components"] > 1].copy()
    plan["target"] = np.minimum(plan["nsteps"], args.cap)
    plan = plan.sort_values("med_step_s")  # cheap first, most progress early

    done = set()
    if os.path.exists(out_csv):
        done = set(pd.read_csv(out_csv)["dataset"].tolist())
    todo = plan[~plan["dataset"].isin(done)]
    if args.limit:
        todo = todo.head(args.limit)
    print("disconnected networks in analysis set: %d | already done: %d | "
          "to run: %d" % (len(plan), len(done), len(todo)))

    base = L.build_base_args(["--data_dir", DATA_DIR])
    for _, row in todo.iterrows():
        ds = row["dataset"]
        path = find_file(ds)
        if path is None:
            print("  MISSING", ds)
            continue
        try:
            graph = R.load_graph(path, base)
        except Exception as exc:
            print("  LOAD-FAIL %-46s %r" % (ds[:46], exc))
            pd.DataFrame([dict(dataset=ds, status="load_fail")]).to_csv(
                out_csv, mode="a", header=not os.path.exists(out_csv),
                index=False)
            continue
        ws = R.CurvatureWorkspace(graph, ALPHA, "auto", 256.0)
        edges = list(ws.edges)
        comp = {}
        for i, c in enumerate(nx.connected_components(graph)):
            for n in c:
                comp[n] = i
        comp_of_edge = np.array([comp[u] for u, v in edges], dtype=int)
        solver = Solver(ws, "sinkhorn")
        t0 = time.time()
        try:
            k = solver.curvature()
            p0, w0, b0 = decompose(k, comp_of_edge)
            step = 0
            stop = "horizon"
            while step < int(row["target"]) and time.time() - t0 < args.budget:
                R.ricci_flow_step(graph, dict(zip(edges, k)), DT)
                k = solver.curvature()
                step += 1
            if time.time() - t0 >= args.budget:
                stop = "budget"
        except Exception as exc:
            k = None
            stop = "error"
            print("  ERROR %-48s %r" % (ds[:46], exc)[:160])
        if k is None:
            pd.DataFrame([dict(dataset=ds, status="error",
                               steps=0)]).to_csv(
                out_csv, mode="a", header=not os.path.exists(out_csv),
                index=False)
            continue
        p1, w1, b1 = decompose(k, comp_of_edge)
        rec = dict(dataset=ds, n_components=int(row["n_components"]),
                   nodes=int(row["nodes"]), edges=int(row["edges"]),
                   steps=step, target=int(row["target"]), stop=stop,
                   V0=p0, V0_within=w0, V_pooled=p1, V_within=w1,
                   V_between=b1,
                   ratio_pooled=p1 / max(p0, 1e-300),
                   ratio_within=w1 / max(w0, 1e-300),
                   between_share=b1 / max(p1, 1e-300),
                   span=weight_span(graph), seconds=time.time() - t0,
                   status="ok")
        pd.DataFrame([rec]).to_csv(
            out_csv, mode="a", header=not os.path.exists(out_csv),
            index=False)
        print("  %-46s comps=%3d steps=%3d | pooled %.3e->%.3e (x%.3f) | "
              "within %.3e->%.3e (x%.3f) | between %.1f%%"
              % (ds[:46], int(row["n_components"]), step, p0, p1,
                 p1 / max(p0, 1e-300), w0, w1, w1 / max(w0, 1e-300),
                 100.0 * b1 / max(p1, 1e-300)))
        sys.stdout.flush()
    print("wrote", out_csv)


if __name__ == "__main__":
    main()
