"""Experiment E: pooled vs within-component curvature variance.

42/110 corpus networks and 6/16 strict-rerun networks are disconnected.  The
flow is component-wise, but the reported V is the variance of the *pooled*
edge-curvature distribution, so between-component differences in the limiting
curvature contribute to it.  If every component reaches its own constant
curvature while those constants differ, V_pooled stalls at a strictly positive
value even though uniformisation succeeded component by component.

We therefore decompose, at the same weight configuration,

    V_pooled = V_within + V_between

and report both.  If V_within -> 0 while V_pooled does not, the stalled
residual is a *component-mismatch* artefact of the pooled estimator, not an
intrinsic obstruction.
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
from fast_sinkhorn import batched_sinkhorn_cost

DATA_DIR = r"D:\100data"
ALPHA = 0.5
DT = 0.05
REG = 0.05

TARGETS = [
    ("scc2034_kilifi_all_contacts_across_households", "stalled"),
    ("scc_rt_oman", "stalled"),
    ("scc_rt_ksa", "converging"),
    ("reptilia-tortoise-network-mc", "converging"),
    ("reptilia-tortoise-network-sg", "converging"),
    ("reptilia-tortoise-network-hw", "converging"),
]


def find_file(dataset_id):
    for path in R.discover_graph_files(DATA_DIR):
        if R.safe_dataset_id(path, DATA_DIR) == dataset_id:
            return path
    return None


def curvatures(ws, graph):
    ws.prepare_step()
    mu, mv, costs, dists = [], [], [], []
    for i in range(len(ws.edges)):
        _, _, a, b, c, d = ws.edge_payload(i)
        mu.append(a)
        mv.append(b)
        costs.append(c)
        dists.append(d)
    vals = batched_sinkhorn_cost(mu, mv, costs, REG, 5000, 1e-9)
    return np.array([1.0 - float(v) / float(d)
                     for v, d in zip(vals, dists)], dtype=float)


def decompose(k, comp_of_edge):
    k = np.asarray(k, dtype=float)
    pooled = float(k.var())
    within, sizes = 0.0, 0
    means = []
    for c in sorted(set(comp_of_edge)):
        kc = k[np.asarray(comp_of_edge) == c]
        if len(kc) == 0:
            continue
        within += float(kc.var()) * len(kc)
        sizes += len(kc)
        means.append(float(kc.mean()))
    within = within / max(sizes, 1)
    between = float(np.var(np.repeat(means, 1))) if means else 0.0
    # exact between-component term, weighted by component edge counts
    counts = np.array([int((np.asarray(comp_of_edge) == c).sum())
                       for c in sorted(set(comp_of_edge))])
    means = np.array(means)
    if counts.sum():
        gm = float(np.average(means, weights=counts))
        between = float(np.average((means - gm) ** 2, weights=counts))
    return pooled, within, between


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="exp_comp")
    ap.add_argument("--steps", type=int, default=250)
    ap.add_argument("--budget", type=float, default=100.0)
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    base = L.build_base_args(["--data_dir", DATA_DIR])
    wanted = [t for t in args.only.split(",") if t]
    rows = []
    for dataset_id, fate in TARGETS:
        if wanted and dataset_id not in wanted:
            continue
        path = find_file(dataset_id)
        if path is None:
            print("MISSING", dataset_id)
            continue
        graph = R.load_graph(path, base)
        ws = R.CurvatureWorkspace(graph, ALPHA, "auto", 256.0)
        edges = list(ws.edges)
        comp = {}
        for i, c in enumerate(nx.connected_components(graph)):
            for n in c:
                comp[n] = i
        comp_of_edge = [comp[u] for u, v in edges]
        n_comp = len(set(comp_of_edge))

        k = curvatures(ws, graph)
        p0, w0, b0 = decompose(k, comp_of_edge)
        t0 = time.time()
        step = 0
        try:
            while step < args.steps and time.time() - t0 < args.budget:
                R.ricci_flow_step(graph, dict(zip(edges, k)), DT)
                R.normalize_weights(graph)
                k = curvatures(ws, graph)
                step += 1
        except Exception as exc:
            print("   [stop] %s: %r" % (dataset_id, exc))
        p1, w1, b1 = decompose(k, comp_of_edge)
        rows.append(dict(dataset=dataset_id, fate=fate, steps=step,
                         n_components=n_comp, n_edges=len(edges),
                         V0_pooled=p0, V0_within=w0, V0_between=b0,
                         V_pooled=p1, V_within=w1, V_between=b1,
                         ratio_pooled=p1 / p0, ratio_within=w1 / w0,
                         between_share_end=b1 / max(p1, 1e-300)))
        print("  %-52s comps=%d steps=%d | pooled %.4e->%.4e (x%.3f) | "
              "within %.4e->%.4e (x%.3f) | between share now %.1f%%"
              % (dataset_id[:52], n_comp, step, p0, p1, p1 / p0,
                 w0, w1, w1 / w0,
                 100.0 * b1 / max(p1, 1e-300)))
    df = pd.DataFrame(rows)
    out = os.path.join(args.out, "components.csv")
    df.to_csv(out, index=False)
    print("\nwrote", out)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
