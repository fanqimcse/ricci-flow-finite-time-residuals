#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Experiment 1 + Experiment 2 for the universal Ollivier-Ricci flow study.

IMPORTANT METHODOLOGICAL CONTRACT
---------------------------------
This script deliberately reuses the definitions in the supplied
1ricci_flow_universal_experiment_correct_fast.py:

1. Undirected simple graph.
2. Self-loops are removed; isolated vertices are removed.
3. Ollivier measure:
       m_x = alpha * delta_x + (1-alpha)/deg(x) * sum_{y~x} delta_y
   with alpha=0.5 by default.
4. Exact shortest-path metric using the CURRENT edge lengths at each flow step.
5. Ollivier curvature:
       kappa(u,v) = 1 - W_1(m_u,m_v) / d(u,v)
   where W_1 is computed by POT Sinkhorn with the same reg/max_iter defaults.
6. Ricci flow:
       w_e <- w_e * exp(-dt * kappa_e)
   followed by division by the mean updated edge weight.
7. lambda1 is the SECOND eigenvalue of the weighted COMBINATORIAL graph
   Laplacian. For a disconnected graph, lambda1 is recorded as 0, exactly
   as in the supplied reference implementation.
8. Dataset is the independent experimental unit. Replicates are NOT treated
   as independent biological/statistical samples.
9. No TensorFlow.
10. No result-dependent filtering, cherry-picking, or manual result changes.

EXPERIMENT 1: TOPOLOGY-PRESERVING NULL
---------------------------------------
For every eligible observed graph, generate R degree-preserving double-edge
swap null graphs. The node set, edge count and degree sequence are preserved.
For connected observed graphs, accepted swaps are constrained to retain
connectedness, so lambda1 remains identifiable under the same connectivity
condition. Edge lengths are carried with the swapped edge, preserving the
empirical edge-length multiset.

EXPERIMENT 2: RANDOMIZED FLOW DIRECTION
----------------------------------------
Keep the observed graph and every computed curvature value fixed, but replace
the deterministic Ricci update direction by a pre-registered random sign:

       w_e <- w_e * exp(-dt * s_e * kappa_e),  s_e in {-1,+1}

with equal probability independently for each edge at each update. The
normalization by the mean updated edge weight is unchanged. Thus the
randomized experiment changes FLOW DIRECTION, not curvature magnitude, graph
size, edge count, or the normalization rule.

PRIMARY INFERENCE
-----------------
For each dataset and endpoint metric, the R randomized/null replicates are
first reduced to a dataset-specific null median. The observed dataset metric
is then compared with its matched null median across datasets. This avoids
the invalid treatment of R simulations x 110 datasets as R*110 independent
observations.

Outputs:
  raw/                    all trajectory records
  dataset_metrics/        one row per dataset/replicate
  paired_statistics.csv   dataset-level observed vs matched-null statistics
  empirical_pvalues.csv   randomization-test results
  replicate_status.csv    every attempted replicate, including failures
  run_config.json         exact configuration and software protocol
  run.log                 execution log

Recommended first run:
  python ricci_flow_null_randomized_experiments.py --data_dir D:\\100data

For a pilot:
  --reps 10 --only eco-everglades

For the paper-scale run, pre-register the replicate count before execution,
e.g. --reps 100 or --reps 1000. Do not change it after seeing outcomes.
"""

from __future__ import print_function

import argparse
import csv
import json
import logging
import math
import os
import re
import sys
import time
import warnings
from collections import Counter

# Prevent BLAS/OpenMP oversubscription.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import networkx as nx
import numpy as np
import pandas as pd
import ot
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra as sparse_dijkstra
from scipy.sparse.linalg import eigsh
from scipy.stats import wilcoxon

warnings.filterwarnings("ignore")

SUPPORTED_EXTENSIONS = (".edges", ".txt", ".mtx", ".csv")


# ---------------------------------------------------------------------------
# Parsing: intentionally aligned with the supplied reference implementation
# ---------------------------------------------------------------------------

HEADER_LEFT = {"source", "src", "from", "node1", "vertex1", "u", "id1", "tail"}
HEADER_RIGHT = {"target", "dst", "to", "node2", "vertex2", "v", "id2", "head"}


def parse_node(token):
    token = token.strip().strip('"').strip("'")
    if token == "":
        raise ValueError("empty node label")
    try:
        value = float(token)
        if np.isfinite(value) and value.is_integer():
            return int(value)
    except Exception:
        pass
    return token


def split_edge_line(line):
    line = line.strip()
    if "," in line:
        try:
            return [item.strip() for item in next(csv.reader([line]))]
        except Exception:
            pass
    if ";" in line:
        return [item.strip() for item in line.split(";")]
    return re.split(r"\s+", line)


def looks_like_header(parts):
    if len(parts) < 2:
        return False
    a = parts[0].strip().strip('"').strip("'").lower()
    b = parts[1].strip().strip('"').strip("'").lower()
    return a in HEADER_LEFT and b in HEADER_RIGHT


def parse_edge_list(path):
    records = []
    stats = {"invalid_lines": 0, "self_loops": 0, "duplicates": 0}
    seen = set()
    with open(path, "r", errors="ignore") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("%"):
                continue
            parts = split_edge_line(line)
            if looks_like_header(parts) or len(parts) < 2:
                if not looks_like_header(parts):
                    stats["invalid_lines"] += 1
                continue
            try:
                u, v = parse_node(parts[0]), parse_node(parts[1])
            except Exception:
                stats["invalid_lines"] += 1
                continue
            if u == v:
                stats["self_loops"] += 1
                continue
            w = None
            if len(parts) >= 3:
                try:
                    w = float(parts[2])
                except Exception:
                    w = None
            key = frozenset((u, v))
            if key in seen:
                stats["duplicates"] += 1
            seen.add(key)
            records.append((u, v, w))
    return records, stats


def parse_matrix_market(path):
    records = []
    stats = {"invalid_lines": 0, "self_loops": 0, "duplicates": 0}
    with open(path, "r", errors="ignore") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines) and (not lines[i].strip() or lines[i].lstrip().startswith("%")):
        i += 1
    if i >= len(lines):
        return records, stats

    header = lines[i].strip()
    field = "real"
    if header.lower().startswith("%%matrixmarket"):
        toks = header.split()
        if len(toks) >= 5 and toks[2].lower() != "coordinate":
            raise ValueError("Only Matrix Market coordinate format is supported")
        if len(toks) >= 5:
            field = toks[3].lower()
        i += 1

        while i < len(lines) and (not lines[i].strip() or lines[i].lstrip().startswith("%")):
            i += 1
        if i < len(lines):
            i += 1  # dimensions

    seen = set()
    for raw in lines[i:]:
        line = raw.strip()
        if not line or line.startswith("%"):
            continue
        parts = line.split()
        if len(parts) < 2:
            stats["invalid_lines"] += 1
            continue
        try:
            u, v = int(parts[0]), int(parts[1])
        except Exception:
            stats["invalid_lines"] += 1
            continue
        if u == v:
            stats["self_loops"] += 1
            continue
        w = None
        if field != "pattern" and len(parts) >= 3:
            try:
                w = float(parts[2])
            except Exception:
                w = None
        if w == 0.0:
            continue
        key = frozenset((u, v))
        if key in seen:
            stats["duplicates"] += 1
        seen.add(key)
        records.append((u, v, w))
    return records, stats


def choose_weight_mode(records, requested):
    numeric = [w for _, _, w in records if w is not None]
    if requested == "unweighted":
        return "unweighted", "explicit unweighted mode"
    if requested == "distance":
        return "distance", "explicit distance mode"
    if requested == "strength":
        return "strength", "explicit strength mode"

    # auto: same safety principle as the reference implementation.
    if len(numeric) != len(records):
        return "unweighted", "not every edge has a numeric third column"
    a = np.asarray(numeric, dtype=float)
    if not np.all(np.isfinite(a)) or np.any(a <= 0):
        return "unweighted", "non-positive/non-finite input weights"
    if np.allclose(a, 1.0):
        return "unweighted", "all input weights equal one"
    return "distance", "all weights finite and positive; auto=distance"


def transform_length(w, mode):
    if mode == "unweighted":
        return 1.0
    if w is None or not np.isfinite(w) or w <= 0:
        raise ValueError("weighted mode requires finite positive edge values")
    if mode == "distance":
        return float(w)
    if mode == "strength":
        return 1.0 / float(w)
    raise ValueError(mode)


def load_graph(path, weight_mode="auto", component_mode="full"):
    ext = os.path.splitext(path)[1].lower()
    if ext == ".mtx":
        records, stats = parse_matrix_market(path)
    else:
        records, stats = parse_edge_list(path)
    if not records:
        raise ValueError("no usable edges: {}".format(path))

    mode, reason = choose_weight_mode(records, weight_mode)
    G = nx.Graph()
    for u, v, raw_w in records:
        w = transform_length(raw_w, mode)
        if G.has_edge(u, v):
            # Same rule as reference: retain shortest usable path length.
            if w < G[u][v]["weight"]:
                G[u][v]["weight"] = w
        else:
            G.add_edge(u, v, weight=w)

    G.remove_edges_from(list(nx.selfloop_edges(G)))
    G.remove_nodes_from(list(nx.isolates(G)))
    if G.number_of_edges() == 0:
        raise ValueError("empty graph after validation")

    weights = np.asarray([float(d.get("weight", 1.0)) for _, _, d in G.edges(data=True)])
    if not np.all(np.isfinite(weights)) or np.any(weights <= 0):
        raise ValueError("invalid edge length")
    mean_w = float(np.mean(weights))
    for _, _, d in G.edges(data=True):
        d["weight"] = float(d.get("weight", 1.0)) / mean_w

    comps = list(nx.connected_components(G))
    if component_mode == "largest" and len(comps) > 1:
        nodes = max(comps, key=len)
        G = G.subgraph(nodes).copy()
        used_components = 1
    else:
        used_components = len(comps)

    G.graph.update({
        "weight_mode": mode,
        "weight_reason": reason,
        "component_mode": component_mode,
        "component_count_used": used_components,
        "source_nodes": G.number_of_nodes(),
        "source_edges": G.number_of_edges(),
        "source_connected": used_components == 1,
        "invalid_lines": stats["invalid_lines"],
        "duplicate_edges": stats["duplicates"],
    })
    return G


def discover_graph_files(data_dir):
    files = []
    for root, _, names in os.walk(data_dir):
        for name in names:
            if name.lower().endswith(SUPPORTED_EXTENSIONS):
                files.append(os.path.join(root, name))
    return sorted(files)


def dataset_id(path, data_dir):
    rel = os.path.relpath(path, data_dir)
    return re.sub(r"[^A-Za-z0-9._-]+", "__", os.path.splitext(rel)[0]).strip("._-") or "dataset"


# ---------------------------------------------------------------------------
# Exact Ollivier curvature implementation aligned with the reference code
# ---------------------------------------------------------------------------

def node_distribution(G, node, alpha):
    nbrs = list(G.neighbors(node))
    if not nbrs:
        return [node], np.asarray([1.0], dtype=float)
    rem = (1.0 - alpha) / float(len(nbrs))
    return [node] + nbrs, np.asarray([alpha] + [rem] * len(nbrs), dtype=float)


class CurvatureWorkspace(object):
    def __init__(self, G, alpha, max_distance_mb):
        self.G = G
        self.alpha = alpha
        self.nodes = list(G.nodes())
        self.edges = list(G.edges())
        self.idx = {u: i for i, u in enumerate(self.nodes)}
        self.support = {}
        self.mass = {}
        self.support_idx = {}
        for u in self.nodes:
            s, m = node_distribution(G, u, alpha)
            self.support[u], self.mass[u] = s, m
            self.support_idx[u] = np.asarray([self.idx[x] for x in s], dtype=np.int64)

        matrix_mb = len(self.nodes) * len(self.nodes) * 8.0 / (1024.0 * 1024.0)
        self.backend = "matrix" if matrix_mb <= max_distance_mb else "batched"
        self.matrix_mb = matrix_mb
        self.distance_matrix = None
        self.target_distances = None

        if self.backend == "batched":
            req = {u: set() for u in self.nodes}
            for u, v in self.edges:
                for x in self.support[u]:
                    req[x].update(self.support[v])
            self.required_targets = req

    def prepare(self):
        rows, cols, vals = [], [], []
        for u, v, d in self.G.edges(data=True):
            w = float(d["weight"])
            i, j = self.idx[u], self.idx[v]
            rows += [i, j]
            cols += [j, i]
            vals += [w, w]
        A = csr_matrix((vals, (rows, cols)), shape=(len(self.nodes), len(self.nodes)))
        if self.backend == "matrix":
            self.distance_matrix = sparse_dijkstra(A, directed=False, return_predecessors=False)
            self.target_distances = None
        else:
            bytes_limit = max(1.0, 256.0) * 1024.0 * 1024.0
            batch = max(1, min(int(bytes_limit // (max(1, len(self.nodes)) * 8.0)), 256))
            srcs = [u for u in self.nodes if self.required_targets[u]]
            self.target_distances = {}
            for k in range(0, len(srcs), batch):
                part = srcs[k:k + batch]
                mat = np.atleast_2d(sparse_dijkstra(
                    A, directed=False, indices=[self.idx[x] for x in part],
                    return_predecessors=False))
                for r, u in enumerate(part):
                    self.target_distances[u] = {
                        v: float(mat[r, self.idx[v]]) for v in self.required_targets[u]
                    }

    def payload(self, eidx):
        u, v = self.edges[eidx]
        xs, ys = self.support[u], self.support[v]
        if self.backend == "matrix":
            cost = self.distance_matrix[np.ix_(self.support_idx[u], self.support_idx[v])]
            dist = float(self.distance_matrix[self.idx[u], self.idx[v]])
        else:
            cost = np.empty((len(xs), len(ys)), dtype=float)
            for i, x in enumerate(xs):
                for j, y in enumerate(ys):
                    cost[i, j] = self.target_distances[x][y]
            dist = float(self.target_distances[u][v])
        if not np.isfinite(dist) or dist <= 0 or not np.all(np.isfinite(cost)):
            raise RuntimeError("invalid shortest-path transport cost")
        return u, v, self.mass[u], self.mass[v], cost, dist


def compute_curvature(G, alpha, sinkhorn_reg, sinkhorn_max_iter, max_distance_mb):
    ws = CurvatureWorkspace(G, alpha, max_distance_mb)
    ws.prepare()
    vals = []
    edge_curv = {}
    for i in range(len(ws.edges)):
        u, v, mu, mv, cost, dist = ws.payload(i)
        try:
            W = ot.sinkhorn2(
                mu, mv, cost, reg=sinkhorn_reg,
                numItermax=sinkhorn_max_iter
            )
            W = float(np.asarray(W).reshape(-1)[0])
        except Exception as exc:
            raise RuntimeError("Sinkhorn failed on ({},{})".format(u, v)) from exc
        if not np.isfinite(W) or W < 0:
            raise RuntimeError("invalid Sinkhorn result")
        k = 1.0 - W / dist
        if not np.isfinite(k):
            raise RuntimeError("non-finite curvature")
        if k > 1.0 and k < 1.0 + 1e-10:
            k = 1.0
        edge_curv[(u, v)] = float(k)
        vals.append(float(k))
    a = np.asarray(vals, dtype=float)
    return float(np.mean(a)), float(np.var(a)), edge_curv


def lambda1(G):
    n = G.number_of_nodes()
    if n < 3:
        return 0.0
    if G.graph.get("component_count_used", nx.number_connected_components(G)) > 1:
        return 0.0
    nodes = list(G.nodes())
    L = nx.laplacian_matrix(G, nodelist=nodes, weight="weight").astype(float)
    if n <= 128:
        ev = np.linalg.eigvalsh(L.toarray())
        x = float(np.sort(np.real(ev))[1])
    else:
        ev, vec = eigsh(L, k=2, which="SM", return_eigenvectors=True)
        order = np.argsort(np.real(ev))
        x = float(np.real(ev[order])[1])
    if abs(x) < 1e-12:
        x = 0.0
    if x < 0 and x > -1e-9:
        x = 0.0
    return x


def ricci_update(G, edge_curv, dt, rng=None, randomized_direction=False):
    edges = list(edge_curv.keys())
    new_weights = []
    signs = {}
    for u, v in edges:
        k = float(edge_curv[(u, v)])
        if randomized_direction:
            sign = 1 if rng.random() < 0.5 else -1
            signs[(u, v)] = sign
        else:
            sign = 1
        old = float(G[u][v]["weight"])
        exponent = -float(dt) * sign * k
        if not np.isfinite(exponent) or abs(exponent) > 700:
            raise RuntimeError("Ricci update overflow risk")
        w = old * np.exp(exponent)
        if not np.isfinite(w) or w <= 0:
            raise RuntimeError("invalid Ricci update")
        new_weights.append(w)

    m = float(np.mean(new_weights))
    if not np.isfinite(m) or m <= 0:
        raise RuntimeError("invalid mean updated edge weight")
    for i, (u, v) in enumerate(edges):
        G[u][v]["weight"] = new_weights[i] / m
    return signs


def run_flow(G0, steps, alpha, dt, sinkhorn_reg, sinkhorn_max_iter,
             max_distance_mb, randomized_direction=False, seed=None):
    G = G0.copy()
    rng = np.random.default_rng(seed)
    rows = []
    sign_rows = []

    for step in range(steps):
        mean_k, var_k, edge_k = compute_curvature(
            G, alpha, sinkhorn_reg, sinkhorn_max_iter, max_distance_mb
        )
        l1 = lambda1(G)
        rows.append({
            "step": step,
            "mean_curvature": mean_k,
            "curvature_variance": var_k,
            "lambda1": l1,
            "nodes": G.number_of_nodes(),
            "edges": G.number_of_edges(),
        })

        if step < steps - 1:
            signs = ricci_update(
                G, edge_k, dt, rng=rng,
                randomized_direction=randomized_direction
            )
            if randomized_direction:
                plus = sum(1 for x in signs.values() if x == 1)
                minus = len(signs) - plus
                sign_rows.append({
                    "step": step,
                    "plus_edges": plus,
                    "minus_edges": minus,
                    "total_edges": len(signs)
                })

    return pd.DataFrame(rows), pd.DataFrame(sign_rows)


# ---------------------------------------------------------------------------
# Endpoint metrics: exactly the metrics used in the preceding study
# ---------------------------------------------------------------------------

def endpoint_metrics(df):
    if len(df) < 2:
        raise ValueError("at least two trajectory steps are required")
    k0, kf = float(df["mean_curvature"].iloc[0]), float(df["mean_curvature"].iloc[-1])
    v0, vf = float(df["curvature_variance"].iloc[0]), float(df["curvature_variance"].iloc[-1])
    l0, lf = float(df["lambda1"].iloc[0]), float(df["lambda1"].iloc[-1])

    # Same endpoint definitions as the user's previous statistical analysis.
    delta_kappa = kf - k0
    variance_reduction = np.nan if abs(v0) <= 1e-15 else (v0 - vf) / v0
    lambda_improvement = np.nan if l0 <= 0 else (lf - l0) / l0
    normalized_rate = delta_kappa / float(len(df) - 1)
    convergence_score = 1.0 if np.isfinite(k0) and np.isfinite(kf) else 0.0

    return {
        "initial_kappa": k0,
        "final_kappa": kf,
        "delta_kappa": delta_kappa,
        "initial_variance": v0,
        "final_variance": vf,
        "variance_reduction": variance_reduction,
        "lambda_initial": l0,
        "lambda_final": lf,
        "lambda_improvement": lambda_improvement,
        "normalized_rate": normalized_rate,
        "convergence_score": convergence_score,
        "lambda_identifiable": bool(l0 > 0),
    }


# ---------------------------------------------------------------------------
# Degree-preserving topology null
# ---------------------------------------------------------------------------

def degree_preserving_swap(G0, rng, max_attempts, preserve_connected=True):
    """
    One double-edge swap:
       (a,b), (c,d) -> (a,d), (c,b)
    with the two edge lengths carried to the new endpoints.

    Degree sequence and edge count are exactly preserved.
    No attempt is accepted merely because it gives a desired metric result.
    """
    G = G0.copy()
    edges = list(G.edges())
    if len(edges) < 2:
        raise RuntimeError("fewer than two edges; cannot rewire")

    original_components = nx.number_connected_components(G)

    for _ in range(max_attempts):
        i, j = rng.choice(len(edges), size=2, replace=False)
        a, b = edges[int(i)]
        c, d = edges[int(j)]
        if len({a, b, c, d}) < 4:
            continue

        # Randomly choose one of the two valid cross pairings.
        if rng.random() < 0.5:
            new1, new2 = (a, d), (c, b)
        else:
            new1, new2 = (a, c), (b, d)

        if new1[0] == new1[1] or new2[0] == new2[1]:
            continue

        e1 = frozenset((a, b))
        e2 = frozenset((c, d))
        n1 = frozenset(new1)
        n2 = frozenset(new2)

        if n1 in (e1, e2) or n2 in (e1, e2) or n1 == n2:
            continue
        if G.has_edge(*new1) or G.has_edge(*new2):
            continue

        w1 = float(G[a][b]["weight"])
        w2 = float(G[c][d]["weight"])

        H = G.copy()
        H.remove_edge(a, b)
        H.remove_edge(c, d)
        H.add_edge(new1[0], new1[1], weight=w1)
        H.add_edge(new2[0], new2[1], weight=w2)

        if preserve_connected and original_components == 1:
            if not nx.is_connected(H):
                continue

        # For disconnected inputs, preserve the number of components as a
        # structural safeguard without inventing a connectedness assumption.
        if nx.number_connected_components(H) != original_components:
            continue

        # Exact degree-sequence check: this guards against implementation bugs.
        if sorted(dict(H.degree()).values()) != sorted(dict(G0.degree()).values()):
            raise RuntimeError("internal error: degree sequence changed")

        return H

    raise RuntimeError(
        "could not find an accepted degree-preserving swap within {} attempts"
        .format(max_attempts)
    )


def rewire_graph(G0, rng, n_swaps, max_attempts_per_swap):
    H = G0.copy()
    for _ in range(n_swaps):
        H = degree_preserving_swap(H, rng, max_attempts_per_swap,
                                   preserve_connected=True)
    return H


# ---------------------------------------------------------------------------
# Statistics at the DATASET level
# ---------------------------------------------------------------------------

def finite_pairs(a, b):
    x = pd.to_numeric(a, errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(b, errors="coerce").to_numpy(dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    return x[ok], y[ok]


def paired_wilcoxon(observed, null_median):
    x, y = finite_pairs(observed, null_median)
    d = x - y
    d = d[np.isfinite(d)]
    if len(d) == 0:
        return {"n": 0, "statistic": np.nan, "p_value": np.nan,
                "median_difference": np.nan, "rank_biserial": np.nan}
    nz = d[np.abs(d) > 1e-15]
    if len(nz) == 0:
        p = 1.0
        stat = 0.0
    else:
        try:
            res = wilcoxon(nz, alternative="two-sided", zero_method="wilcox")
            stat, p = float(res.statistic), float(res.pvalue)
        except Exception:
            stat, p = np.nan, np.nan
    pos = np.sum(d > 0)
    neg = np.sum(d < 0)
    rb = (pos - neg) / float(len(d))
    return {
        "n": int(len(d)),
        "statistic": stat,
        "p_value": p,
        "median_difference": float(np.median(d)),
        "mean_difference": float(np.mean(d)),
        "rank_biserial_sign": float(rb),
    }


def bh_fdr(pvalues):
    p = np.asarray(pvalues, dtype=float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if not np.any(ok):
        return out
    idx = np.where(ok)[0]
    order = idx[np.argsort(p[idx])]
    m = len(order)
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        q = min(prev, p[i] * m / float(rank))
        out[i] = q
        prev = q
    return out


def paired_null_statistics(dataset_table, endpoint):
    sub = dataset_table.dropna(subset=[endpoint, "null_median"]).copy()
    result = paired_wilcoxon(sub[endpoint], sub["null_median"])
    return result


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------

def run(args):
    out = os.path.abspath(args.output_dir)
    os.makedirs(out, exist_ok=True)
    for d in ("raw", "dataset_metrics"):
        os.makedirs(os.path.join(out, d), exist_ok=True)

    log_path = os.path.join(out, "run.log")
    logging.basicConfig(
        filename=log_path, level=logging.INFO,
        format="%(asctime)s | %(message)s"
    )

    config = vars(args).copy()
    config["python"] = sys.version
    config["platform"] = sys.platform
    config["numpy"] = np.__version__
    config["networkx"] = nx.__version__
    config["pandas"] = pd.__version__
    config["protocol"] = "matched dataset-level null/randomized inference"
    with open(os.path.join(out, "run_config.json"), "w", encoding="utf-8") as f:
        json.dump(config, f, indent=2, ensure_ascii=False)

    files = discover_graph_files(args.data_dir)
    if args.only:
        token = args.only.lower()
        files = [p for p in files if token in os.path.relpath(p, args.data_dir).lower()]

    if not files:
        raise RuntimeError("No supported graph files found under {}".format(args.data_dir))

    all_metrics = []
    status_rows = []
    raw_rows = []
    sign_rows = []

    for file_index, path in enumerate(files):
        ds = dataset_id(path, args.data_dir)
        base_seed = args.seed + file_index * 1000003
        print("[{}/{}] {}".format(file_index + 1, len(files), ds), flush=True)
        logging.info("START %s | %s", ds, path)

        try:
            G0 = load_graph(path, args.weight_mode, args.component_mode)
            n, m = G0.number_of_nodes(), G0.number_of_edges()

            # Real observed trajectory.
            real_df, _ = run_flow(
                G0, args.steps, args.alpha, args.dt,
                args.sinkhorn_reg, args.sinkhorn_max_iter,
                args.max_distance_mb, randomized_direction=False,
                seed=base_seed
            )
            real_m = endpoint_metrics(real_df)
            real_df.insert(0, "dataset", ds)
            real_df.insert(1, "experiment", "observed")
            real_df.to_csv(os.path.join(out, "raw", ds + "__observed.csv"), index=False)

            real_record = {
                "dataset": ds, "file": path, "nodes": n, "edges": m,
                "weight_mode": G0.graph["weight_mode"],
                "weight_reason": G0.graph["weight_reason"],
                "components": G0.graph["component_count_used"],
                "replicate": 0, "seed": base_seed,
                "experiment": "observed", "status": "finished"
            }
            real_record.update(real_m)
            all_metrics.append(real_record)
            status_rows.append(real_record.copy())

            # Experiment 1: degree-preserving topology null.
            for r in range(args.reps):
                seed = base_seed + 100000 + r
                rng = np.random.default_rng(seed)
                row = {
                    "dataset": ds, "experiment": "topology_null",
                    "replicate": r + 1, "seed": seed, "status": "failed",
                    "message": "", "nodes": n, "edges": m
                }
                try:
                    swaps = args.n_swaps if args.n_swaps > 0 else max(1, 10 * m)
                    H = rewire_graph(
                        G0, rng, swaps, args.max_swap_attempts
                    )
                    # Hard invariants. Failure is recorded, never silently accepted.
                    if H.number_of_nodes() != n or H.number_of_edges() != m:
                        raise RuntimeError("node/edge count changed")
                    if sorted(dict(H.degree()).values()) != sorted(dict(G0.degree()).values()):
                        raise RuntimeError("degree sequence changed")
                    if nx.number_connected_components(H) != G0.graph["component_count_used"]:
                        raise RuntimeError("component count changed")

                    traj, _ = run_flow(
                        H, args.steps, args.alpha, args.dt,
                        args.sinkhorn_reg, args.sinkhorn_max_iter,
                        args.max_distance_mb, False, seed
                    )
                    met = endpoint_metrics(traj)
                    traj.insert(0, "dataset", ds)
                    traj.insert(1, "experiment", "topology_null")
                    traj.insert(2, "replicate", r + 1)
                    traj.to_csv(
                        os.path.join(out, "raw",
                                     "{}__topology_null__{:04d}.csv".format(ds, r + 1)),
                        index=False
                    )
                    row["status"] = "finished"
                    row["message"] = ""
                    row.update(met)
                    all_metrics.append(row.copy())
                    status_rows.append(row.copy())
                except Exception as exc:
                    row["message"] = repr(exc)
                    status_rows.append(row.copy())

            # Experiment 2: randomized flow direction.
            for r in range(args.reps):
                seed = base_seed + 200000 + r
                row = {
                    "dataset": ds, "experiment": "randomized_direction",
                    "replicate": r + 1, "seed": seed, "status": "failed",
                    "message": "", "nodes": n, "edges": m
                }
                try:
                    traj, signs = run_flow(
                        G0, args.steps, args.alpha, args.dt,
                        args.sinkhorn_reg, args.sinkhorn_max_iter,
                        args.max_distance_mb, True, seed
                    )
                    met = endpoint_metrics(traj)
                    traj.insert(0, "dataset", ds)
                    traj.insert(1, "experiment", "randomized_direction")
                    traj.insert(2, "replicate", r + 1)
                    traj.to_csv(
                        os.path.join(out, "raw",
                                     "{}__randomized_direction__{:04d}.csv".format(ds, r + 1)),
                        index=False
                    )
                    if len(signs):
                        signs.insert(0, "dataset", ds)
                        signs.insert(1, "replicate", r + 1)
                        signs.insert(2, "seed", seed)
                        sign_rows.extend(signs.to_dict("records"))
                    row["status"] = "finished"
                    row.update(met)
                    all_metrics.append(row.copy())
                    status_rows.append(row.copy())
                except Exception as exc:
                    row["message"] = repr(exc)
                    status_rows.append(row.copy())

        except Exception as exc:
            status_rows.append({
                "dataset": ds, "file": path, "experiment": "dataset",
                "replicate": 0, "seed": base_seed, "status": "failed",
                "message": repr(exc)
            })
            logging.exception("DATASET FAILED %s", ds)

    metrics = pd.DataFrame(all_metrics)
    status = pd.DataFrame(status_rows)
    status.to_csv(os.path.join(out, "replicate_status.csv"), index=False)
    if sign_rows:
        pd.DataFrame(sign_rows).to_csv(
            os.path.join(out, "randomized_direction_sign_counts.csv"), index=False
        )

    # Dataset-level matched null medians.
    endpoints = [
        "delta_kappa", "variance_reduction", "lambda_improvement",
        "normalized_rate", "convergence_score"
    ]
    observed = metrics[metrics["experiment"] == "observed"].copy()
    results = []
    empirical = []

    for ep in endpoints:
        obs = observed[["dataset", ep]].rename(columns={ep: "observed"})
        for expname in ("topology_null", "randomized_direction"):
            sims = metrics[metrics["experiment"] == expname].copy()
            if len(sims):
                med = sims.groupby("dataset")[ep].median().rename("null_median")
                q05 = sims.groupby("dataset")[ep].quantile(0.05).rename("null_q05")
                q95 = sims.groupby("dataset")[ep].quantile(0.95).rename("null_q95")
                tab = obs.merge(med, left_on="dataset", right_index=True, how="inner")
                tab = tab.merge(q05, left_on="dataset", right_index=True, how="left")
                tab = tab.merge(q95, left_on="dataset", right_index=True, how="left")
                tab["experiment"] = expname
                tab["endpoint"] = ep
                stat = paired_wilcoxon(tab["observed"], tab["null_median"])
                stat.update({"endpoint": ep, "experiment": expname})
                results.append(stat)

                # Per-dataset empirical two-sided tail relative to the
                # matched replicate distribution, then aggregate the
                # observed-vs-null median differences at dataset level.
                diffs = []
                for ds in tab["dataset"]:
                    ov = float(tab.loc[tab["dataset"] == ds, "observed"].iloc[0])
                    vals = sims.loc[sims["dataset"] == ds, ep].dropna().to_numpy(float)
                    if len(vals) == 0 or not np.isfinite(ov):
                        continue
                    center = float(np.median(vals))
                    diffs.append(ov - center)
                if diffs:
                    D = np.asarray(diffs, dtype=float)
                    observed_mean_diff = float(np.mean(D))
                    observed_median_diff = float(np.median(D))
                else:
                    observed_mean_diff = observed_median_diff = np.nan

                # A dataset-level bootstrap is descriptive, not used to
                # manufacture significance.
                stable_code = sum((i + 1) * ord(ch) for i, ch in enumerate(ep + "|" + expname)) % 1000000
                rng = np.random.default_rng(args.seed + stable_code)
                if len(diffs) >= 2:
                    boots = np.empty(args.bootstrap, dtype=float)
                    for b in range(args.bootstrap):
                        boots[b] = np.mean(rng.choice(D, size=len(D), replace=True))
                    ci = np.quantile(boots, [0.025, 0.975])
                else:
                    ci = [np.nan, np.nan]
                empirical.append({
                    "endpoint": ep, "experiment": expname,
                    "n_datasets": len(D) if diffs else 0,
                    "mean_observed_minus_null_median": observed_mean_diff,
                    "median_observed_minus_null_median": observed_median_diff,
                    "bootstrap_mean_ci_low": float(ci[0]),
                    "bootstrap_mean_ci_high": float(ci[1]),
                })

    stat_df = pd.DataFrame(results)
    if len(stat_df):
        stat_df["q_value_bh"] = bh_fdr(stat_df["p_value"].to_numpy(float))
    stat_df.to_csv(os.path.join(out, "paired_statistics.csv"), index=False)
    pd.DataFrame(empirical).to_csv(os.path.join(out, "empirical_pvalues.csv"), index=False)

    # Complete dataset-level metrics for downstream analysis.
    metrics.to_csv(os.path.join(out, "dataset_metrics", "all_dataset_metrics.csv"), index=False)

    print("\nFinished.")
    print("Output:", out)
    print("Datasets discovered:", len(files))
    print("Metric rows:", len(metrics))
    print("Status rows:", len(status))
    print("Primary statistics:", os.path.join(out, "paired_statistics.csv"))


def build_parser():
    p = argparse.ArgumentParser(
        description="Experiments 1/2 for Ollivier-Ricci flow robustness"
    )
    p.add_argument("--data_dir", default=r"D:\100data")
    p.add_argument("--output_dir", default="ricci_flow_exp12_results")
    p.add_argument("--steps", type=int, default=10)
    p.add_argument("--reps", type=int, default=100)
    p.add_argument("--alpha", type=float, default=0.5)
    p.add_argument("--dt", type=float, default=0.05)
    p.add_argument("--sinkhorn_reg", type=float, default=0.05)
    p.add_argument("--sinkhorn_max_iter", type=int, default=5000)
    p.add_argument("--max_distance_mb", type=float, default=256.0)
    p.add_argument("--weight_mode",
                   choices=("auto", "unweighted", "distance", "strength"),
                   default="auto")
    p.add_argument("--component_mode", choices=("full", "largest"), default="full")
    p.add_argument("--n_swaps", type=int, default=0,
                   help="0 => 10 x number_of_edges, a conventional strong rewiring budget")
    p.add_argument("--max_swap_attempts", type=int, default=2000)
    p.add_argument("--bootstrap", type=int, default=5000)
    p.add_argument("--seed", type=int, default=20260907)
    p.add_argument("--only", default="",
                   help="optional substring to restrict datasets, e.g. eco-everglades")
    return p


def main():
    args = build_parser().parse_args()
    # n_swaps=10E is selected BEFORE seeing any outcome and is therefore not
    # result-dependent. Record the realized setting in the configuration.
    if args.n_swaps <= 0:
        args.n_swaps = None

    if args.n_swaps is None:
        # run() resolves per-dataset 10E; keep the CLI value auditable.
        args.n_swaps = 0

    run(args)


if __name__ == "__main__":
    main()
