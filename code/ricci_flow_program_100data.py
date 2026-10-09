#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Robust and accelerated universal Ollivier-Ricci flow experiment.

Designed for Python 3.6+ and heterogeneous graph datasets stored under
``data/<dataset_name>/``.  The program recursively discovers graph files,
reads whitespace/comma separated edge lists and Matrix Market files, validates
edge lengths, computes regularized Ollivier curvature, performs weighted Ricci
flow, and writes resumable per-dataset results.

Important defaults
------------------
* Undirected simple graph, matching the original experiment.
* Full graph after removing self-loops and isolated vertices.
* Weight mode ``auto``:
    - README says unweighted -> every edge length is 1;
    - no usable third column -> every edge length is 1;
    - non-positive/signed third column -> every edge length is 1 and a warning
      is recorded, because Dijkstra requires a non-negative metric;
    - all third-column values are finite and positive -> preserve the original
      interpretation and use them as edge lengths.
* Exact shortest-path distances are reused once per flow step.
* Sinkhorn errors are never silently converted to zero.
* Very dense/large graphs are skipped by default instead of hanging for days.
  Use ``--force_large`` only when you intentionally want to force the exact
  calculation.

Direct use
----------
    python3 ricci_flow_universal_experiment_correct_fast.py

The default command scans ``D:\100data`` and stores outputs in
``results_correct_fast``.
"""

from __future__ import print_function

# Avoid nested BLAS/OpenMP oversubscription when several worker processes run.
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import argparse
import csv
import html
import logging
import multiprocessing as mp
import re
import sys
import time
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import ot

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra as sparse_dijkstra
from scipy.sparse.linalg import eigsh

warnings.filterwarnings("ignore")

SUPPORTED_EXTENSIONS = (".edges", ".txt", ".mtx", ".csv")
DEFAULT_RESULT_DIR = "results_correct_fast"

HEADER_LEFT = {
    "source", "src", "from", "node1", "vertex1", "u", "id1", "tail"
}
HEADER_RIGHT = {
    "target", "dst", "to", "node2", "vertex2", "v", "id2", "head"
}


class DatasetSkipped(Exception):
    """Raised when a graph is valid but intentionally skipped for practicality."""


class DatasetFormatError(Exception):
    """Raised when a graph file cannot be interpreted safely."""


# ---------------------------------------------------------------------------
# Logging and output directories
# ---------------------------------------------------------------------------

RESULT_DIR = DEFAULT_RESULT_DIR


def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def prepare_dirs():
    for path in (
        RESULT_DIR,
        os.path.join(RESULT_DIR, "figures"),
        os.path.join(RESULT_DIR, "npz"),
        os.path.join(RESULT_DIR, "per_dataset"),
        os.path.join(RESULT_DIR, "logs")
    ):
        ensure_dir(path)


def setup_logger():
    logging.basicConfig(
        filename=os.path.join(RESULT_DIR, "logs", "experiment.log"),
        level=logging.INFO,
        format="%(asctime)s | %(message)s"
    )


def log(message):
    print(message, flush=True)
    logging.info(message)


# ---------------------------------------------------------------------------
# README metadata
# ---------------------------------------------------------------------------


def _html_to_text(value):
    value = html.unescape(value)
    value = re.sub(r"<script[^>]*>.*?</script>", " ", value,
                   flags=re.IGNORECASE | re.DOTALL)
    value = re.sub(r"<style[^>]*>.*?</style>", " ", value,
                   flags=re.IGNORECASE | re.DOTALL)
    value = re.sub(r"<[^>]+>", " ", value)
    return re.sub(r"\s+", " ", value).strip()


def read_local_metadata(graph_path):
    """Read lightweight metadata from README files in the dataset directory."""
    directory = os.path.dirname(graph_path)
    metadata = {
        "readme": "",
        "edge_weights": "unknown",
        "direction": "unknown"
    }

    try:
        candidates = []
        for name in os.listdir(directory):
            lower = name.lower()
            if lower.startswith("readme") and lower.endswith((".html", ".htm", ".txt", ".md")):
                candidates.append(os.path.join(directory, name))
        candidates.sort()
    except Exception:
        candidates = []

    for path in candidates:
        try:
            with open(path, "r", errors="ignore") as handle:
                raw = handle.read(2 * 1024 * 1024)
        except Exception:
            continue

        text = _html_to_text(raw)
        lower = text.lower()
        metadata["readme"] = os.path.basename(path)

        # Check "unweighted" before "weighted" because it contains the latter.
        weight_pos = lower.find("edge weights")
        if weight_pos < 0:
            weight_pos = lower.find("edge weight")
        weight_window = lower[weight_pos:weight_pos + 300] if weight_pos >= 0 else lower

        if "unweighted" in weight_window:
            metadata["edge_weights"] = "unweighted"
        elif "weighted" in weight_window:
            metadata["edge_weights"] = "weighted"

        format_pos = lower.find("format")
        format_window = lower[format_pos:format_pos + 200] if format_pos >= 0 else lower
        if "directed" in format_window and "undirected" not in format_window:
            metadata["direction"] = "directed"
        elif "undirected" in format_window:
            metadata["direction"] = "undirected"
        break

    return metadata


# ---------------------------------------------------------------------------
# Graph parsing
# ---------------------------------------------------------------------------


def detect_format(path):
    return os.path.splitext(path)[1].lower()


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
    left = parts[0].strip().strip('"').strip("'").lower()
    right = parts[1].strip().strip('"').strip("'").lower()
    return left in HEADER_LEFT and right in HEADER_RIGHT


def parse_optional_weight(parts):
    if len(parts) < 3:
        return None
    try:
        return float(parts[2])
    except Exception:
        return None


def parse_edge_list(path):
    records = []
    stats = {
        "total_lines": 0,
        "comment_or_blank_lines": 0,
        "header_lines": 0,
        "invalid_lines": 0,
        "self_loops": 0,
        "numeric_weight_rows": 0,
        "missing_weight_rows": 0,
        "format": "edge_list",
        "matrix_rows": "",
        "matrix_cols": ""
    }
    invalid_examples = []

    with open(path, "r", errors="ignore") as handle:
        for raw_line in handle:
            stats["total_lines"] += 1
            line = raw_line.strip()
            if not line or line.startswith("#") or line.startswith("%"):
                stats["comment_or_blank_lines"] += 1
                continue

            parts = split_edge_line(line)
            if looks_like_header(parts):
                stats["header_lines"] += 1
                continue
            if len(parts) < 2:
                stats["invalid_lines"] += 1
                if len(invalid_examples) < 5:
                    invalid_examples.append(line)
                continue

            try:
                u = parse_node(parts[0])
                v = parse_node(parts[1])
            except Exception:
                stats["invalid_lines"] += 1
                if len(invalid_examples) < 5:
                    invalid_examples.append(line)
                continue

            if u == v:
                stats["self_loops"] += 1
                continue

            raw_weight = parse_optional_weight(parts)
            if raw_weight is None:
                stats["missing_weight_rows"] += 1
            else:
                stats["numeric_weight_rows"] += 1

            records.append((u, v, raw_weight))

    stats["invalid_examples"] = " | ".join(invalid_examples)
    return records, stats


def parse_matrix_market(path):
    records = []
    stats = {
        "total_lines": 0,
        "comment_or_blank_lines": 0,
        "header_lines": 0,
        "invalid_lines": 0,
        "self_loops": 0,
        "numeric_weight_rows": 0,
        "missing_weight_rows": 0,
        "format": "matrix_market",
        "matrix_rows": "",
        "matrix_cols": ""
    }

    with open(path, "r", errors="ignore") as handle:
        lines = handle.readlines()

    stats["total_lines"] = len(lines)
    header = ""
    first_nonempty = None
    for index, raw_line in enumerate(lines):
        stripped = raw_line.strip()
        if stripped:
            first_nonempty = index
            header = stripped
            break

    if first_nonempty is None:
        stats["invalid_examples"] = "empty file"
        return records, stats

    field = "real"
    matrix_format = "coordinate"
    if header.lower().startswith("%%matrixmarket"):
        tokens = header.split()
        if len(tokens) >= 5:
            matrix_format = tokens[2].lower()
            field = tokens[3].lower()
        if matrix_format != "coordinate":
            raise DatasetFormatError(
                "Only Matrix Market coordinate format is supported: {}".format(path)
            )
        cursor = first_nonempty + 1
    else:
        # Some repository files use the .mtx suffix without a formal banner.
        cursor = first_nonempty

    while cursor < len(lines):
        stripped = lines[cursor].strip()
        if not stripped or stripped.startswith("%"):
            stats["comment_or_blank_lines"] += 1
            cursor += 1
            continue
        break

    # The first non-comment row in a proper Matrix Market file is dimensions.
    if header.lower().startswith("%%matrixmarket") and cursor < len(lines):
        dims = lines[cursor].split()
        if len(dims) < 3:
            raise DatasetFormatError("Invalid Matrix Market dimension line in {}".format(path))
        try:
            stats["matrix_rows"] = int(dims[0])
            stats["matrix_cols"] = int(dims[1])
        except Exception:
            raise DatasetFormatError("Invalid Matrix Market dimensions in {}".format(path))
        stats["header_lines"] += 1
        cursor += 1

    invalid_examples = []
    for raw_line in lines[cursor:]:
        line = raw_line.strip()
        if not line or line.startswith("%"):
            stats["comment_or_blank_lines"] += 1
            continue

        parts = line.split()
        if len(parts) < 2:
            stats["invalid_lines"] += 1
            if len(invalid_examples) < 5:
                invalid_examples.append(line)
            continue

        try:
            u = int(parts[0])
            v = int(parts[1])
        except Exception:
            stats["invalid_lines"] += 1
            if len(invalid_examples) < 5:
                invalid_examples.append(line)
            continue

        if u == v:
            stats["self_loops"] += 1
            continue

        raw_weight = None
        if field != "pattern" and len(parts) >= 3:
            try:
                raw_weight = float(parts[2])
            except Exception:
                raw_weight = None

        # An explicitly stored zero is not an adjacency edge.
        if raw_weight == 0.0:
            continue

        if raw_weight is None:
            stats["missing_weight_rows"] += 1
        else:
            stats["numeric_weight_rows"] += 1

        records.append((u, v, raw_weight))

    stats["invalid_examples"] = " | ".join(invalid_examples)
    return records, stats


def choose_weight_mode(records, metadata, requested_mode,
                       positive_weight_interpretation):
    numeric = [weight for _, _, weight in records if weight is not None]
    all_rows_have_numeric_weight = len(numeric) == len(records) and len(records) > 0

    if requested_mode != "auto":
        return requested_mode, "explicit --weight_mode {}".format(requested_mode)

    if metadata.get("edge_weights") == "unweighted":
        return "unweighted", "README metadata says unweighted"

    if not all_rows_have_numeric_weight:
        return "unweighted", "not every edge has a numeric third column"

    numeric_array = np.asarray(numeric, dtype=float)
    if not np.all(np.isfinite(numeric_array)):
        return "unweighted", "non-finite input weights cannot define path lengths"

    if np.any(numeric_array <= 0.0):
        return "unweighted", (
            "signed/zero input weights cannot be used by Dijkstra; "
            "using topology-only unit lengths"
        )

    if np.allclose(numeric_array, 1.0):
        return "unweighted", "all input weights equal one"

    return positive_weight_interpretation, (
        "all input weights are finite and positive; auto interpretation={}".format(
            positive_weight_interpretation
        )
    )


def transform_length(raw_weight, mode):
    if mode == "unweighted":
        return 1.0
    if raw_weight is None or not np.isfinite(raw_weight) or raw_weight <= 0.0:
        raise DatasetFormatError(
            "Weight mode '{}' requires every edge value to be finite and positive".format(mode)
        )
    if mode == "distance":
        return float(raw_weight)
    if mode == "strength":
        return 1.0 / float(raw_weight)
    raise ValueError("Unknown weight mode: {}".format(mode))


def build_graph(records, weight_mode):
    graph = nx.Graph()
    duplicate_count = 0

    for u, v, raw_weight in records:
        length = transform_length(raw_weight, weight_mode)
        if graph.has_edge(u, v):
            duplicate_count += 1
            # For parallel observations, retain the shortest usable route.
            if length < graph[u][v]["weight"]:
                graph[u][v]["weight"] = length
        else:
            graph.add_edge(u, v, weight=length)

    graph.graph["duplicate_edges"] = duplicate_count
    return graph


def normalize_weights(graph):
    weights = np.asarray(
        [float(data.get("weight", 1.0)) for _, _, data in graph.edges(data=True)],
        dtype=float
    )
    if weights.size == 0:
        raise DatasetFormatError("The parsed graph contains no usable edges")
    if not np.all(np.isfinite(weights)) or np.any(weights <= 0.0):
        raise DatasetFormatError("All path lengths must be finite and strictly positive")

    mean_weight = float(np.mean(weights))
    if not np.isfinite(mean_weight) or mean_weight <= 0.0:
        raise DatasetFormatError("Mean path length is invalid")

    for _, _, data in graph.edges(data=True):
        data["weight"] = float(data.get("weight", 1.0)) / mean_weight


def load_graph(path, args):
    metadata = read_local_metadata(path)
    if detect_format(path) == ".mtx":
        records, parse_stats = parse_matrix_market(path)
    else:
        records, parse_stats = parse_edge_list(path)

    if not records:
        message = "No valid edges were parsed from {}".format(path)
        examples = parse_stats.get("invalid_examples", "")
        if examples:
            message += "; examples: " + examples
        raise DatasetFormatError(message)

    weight_mode, weight_reason = choose_weight_mode(
        records,
        metadata,
        args.weight_mode,
        args.positive_weight_interpretation
    )
    graph = build_graph(records, weight_mode)
    graph.remove_edges_from(list(nx.selfloop_edges(graph)))
    graph.remove_nodes_from(list(nx.isolates(graph)))

    if graph.number_of_edges() == 0:
        raise DatasetFormatError("The graph is empty after validation")

    normalize_weights(graph)

    components = list(nx.connected_components(graph))
    component_sizes = sorted([len(component) for component in components], reverse=True)

    graph.graph.update({
        "source_path": path,
        "readme": metadata.get("readme", ""),
        "readme_edge_weights": metadata.get("edge_weights", "unknown"),
        "readme_direction": metadata.get("direction", "unknown"),
        "weight_mode": weight_mode,
        "weight_reason": weight_reason,
        "parse_stats": parse_stats,
        "original_nodes": graph.number_of_nodes(),
        "original_edges": graph.number_of_edges(),
        "component_count": len(components),
        "component_sizes": component_sizes,
        "largest_component_nodes": component_sizes[0] if component_sizes else 0
    })

    if args.component_mode == "largest" and len(components) > 1:
        largest_nodes = max(components, key=len)
        reduced = graph.subgraph(largest_nodes).copy()
        reduced.graph.update(graph.graph)
        reduced.graph["component_mode"] = "largest"
        reduced.graph["component_count_used"] = 1
        graph = reduced
    else:
        graph.graph["component_mode"] = "full"
        graph.graph["component_count_used"] = len(components)

    return graph


# ---------------------------------------------------------------------------
# Complexity guards
# ---------------------------------------------------------------------------


def graph_work_estimate(graph):
    degree = dict(graph.degree())
    transport_cells = 0
    max_edge_cells = 0
    for u, v in graph.edges():
        cells = (degree[u] + 1) * (degree[v] + 1)
        transport_cells += cells
        if cells > max_edge_cells:
            max_edge_cells = cells

    shortest_path_work = graph.number_of_nodes() * graph.number_of_edges()
    return {
        "transport_cells_per_step": int(transport_cells),
        "max_edge_transport_cells": int(max_edge_cells),
        "shortest_path_work": int(shortest_path_work)
    }


def enforce_work_limits(graph, args):
    estimate = graph_work_estimate(graph)
    graph.graph.update(estimate)

    reasons = []
    if args.max_transport_cells > 0 and (
            estimate["transport_cells_per_step"] > args.max_transport_cells):
        reasons.append(
            "transport cells/step {} > limit {}".format(
                estimate["transport_cells_per_step"], args.max_transport_cells
            )
        )
    if args.max_shortest_path_work > 0 and (
            estimate["shortest_path_work"] > args.max_shortest_path_work):
        reasons.append(
            "nodes*edges {} > limit {}".format(
                estimate["shortest_path_work"], args.max_shortest_path_work
            )
        )

    if reasons and not args.force_large:
        raise DatasetSkipped(
            "; ".join(reasons) +
            ". This graph is valid but exact per-edge Sinkhorn Ricci flow is "
            "impractical under the default safety limits. Use --force_large "
            "only when you accept a potentially very long run."
        )

    return estimate


# ---------------------------------------------------------------------------
# Exact reusable shortest paths and curvature
# ---------------------------------------------------------------------------


def node_distribution(graph, node, alpha):
    neighbors = list(graph.neighbors(node))
    if not neighbors:
        return [node], np.asarray([1.0], dtype=float)
    remaining = (1.0 - alpha) / float(len(neighbors))
    support = [node] + neighbors
    mass = np.asarray([alpha] + [remaining] * len(neighbors), dtype=float)
    return support, mass


class CurvatureWorkspace(object):
    """Topology cache plus exact per-step distance workspace."""

    def __init__(self, graph, alpha, backend, max_distance_mb):
        self.graph = graph
        self.alpha = alpha
        self.nodes = list(graph.nodes())
        self.edges = list(graph.edges())
        self.node_index = dict((node, index) for index, node in enumerate(self.nodes))
        self.support = {}
        self.mass = {}
        self.support_index = {}
        self.distance_matrix = None
        self.target_distances = None
        self.max_distance_mb = float(max_distance_mb)

        for node in self.nodes:
            support, mass = node_distribution(graph, node, alpha)
            self.support[node] = support
            self.mass[node] = mass
            self.support_index[node] = np.asarray(
                [self.node_index[item] for item in support], dtype=np.int64
            )

        matrix_mb = len(self.nodes) * len(self.nodes) * 8.0 / (1024.0 * 1024.0)
        self.estimated_matrix_mb = matrix_mb
        if backend == "auto":
            self.backend = "matrix" if matrix_mb <= max_distance_mb else "batched"
        else:
            self.backend = backend

        self.required_targets = None
        if self.backend == "batched":
            required = dict((node, set()) for node in self.nodes)
            for u, v in self.edges:
                ys = self.support[v]
                for source in self.support[u]:
                    required[source].update(ys)
            self.required_targets = required

    def build_sparse_adjacency(self):
        rows = []
        cols = []
        values = []
        for u, v, data in self.graph.edges(data=True):
            weight = float(data.get("weight", 1.0))
            if not np.isfinite(weight) or weight <= 0.0:
                raise ValueError("Dijkstra requires finite positive edge lengths")
            i = self.node_index[u]
            j = self.node_index[v]
            rows.extend((i, j))
            cols.extend((j, i))
            values.extend((weight, weight))
        return csr_matrix(
            (values, (rows, cols)),
            shape=(len(self.nodes), len(self.nodes)),
            dtype=float
        )

    def prepare_step(self):
        start = time.time()
        adjacency = self.build_sparse_adjacency()

        if self.backend == "matrix":
            self.distance_matrix = sparse_dijkstra(
                csgraph=adjacency,
                directed=False,
                return_predecessors=False
            )
            self.target_distances = None

        elif self.backend == "batched":
            node_count = max(1, len(self.nodes))
            bytes_limit = max(1.0, self.max_distance_mb) * 1024.0 * 1024.0
            batch_size = int(bytes_limit // (node_count * 8.0))
            batch_size = max(1, min(batch_size, 256))

            distances = {}
            source_nodes = [node for node in self.nodes if self.required_targets[node]]
            for start_index in range(0, len(source_nodes), batch_size):
                batch_nodes = source_nodes[start_index:start_index + batch_size]
                batch_indices = [self.node_index[node] for node in batch_nodes]
                batch_matrix = sparse_dijkstra(
                    csgraph=adjacency,
                    directed=False,
                    indices=batch_indices,
                    return_predecessors=False
                )
                batch_matrix = np.atleast_2d(batch_matrix)

                for row_index, source in enumerate(batch_nodes):
                    target_map = {}
                    for target in self.required_targets[source]:
                        value = float(batch_matrix[row_index, self.node_index[target]])
                        target_map[target] = value
                    distances[source] = target_map

            self.target_distances = distances
            self.distance_matrix = None
        else:
            raise ValueError("Unknown distance backend: {}".format(self.backend))

        return time.time() - start

    def edge_payload(self, index):
        u, v = self.edges[index]
        xs = self.support[u]
        ys = self.support[v]

        if self.backend == "matrix":
            x_index = self.support_index[u]
            y_index = self.support_index[v]
            cost = self.distance_matrix[np.ix_(x_index, y_index)]
            edge_distance = float(
                self.distance_matrix[self.node_index[u], self.node_index[v]]
            )
        else:
            cost = np.empty((len(xs), len(ys)), dtype=float)
            for i, source in enumerate(xs):
                source_map = self.target_distances[source]
                for j, target in enumerate(ys):
                    cost[i, j] = source_map[target]
            edge_distance = float(self.target_distances[u][v])

        if not np.isfinite(edge_distance) or edge_distance <= 0.0:
            raise ValueError("Invalid edge distance for ({}, {})".format(u, v))
        if not np.all(np.isfinite(cost)):
            raise ValueError("Non-finite local transport cost for ({}, {})".format(u, v))

        return u, v, self.mass[u], self.mass[v], cost, edge_distance


_WORKSPACE = None
_WORKER_REG = 0.05
_WORKER_MAX_ITER = 5000


def _edge_curvature_worker(index):
    u, v, mass_u, mass_v, cost, edge_distance = _WORKSPACE.edge_payload(index)
    try:
        value = ot.sinkhorn2(
            mass_u,
            mass_v,
            cost,
            reg=_WORKER_REG,
            numItermax=_WORKER_MAX_ITER
        )
        value = float(np.asarray(value).reshape(-1)[0])
    except Exception as exc:
        raise RuntimeError(
            "Sinkhorn failed for edge ({}, {}): {}".format(u, v, repr(exc))
        )

    if not np.isfinite(value) or value < 0.0:
        raise RuntimeError(
            "Sinkhorn returned invalid value {} for edge ({}, {})".format(value, u, v)
        )

    curvature = 1.0 - value / edge_distance
    if not np.isfinite(curvature):
        raise RuntimeError("Non-finite curvature for edge ({}, {})".format(u, v))
    if curvature > 1.0 and curvature < 1.0 + 1e-10:
        curvature = 1.0
    return float(curvature)


def compute_curvature(workspace, workers, reg, max_iter):
    distance_time = workspace.prepare_step()

    global _WORKSPACE
    global _WORKER_REG
    global _WORKER_MAX_ITER
    _WORKSPACE = workspace
    _WORKER_REG = float(reg)
    _WORKER_MAX_ITER = int(max_iter)

    indices = list(range(len(workspace.edges)))
    start = time.time()

    use_parallel = workers > 1 and sys.platform.startswith("linux")
    if use_parallel:
        context = mp.get_context("fork")
        chunksize = max(1, len(indices) // max(1, workers * 8))
        pool = context.Pool(processes=workers)
        try:
            curvatures = pool.map(
                _edge_curvature_worker,
                indices,
                chunksize=chunksize
            )
        finally:
            pool.close()
            pool.join()
    else:
        curvatures = [_edge_curvature_worker(index) for index in indices]

    sinkhorn_time = time.time() - start
    values = np.asarray(curvatures, dtype=float)
    if values.size == 0 or not np.all(np.isfinite(values)):
        raise RuntimeError("Curvature computation returned no finite values")

    return (
        float(np.mean(values)),
        float(np.var(values)),
        dict(zip(workspace.edges, curvatures)),
        distance_time,
        sinkhorn_time
    )


# ---------------------------------------------------------------------------
# Algebraic connectivity and Ricci flow update
# ---------------------------------------------------------------------------


class LambdaWorkspace(object):
    def __init__(self, graph):
        self.graph = graph
        self.nodes = list(graph.nodes())
        self.v0 = None

    def compute(self):
        node_count = self.graph.number_of_nodes()
        if node_count < 3:
            return 0.0

        component_count = self.graph.graph.get("component_count_used")
        if component_count is None:
            component_count = nx.number_connected_components(self.graph)
        if component_count > 1:
            return 0.0

        laplacian = nx.laplacian_matrix(
            self.graph,
            nodelist=self.nodes,
            weight="weight"
        ).astype(float)

        try:
            if node_count <= 128:
                eigenvalues = np.linalg.eigvalsh(laplacian.toarray())
                value = float(np.sort(np.real(eigenvalues))[1])
            else:
                eigenvalues, eigenvectors = eigsh(
                    laplacian,
                    k=2,
                    which="SM",
                    return_eigenvectors=True,
                    v0=self.v0
                )
                order = np.argsort(np.real(eigenvalues))
                eigenvalues = np.real(eigenvalues[order])
                eigenvectors = np.real(eigenvectors[:, order])
                value = float(eigenvalues[1])
                self.v0 = eigenvectors[:, 1]
        except Exception as exc:
            raise RuntimeError("lambda1 computation failed: {}".format(repr(exc)))

        if abs(value) < 1e-12:
            value = 0.0
        if value < 0.0 and value > -1e-9:
            value = 0.0
        return value


def ricci_flow_step(graph, curvature, dt):
    edges = list(curvature.keys())
    candidate_weights = []

    for u, v in edges:
        old_weight = float(graph[u][v].get("weight", 1.0))
        exponent = -float(dt) * float(curvature[(u, v)])
        if not np.isfinite(exponent) or abs(exponent) > 700.0:
            raise RuntimeError(
                "Ricci update overflow risk on edge ({}, {}), exponent={}".format(
                    u, v, exponent
                )
            )
        candidate = old_weight * np.exp(exponent)
        if not np.isfinite(candidate) or candidate <= 0.0:
            raise RuntimeError("Invalid updated edge weight on ({}, {})".format(u, v))
        candidate_weights.append(float(candidate))

    mean_weight = float(np.mean(candidate_weights))
    if not np.isfinite(mean_weight) or mean_weight <= 0.0:
        raise RuntimeError("Ricci flow produced an invalid mean edge weight")

    for index, (u, v) in enumerate(edges):
        graph[u][v]["weight"] = candidate_weights[index] / mean_weight


def run_ricci_flow(graph, args):
    records = []
    workspace = CurvatureWorkspace(
        graph,
        alpha=args.alpha,
        backend=args.distance_backend,
        max_distance_mb=args.max_distance_mb
    )
    lambda_workspace = LambdaWorkspace(graph)

    log(
        "Distance backend: {} (full distance matrix estimate: {:.2f} MB)".format(
            workspace.backend, workspace.estimated_matrix_mb
        )
    )
    log("Curvature workers: {}".format(args.workers))

    for step in range(args.steps):
        step_start = time.time()
        mean_k, var_k, edge_k, distance_time, sinkhorn_time = compute_curvature(
            workspace,
            workers=args.workers,
            reg=args.sinkhorn_reg,
            max_iter=args.sinkhorn_max_iter
        )

        lambda_start = time.time()
        lambda1 = lambda_workspace.compute()
        lambda_time = time.time() - lambda_start
        elapsed = time.time() - step_start

        records.append({
            "step": step,
            "mean_curvature": mean_k,
            "curvature_variance": var_k,
            "lambda1": lambda1,
            "nodes": graph.number_of_nodes(),
            "edges": graph.number_of_edges(),
            "time": elapsed,
            "distance_time": distance_time,
            "sinkhorn_time": sinkhorn_time,
            "lambda_time": lambda_time
        })

        log(
            "\nStep {}\n"
            "k={:.6f}\n"
            "var={:.6f}\n"
            "lambda1={:.12e}\n"
            "distance_time={:.3f}s\n"
            "sinkhorn_time={:.3f}s\n"
            "lambda_time={:.3f}s\n"
            "total_time={:.3f}s\n".format(
                step,
                mean_k,
                var_k,
                lambda1,
                distance_time,
                sinkhorn_time,
                lambda_time,
                elapsed
            )
        )

        # The final update is not observed by any output, so skipping it saves
        # time without changing the recorded experiment results.
        if step < args.steps - 1:
            ricci_flow_step(graph, edge_k, dt=args.dt)

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Result handling
# ---------------------------------------------------------------------------


def discover_graph_files(data_dir):
    files = []
    for root, _, filenames in os.walk(data_dir):
        for filename in filenames:
            if filename.lower().endswith(SUPPORTED_EXTENSIONS):
                files.append(os.path.join(root, filename))
    files.sort()
    return files


def safe_dataset_id(filepath, data_dir):
    relative = os.path.relpath(filepath, data_dir)
    relative = os.path.splitext(relative)[0]
    safe = re.sub(r"[^A-Za-z0-9._-]+", "__", relative)
    return safe.strip("._-") or "dataset"


def calculate_slope(dataframe):
    clean = dataframe[["mean_curvature", "lambda1"]].replace(
        [np.inf, -np.inf], np.nan
    ).dropna()
    if len(clean) < 2:
        return 0.0, 0.0

    x = clean["mean_curvature"].values.astype(float)
    y = clean["lambda1"].values.astype(float)
    design = np.column_stack((x, np.ones_like(x)))
    coefficient, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    return float(coefficient[0]), float(coefficient[1])


def save_dataset_result(dataframe, dataset_id):
    csv_path = os.path.join(RESULT_DIR, "per_dataset", dataset_id + ".csv")
    npz_path = os.path.join(RESULT_DIR, "npz", dataset_id + ".npz")
    dataframe.to_csv(csv_path, index=False)
    np.savez(
        npz_path,
        step=dataframe["step"].values,
        curvature=dataframe["mean_curvature"].values,
        variance=dataframe["curvature_variance"].values,
        lambda1=dataframe["lambda1"].values
    )
    return csv_path


def normalize_curve(dataframe):
    x = dataframe["step"].values.astype(float)
    y = dataframe["lambda1"].values.astype(float)

    if len(x) < 2 or x[-1] == x[0]:
        xn = np.zeros_like(x)
    else:
        xn = (x - x[0]) / (x[-1] - x[0])

    if len(y) < 2 or not np.all(np.isfinite(y)) or abs(y[-1] - y[0]) < 1e-12:
        yn = np.zeros_like(y)
    else:
        yn = (y - y[0]) / (y[-1] - y[0])
    return xn, yn


def save_plots(all_results):
    figure_specs = [
        (
            "ricci_flow_collapse.png",
            "Normalized Ricci Flow Time",
            "Normalized Algebraic Connectivity",
            "Universal Ricci Flow Collapse",
            True,
            "step",
            "lambda1"
        ),
        (
            "curvature_evolution.png",
            "Flow step",
            "Mean Ollivier Curvature",
            "Curvature Evolution",
            False,
            "step",
            "mean_curvature"
        ),
        (
            "lambda1_evolution.png",
            "Flow step",
            "lambda1",
            "Spectral Evolution",
            False,
            "step",
            "lambda1"
        )
    ]

    for filename, xlabel, ylabel, title, normalize, x_column, y_column in figure_specs:
        plt.figure(figsize=(9, 6))
        for name, dataframe in all_results.items():
            if normalize:
                x, y = normalize_curve(dataframe)
            else:
                x = dataframe[x_column].values
                y = dataframe[y_column].values
            plt.plot(x, y, marker="o", linewidth=1.8, label=name)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        if all_results:
            plt.legend(fontsize=7)
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(RESULT_DIR, "figures", filename), dpi=300)
        plt.close()


def summary_from_graph(graph, filepath, dataset_id, status, message, slope, intercept):
    parse_stats = graph.graph.get("parse_stats", {}) if graph is not None else {}
    return {
        "dataset": dataset_id,
        "file": filepath,
        "status": status,
        "message": message,
        "nodes": graph.number_of_nodes() if graph is not None else "",
        "edges": graph.number_of_edges() if graph is not None else "",
        "components": graph.graph.get("component_count_used", "") if graph is not None else "",
        "weight_mode": graph.graph.get("weight_mode", "") if graph is not None else "",
        "weight_reason": graph.graph.get("weight_reason", "") if graph is not None else "",
        "readme": graph.graph.get("readme", "") if graph is not None else "",
        "readme_direction": graph.graph.get("readme_direction", "") if graph is not None else "",
        "duplicates": graph.graph.get("duplicate_edges", "") if graph is not None else "",
        "invalid_lines": parse_stats.get("invalid_lines", ""),
        "transport_cells_per_step": graph.graph.get("transport_cells_per_step", "") if graph is not None else "",
        "shortest_path_work": graph.graph.get("shortest_path_work", "") if graph is not None else "",
        "slope": slope,
        "intercept": intercept
    }


def analyze_single_graph(filepath, args):
    log("\nProcessing: {}".format(filepath))
    graph = load_graph(filepath, args)
    estimate = enforce_work_limits(graph, args)
    parse_stats = graph.graph.get("parse_stats", {})

    log(
        "\nNodes: {}\n"
        "Edges: {}\n"
        "Connected components (used graph): {}\n"
        "Largest component nodes: {}\n"
        "Weight mode: {}\n"
        "Weight decision: {}\n"
        "README: {}\n"
        "Invalid input lines skipped: {}\n"
        "Transport cells per step: {}\n"
        "Shortest-path work estimate (nodes*edges): {}\n".format(
            graph.number_of_nodes(),
            graph.number_of_edges(),
            graph.graph.get("component_count_used", 0),
            graph.graph.get("largest_component_nodes", 0),
            graph.graph.get("weight_mode", ""),
            graph.graph.get("weight_reason", ""),
            graph.graph.get("readme", "none"),
            parse_stats.get("invalid_lines", 0),
            estimate["transport_cells_per_step"],
            estimate["shortest_path_work"]
        )
    )

    if args.inspect_only:
        return graph, None

    dataframe = run_ricci_flow(graph, args)
    return graph, dataframe


def run_all(args):
    files = discover_graph_files(args.data_dir)
    if args.only:
        token = args.only.lower()
        files = [path for path in files if token in os.path.relpath(path, args.data_dir).lower()]

    log("\nTotal datasets: {}\n".format(len(files)))
    if not files:
        raise DatasetFormatError(
            "No graph files found under '{}'. Supported extensions: {}".format(
                args.data_dir, ", ".join(SUPPORTED_EXTENSIONS)
            )
        )

    all_results = {}
    summaries = []

    for filepath in files:
        dataset_id = safe_dataset_id(filepath, args.data_dir)
        cached_path = os.path.join(RESULT_DIR, "per_dataset", dataset_id + ".csv")

        if os.path.exists(cached_path) and not args.overwrite and not args.inspect_only:
            try:
                dataframe = pd.read_csv(cached_path)
                all_results[dataset_id] = dataframe
                slope, intercept = calculate_slope(dataframe)
                summaries.append({
                    "dataset": dataset_id,
                    "file": filepath,
                    "status": "cached",
                    "message": "existing per-dataset result reused",
                    "nodes": int(dataframe["nodes"].iloc[0]) if len(dataframe) else "",
                    "edges": int(dataframe["edges"].iloc[0]) if len(dataframe) else "",
                    "components": "",
                    "weight_mode": "",
                    "weight_reason": "",
                    "readme": "",
                    "readme_direction": "",
                    "duplicates": "",
                    "invalid_lines": "",
                    "transport_cells_per_step": "",
                    "shortest_path_work": "",
                    "slope": slope,
                    "intercept": intercept
                })
                log("\nCached result reused: {}\n".format(dataset_id))
                continue
            except Exception as exc:
                log("Could not reuse cached result {}: {}".format(cached_path, repr(exc)))

        graph = None
        try:
            graph, dataframe = analyze_single_graph(filepath, args)
            if args.inspect_only:
                summaries.append(summary_from_graph(
                    graph, filepath, dataset_id, "inspected", "", "", ""
                ))
                continue

            save_dataset_result(dataframe, dataset_id)
            slope, intercept = calculate_slope(dataframe)
            all_results[dataset_id] = dataframe
            summaries.append(summary_from_graph(
                graph, filepath, dataset_id, "finished", "", slope, intercept
            ))
            log("\nFinished: {}\nslope={:.6f}\n".format(dataset_id, slope))

        except DatasetSkipped as exc:
            message = str(exc)
            summaries.append(summary_from_graph(
                graph, filepath, dataset_id, "skipped_large", message, "", ""
            ))
            log("\nSKIPPED: {}\n{}\n".format(dataset_id, message))

        except KeyboardInterrupt:
            log("\nInterrupted by user. Completed results have already been saved.\n")
            raise

        except Exception as exc:
            message = repr(exc)
            summaries.append(summary_from_graph(
                graph, filepath, dataset_id, "failed", message, "", ""
            ))
            log("\nFAILED: {}\n{}\n".format(dataset_id, message))

    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(os.path.join(RESULT_DIR, "dataset_summary.csv"), index=False)

    if all_results:
        combined = []
        for name, dataframe in all_results.items():
            item = dataframe.copy()
            item.insert(0, "dataset", name)
            combined.append(item)
        pd.concat(combined, ignore_index=True).to_csv(
            os.path.join(RESULT_DIR, "raw_results.csv"), index=False
        )
        save_plots(all_results)

    log(
        "\n================================\n"
        "ALL AVAILABLE EXPERIMENTS FINISHED\n"
        "Summary: {}\n"
        "================================\n".format(
            os.path.join(RESULT_DIR, "dataset_summary.csv")
        )
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser():
    parser = argparse.ArgumentParser(
        description="Robust accelerated universal Ollivier-Ricci flow experiment"
    )
    parser.add_argument("--data_dir", type=str, default=r"D:\100data")
    parser.add_argument("--result_dir", type=str, default=DEFAULT_RESULT_DIR)
    parser.add_argument("--steps", type=int, default=10)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--dt", type=float, default=0.05)

    parser.add_argument(
        "--weight_mode",
        choices=("auto", "unweighted", "distance", "strength"),
        default="auto",
        help=(
            "auto validates each dataset; distance uses positive values directly; "
            "strength converts positive values to 1/value"
        )
    )
    parser.add_argument(
        "--positive_weight_interpretation",
        choices=("distance", "strength", "unweighted"),
        default="distance",
        help="interpretation used by auto mode when every third-column value is positive"
    )
    parser.add_argument(
        "--component_mode",
        choices=("full", "largest"),
        default="full",
        help="full preserves the original non-isolated graph; largest changes the experiment"
    )
    parser.add_argument(
        "--distance_backend",
        choices=("auto", "matrix", "batched"),
        default="auto",
        help="both matrix and batched backends compute exact shortest-path distances"
    )
    parser.add_argument(
        "--max_distance_mb",
        type=float,
        default=256.0,
        help="memory budget for a full or temporary distance matrix"
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=0,
        help="0 chooses min(8, available CPU cores); 1 disables multiprocessing"
    )
    parser.add_argument("--sinkhorn_reg", type=float, default=0.05)
    parser.add_argument("--sinkhorn_max_iter", type=int, default=5000)

    parser.add_argument(
        "--max_transport_cells",
        type=int,
        default=100000000,
        help="skip exact runs exceeding this local OT matrix-cell count per step; 0 disables"
    )
    parser.add_argument(
        "--max_shortest_path_work",
        type=int,
        default=500000000,
        help="skip exact runs when nodes*edges exceeds this value; 0 disables"
    )
    parser.add_argument(
        "--force_large",
        action="store_true",
        help="ignore large-graph safety limits; may take days or exhaust memory"
    )
    parser.add_argument(
        "--inspect_only",
        action="store_true",
        help="only parse and validate datasets; do not run Ricci flow"
    )
    parser.add_argument(
        "--only",
        type=str,
        default="",
        help="run only dataset paths containing this substring"
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="recompute datasets even when a saved per-dataset CSV exists"
    )
    return parser


def main():
    global RESULT_DIR

    parser = build_parser()
    args = parser.parse_args()
    RESULT_DIR = args.result_dir

    if args.steps < 1:
        parser.error("--steps must be at least 1")
    if not (0.0 <= args.alpha <= 1.0):
        parser.error("--alpha must be in [0, 1]")
    if args.dt <= 0.0:
        parser.error("--dt must be positive")
    if args.max_distance_mb <= 0.0:
        parser.error("--max_distance_mb must be positive")
    if args.sinkhorn_reg <= 0.0:
        parser.error("--sinkhorn_reg must be positive")
    if args.sinkhorn_max_iter < 1:
        parser.error("--sinkhorn_max_iter must be positive")

    if args.workers == 0:
        args.workers = max(1, min(8, mp.cpu_count()))
    elif args.workers < 1:
        parser.error("--workers must be zero or a positive integer")

    prepare_dirs()
    setup_logger()

    log(
        "Configuration: data_dir={}, result_dir={}, steps={}, workers={}, "
        "weight_mode={}, component_mode={}\n".format(
            args.data_dir,
            args.result_dir,
            args.steps,
            args.workers,
            args.weight_mode,
            args.component_mode
        )
    )

    run_all(args)


if __name__ == "__main__":
    main()


