"""Diagnose why surgery never triggers.

Checks two suspects:
  1. nx.bridges raising (silently swallowed -> zero removals forever)
  2. the weight attribute not being read correctly
"""

import os
import sys
import traceback

import networkx as nx
import numpy as np

SOURCE_DIR = os.path.dirname(os.path.abspath(__file__))
if SOURCE_DIR not in sys.path:
    sys.path.insert(0, SOURCE_DIR)

import ricci_flow_program_100data as R                    # noqa: E402
from ricci_flow_longrun import compute_curvature, build_base_args  # noqa: E402
from surgery_check import CFG, weights_of                 # noqa: E402

OUT = os.path.join(SOURCE_DIR, "_surg_diag.txt")
DATA_DIR = r"D:\100data"


def main():
    lines = []

    # 1. does nx.bridges survive a disconnected graph?
    g = nx.Graph()
    g.add_edges_from([(1, 2), (2, 3), (3, 1), (4, 5)])
    try:
        lines.append("nx.bridges(disconnected) -> {}".format(
            sorted(map(sorted, nx.bridges(g)))))
    except Exception as exc:
        lines.append("nx.bridges(disconnected) RAISED {}: {}".format(
            type(exc).__name__, exc))

    # 2. real flow: what does surgery_pass actually see?
    base = build_base_args(["--data_dir", DATA_DIR])
    files = R.discover_graph_files(DATA_DIR)
    for name in ["reptilia-tortoise-network-pv", "ENZYMES_g103"]:
        try:
            path = [f for f in files
                    if R.safe_dataset_id(f, DATA_DIR) == name][0]
        except IndexError:
            lines.append("{} not found".format(name))
            continue
        graph = R.load_graph(path, base)
        lines.append("{}: nodes={} edges={} connected_components={}".format(
            name, graph.number_of_nodes(), graph.number_of_edges(),
            nx.number_connected_components(graph)))
        workspace = R.CurvatureWorkspace(
            graph, alpha=CFG["alpha"], backend=base.distance_backend,
            max_distance_mb=base.max_distance_mb)
        for step in range(30):
            mean_k, var_k, edge_k, _, _ = compute_curvature(
                workspace, CFG, CFG["solver"])
            R.ricci_flow_step(graph, edge_k, dt=CFG["dt"])
        w = weights_of(graph)
        lines.append("{}: after 30 steps weight min={:.4f} mean={:.4f} "
                     "max={:.4f}; threshold(0.05)={:.4f}; would_cut={}".format(
                         name, w.min(), w.mean(), w.max(), 0.05 * w.mean(),
                         int((w < 0.05 * w.mean()).sum())))
        try:
            b = list(nx.bridges(graph))
            lines.append("{}: nx.bridges ok, {} bridges".format(name, len(b)))
        except Exception as exc:
            lines.append("{}: nx.bridges RAISED {}: {}".format(
                name, type(exc).__name__, exc))

    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    try:
        main()
    except Exception:
        with open(OUT, "w", encoding="utf-8") as fh:
            fh.write(traceback.format_exc())
