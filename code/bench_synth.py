"""Quick throughput check for the synthetic sweep.

Answers one question: how long does one Ricci-flow step take on a small
configuration-model graph, so we can tell whether the running sweep is
progressing normally or hung.
"""

import os
import sys
import time

import numpy as np
import networkx as nx

SOURCE_DIR = os.path.dirname(os.path.abspath(__file__))
if SOURCE_DIR not in sys.path:
    sys.path.insert(0, SOURCE_DIR)

from synthetic_assortativity import base_graph, run_flow  # noqa: E402


def main():
    for n, mean_degree, steps in [(40, 5, 6), (80, 6, 6)]:
        graph = base_graph(n, mean_degree, seed=555)
        t0 = time.time()
        res = run_flow(graph.copy(), max_steps=steps)
        dt = time.time() - t0
        print("n={} mean_deg={} edges={} | {} steps in {:.1f}s "
              "-> {:.2f} s/step | status={} ratio={:.4f}".format(
                  graph.number_of_nodes(), mean_degree,
                  graph.number_of_edges(), steps, dt, dt / max(steps, 1),
                  res["status"], res["var_ratio"]), flush=True)


if __name__ == "__main__":
    main()
