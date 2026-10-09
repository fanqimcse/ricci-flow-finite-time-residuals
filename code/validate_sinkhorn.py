"""Validate the batched Sinkhorn against POT before using it in production runs."""

from __future__ import print_function

import os
import sys
import time

import numpy as np
import ot

BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE)

import ricci_flow_program_100data as R
from fast_sinkhorn import batched_sinkhorn_cost

REG = 0.05
MAX_ITER = 5000
STOP_THR = 1e-9


def check_random():
    rng = np.random.default_rng(20261001)
    problems = []
    for _ in range(120):
        n = int(rng.integers(2, 9))
        m = int(rng.integers(2, 9))
        a = rng.random(n) + 0.1
        a /= a.sum()
        b = rng.random(m) + 0.1
        b /= b.sum()
        c = rng.random((n, m)) * float(rng.integers(1, 6))
        problems.append((a, b, c))

    reference = np.array([ot.sinkhorn2(a, b, c, reg=REG, numItermax=MAX_ITER)
                          for a, b, c in problems])
    fast = batched_sinkhorn_cost([p[0] for p in problems],
                                 [p[1] for p in problems],
                                 [p[2] for p in problems],
                                 REG, MAX_ITER, STOP_THR)
    diff = np.abs(reference - fast)
    print("random problems: n={}".format(len(problems)))
    print("  max abs diff = {:.3e}".format(diff.max()))
    print("  max rel diff = {:.3e}".format(
        np.max(diff / np.maximum(np.abs(reference), 1e-12))))
    return diff.max()


def check_graphs(names, steps=10):
    parser = R.build_parser()
    args = parser.parse_args(["--data_dir", r"D:\100data"])
    args.alpha = 0.5
    args.dt = 0.05
    args.sinkhorn_reg = REG
    args.sinkhorn_max_iter = MAX_ITER
    args.weight_mode = "auto"
    args.positive_weight_interpretation = "distance"
    args.component_mode = "full"
    args.distance_backend = "auto"
    args.max_distance_mb = 256.0
    args.workers = 1
    args.inspect_only = False
    args.force_large = False

    files = R.discover_graph_files(args.data_dir)
    selected = []
    for name in names:
        for path in files:
            if R.safe_dataset_id(path, args.data_dir) == name:
                selected.append(path)
                break

    worst = 0.0
    for path in selected:
        dataset_id = R.safe_dataset_id(path, args.data_dir)
        graph = R.load_graph(path, args)
        R.enforce_work_limits(graph, args)
        workspace = R.CurvatureWorkspace(graph, alpha=0.5, backend="auto",
                                         max_distance_mb=256.0)
        t_pot = 0.0
        t_fast = 0.0
        for step in range(steps):
            workspace.prepare_step()
            payloads = [workspace.edge_payload(i)
                        for i in range(len(workspace.edges))]

            start = time.time()
            reference = []
            for _, _, mu, mv, cost, dist in payloads:
                value = float(np.asarray(
                    ot.sinkhorn2(mu, mv, cost, reg=REG, numItermax=MAX_ITER)
                ).reshape(-1)[0])
                reference.append(1.0 - value / dist)
            t_pot += time.time() - start

            start = time.time()
            values = batched_sinkhorn_cost([p[2] for p in payloads],
                                           [p[3] for p in payloads],
                                           [p[4] for p in payloads],
                                           REG, MAX_ITER, STOP_THR)
            fast = [1.0 - v / p[5] for v, p in zip(values, payloads)]
            t_fast += time.time() - start

            diff = np.max(np.abs(np.array(reference) - np.array(fast)))
            worst = max(worst, diff)

            edge_k = dict(zip(workspace.edges, reference))
            R.ricci_flow_step(graph, edge_k, dt=0.05)

        print("  {:<28s} steps={}  max|dkappa|={:.3e}  POT={:.2f}s fast={:.2f}s "
              "speedup={:.1f}x".format(dataset_id, steps, worst, t_pot, t_fast,
                                       t_pot / max(t_fast, 1e-9)))
    return worst


if __name__ == "__main__":
    print("POT version: {}".format(ot.__version__))
    worst_random = check_random()
    print()
    print("real graphs, 10 Ricci-flow steps each:")
    worst_graph = check_graphs([
        "ENZYMES_g1", "ENZYMES_g10", "ENZYMES_g100",
        "aves-weaver-social", "mammalia-raccoon-proximity", "scc_rt_obama"
    ])
    print()
    print("OVERALL: random={:.3e} graphs={:.3e}".format(worst_random, worst_graph))
