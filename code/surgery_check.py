"""Surgery check: does the curvature plateau survive if vanishing edges are removed?

Why this is decisive
--------------------
The literature (Ma-Yang 2025 piecewise-linear flow + surgery, Naama et al. 2024
dRfge) proves convergence to *constant* curvature for flows that DELETE edges
whose weight collapses. Our plateau was measured with a flow that never deletes
anything. A referee will ask, immediately:

    "Isn't your residual simply the cost of not doing surgery?"

Two possible answers, both publishable but very different:

  A) residual survives surgery  -> our obstruction is NOT the surgery artefact;
     genuinely new, sits against the existing convergence literature.
  B) residual vanishes with surgery -> the residual measures *how far the graph
     is from being uniformisable without changing its topology*. Still a real
     invariant, but the claim must be rewritten that way.

Design
------
Same numerics as `ricci_flow_longrun.py`. After each update, delete every edge
whose weight falls below `tol` x mean weight, EXCEPT bridges (so connectivity
and component structure are preserved). Rebuild the workspace when the edge set
changes. Record the variance trajectory so the no-surgery and surgery arms can
be compared step by step.
"""

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

import ricci_flow_program_100data as R                    # noqa: E402
from ricci_flow_longrun import compute_curvature, build_base_args  # noqa: E402

CFG = {
    "alpha": 0.5,
    "dt": 0.05,
    "sinkhorn_reg": 0.05,
    "sinkhorn_max_iter": 5000,
    "stop_thr": 1e-9,
    "solver": "batched",
}


def weights_of(graph):
    return np.asarray([float(d.get("weight", 1.0))
                       for _, _, d in graph.edges(data=True)], dtype=float)


def surgery_pass(graph, tol):
    """Delete low-weight non-bridge edges. Returns (removed_count, kept_graph)."""
    w = weights_of(graph)
    if w.size == 0:
        return 0
    threshold = tol * float(np.mean(w))
    doomed = [(u, v) for (u, v), weight in zip(list(graph.edges()), w)
              if weight < threshold]
    if not doomed:
        return 0
    # never disconnect the graph: bridges are protected.
    # nx.bridges is undefined on trivial graphs and can fail on exotic inputs,
    # so fall back to "do not cut anything" rather than crashing the run.
    try:
        bridges = set(map(frozenset, nx.bridges(graph)))
    except Exception:
        return 0
    kept = [(u, v) for (u, v) in doomed if frozenset((u, v)) not in bridges]
    if not kept:
        return 0
    graph.remove_edges_from(kept)
    return len(kept)


def run_one(filepath, dataset_id, base_args, tol, max_steps, budget_s):
    graph = R.load_graph(filepath, base_args)
    try:
        R.enforce_work_limits(graph, base_args)
    except Exception:
        return None

    workspace = R.CurvatureWorkspace(
        graph, alpha=CFG["alpha"], backend=base_args.distance_backend,
        max_distance_mb=base_args.max_distance_mb)

    records = []
    var0 = None
    prev = None
    stable = 0
    status = "horizon_reached"
    removed_total = 0
    started = time.time()

    for step in range(max_steps):
        try:
            mean_k, var_k, edge_k, _, _ = compute_curvature(
                workspace, CFG, CFG["solver"])
        except Exception as exc:
            status = "failed"
            records.append({"step": step, "var": np.nan,
                            "message": repr(exc)[:90]})
            break

        if step == 0:
            var0 = float(var_k)
        ratio = float(var_k) / var0 if var0 and var0 > 0 else np.nan
        records.append({
            "step": step, "var": float(var_k), "ratio": ratio,
            "edges": graph.number_of_edges(), "removed": removed_total,
            "mean_k": float(mean_k),
        })

        if prev is not None and var0 and var0 > 0:
            move = abs(float(var_k) - prev) / var0
            stable = stable + 1 if move < 1e-4 else 0
            if stable >= 5:
                status = "saturated"
                break
        else:
            # first step, or the step right after an edge removal: the variance
            # level is not comparable with the previous one.
            stable = 0
        if var0 and var0 > 0 and float(var_k) <= 1e-6 * var0:
            status = "constant_curvature"
            break
        prev = float(var_k)

        try:
            R.ricci_flow_step(graph, edge_k, dt=CFG["dt"])
        except Exception as exc:
            status = "diverged"
            records.append({"step": step + 1, "var": np.nan,
                            "message": repr(exc)[:90]})
            break

        try:
            n_removed = surgery_pass(graph, tol)
        except Exception as exc:
            print("  surgery skipped on step {}: {}".format(step, repr(exc)[:60]),
                  flush=True)
            n_removed = 0
        if n_removed:
            removed_total += n_removed
            workspace = R.CurvatureWorkspace(
                graph, alpha=CFG["alpha"],
                backend=base_args.distance_backend,
                max_distance_mb=base_args.max_distance_mb)
            prev = None          # variance is not comparable across edge sets
            stable = 0

        if time.time() - started > budget_s:
            status = "budget_exceeded"
            break

    frame = pd.DataFrame(records)
    frame["dataset"] = dataset_id
    frame["surgery_tol"] = tol
    frame["status"] = status
    return frame


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir", default=r"D:\100data")
    parser.add_argument("--out_dir", default=os.path.join(
        SOURCE_DIR, "surgery_results"))
    parser.add_argument("--datasets", default="",
                        help="comma separated dataset ids (default: all)")
    parser.add_argument("--tol", type=float, default=0.05)
    parser.add_argument("--max_steps", type=int, default=300)
    parser.add_argument("--budget", type=float, default=1800.0)
    args = parser.parse_args()

    os.makedirs(os.path.join(args.out_dir, "per_dataset"), exist_ok=True)
    base_args = build_base_args(["--data_dir", args.data_dir])
    files = R.discover_graph_files(args.data_dir)

    wanted = [s for s in args.datasets.split(",") if s.strip()]
    if wanted:
        files = [p for p in files
                 if R.safe_dataset_id(p, args.data_dir) in wanted]
    else:
        list_file = os.path.join(SOURCE_DIR, "converged_clean_list.csv")
        if os.path.exists(list_file):
            ids = pd.read_csv(list_file)
            col = "dataset" if "dataset" in ids.columns else ids.columns[0]
            keep = set(str(v) for v in ids[col])
            files = [p for p in files
                     if R.safe_dataset_id(p, args.data_dir) in keep]

    print("surgery check: {} datasets, tol={}".format(len(files), args.tol),
          flush=True)
    log_path = os.path.join(args.out_dir, "progress.log")
    log_fh = open(log_path, "a", encoding="utf-8")

    def log(msg):
        print(msg, flush=True)
        log_fh.write(msg + "\n")
        log_fh.flush()

    log("=== run start tol={} max_steps={} budget={} ===".format(
        args.tol, args.max_steps, args.budget))
    summary = []
    for index, filepath in enumerate(files, 1):
        dataset_id = R.safe_dataset_id(filepath, args.data_dir)
        target = os.path.join(args.out_dir, "per_dataset", dataset_id + ".csv")
        if os.path.exists(target):
            log("[{}/{}] cached {}".format(index, len(files), dataset_id))
            continue
        log("[{}/{}] starting {}".format(index, len(files), dataset_id))
        t0 = time.time()
        try:
            frame = run_one(filepath, dataset_id, base_args, args.tol,
                            args.max_steps, args.budget)
        except Exception as exc:
            log("[{}/{}] ERROR {} {}".format(index, len(files), dataset_id,
                                             repr(exc)[:80]))
            continue
        if frame is None or frame.empty:
            log("[{}/{}] skipped {}".format(index, len(files), dataset_id))
            continue
        frame.to_csv(target, index=False)
        final = frame["ratio"].dropna()
        summary.append({
            "dataset": dataset_id,
            "status": frame["status"].iloc[-1],
            "steps": int(frame["step"].max()) + 1,
            "final_ratio": float(final.iloc[-1]) if len(final) else np.nan,
            "min_ratio": float(final.min()) if len(final) else np.nan,
            "edges_removed": int(frame["removed"].max()),
            "wall": round(time.time() - t0, 1),
        })
        log("[{}/{}] {} -> status={} steps={} final_ratio={:.5f} removed={} "
            "({:.0f}s)".format(
                index, len(files), dataset_id,
                summary[-1]["status"], summary[-1]["steps"],
                summary[-1]["final_ratio"], summary[-1]["edges_removed"],
                summary[-1]["wall"]))
        pd.DataFrame(summary).to_csv(
            os.path.join(args.out_dir, "surgery_summary.csv"), index=False)

    pd.DataFrame(summary).to_csv(
        os.path.join(args.out_dir, "surgery_summary.csv"), index=False)
    log("DONE -> {}".format(os.path.join(args.out_dir,
                                         "surgery_summary.csv")))
    log_fh.close()


if __name__ == "__main__":
    main()
