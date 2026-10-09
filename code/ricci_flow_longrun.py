"""Long-horizon discrete Ollivier-Ricci flow: run every dataset to saturation.

This driver reuses the numerics of ``ricci_flow_program_100data.py`` verbatim
(same alpha-Ricci curvature, same Sinkhorn solver, same weight update and
per-step mean-weight renormalisation) and changes only the *control logic*:

1. adaptive horizon - the number of steps is chosen from the measured cost of
   the first step so that each dataset stays inside a wall-clock budget;
2. saturation stop - the run ends when the edge-curvature variance stops
   moving relative to its initial value;
3. graceful divergence - a Ricci update that would overflow ends that dataset
   cleanly and keeps the trajectory recorded up to that point;
4. checkpoints - one CSV per dataset is written immediately, so a long run can
   be resumed and partially inspected.

Nothing about the flow itself is changed, so the first 10 steps reproduce the
original experiment.
"""

from __future__ import print_function

import argparse
import json
import logging
import os
import sys
import time

import numpy as np
import pandas as pd

SOURCE_DIR = os.path.dirname(os.path.abspath(__file__))
if SOURCE_DIR not in sys.path:
    sys.path.insert(0, SOURCE_DIR)

import ricci_flow_program_100data as R
from fast_sinkhorn import batched_sinkhorn_cost


def compute_curvature(workspace, cfg, solver):
    """Curvature for every edge. 'pot' reproduces the original run exactly,
    'batched' is a validated vectorised re-implementation (agreement 8.9e-16)."""
    if solver == "pot":
        return R.compute_curvature(
            workspace,
            workers=1,
            reg=cfg["sinkhorn_reg"],
            max_iter=cfg["sinkhorn_max_iter"]
        )

    distance_time = workspace.prepare_step()
    start = time.time()
    edges = workspace.edges
    mass_u = []
    mass_v = []
    costs = []
    distances = []
    for index in range(len(edges)):
        _, _, mu, mv, cost, edge_distance = workspace.edge_payload(index)
        mass_u.append(mu)
        mass_v.append(mv)
        costs.append(cost)
        distances.append(edge_distance)

    values = batched_sinkhorn_cost(
        mass_u, mass_v, costs,
        cfg["sinkhorn_reg"], cfg["sinkhorn_max_iter"], cfg["stop_thr"]
    )

    curvatures = []
    for value, edge_distance in zip(values, distances):
        curvature = 1.0 - float(value) / float(edge_distance)
        if not np.isfinite(curvature):
            raise RuntimeError("Non-finite curvature")
        if curvature > 1.0 and curvature < 1.0 + 1e-10:
            curvature = 1.0
        curvatures.append(float(curvature))

    sinkhorn_time = time.time() - start
    array = np.asarray(curvatures, dtype=float)
    if array.size == 0 or not np.all(np.isfinite(array)):
        raise RuntimeError("Curvature computation returned no finite values")
    return (float(array.mean()), float(array.var()),
            dict(zip(edges, curvatures)), distance_time, sinkhorn_time)


def build_base_args(argv=None):
    parser = R.build_parser()
    args = parser.parse_args(argv or [])
    # keep the original numerics exactly as they were
    args.alpha = 0.5
    args.dt = 0.05
    args.sinkhorn_reg = 0.05
    args.sinkhorn_max_iter = 5000
    args.weight_mode = "auto"
    args.positive_weight_interpretation = "distance"
    args.component_mode = "full"
    args.distance_backend = "auto"
    args.max_distance_mb = 256.0
    args.workers = 1
    args.inspect_only = False
    args.force_large = False
    return args


def edge_weight_stats(graph):
    weights = np.asarray(
        [float(data.get("weight", 1.0)) for _, _, data in graph.edges(data=True)],
        dtype=float
    )
    if weights.size == 0:
        return 0.0, 0.0, 0.0
    return (
        float(np.min(weights)),
        float(np.median(weights)),
        float(np.max(weights))
    )


def run_one_dataset(filepath, dataset_id, base_args, cfg, out_dir):
    per_dataset_dir = os.path.join(out_dir, "per_dataset")
    if not os.path.isdir(per_dataset_dir):
        os.makedirs(per_dataset_dir)
    target = os.path.join(per_dataset_dir, dataset_id + ".csv")

    graph = None
    try:
        graph = R.load_graph(filepath, base_args)
        estimate = R.enforce_work_limits(graph, base_args)
    except R.DatasetSkipped as exc:
        return {"dataset": dataset_id, "file": filepath, "status": "skipped_large",
                "message": str(exc), "steps": 0, "stop_reason": "work_limit"}
    except Exception as exc:
        return {"dataset": dataset_id, "file": filepath, "status": "failed",
                "message": repr(exc), "steps": 0, "stop_reason": "load_error"}

    workspace = R.CurvatureWorkspace(
        graph,
        alpha=cfg["alpha"],
        backend=base_args.distance_backend,
        max_distance_mb=base_args.max_distance_mb
    )
    lambda_workspace = R.LambdaWorkspace(graph)

    nodes = graph.number_of_nodes()
    edges = graph.number_of_edges()
    components = graph.graph.get("component_count_used", "")

    records = []
    var0 = None
    prev_var = None
    stable = 0
    horizon = cfg["max_steps"]
    started = time.time()
    status = "horizon_reached"
    message = ""

    for step in range(horizon):
        step_start = time.time()
        try:
            mean_k, var_k, edge_k, distance_time, sinkhorn_time = compute_curvature(
                workspace, cfg, cfg["solver"]
            )
        except Exception as exc:
            status = "failed"
            message = "curvature computation: " + repr(exc)
            break

        if step == 0:
            var0 = float(var_k)

        lambda1 = float("nan")
        lambda_time = 0.0
        if step % cfg["lambda_every"] == 0 or step == horizon - 1:
            lambda_start = time.time()
            try:
                lambda1 = lambda_workspace.compute()
            except Exception:
                lambda1 = float("nan")
            lambda_time = time.time() - lambda_start

        min_w, med_w, max_w = edge_weight_stats(graph)
        kappa_abs = np.abs(np.asarray(list(edge_k.values()), dtype=float))
        elapsed = time.time() - step_start

        records.append({
            "step": step,
            "t": step * cfg["dt"],
            "mean_curvature": float(mean_k),
            "curvature_variance": float(var_k),
            "var_ratio": (float(var_k) / var0) if var0 and var0 > 0 else float("nan"),
            "lambda1": lambda1,
            "min_weight": min_w,
            "median_weight": med_w,
            "max_weight": max_w,
            "max_abs_kappa": float(np.max(kappa_abs)) if kappa_abs.size else float("nan"),
            "nodes": nodes,
            "edges": edges,
            "step_time": elapsed,
            "distance_time": distance_time,
            "sinkhorn_time": sinkhorn_time,
            "lambda_time": lambda_time
        })

        if step == 0:
            cost = max(elapsed, 1e-6)
            horizon = min(
                cfg["max_steps"],
                max(cfg["min_steps"], int(cfg["budget_seconds"] / cost))
            )
            if horizon * cost > cfg["budget_seconds"] * cfg["budget_slack"]:
                status = "budget_insufficient"
                message = ("one step costs {:.1f}s, the minimum horizon of {} steps "
                           "would need {:.0f}s > budget {:.0f}s; the original "
                           "10-step result stands".format(
                               cost, horizon, horizon * cost, cfg["budget_seconds"]))
                break

        if step > 0:
            # "initial" (the historical behaviour) divides by the variance at
            # step 0; because that number is fixed while the variance decays,
            # ANY slowly decaying run eventually satisfies the test, which is
            # how 28/28 runs came to be labelled saturated while their tails
            # were still falling.  "current" divides by the variance at the
            # current step, so tol is a bound on the per-step exponential rate:
            # tol = 1e-5 means "the variance moves by <1% over the next 1000
            # steps".
            if cfg.get("tol_mode") == "current":
                denom = max(float(var_k), 1e-300)
            else:
                denom = max(var0, 1e-30) if var0 else 1e-30
            relative_move = abs(float(var_k) - float(prev_var)) / denom
            if relative_move < cfg["tol"]:
                stable += 1
            else:
                stable = 0
            if stable >= cfg["patience"]:
                status = "saturated_current" if cfg.get("tol_mode") == "current" else "saturated"
                break

        if var0 and var0 > 0 and float(var_k) <= cfg["var_floor"] * var0:
            status = "constant_curvature"
            break

        if (step + 1 >= cfg["min_steps"]
                and time.time() - started > cfg["budget_seconds"] * cfg["budget_slack"]):
            status = "budget_exceeded"
            break

        prev_var = float(var_k)

        try:
            R.ricci_flow_step(graph, edge_k, dt=cfg["dt"])
        except Exception as exc:
            status = "diverged"
            message = repr(exc)
            break

    frame = pd.DataFrame(records)
    frame.insert(0, "dataset", dataset_id)
    frame.to_csv(target, index=False)

    summary = {
        "dataset": dataset_id,
        "file": filepath,
        "status": status,
        "message": message,
        "steps": len(records),
        "stop_reason": status,
        "nodes": nodes,
        "edges": edges,
        "components": components,
        "weight_mode": graph.graph.get("weight_mode", ""),
        "weight_reason": graph.graph.get("weight_reason", ""),
        "transport_cells_per_step": graph.graph.get("transport_cells_per_step", ""),
        "horizon_planned": horizon,
        "initial_variance": var0,
        "final_variance": float(frame["curvature_variance"].iloc[-1]) if len(frame) else None,
        "var_ratio_final": float(frame["var_ratio"].iloc[-1]) if len(frame) else None,
        "wall_seconds": round(time.time() - started, 2)
    }
    return summary


def load_prior_status(out_dir):
    """dataset -> stopping reason, collected from every existing summary file.

    Without this, a restarted worker re-labels finished datasets as ``cached``
    and the original stopping reasons are lost.
    """
    import glob
    status = {}
    for path in glob.glob(os.path.join(out_dir, "longrun_summary*.csv")):
        try:
            frame = pd.read_csv(path)
        except Exception:
            continue
        if not {"dataset", "status"} <= set(frame.columns):
            continue
        for _, row in frame.iterrows():
            label = str(row["status"])
            dataset = str(row["dataset"])
            if label == "cached":
                continue
            # a real label always wins over an earlier one
            if dataset not in status or status[dataset] == "cached":
                status[dataset] = label
    return status


def main():
    parser = argparse.ArgumentParser(
        description="Run the discrete Ollivier-Ricci flow to saturation on every dataset"
    )
    parser.add_argument("--data_dir", type=str, default=r"D:\100data")
    parser.add_argument("--out_dir", type=str,
                        default=os.path.join(SOURCE_DIR, "ricci_flow_longrun_results"))
    parser.add_argument("--max_steps", type=int, default=400)
    parser.add_argument("--min_steps", type=int, default=10)
    parser.add_argument("--budget_seconds", type=float, default=900.0)
    parser.add_argument("--budget_slack", type=float, default=1.05)
    parser.add_argument("--tol", type=float, default=1e-4)
    parser.add_argument("--tol_mode", choices=("initial", "current"),
                        default="initial",
                        help="denominator of the saturation test: the variance "
                             "at step 0 (historical, unsafe) or at the current "
                             "step (a bound on the per-step decay rate)")
    parser.add_argument("--patience", type=int, default=5)
    parser.add_argument("--var_floor", type=float, default=1e-6)
    parser.add_argument("--lambda_every", type=int, default=10)
    parser.add_argument("--solver", choices=("batched", "pot"), default="batched")
    parser.add_argument("--only", type=str, default="")
    parser.add_argument("--only_exact", type=str, default="",
                        help="comma separated dataset ids, matched exactly "
                             "(the --only flag is a substring match and would "
                             "make ENZYMES_g11 also select g110/g111)")
    parser.add_argument("--exclude", type=str, default="",
                        help="comma separated tokens; a dataset whose relative "
                             "path contains any of them is skipped")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--reverse", action="store_true",
                        help="walk the dataset list backwards; lets a second "
                             "worker converge on the first one from the tail")
    parser.add_argument("--worker_tag", type=str, default="",
                        help="if set, this worker writes longrun_summary_<tag>.csv "
                             "instead of the shared longrun_summary.csv")
    args = parser.parse_args()

    cfg = {
        "alpha": 0.5,
        "dt": 0.05,
        "sinkhorn_reg": 0.05,
        "sinkhorn_max_iter": 5000,
        "max_steps": args.max_steps,
        "min_steps": args.min_steps,
        "budget_seconds": args.budget_seconds,
        "budget_slack": args.budget_slack,
        "tol": args.tol,
        "tol_mode": args.tol_mode,
        "patience": args.patience,
        "var_floor": args.var_floor,
        "lambda_every": args.lambda_every,
        "stop_thr": 1e-9,
        "solver": args.solver
    }

    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)
    logs_dir = os.path.join(args.out_dir, "logs")
    if not os.path.isdir(logs_dir):
        os.makedirs(logs_dir)
    logging.basicConfig(
        filename=os.path.join(logs_dir, "longrun.log"),
        level=logging.INFO,
        format="%(asctime)s | %(message)s"
    )

    base_args = build_base_args(["--data_dir", args.data_dir])

    files = R.discover_graph_files(args.data_dir)
    if args.only_exact:
        wanted = {t.strip() for t in args.only_exact.split(",") if t.strip()}
        files = [p for p in files
                 if R.safe_dataset_id(p, args.data_dir) in wanted]
    elif args.only:
        token = args.only.lower()
        files = [p for p in files if token in os.path.relpath(p, args.data_dir).lower()]
    if args.exclude:
        tokens = [t.strip().lower() for t in args.exclude.split(",") if t.strip()]
        files = [p for p in files
                 if not any(t in os.path.relpath(p, args.data_dir).lower()
                            for t in tokens)]
    files.sort()
    if args.reverse:
        files = files[::-1]
    if args.limit:
        files = files[:args.limit]

    summary_name = ("longrun_summary.csv" if not args.worker_tag
                    else "longrun_summary_{}.csv".format(args.worker_tag))
    summary_path = os.path.join(args.out_dir, summary_name)

    with open(os.path.join(args.out_dir, "run_config.json"), "w", encoding="utf-8") as handle:
        json.dump({"cfg": cfg, "cli": vars(args), "n_files": len(files),
                   "python": sys.version}, handle, indent=2)

    print("datasets to process: {}".format(len(files)), flush=True)
    logging.info("datasets to process: {}".format(len(files)))

    summaries = []
    prior_status = load_prior_status(args.out_dir)
    for index, filepath in enumerate(files, 1):
        dataset_id = R.safe_dataset_id(filepath, args.data_dir)
        cached = os.path.join(args.out_dir, "per_dataset", dataset_id + ".csv")
        if os.path.exists(cached) and not args.overwrite:
            frame = pd.read_csv(cached)
            # reuse the previously recorded stopping reason: overwriting it with
            # "cached" would destroy the labels of every finished dataset
            prior = prior_status.get(dataset_id, "cached")
            summaries.append({
                "dataset": dataset_id, "file": filepath, "status": prior,
                "message": "reused existing trajectory", "steps": len(frame),
                "stop_reason": "cached",
                "final_variance": float(frame["curvature_variance"].iloc[-1]) if len(frame) else None,
                "var_ratio_final": float(frame["var_ratio"].iloc[-1]) if len(frame) else None
            })
            print("[{}/{}] cached {}".format(index, len(files), dataset_id), flush=True)
            continue

        print("[{}/{}] running {}".format(index, len(files), dataset_id), flush=True)
        summary = run_one_dataset(filepath, dataset_id, base_args, cfg, args.out_dir)
        summaries.append(summary)
        logging.info("{} -> {} steps={} final_ratio={}".format(
            dataset_id, summary["status"], summary.get("steps"),
            summary.get("var_ratio_final")))
        print("    -> {}: steps={} var_ratio={} wall={}s".format(
            summary["status"], summary.get("steps"),
            summary.get("var_ratio_final"), summary.get("wall_seconds")), flush=True)

        pd.DataFrame(summaries).to_csv(summary_path, index=False)

    pd.DataFrame(summaries).to_csv(summary_path, index=False)

    if not args.worker_tag:
        frames = []
        per_dir = os.path.join(args.out_dir, "per_dataset")
        for name in sorted(os.listdir(per_dir)):
            if name.endswith(".csv"):
                frames.append(pd.read_csv(os.path.join(per_dir, name)))
        if frames:
            pd.concat(frames, ignore_index=True).to_csv(
                os.path.join(args.out_dir, "longrun_raw.csv"), index=False)

    print("\nDONE. summary: {}".format(summary_path), flush=True)


if __name__ == "__main__":
    main()
