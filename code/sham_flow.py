"""Sham-flow control: does homogenisation need the curvature-to-edge mapping?

The question this settles
-------------------------
Long-horizon runs show curvature variance dropping to ~20% of its initial
value and then stopping. Two readings are compatible with that observation:

  H1 (structural)  the flow reads real geometry; the residual is what the
                   topology refuses to give up, so it should be predicted by
                   structural covariates (assortativity etc.)
  H2 (mechanical)  the update w <- w*exp(-dt*k) followed by mean-weight
                   renormalisation compresses the weight spread no matter
                   *which* edge gets which curvature. The drop is then a
                   property of the update rule, and the residual is just
                   "where the clock ran out".

The sham flow destroys only the curvature->edge correspondence: each step,
the curvature values are randomly permuted across edges before the weight
update. Topology, the curvature *distribution*, the update rule and the
renormalisation are all untouched. The trajectory that is recorded is always
the TRUE curvature variance, so the two arms are directly comparable.

If the sham arm reaches the same plateau as the real arm, H2 wins and the
"residual is structural" claim is dead. If the sham arm plateaus much higher
(or fails to drop), H1 survives.

Design choices
--------------
* Step count is fixed to the number of steps the real run actually took, so
  the comparison is exactly paired - sham cannot "lose" merely by being
  stopped earlier.
* Permutation is redrawn independently at every step.
* Both arms are recomputed here (mode=real) for a reproducibility check
  against the stored long-run trajectories.

Outputs: sham_results/per_dataset/<dataset>__s<seed>.csv and __real.csv
"""

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

SOURCE_DIR = os.path.dirname(os.path.abspath(__file__))
if SOURCE_DIR not in sys.path:
    sys.path.insert(0, SOURCE_DIR)

import ricci_flow_program_100data as R                      # noqa: E402
from ricci_flow_longrun import (                            # noqa: E402
    compute_curvature, build_base_args, edge_weight_stats)

BASE_ARGS = None
CFG = {
    "alpha": 0.5,
    "dt": 0.05,
    "sinkhorn_reg": 0.05,
    "sinkhorn_max_iter": 5000,
    "stop_thr": 1e-9,
    "solver": "batched",
}


def shuffle_curvature(edge_k, rng):
    keys = list(edge_k.keys())
    values = np.asarray([edge_k[k] for k in keys], dtype=float)
    order = rng.permutation(values.size)
    return {keys[i]: float(values[order[i]]) for i in range(len(keys))}


def run_one(filepath, dataset_id, n_steps, mode, seed, out_dir, base_args):
    per_dir = os.path.join(out_dir, "per_dataset")
    if not os.path.isdir(per_dir):
        os.makedirs(per_dir)
    suffix = "__real.csv" if mode == "real" else "__s{}.csv".format(seed)
    target = os.path.join(per_dir, dataset_id + suffix)

    try:
        graph = R.load_graph(filepath, base_args)
        R.enforce_work_limits(graph, base_args)
    except R.DatasetSkipped as exc:
        return {"dataset": dataset_id, "mode": mode, "seed": seed,
                "status": "skipped_large", "message": str(exc)[:120], "steps": 0}
    except Exception as exc:
        return {"dataset": dataset_id, "mode": mode, "seed": seed,
                "status": "failed", "message": repr(exc)[:120], "steps": 0}

    workspace = R.CurvatureWorkspace(
        graph, alpha=CFG["alpha"], backend=base_args.distance_backend,
        max_distance_mb=base_args.max_distance_mb)
    rng = np.random.default_rng(seed)

    records = []
    var0 = None
    status = "completed"
    message = ""
    started = time.time()

    for step in range(n_steps):
        step_start = time.time()
        try:
            mean_k, var_k, edge_k, distance_time, sinkhorn_time = compute_curvature(
                workspace, CFG, CFG["solver"])
        except Exception as exc:
            status = "failed"
            message = "curvature: " + repr(exc)[:120]
            break

        if step == 0:
            var0 = float(var_k)

        min_w, med_w, max_w = edge_weight_stats(graph)
        records.append({
            "step": step,
            "t": step * CFG["dt"],
            "mean_curvature": float(mean_k),
            "curvature_variance": float(var_k),
            "var_ratio": (float(var_k) / var0) if var0 and var0 > 0 else float("nan"),
            "min_weight": min_w,
            "median_weight": med_w,
            "max_weight": max_w,
            "step_time": time.time() - step_start,
            "distance_time": distance_time,
            "sinkhorn_time": sinkhorn_time,
        })

        if step + 1 >= n_steps:
            break

        try:
            update = edge_k if mode == "real" else shuffle_curvature(edge_k, rng)
            R.ricci_flow_step(graph, update, dt=CFG["dt"])
        except Exception as exc:
            status = "diverged"
            message = repr(exc)[:120]
            break

    if not records:
        return {"dataset": dataset_id, "mode": mode, "seed": seed,
                "status": "failed", "message": "no records", "steps": 0}

    frame = pd.DataFrame(records)
    frame.insert(0, "dataset", dataset_id)
    frame.insert(1, "mode", mode)
    frame.insert(2, "seed", seed)
    frame.to_csv(target, index=False)

    return {
        "dataset": dataset_id, "mode": mode, "seed": seed, "status": status,
        "message": message,
        "steps": len(records),
        "initial_variance": var0,
        "final_variance": float(frame["curvature_variance"].iloc[-1]),
        "var_ratio_final": float(frame["var_ratio"].iloc[-1]),
        "wall_seconds": round(time.time() - started, 2),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir", default=r"D:\100data")
    parser.add_argument("--out_dir", default=os.path.join(
        SOURCE_DIR, "sham_results"))
    parser.add_argument("--reference_dir", default=os.path.join(
        SOURCE_DIR, "ricci_flow_longrun_results", "per_dataset"))
    parser.add_argument("--mode", choices=("sham", "real", "both"),
                        default="sham")
    parser.add_argument("--seeds", default="1,2,3")
    parser.add_argument("--list_file", default="",
                        help="csv with a 'dataset' column; default = all "
                             "datasets that have a reference trajectory")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--max_steps", type=int, default=100000)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    base_args = build_base_args([])
    base_args.data_dir = args.data_dir

    files = R.discover_graph_files(args.data_dir)
    files.sort()
    by_id = {R.safe_dataset_id(p, args.data_dir): p for p in files}

    if args.list_file and os.path.exists(args.list_file):
        wanted = pd.read_csv(args.list_file)["dataset"].astype(str).tolist()
    else:
        wanted = sorted(by_id.keys())

    seeds = [int(s) for s in args.seeds.split(",") if s.strip()]

    jobs = []
    for dataset in wanted:
        if dataset not in by_id:
            continue
        ref = os.path.join(args.reference_dir, dataset + ".csv")
        if not os.path.exists(ref):
            continue
        try:
            n_steps = len(pd.read_csv(ref))
        except Exception:
            continue
        n_steps = min(n_steps, args.max_steps)
        if n_steps < 2:
            continue
        jobs.append((dataset, by_id[dataset], n_steps))
    if args.limit:
        jobs = jobs[:args.limit]

    per_dir = os.path.join(args.out_dir, "per_dataset")
    if not os.path.isdir(per_dir):
        os.makedirs(per_dir)

    print("sham-flow control: {} datasets, mode={}, seeds={}".format(
        len(jobs), args.mode, seeds), flush=True)

    summaries = []
    for index, (dataset, path, n_steps) in enumerate(jobs, 1):
        modes = []
        if args.mode in ("real", "both"):
            modes.append(("real", 0))
        if args.mode in ("sham", "both"):
            modes.extend(("sham", s) for s in seeds)

        for mode, seed in modes:
            suffix = "__real.csv" if mode == "real" else "__s{}.csv".format(seed)
            target = os.path.join(per_dir, dataset + suffix)
            if os.path.exists(target) and not args.overwrite:
                continue
            summary = run_one(path, dataset, n_steps, mode, seed,
                              args.out_dir, base_args)
            summaries.append(summary)
            print("[{}/{}] {:<46s} {:>5s} seed={} steps={} ratio={} ({:.1f}s)".format(
                index, len(jobs), dataset[:46], mode, seed, summary["steps"],
                ("{:.4f}".format(summary["var_ratio_final"])
                 if summary.get("var_ratio_final") is not None else "NA"),
                summary.get("wall_seconds", 0)), flush=True)

    if summaries:
        frame = pd.DataFrame(summaries)
        out_csv = os.path.join(args.out_dir, "sham_summary.csv")
        if os.path.exists(out_csv) and not args.overwrite:
            try:
                old = pd.read_csv(out_csv)
                frame = pd.concat([old, frame], ignore_index=True)
                frame = frame.drop_duplicates(
                    subset=["dataset", "mode", "seed"], keep="last")
            except Exception:
                pass
        frame.to_csv(out_csv, index=False)
        print("\nsummary written: {}".format(out_csv), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
