#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Batch scheduler for the strict-criterion re-run (paper A, P0 fix).

Why this exists
---------------
The stopping rule inside ``ricci_flow_longrun.py`` compares the per-step move
of the curvature variance with the *initial* variance::

    relative_move = abs(var_k - prev_var) / var0

Any slowly decaying variance therefore eventually satisfies this test, so
"saturated" labels are unreliable (28/28 previously labelled saturated were
still falling in their tail).  Instead of trusting the label we now

  1. run each network as far as the wall-clock budget allows with the
     saturation test effectively disabled (``--tol 1e-12 --patience 100000``),
  2. classify the plateau *after the fact* from the trajectory tail
     (``strict_plateau.py``).

Because a single step costs 0.02-20 s and the corpus would need ~80 h to give
every network 400 steps, the re-run is stratified by cost: the cheap networks
(which are also the small ones) are run to a verified plateau, the expensive
ones keep their old trajectory and contribute an *upper bound* only (the
variance is monotone non-increasing, so the last observed value bounds the
true plateau from above).

Progress is cached per dataset in ``--out_dir/per_dataset`` so an interrupted
batch can simply be relaunched.
"""

from __future__ import print_function

import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd

SOURCE_DIR = os.path.dirname(os.path.abspath(__file__))
MASTER = os.path.join(SOURCE_DIR, "ricci_flow_longrun_results", "analysis",
                      "07_merged_with_structure.csv")
OLD_PER = os.path.join(SOURCE_DIR, "ricci_flow_longrun_results", "per_dataset")
DRIVER = os.path.join(SOURCE_DIR, "ricci_flow_longrun.py")
HIST_COST_FACTOR = 2.0   # measured cost / cost recorded in the original long run
PYTHON = sys.executable


def cost_table(strict_dirs=None):
    """Per-step wall-clock cost + degenerate-case flags.

    The cost recorded during the original long run underestimates the current
    cost by roughly a factor of two: as the weights collapse the entropic OT
    problem becomes ill conditioned and Sinkhorn needs more iterations.  Where
    a strict run already exists we therefore use *its own* measured cost.
    """
    m = pd.read_csv(MASTER)
    measured = {}
    for strict_dir in (strict_dirs or []):
        per = os.path.join(strict_dir, "per_dataset")
        if os.path.isdir(per):
            for fn in os.listdir(per):
                if not fn.endswith(".csv"):
                    continue
                try:
                    tr = pd.read_csv(os.path.join(per, fn))
                except Exception:
                    continue
                if len(tr) >= 5 and "step_time" in tr.columns:
                    measured[fn[:-4]] = float(np.median(
                        tr["step_time"].astype(float).values))
    rows = []
    for _, r in m.iterrows():
        p = os.path.join(OLD_PER, str(r["dataset"]) + ".csv")
        if not os.path.exists(p):
            continue
        try:
            t = pd.read_csv(p)
        except Exception:
            continue
        if len(t) < 5:
            continue
        sps = float(np.median(t["step_time"].astype(float).values))
        var0 = float(t["curvature_variance"].astype(float).values[0])
        rows.append({
            "dataset": r["dataset"],
            "nodes": float(r["nodes"]),
            "sps": sps,
            "sps_est": measured.get(r["dataset"], HIST_COST_FACTOR * sps),
            "old_steps": len(t),
            "old_vr": float(t["var_ratio"].astype(float).values[-1]),
            "var0": var0,
        })
    d = pd.DataFrame(rows)
    # degenerate: zero initial variance -> var_ratio is a 0/0 artefact
    d["degenerate"] = (~np.isfinite(d["var0"])) | (d["var0"] <= 1e-30)
    return d.sort_values("sps").reset_index(drop=True)


def progress_path(out_dir):
    return os.path.join(out_dir, "batch_progress.json")


def load_progress(out_dir):
    p = progress_path(out_dir)
    if not os.path.exists(p):
        return {}
    try:
        with open(p, "r", encoding="utf-8") as fh:
            return json.load(fh)
    except Exception:
        return {}


def already_done(out_dir, dataset, target_steps, need):
    """A dataset is finished if it reached the step target, or if it already
    received at least as much wall-clock as we are about to give it (otherwise
    a cheap-but-slow network whose cost estimate is too low would be re-run
    on every single call)."""
    p = os.path.join(out_dir, "per_dataset", dataset + ".csv")
    if not os.path.exists(p):
        return False
    try:
        n = len(pd.read_csv(p))
    except Exception:
        return False
    if n >= 0.95 * target_steps:
        return True
    rec = load_progress(out_dir).get(dataset)
    if not rec:
        return False
    # A run that stopped on its own (saturation criterion, or any other early
    # exit) used less wall-clock than it was given; a run that was cut off by
    # its budget used (nearly) all of it and deserves a bigger budget instead.
    stopped_early = (rec.get("wall") is not None
                     and rec["wall"] < 0.9 * rec.get("budget", 1e9)
                     and rec.get("steps", 0) < 0.95 * target_steps)
    if rec.get("budget", 0) >= 0.95 * need and stopped_early:
        return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out_dir", default=os.path.join(SOURCE_DIR, "ricci_flow_v2"))
    ap.add_argument("--max_steps", type=int, default=400)
    ap.add_argument("--min_steps", type=int, default=60)
    ap.add_argument("--tol_mode", default="current")
    ap.add_argument("--tol", type=float, default=1e-5)
    ap.add_argument("--patience", type=int, default=50)
    ap.add_argument("--sps_max", type=float, default=1.6,
                    help="only run networks whose step costs at most this many seconds")
    ap.add_argument("--sps_min", type=float, default=0.0)
    ap.add_argument("--budget_cap", type=float, default=700.0)
    ap.add_argument("--retry_mult", type=float, default=1.6,
                    help="budget multiplier applied when the previous attempt "
                         "was cut off by its own budget")
    ap.add_argument("--budget_mult", type=float, default=1.35,
                    help="budget = budget_mult * max_steps * seconds-per-step")
    ap.add_argument("--budget_floor", type=float, default=60.0)
    ap.add_argument("--launch_frac", type=float, default=0.55,
                    help="launch the next item while at least this fraction of "
                         "its budget remains in the current call")
    ap.add_argument("--wall", type=float, default=480.0,
                    help="wall-clock budget for THIS batch call (seconds)")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--tag", default="",
                    help="worker tag forwarded to the driver so that parallel "
                         "workers do not fight over longrun_summary.csv")
    ap.add_argument("--shard_n", type=int, default=1)
    ap.add_argument("--shard_k", type=int, default=0,
                    help="keep only items whose position modulo shard_n "
                         "equals shard_k (lets two workers share one tier)")
    args = ap.parse_args()

    d = cost_table([args.out_dir, os.path.join(SOURCE_DIR, "ricci_flow_strict")])
    sel = d[(~d["degenerate"]) & (d["sps"] <= args.sps_max) & (d["sps"] >= args.sps_min)]
    prog = load_progress(args.out_dir)
    todo = []
    budgets = {}
    for _, r in sel.iterrows():
        need = min(args.budget_cap,
                   max(args.budget_floor, args.budget_mult * args.max_steps * r["sps_est"]))
        if already_done(args.out_dir, r["dataset"], args.max_steps, need):
            continue
        rec = prog.get(r["dataset"])
        if (rec and rec.get("budget", 0) >= 0.95 * need
                and rec.get("wall") is not None
                and rec["wall"] >= 0.9 * rec["budget"]):
            # last attempt ran out of time before reaching the step target
            need = min(args.budget_cap * 1.5, need * args.retry_mult)
        budgets[r["dataset"]] = need
        todo.append(r)
    if args.shard_n > 1:
        todo = [r for i, r in enumerate(todo) if i % args.shard_n == args.shard_k]
    if args.limit:
        todo = todo[:args.limit]

    print("candidates: %d (sps<=%.2f)   still to run: %d" % (len(sel), args.sps_max, len(todo)))
    print("estimated cost: %.0f s" % sum(budgets.values()))
    print("wall budget for this call: %.0f s" % args.wall)
    print("-" * 60)

    env = dict(os.environ)
    env.update({"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                "MKL_NUM_THREADS": "1", "PYTHONUNBUFFERED": "1"})

    started = time.time()
    done, skipped = [], []
    for r in todo:
        name = r["dataset"]
        need = budgets[name]
        left = args.wall - (time.time() - started)
        if left < need * args.launch_frac:
            print("  stop: %.0f s left, next item needs %.0f s" % (left, need))
            break
        cmd = [PYTHON, DRIVER,
               "--only_exact", name,
               "--out_dir", args.out_dir,
               "--max_steps", str(args.max_steps),               "--min_steps", str(args.min_steps),
               "--budget_seconds", "%.1f" % need,
               "--tol_mode", args.tol_mode,
               "--tol", "%g" % args.tol, "--patience", str(args.patience),
               "--lambda_every", "500",
               "--overwrite"]
        if args.tag:
            cmd += ["--worker_tag", args.tag]
        t0 = time.time()
        proc = subprocess.run(cmd, env=env, cwd=SOURCE_DIR,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        dt = time.time() - t0
        n_steps = 0
        p_out = os.path.join(args.out_dir, "per_dataset", name + ".csv")
        if os.path.exists(p_out):
            try:
                n_steps = len(pd.read_csv(p_out))
            except Exception:
                n_steps = 0
        tail = (proc.stdout or b"").decode("utf-8", "replace").strip().splitlines()
        last = [x for x in tail if "->" in x or "DONE" in x]
        print("[%5.0fs/%.0fs] %-42s steps=%-4d %s"
              % (dt, time.time() - started, name, n_steps,
                 (last[-1].strip()[:46] if last else "rc=%d" % proc.returncode)))
        sys.stdout.flush()
        prog[name] = {"budget": need, "steps": n_steps, "wall": round(dt, 1),
                      "rc": proc.returncode}
        with open(progress_path(args.out_dir), "w", encoding="utf-8") as fh:
            json.dump(prog, fh, indent=1)
        if proc.returncode != 0:
            skipped.append(name)
        else:
            done.append(name)

    print("-" * 60)
    print("completed in this call: %d ; failed: %d ; remaining in tier: %d"
          % (len(done), len(skipped),
             len(todo) - len(done) - len(skipped)))
    if skipped:
        print("failed:", skipped)


if __name__ == "__main__":
    main()
