#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Post-hoc plateau classification (replaces the unreliable stopping label).

The stopping rule used while generating the trajectories compares the per-step
move of the curvature variance with the *initial* variance, so a slowly
decaying run is eventually declared "saturated".  Here the plateau is decided
from the trajectory tail instead, and three numbers are reported per network:

  ``v_end``     last observed normalised variance  -- an UPPER bound on the
                plateau, valid for every run because the variance is monotone
                non-increasing in this flow;
  ``v_inf_fit`` asymptote of a three-parameter fit
                V(t) = Vinf + (V0 - Vinf) exp(-r t)  -- an ESTIMATE;
  ``plateaued`` True only if the run is long enough (>= ``--min_steps``) and
                the variance moved by less than ``--tol`` (relative) across the
                last ``--tail`` fraction of the trajectory.

Networks that never plateau contribute an upper bound only; they are reported
separately so that no headline number is computed by mixing the two.
"""

from __future__ import print_function

import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit

SOURCE_DIR = os.path.dirname(os.path.abspath(__file__))
MASTER = os.path.join(SOURCE_DIR, "ricci_flow_longrun_results", "analysis",
                      "07_merged_with_structure.csv")


def fam_of(d):
    for f in ("scc_rt", "scc2034", "reptilia", "mammalia", "aves", "insecta",
              "bn", "ENZYMES", "bio", "ia", "power", "eco", "social"):
        if d.startswith(f):
            return f
    return d.split("_")[0]


def three_param(t, v_inf, v0, r):
    return v_inf + (v0 - v_inf) * np.exp(-r * t)


def analyse(tr, tail=0.30, tol=0.02, min_steps=150):
    """Return a dict of plateau statistics for one trajectory frame."""
    step = tr["step"].astype(float).values
    v = tr["curvature_variance"].astype(float).values
    v0 = v[0]
    if not np.isfinite(v0) or v0 <= 1e-30:
        return None
    vn = v / v0
    n = len(step)
    k = max(1, int(np.floor(n * (1.0 - tail))))
    tail_net = float(vn[-1] / vn[k] - 1.0) if vn[k] > 0 else np.nan

    win = slice(k, n)
    if n - k >= 10:
        slope = stats.linregress(step[win], np.log(np.clip(vn[win], 1e-300, None))).slope
    else:
        slope = np.nan
    extrap = float(np.exp(slope * (n - k)) - 1.0) if np.isfinite(slope) else np.nan

    v_inf = np.nan
    try:
        popt, _ = curve_fit(three_param, step, vn,
                            p0=[max(vn[-1] * 0.5, 1e-6), 1.0, 0.05],
                            bounds=([0.0, 0.5, 1e-6], [1.0, 1.5, 5.0]),
                            maxfev=20000)
        v_inf = float(popt[0])
    except Exception:
        pass

    plateaued = bool(n >= min_steps and abs(tail_net) < tol)
    return {
        "steps": n,
        "v_end": float(vn[-1]),
        "v_inf_fit": v_inf,
        "tail_net": tail_net,
        "tail_slope": float(slope) if np.isfinite(slope) else np.nan,
        "tail_extrap": extrap,
        "plateaued": plateaued,
        "mean_k_end": float(tr["mean_curvature"].astype(float).values[-1]),
        "min_w_end": float(tr["min_weight"].astype(float).values[-1]),
        "max_w_end": float(tr["max_weight"].astype(float).values[-1]),
    }


def load_best(name, dirs):
    """Return (frame, dir) for the *longest* available trajectory of a network.

    The same network may exist in several output directories (a 400-step
    screening run and a much longer "deep" run); the longest one carries the
    most information about the tail, so it wins.
    """
    best, src = None, None
    for base in dirs:
        p = os.path.join(base, "per_dataset", name + ".csv")
        if not os.path.exists(p):
            continue
        try:
            tr = pd.read_csv(p)
        except Exception:
            continue
        if len(tr) < 5:
            continue
        if best is None or len(tr) > len(best):
            best, src = tr, base
    return best, src


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--strict_dir", nargs="*",
                    default=[os.path.join(SOURCE_DIR, "ricci_flow_v2"),
                             os.path.join(SOURCE_DIR, "ricci_flow_deep"),
                             os.path.join(SOURCE_DIR, "ricci_flow_strict")])
    ap.add_argument("--old_dir",
                    default=os.path.join(SOURCE_DIR, "ricci_flow_longrun_results"))
    ap.add_argument("--out_dir", default=os.path.join(SOURCE_DIR, "strict_results"))
    ap.add_argument("--tail", type=float, default=0.30)
    ap.add_argument("--tol", type=float, default=0.02)
    ap.add_argument("--min_steps", type=int, default=150)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    master = pd.read_csv(MASTER)
    meta = {r["dataset"]: r for _, r in master.iterrows()}

    rows = []
    for name in sorted(meta):
        rec = {"dataset": name, "family": fam_of(name), "source": "old",
               "nodes": meta[name].get("nodes"), "edges": meta[name].get("edges")}
        src_dirs = list(args.strict_dir)
        tr, src = load_best(name, src_dirs)
        if tr is None:
            tr, src = load_best(name, [args.old_dir])
            if tr is None:
                continue
            src = "old"
        st = analyse(tr, args.tail, args.tol, args.min_steps)
        if st is None:
            continue
        rec.update(st)
        rec["source"] = "old" if src == "old" else os.path.basename(os.path.dirname(src))
        rec["file_steps"] = len(tr)
        if "steps" in rec:
            rows.append(rec)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(args.out_dir, "plateau_table.csv"), index=False)

    lines = []
    w = lines.append
    w("PLATEAU CLASSIFICATION (post-hoc, strict tail test)")
    w("=" * 68)
    w("criterion: run length >= %d steps AND |relative move over the last "
      "%.0f%%| < %.1f%%" % (args.min_steps, args.tail * 100, args.tol * 100))
    w("")
    w("networks analysed: %d (from strict re-run: %d ; from the original "
      "long run: %d)" % (len(df), (df["source"] == "strict").sum(),
                         (df["source"] == "old").sum()))
    w("")

    st = df[~df["source"].isin(["old"])]
    if len(st):
        w("--- STRICT RE-RUN SUBSET (n=%d) ---" % len(st))
        pl = st[st["plateaued"]]
        w("plateaued: %d / %d   (families: %s)"
          % (len(pl), len(st), sorted(set(pl["family"])) if len(pl) else "-"))
        if len(pl):
            w("  median v_end (plateaued)      = %.4f" % pl["v_end"].median())
            w("  IQR                           = %.4f .. %.4f"
              % (pl["v_end"].quantile(0.25), pl["v_end"].quantile(0.75)))
            w("  median fitted asymptote Vinf  = %.4f" % pl["v_inf_fit"].median())
            w("  median |tail move|            = %.4f%%" % (100 * pl["tail_net"].abs().median()))
        nf = st[~st["plateaued"]]
        if len(nf):
            w("still falling: %d" % len(nf))
            w("  median v_end (UPPER bound)    = %.4f" % nf["v_end"].median())
            w("  median tail move              = %+.1f%%" % (100 * nf["tail_net"].median()))
            w("  median fitted asymptote Vinf  = %.4f" % nf["v_inf_fit"].median())
        w("")
        w("  truncation bias check (same networks, old vs strict end value):")
        both = []
        for name in st["dataset"]:
            p_old = os.path.join(args.old_dir, "per_dataset", name + ".csv")
            if not os.path.exists(p_old):
                continue
            vo = pd.read_csv(p_old)["curvature_variance"].astype(float).values
            vs_row = st[st["dataset"] == name].iloc[0]
            if vo[0] > 0:
                both.append((name, vo[-1] / vo[0], vs_row["v_end"],
                             vs_row["v_end"] / (vo[-1] / vo[0])))
        if both:
            b = pd.DataFrame(both, columns=["dataset", "old", "strict", "ratio"])
            w("    n=%d ; old median %.4f -> strict median %.4f ; median ratio "
              "%.3f" % (len(b), b["old"].median(), b["strict"].median(),
                        b["ratio"].median()))
        w("")

    w("--- ALL NETWORKS (v_end is an upper bound on the plateau) ---")
    w("median v_end over all %d networks = %.4f" % (len(df), df["v_end"].median()))
    conv = df[df["dataset"].isin(master.loc[master["converged"] == True, "dataset"])] \
        if "converged" in master.columns else df
    w("median v_end over the previously 'converged' subset (n=%d) = %.4f"
      % (len(conv), conv["v_end"].median()))
    w("")
    w("--- weight dynamics at the end of the strict runs ---")
    if len(st):
        w("min_weight < 1e-3 in %d / %d strict runs"
          % ((st["min_w_end"] < 1e-3).sum(), len(st)))
        w("median weight span max/min at the end: %.3g"
          % float(np.median(st["max_w_end"] / np.maximum(st["min_w_end"], 1e-300))))
    w("")
    w("--- per-network table ---")
    cols = ["dataset", "family", "source", "steps", "v_end", "v_inf_fit",
            "tail_net", "plateaued"]
    w(df[cols].sort_values(["plateaued", "v_end"], ascending=[False, False])
      .to_string(index=False, float_format=lambda x: "%.4f" % x))

    out = os.path.join(args.out_dir, "PLATEAU_REPORT.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    print("\n".join(lines[:40]))
    print("\nwrote", out)


if __name__ == "__main__":
    main()
