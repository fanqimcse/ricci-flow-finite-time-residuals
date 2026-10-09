# -*- coding: utf-8 -*-
"""Experiment D: the definitive reading, using every horizon available.

Sources
-------
  clustering_merged.csv          first pass, unsafe rule   (V at index steps-1)
  clust_strict_b*.csv            corrected rule re-runs    (V at index steps-1)
  ckpt_base*_t*/hist.csv         per-step history of the chunked continuations
                                 (clust_continue.py), which reproduces the
                                 original trajectory bit-for-bit (verified:
                                 V(249) = 0.047023126080126196 in both)

All reported var_ratio values sit at index (steps - 1); the history gives the
index explicitly, so the two are merged on the index.

What is computed
----------------
  * V at fixed horizons (125, 250, 400, 700, 900 steps)
  * the exponential rate over consecutive windows between those horizons
  * the terminal state (constant_curvature / saturated / horizon_reached)
"""

import os
import glob
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SYN = os.path.join(HERE, "synthetic_results")
OUT_TXT = os.path.join(SYN, "CLUSTERING_STRICT_REPORT.txt")
OUT_CSV = os.path.join(SYN, "clustering_strict_merged.csv")
HORIZONS = [125, 250, 400, 700, 900]

lines = []


def note(x=""):
    print(x)
    lines.append(str(x))


def collect():
    """cell key -> {index: V}, plus terminal info."""
    cells = {}

    def put(base, level, idx, v, term=None):
        key = (int(base), round(float(level), 3))
        c = cells.setdefault(key, {"V": {}, "term": None, "steps": 0})
        c["V"][int(idx)] = float(v)
        if term:
            c["term"] = term
            c["steps"] = max(c["steps"], int(term["steps"]))

    # 1. per-step histories (authoritative when present)
    for d in sorted(glob.glob(os.path.join(SYN, "ckpt_base*_t*"))):
        h = os.path.join(d, "hist.csv")
        meta = os.path.join(d, "meta.json")
        if not os.path.exists(h):
            continue
        name = os.path.basename(d)
        base = int(name.split("base")[1].split("_")[0])
        level = int(name.split("_t")[1]) / 100.0
        df = pd.read_csv(h)
        for _, r in df.iterrows():
            put(base, level, int(r["step"]), float(r["var_ratio"]))
        if os.path.exists(meta):
            import json
            m = json.load(open(meta, encoding="utf-8"))
            put(base, level, int(m["step"]) - 1, np.nan,
                {"status": m["status"], "steps": int(m["step"])})

    # 2. single-process result CSVs (V sits at index steps-1)
    def add_csv(path):
        try:
            df = pd.read_csv(path)
        except Exception:
            return
        for _, r in df.iterrows():
            v = float(r["var_ratio"])
            if not np.isfinite(v) or v <= 0:
                continue
            st = str(r["status"])
            put(r["base"], r["target_T"], int(r["steps"]) - 1, v,
                None if st == "running" else {"status": st,
                                              "steps": int(r["steps"])})

    p = os.path.join(SYN, "clustering_merged.csv")
    if os.path.exists(p):
        add_csv(p)
    for f in sorted(glob.glob(os.path.join(SYN, "clust_strict_b*.csv"))):
        add_csv(f)
    # drop the nan placeholder written alongside the terminal meta
    for c in cells.values():
        c["V"] = {k: v for k, v in c["V"].items() if np.isfinite(v) and v > 0}
    return cells


def v_at(c, horizon):
    """V at index horizon-1 (largest available index <= horizon-1)."""
    idx = horizon - 1
    if idx in c["V"]:
        return c["V"][idx], idx
    cand = [i for i in c["V"] if i <= idx]
    if not cand:
        return np.nan, None
    i = max(cand)
    return c["V"][i], i


def main():
    cells = collect()
    note("EXPERIMENT D: TRANSITIVITY, ALL HORIZONS")
    note("=" * 72)

    rows = []
    for (b, lvl) in sorted(cells):
        c = cells[(b, lvl)]
        if not c["V"]:
            continue
        vals = {h: v_at(c, h) for h in HORIZONS}
        windows = []
        hs = [h for h in HORIZONS if vals[h][1] is not None]
        for a, z in zip(hs[:-1], hs[1:]):
            v1, i1 = vals[a]
            v2, i2 = vals[z]
            if i2 > i1:
                windows.append("%d-%d:%+.2e" % (i1 + 1, i2 + 1,
                                                np.log(v2 / v1) / (i2 - i1)))
        rows.append({
            "base": b, "target_T": lvl,
            "V_125": vals[125][0], "V_250": vals[250][0], "V_400": vals[400][0],
            "V_700": vals[700][0], "V_900": vals[900][0],
            "steps_terminal": c["steps"],
            "status": (c["term"] or {}).get("status", ""),
            "windows": "; ".join(windows),
        })
    res = pd.DataFrame(rows)
    res.to_csv(OUT_CSV, index=False)

    note("")
    note("--- 1. residual at matched horizons (V = sigma^2(t)/sigma^2(0)) ---")
    note("  base  level    V(125)     V(250)     V(400)     V(700)     V(900)   terminal")
    for _, r in res.iterrows():
        f = lambda x: ("%9.3e" % x) if np.isfinite(x) else "        -"
        note("  %4d  %.2f  %s %s %s %s %s   %s @%d"
             % (r["base"], r["target_T"], f(r["V_125"]), f(r["V_250"]),
                f(r["V_400"]), f(r["V_700"]), f(r["V_900"]),
                r["status"], r["steps_terminal"]))

    note("")
    note("--- 2. per-window exponential rate (per step) ---")
    for _, r in res.iterrows():
        note("  base %d T~%.2f: %s" % (r["base"], r["target_T"], r["windows"]))

    note("")
    note("--- 3. the fold at matched horizons ---")
    for b, g in res.groupby("base"):
        g = g.sort_values("target_T")
        if len(g) < 2:
            continue
        lo, hi = g.iloc[0], g.iloc[-1]
        note("  base %d (T %.2f vs %.2f):" % (b, lo["target_T"], hi["target_T"]))
        for h in HORIZONS:
            a, c_ = lo["V_%d" % h], hi["V_%d" % h]
            if np.isfinite(a) and np.isfinite(c_) and a > 0:
                note("     at %3d steps: %.3e vs %.3e  fold = %.3gx"
                     % (h, a, c_, c_ / a))

    note("")
    note("--- 4. do the arms reach the SAME terminal state? ---")
    for b, g in res.groupby("base"):
        g = g.sort_values("target_T")
        if len(g) < 2:
            continue
        lo, hi = g.iloc[0], g.iloc[-1]
        note("  base %d: low-T %s @%d ; high-T %s @%d"
             % (b, lo["status"], lo["steps_terminal"],
                hi["status"], hi["steps_terminal"]))
        if "constant_curvature" in str(lo["status"]) and \
           "constant_curvature" in str(hi["status"]):
            note("     -> BOTH reach the constant-curvature floor: transitivity")
            note("        delays uniformisation by %.1fx in steps, it does not"
                 % (hi["steps_terminal"] / lo["steps_terminal"]))
            note("        prevent it. Any fold measured before both arms land is")
            note("        a TIME-SCALE quantity, not a plateau ratio.")

    note("")
    note("--- 5. what this licenses ---")
    note("  LICENSED: (a) at every matched horizon up to 400 steps the high-T")
    note("            residual exceeds the low-T one by ~10^3; (b) on base 0")
    note("            both arms reach the constant-curvature floor; (c) the")
    note("            delay in steps is a factor ~2.3 (381 vs 884).")
    note("  NOT LICENSED: a plateau ratio (it is 1 on the one graph where both")
    note("            arms finished), or 'transitivity creates the stalled branch'.")
    note("  WITHDRAWN: the first-pass 25-37x plateau reading AND our own")
    note("            intermediate reading ('moves the graph towards the stalled")
    note("            branch') -- both were horizon artefacts.")

    with open(OUT_TXT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\nwrote %s" % OUT_TXT)


if __name__ == "__main__":
    main()
