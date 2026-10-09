# -*- coding: utf-8 -*-
"""Rebuild Supplementary Table S9 with the MATCHED-HORIZON epsilon comparison.

eps_sweep.csv compares V at each run's own stopping step, which is itself a
function of epsilon -- a confound. Here we add the values read at the largest
step index common to all six epsilon values (239 for g105, 400 for g106),
together with the exact-transport re-scoring at the nearest checkpoint.
"""
import csv, io, os, glob, re, statistics as st

D = r"D:\2026codes\2026821100data\exp_round4"
T = os.path.join(D, "traj_eps")
OUT = r"D:\2026codes\2026821100data\paper_a\manuscript\Supplementary_Table_S9_epsilon_trajectory_sweep.csv"

traj = {}
for f in sorted(glob.glob(os.path.join(T, "*.csv"))):
    m = re.match(r"(ENZYMES_g\d+)__reg([\d.]+)\.csv", os.path.basename(f))
    if not m:
        continue
    traj.setdefault(m.group(1), {})[float(m.group(2))] = list(
        csv.DictReader(io.open(f, encoding="utf-8")))

sweep = list(csv.DictReader(io.open(os.path.join(D, "eps_sweep.csv"), encoding="utf-8")))


def at(rows, step):
    for r in rows:
        if int(r["step"]) == step:
            return r
    return None


def last_exact(rows, step):
    """nearest exact-OT checkpoint at or before `step`"""
    best = None
    for r in rows:
        if int(r["step"]) > step:
            break
        if r.get("V_exact"):
            best = r
    return best


matched = {g: min(len(traj[g][r]) - 1 for r in traj[g]) for g in traj}

cols = ["dataset", "epsilon", "n_steps_run", "t_end", "matched_step", "t_matched",
        "checkpoint_step_used", "V0", "V_end_own_horizon", "V_at_matched",
        "V_exact_at_matched", "V_exact_over_V_sinkhorn", "max_dkappa_vs_exact",
        "weight_span_at_matched", "included", "note"]

records = []
for r in sweep:
    g, eps = r["dataset"], float(r["reg"])
    ms = matched[g]
    rr = at(traj[g][eps], ms)
    le = last_exact(traj[g][eps], ms)
    vx = float(le["V_exact"]) if le else float("nan")
    md = float(le["max_dkappa_vs_exact"]) if le and le.get("max_dkappa_vs_exact") else float("nan")
    v = float(rr["V"]) if rr else float("nan")
    incl = "no" if md > 1.0 else "yes"
    note = ""
    if incl == "no":
        note = "excluded: Sinkhorn solve degraded (max|dkappa|>1); V_exact out of line"
    elif vx / v > 2.0:
        note = "Sinkhorn understates V by %.1fx vs exact on the same weights" % (vx / v)
    elif vx / v > 1.1:
        note = "Sinkhorn mildly understates V (%.2fx)" % (vx / v)
    records.append(dict(
        dataset=g, epsilon=eps, n_steps_run=r["steps"], t_end=float(r["T_phys"]),
        matched_step=ms, t_matched=0.05 * ms,
        checkpoint_step_used=(int(le["step"]) if le else ""),
        V0=float(r["V0"]), V_end_own_horizon=float(r["V_end"]),
        V_at_matched=v, V_exact_at_matched=vx, V_exact_over_V_sinkhorn=vx / v,
        max_dkappa_vs_exact=md,
        weight_span_at_matched=(float(rr["span"]) if rr and rr.get("span") else ""),
        included=incl, note=note))

with io.open(OUT, "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    for rec in records:
        w.writerow(rec)

print("wrote", OUT, len(records), "rows\n")
hdr = "%-13s %6s %5s %8s %8s %8s %7s %6s" % ("dataset", "eps", "step", "V", "V_exact", "Vx/V", "maxdk", "incl")
print(hdr)
for rec in records:
    print("%-13s %6g %5d %8.4g %8.4g %8.2f %7.3g %6s" % (
        rec["dataset"], rec["epsilon"], rec["matched_step"], rec["V_at_matched"],
        rec["V_exact_at_matched"], rec["V_exact_over_V_sinkhorn"],
        rec["max_dkappa_vs_exact"], rec["included"]))

print("\n== spread of the INCLUDED runs, at the matched horizon ==")
for g in sorted(traj):
    inc = [r for r in records if r["dataset"] == g and r["included"] == "yes"]
    vs = [r["V_at_matched"] for r in inc]
    vxs = [r["V_exact_at_matched"] for r in inc]
    eps = [r["epsilon"] for r in inc]
    print("%s: eps %g..%g (n=%d, %g-fold)" % (g, min(eps), max(eps), len(inc), max(eps) / min(eps)))
    print("   V (Sinkhorn) : %.4g .. %.4g   spread %.1f%%" %
          (min(vs), max(vs), 100 * (max(vs) - min(vs)) / min(vs)))
    print("   V (exact OT) : %.4g .. %.4g   spread %.1f%%" %
          (min(vxs), max(vxs), 100 * (max(vxs) - min(vxs)) / min(vxs)))
