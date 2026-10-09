# -*- coding: utf-8 -*-
"""Analyse experiment C: community structure and the 92x magnitude gap.

Questions, fixed in advance
---------------------------
1. Manipulation check: did rewiring hit the target cross-group fraction, and
   does planted modularity Q rise as mu falls?
2. Does the plateau RISE as community structure strengthens (mu falls, Q rises)?
3. Within base graphs (degree sequence fixed), is the effect consistent?
4. With base fixed effects, do modularity and assortativity survive together,
   or does one absorb the other?
5. Magnitude: does strong community structure close the 92x gap to real
   networks (synthetic mu=0.15 plateau vs real median 0.19668)?

Expected signs: Spearman(plateau, mu) < 0, Spearman(plateau, Q) > 0.
"""

import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

BASE = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(BASE, "synthetic_results")
REAL_MEDIAN = 0.19668
ASSORT_SYNTH_MEDIAN = 0.00214


def rho_p(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3 or np.std(x[m]) == 0 or np.std(y[m]) == 0:
        return np.nan, np.nan, int(m.sum())
    r, p = spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


def ols(y, X, names):
    X = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    xtx_inv = np.linalg.pinv(X.T @ X)
    sigma2 = float(resid @ resid) / dof
    se = np.sqrt(np.maximum(np.diag(xtx_inv) * sigma2, 1e-300))
    from scipy.stats import t as tdist
    tvals = beta / se
    pvals = 2 * (1 - tdist.cdf(np.abs(tvals), dof))
    r2 = 1 - float(resid @ resid) / float(((y - y.mean()) ** 2).sum())
    return {"names": ["intercept"] + list(names), "coef": beta, "se": se,
            "t": tvals, "p": pvals, "r2": r2, "dof": dof}


def show(model):
    out = ["    {:12s} {:>10s} {:>9s} {:>8s} {:>10s}".format(
        "term", "coef", "se", "t", "p")]
    for i, name in enumerate(model["names"]):
        out.append("    {:12s} {:+10.4f} {:9.4f} {:+8.2f} {:10.3g}".format(
            name, model["coef"][i], model["se"][i],
            model["t"][i], model["p"][i]))
    out.append("    R2 = {:.3f}  dof = {}".format(model["r2"], model["dof"]))
    return out


def main():
    frames = []
    for name in sorted(os.listdir(SRC_DIR)):
        if name.startswith("community_b") and name.endswith(".csv"):
            try:
                frames.append(pd.read_csv(os.path.join(SRC_DIR, name)))
            except Exception:
                pass
    if not frames:
        raise SystemExit("no community sweep output yet")
    d = pd.concat(frames, ignore_index=True)
    d = d.drop_duplicates(subset=["base", "target_mu", "rep"], keep="last")

    L = []
    A = L.append
    A("EXPERIMENT C: COMMUNITY STRUCTURE AND THE MAGNITUDE GAP")
    A("=" * 62)
    A("graphs: {}".format(len(d)))
    A("status: {}".format(d["status"].value_counts().to_dict()))
    A("")

    ok = d[np.isfinite(d["var_ratio"].astype(float))].copy()
    ok["y"] = np.log10(ok["var_ratio"].astype(float))

    A("--- 1. manipulation check ---")
    mc = ok.groupby("target_mu").agg(
        n=("achieved_mu", "size"),
        achieved_mu=("achieved_mu", "mean"),
        modularity=("modularity", "mean"),
        assortativity=("final_r", "mean"),
        transitivity=("transitivity", "mean"),
        plateau=("var_ratio", "median"),
    )
    A(mc.to_string(float_format=lambda v: "{:+.4f}".format(v)))
    A("  median |achieved - target| = {:.3f}".format(
        float((ok["achieved_mu"] - ok["target_mu"]).abs().median())))
    A("")

    A("--- 2. plateau vs community strength ---")
    for cov, expect in (("target_mu", "negative"), ("modularity", "positive"),
                        ("final_r", "?"), ("transitivity", "?"),
                        ("edges", "null")):
        r, p, n = rho_p(ok[cov], ok["y"])
        A("  Spearman(log10 plateau, {:14s}) rho={:+.3f} p={:.3g} n={}  "
          "(expect {})".format(cov, r, p, n, expect))
    A("")

    A("--- 3. within-base (degree sequence fixed) ---")
    signs, monos = [], []
    for base, sub in ok.groupby("base"):
        rm, pm, nm = rho_p(sub["modularity"], sub["y"])
        ru, pu, nu = rho_p(sub["target_mu"], sub["y"])
        srt = sub.sort_values("modularity")
        mono = bool(srt["var_ratio"].is_monotonic_increasing) \
            if len(srt) > 1 else None
        A("  base {} n={:<3d} rho(modularity)={:+.3f} (p={:.3f}) | "
          "rho(mu)={:+.3f} (p={:.3f}) | monotone_up={}".format(
              base, len(sub), rm, pm, ru, pu, mono))
        if np.isfinite(rm):
            signs.append(rm)
        if mono is not None:
            monos.append(mono)
    A("")
    if signs:
        A("  bases with rho(modularity) > 0: {}/{}".format(
            sum(1 for s in signs if s > 0), len(signs)))
        A("  bases with strictly monotone increase: {}/{}".format(
            sum(1 for m in monos if m), len(monos)))
        A("  median within-base rho = {:+.3f}".format(float(np.median(signs))))
        if len(signs) >= 5:
            try:
                from scipy.stats import wilcoxon
                w = wilcoxon(np.asarray(signs))
                A("  Wilcoxon(within-base rho > 0): p={:.3g}".format(w.pvalue))
            except Exception:
                pass
    A("")
    A("  NOTE: cells are matched-horizon (horizon_reached), not saturated.")
    A("  A higher residual may mean slower convergence, not a higher fixed")
    A("  point. See the saturation check before claiming the latter.")
    A("")

    A("--- 4. joint model with base fixed effects ---")
    use = ok.dropna(subset=["y", "modularity", "final_r", "transitivity"])
    if len(use) >= 10:
        y = use["y"].values
        X = use[["modularity", "final_r", "transitivity"]].astype(float).values
        names_cov = ["modularity", "final_r", "transitivity"]
        dummies = np.zeros((len(use), max(use["base"].nunique() - 1, 0)))
        names = []
        bases = sorted(use["base"].unique().tolist())
        for j, b in enumerate(bases[1:]):
            dummies[:, j] = (use["base"].values == b).astype(float)
            names.append("base{}".format(b))
        X_full = np.column_stack([X, dummies])
        model = ols(y, X_full, names_cov + names)
        L.extend(show(model))
        i_q = model["names"].index("modularity")
        i_r = model["names"].index("final_r")
        i_t = model["names"].index("transitivity")
        A("  -> modularity beta={:+.3f} (p={:.3g}); assortativity beta={:+.3f} "
          "(p={:.3g}); transitivity beta={:+.3f} (p={:.3g})".format(
              model["coef"][i_q], model["p"][i_q],
              model["coef"][i_r], model["p"][i_r],
              model["coef"][i_t], model["p"][i_t]))
    A("")

    A("--- 5. magnitude: does community structure close the 92x gap? ---")
    weak = ok[np.isclose(ok["target_mu"], ok["target_mu"].max())]
    strong = ok[np.isclose(ok["target_mu"], ok["target_mu"].min())]
    A("  weakest community (mu={:.2f}) plateau median = {:.5f}  (n={})".format(
        float(weak["target_mu"].iloc[0]), float(weak["var_ratio"].median()),
        len(weak)))
    A("  strongest community (mu={:.2f}) plateau median = {:.5f}  (n={})".format(
        float(strong["target_mu"].iloc[0]), float(strong["var_ratio"].median()),
        len(strong)))
    if len(weak) and len(strong):
        fold = float(strong["var_ratio"].median()) / \
            max(float(weak["var_ratio"].median()), 1e-12)
        A("  fold change from weakest to strongest = {:.1f}x".format(fold))
        A("  real-network median plateau           = {:.5f}".format(REAL_MEDIAN))
        A("  remaining gap to real networks        = {:.1f}x".format(
            REAL_MEDIAN / max(float(strong["var_ratio"].median()), 1e-12)))
        A("  (assortativity sweep synthetic median = {:.5f})".format(
            ASSORT_SYNTH_MEDIAN))
    A("")

    A("--- verdict ---")
    rq, pq, _ = rho_p(ok["modularity"], ok["y"])
    if np.isfinite(rq):
        A("  Spearman(plateau, modularity) = {:+.3f} (p={:.3g})".format(rq, pq))
        A("  If positive and the magnitude gap closes, community structure is")
        A("  the dominant driver of the residual SIZE and the paper gains a")
        A("  two-factor law: plateau = f(assortativity, community strength).")
    with open(os.path.join(SRC_DIR, "COMMUNITY_REPORT.txt"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(L))
    ok.to_csv(os.path.join(SRC_DIR, "community_merged.csv"), index=False)
    print("\n".join(L))


if __name__ == "__main__":
    main()
