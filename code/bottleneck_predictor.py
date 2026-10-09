# -*- coding: utf-8 -*-
"""L3: is the residual set by the heterogeneity of edge bottleneck-ness?

Why this hypothesis
-------------------
Experiment B and C established that, *within* a controlled family, both
assortativity and planted community strength raise the residual plateau.
But the run-to-saturation check showed that the apparent community effect at
a matched horizon is largely a convergence-speed artefact, and even saturated
synthetic graphs sit near 0.002 while real networks sit near 0.197 -- a gap of
roughly 90x that neither factor explains.

Ollivier curvature is known to be strongly negative on bottleneck edges and
positive inside well-connected clusters. So the residual variance of the
curvature vector should be driven by how *heterogeneous* the bottleneck-ness
of edges is. Configuration-model graphs have nearly homogeneous edge
betweenness; real networks have heavy-tailed edge betweenness. If the residual
tracks betweenness heterogeneity, the 90x gap stops being mysterious.

Design
------
  * one pooled sample: real converged networks + synthetic graphs (B and C)
  * predictor computed on the INITIAL graph (before any flow), so it cannot be
    contaminated by the flow itself
  * primary statistic: Spearman(log10 plateau, CV of edge betweenness)
  * secondary: does the synthetic/real gap in the predictor match the gap in
    the plateau?

Caveat: betweenness is computed on the unweighted skeleton for comparability
across weighted and unweighted datasets.
"""

import glob
import os
import sys

import numpy as np
import pandas as pd
import networkx as nx
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import ricci_flow_program_100data as R                     # noqa: E402
from ricci_flow_longrun import build_base_args             # noqa: E402
from synthetic_assortativity import (base_graph,           # noqa: E402
                                     rewire_to_target)
from synthetic_community2 import rewire_to_mu, cross_fraction  # noqa: E402

DATA_DIR = r"D:\100data"
ANALYSIS = os.path.join(HERE, "ricci_flow_longrun_results", "analysis")
SYNTH = os.path.join(HERE, "synthetic_results")
OUT = os.path.join(ANALYSIS, "BOTTLENECK_REPORT.txt")


def bottleneck_features(graph):
    """Edge-betweenness heterogeneity and friends, on the skeleton."""
    skeleton = nx.Graph(graph)
    feats = {"nodes": skeleton.number_of_nodes(),
             "edges": skeleton.number_of_edges()}
    try:
        eb = nx.edge_betweenness_centrality(skeleton, normalized=True)
        vals = np.asarray(list(eb.values()), dtype=float)
    except Exception:
        return None
    if len(vals) < 3:
        return None
    mu = float(vals.mean())
    sd = float(vals.std(ddof=1))
    feats["eb_mean"] = mu
    feats["eb_sd"] = sd
    feats["eb_cv"] = sd / mu if mu > 0 else np.nan
    top = np.sort(vals)[::-1]
    k = max(1, int(round(0.05 * len(vals))))
    feats["eb_top5_share"] = float(top[:k].sum()) / float(vals.sum()) \
        if vals.sum() > 0 else np.nan
    feats["eb_gini"] = float(
        np.mean(np.abs(vals[:, None] - vals[None, :])) / (2 * mu)) \
        if mu > 0 else np.nan
    try:
        feats["transitivity"] = float(nx.transitivity(skeleton))
        feats["assortativity"] = float(
            nx.degree_assortativity_coefficient(skeleton))
        deg = np.asarray([d for _, d in skeleton.degree()], dtype=float)
        feats["degree_cv"] = float(deg.std(ddof=1) / deg.mean())
    except Exception:
        pass
    return feats


def real_rows():
    path = os.path.join(ANALYSIS, "07_merged_with_structure.csv")
    meta = pd.read_csv(path)
    keep = meta[(meta["converged"] == True)                      # noqa: E712
                & (meta["initial_variance"] > 1e-20)
                & (meta["nodes"] >= 20) & (meta["edges"] >= 20)
                & np.isfinite(meta["var_ratio_final"])
                & (meta["var_ratio_final"] > 0)].copy()
    base_args = build_base_args(["--data_dir", DATA_DIR])
    files = R.discover_graph_files(DATA_DIR)
    index = {}
    for p in files:
        try:
            index[R.safe_dataset_id(p, DATA_DIR)] = p
        except Exception:
            continue
    rows = []
    for _, rec in keep.iterrows():
        name = rec["dataset"]
        if name not in index:
            continue
        try:
            graph = R.load_graph(index[name], base_args)
        except Exception:
            continue
        feats = bottleneck_features(graph)
        if feats is None:
            continue
        feats.update({"dataset": name, "source": "real",
                      "family": rec.get("family", ""),
                      "plateau": float(rec["var_ratio_final"])})
        rows.append(feats)
    return rows


def synth_rows():
    rows = []
    # experiment B: assortativity sweep
    for path in sorted(glob.glob(os.path.join(SYNTH, "assort_base*.csv"))):
        try:
            df = pd.read_csv(path)
        except Exception:
            continue
        for _, rec in df.iterrows():
            base = int(rec["base"])
            level = float(rec["target_r"])
            rep = int(rec["rep"])
            nodes = int(rec["nodes"])
            graph = base_graph(nodes, 6, 1000 + base)
            rng = np.random.default_rng(
                1000 + 97 * base + 13 * rep + int(level * 10) + 1000)
            rewired, _ = rewire_to_target(graph, level, rng)
            feats = bottleneck_features(rewired)
            if feats is None:
                continue
            feats.update({"dataset": "synB_b{}_r{}_rep{}".format(
                             base, level, rep), "source": "synthB",
                          "family": "synthB",
                          "plateau": float(rec["var_ratio"])})
            rows.append(feats)
    # experiment C: community sweep
    seen = set()
    dfs = []
    for path in sorted(glob.glob(os.path.join(SYNTH, "community_*.csv"))):
        if os.path.basename(path) == "community_merged.csv":
            continue
        try:
            dfs.append(pd.read_csv(path))
        except Exception:
            continue
    if dfs:
        df = pd.concat(dfs, ignore_index=True)
        df = df.drop_duplicates(subset=["base", "target_mu", "rep"])
        for _, rec in df.iterrows():
            key = (int(rec["base"]), float(rec["target_mu"]), int(rec["rep"]))
            if key in seen:
                continue
            seen.add(key)
            base = int(rec["base"])
            level = float(rec["target_mu"])
            rep = int(rec["rep"])
            nodes = int(rec["nodes"])
            graph = base_graph(nodes, 6, 2000 + base)
            rng0 = np.random.default_rng(2000 + 31 * base)
            order = list(graph.nodes())
            rng0.shuffle(order)
            groups = {node: pos % 4 for pos, node in enumerate(order)}
            rng = np.random.default_rng(
                2000 + 97 * base + 13 * rep + int(level * 100) + 500)
            rewired, _ = rewire_to_mu(graph, groups, level, rng)
            feats = bottleneck_features(rewired)
            if feats is None:
                continue
            feats.update({"dataset": "synC_b{}_mu{}_rep{}".format(
                             base, level, rep), "source": "synthC",
                          "family": "synthC",
                          "plateau": float(rec["var_ratio"])})
            rows.append(feats)
    return rows


def spear(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 5 or np.std(x[m]) == 0 or np.std(y[m]) == 0:
        return np.nan, np.nan, int(m.sum())
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


def main():
    rows = real_rows() + synth_rows()
    d = pd.DataFrame(rows)
    d = d[np.isfinite(d["plateau"]) & (d["plateau"] > 0)].copy()
    d["y"] = np.log10(d["plateau"].astype(float))
    d.to_csv(os.path.join(ANALYSIS, "bottleneck_merged.csv"), index=False)

    L = []
    A = L.append
    A("BOTTLENECK HETEROGENEITY AND THE RESIDUAL PLATEAU")
    A("=" * 68)
    A("sample: {} graphs ({})".format(len(d), d["source"].value_counts().to_dict()))
    A("")
    A("--- 1. does the synthetic/real plateau gap match a predictor gap? ---")
    grp = d.groupby("source").agg(
        n=("plateau", "size"),
        plateau_median=("plateau", "median"),
        eb_cv=("eb_cv", "median"),
        eb_gini=("eb_gini", "median"),
        eb_top5=("eb_top5_share", "median"),
        transitivity=("transitivity", "median"),
        degree_cv=("degree_cv", "median"),
    )
    A(grp.to_string(float_format=lambda v: "{:+.4f}".format(v)))
    real_med = float(d[d["source"] == "real"]["plateau"].median())
    syn_med = float(d[d["source"] != "real"]["plateau"].median())
    A("")
    A("plateau gap real/synthetic = {:.1f}x".format(real_med / syn_med))
    for cov in ("eb_cv", "eb_gini", "eb_top5_share", "transitivity",
                "degree_cv"):
        r_med = float(d[d["source"] == "real"][cov].median())
        s_med = float(d[d["source"] != "real"][cov].median())
        A("  {:14s} real={:+.4f} synth={:+.4f} ratio={:.2f}".format(
            cov, r_med, s_med, r_med / s_med if s_med else np.nan))
    A("")

    A("--- 2. pooled Spearman(log10 plateau, predictor) ---")
    for cov in ("eb_cv", "eb_gini", "eb_top5_share", "eb_mean", "eb_sd",
                "transitivity", "degree_cv", "assortativity", "nodes",
                "edges"):
        if cov not in d.columns:
            continue
        rho, p, n = spear(d[cov], d["y"])
        A("  {:14s} rho={:+.3f} p={:.3g} n={}".format(cov, rho, p, n))
    A("")

    A("--- 3. within real networks only (n={}) ---".format(
        int((d["source"] == "real").sum())))
    real = d[d["source"] == "real"]
    for cov in ("eb_cv", "eb_gini", "eb_top5_share", "transitivity",
                "degree_cv", "assortativity"):
        rho, p, n = spear(real[cov], real["y"])
        A("  {:14s} rho={:+.3f} p={:.3g} n={}".format(cov, rho, p, n))
    A("")

    A("--- 4. within synthetic only (the controlled side) ---")
    syn = d[d["source"] != "real"]
    for cov in ("eb_cv", "eb_gini", "eb_top5_share", "transitivity",
                "degree_cv", "assortativity"):
        rho, p, n = spear(syn[cov], syn["y"])
        A("  {:14s} rho={:+.3f} p={:.3g} n={}".format(cov, rho, p, n))
    A("")

    A("--- 5. does the predictor separate real from synthetic? ---")
    for cov in ("eb_cv", "eb_gini", "eb_top5_share", "transitivity",
                "degree_cv"):
        a = d[d["source"] == "real"][cov].astype(float).dropna().values
        b = d[d["source"] != "real"][cov].astype(float).dropna().values
        if len(a) < 5 or len(b) < 5:
            continue
        try:
            u = stats.mannwhitneyu(a, b, alternative="two-sided")
            A("  {:14s} MWU p={:.3g}  (real median {:.4f} vs synth {:.4f})"
              .format(cov, u.pvalue, float(np.median(a)), float(np.median(b))))
        except Exception:
            pass
    A("")

    A("--- reading ---")
    A("If eb_cv (or eb_gini) is both (a) far larger in real networks and")
    A("(b) strongly rank-correlated with the plateau POOLED across real and")
    A("synthetic graphs, then the 90x gap is a betweenness-heterogeneity")
    A("effect and the paper gains a single mechanism that covers both the")
    A("within-family causal results and the cross-family magnitude gap.")
    A("If it is large but NOT correlated, report it as a rejected hypothesis.")

    text = "\n".join(L)
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
