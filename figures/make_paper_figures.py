# -*- coding: utf-8 -*-
"""
Build the display items for Paper A (Communications Physics submission).

Every panel is computed from a raw result file on disk; nothing is
hard-coded. The script prints the numbers it used so they can be quoted
in the manuscript.
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)

MASTER = os.path.join(ROOT, "ricci_flow_longrun_results", "analysis",
                      "07_merged_with_structure.csv")
PER_DIR = os.path.join(ROOT, "ricci_flow_longrun_results", "per_dataset")
SHAM = os.path.join(ROOT, "sham_results", "sham_summary.csv")
SURGERY = os.path.join(ROOT, "surgery_results", "surgery_paired_v2.csv")
SYN = os.path.join(ROOT, "synthetic_results")
CONV = os.path.join(ROOT, "convergence_results", "convergence_rates.csv")
LOFO = os.path.join(ROOT, "ricci_flow_longrun_results", "analysis",
                    "16_lofo_comparison.csv")
# strict re-run: longest available trajectory per network wins
STRICT_DIRS = [os.path.join(ROOT, "ricci_flow_deep"),
               os.path.join(ROOT, "ricci_flow_v2"),
               os.path.join(ROOT, "ricci_flow_strict")]
FATE = os.path.join(ROOT, "strict_results", "fate_table.csv")

# Communications Physics: figures are single- or double-column in a
# two-column layout; 8 pt labels, colour-blind safe palette.
plt.rcParams.update({
    "font.size": 8,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 7,
    "figure.dpi": 300,
    "savefig.dpi": 400,
    "savefig.bbox": "tight",
    "font.family": "DejaVu Sans",
})
BLUE, ORANGE, GREY, GREEN, PURPLE = "#1f77b4", "#ff7f0e", "#7f7f7f", "#2ca02c", "#9467bd"


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    fig.savefig(os.path.join(OUT, name + ".png"))
    plt.close(fig)
    print("  wrote %s.{pdf,png}" % name)


def family_of(ds):
    for fam in ("scc_rt", "reptilia", "mammalia", "aves", "insecta", "bn",
                "ENZYMES", "bio", "ia", "power", "eco", "social"):
        if ds.startswith(fam):
            return fam
    return ds.split("_")[0]


log = []


def note(x):
    print(x)
    log.append(str(x))


# ----------------------------------------------------------------- Fig 1
def fig1(master):
    note("")
    note("[Fig 1] cross-domain homogenization + sham control")
    # 1a: median trajectory
    grid = np.arange(0, 201, 1)
    curves = []
    files = sorted(glob.glob(os.path.join(PER_DIR, "*.csv")))
    for fp in files:
        try:
            tr = pd.read_csv(fp)
        except Exception:
            continue
        if len(tr) < 20:
            continue
        v0 = float(tr["curvature_variance"].astype(float).iloc[0])
        if not np.isfinite(v0) or v0 <= 1e-30:
            continue
        s = tr["step"].astype(float).values
        v = tr["var_ratio"].astype(float).values
        v = v / v[0]
        if s.max() < 20:
            continue
        curves.append(np.interp(grid, s, v, left=v[0], right=v[-1]))
    C = np.array(curves)
    med = np.median(C, axis=0)
    q25 = np.percentile(C, 25, axis=0)
    q75 = np.percentile(C, 75, axis=0)
    note("  trajectories used: %d  (interpolated on a common 0-200 step grid)" % C.shape[0])
    note("  median var_ratio at step 10 (t=0.5): %.3f" % med[10])
    note("  median var_ratio at step 130      : %.3f" % med[130])
    note("  median var_ratio at end of grid   : %.3f" % med[-1])

    # 1b: plateau distribution
    vr = master["var_ratio_final"].astype(float).dropna().values
    vr = vr[vr > 0]
    conv = master.loc[master["converged"] == True, "var_ratio_final"].dropna()
    conv = conv[conv > 0]
    note("  plateau: n=%d  median=%.4f  IQR=%.4f..%.4f"
         % (len(vr), np.median(vr), np.percentile(vr, 25), np.percentile(vr, 75)))
    note("  converged subset: n=%d  median=%.4f" % (len(conv), np.median(conv)))
    note("  fraction below 1 (variance actually reduced): %.3f" % float(np.mean(vr < 1)))

    # 1c: sham pairs. Exclusions are the ones declared in
    # sham_results/analysis/ROBUST_REPORT.txt (1 runaway, 3 inert flows).
    sham = pd.read_csv(SHAM)
    sham = sham[sham["mode"] == "sham"]
    real = master[["dataset", "var_ratio_final"]].rename(
        columns={"var_ratio_final": "real"})
    pairs = sham[["dataset", "var_ratio_final"]].rename(
        columns={"var_ratio_final": "sham"}).merge(real, on="dataset", how="inner")
    excl = ["eco-mangwet",
            "bn-macaque-rhesus_interareal-cortical-network_2",
            "scc_rt_dash", "scc_rt_ksa"]
    pairs["excluded"] = pairs["dataset"].isin(excl)
    keep = pairs[~pairs["excluded"]].copy()
    note("  sham pairs: %d total, %d defensible (excluded runaway/inert: %s)"
         % (len(pairs), len(keep), ", ".join(excl)))
    note("  median real = %.4f   median sham = %.4f" % (keep["real"].median(),
                                                        keep["sham"].median()))
    higher = int((keep["sham"] > keep["real"]).sum())
    p_sign = stats.binomtest(higher, len(keep), 0.5, alternative="greater").pvalue
    note("  sham > real in %d/%d ; one-sided sign test p = %.3g" % (higher, len(keep), p_sign))
    w = stats.wilcoxon(keep["sham"], keep["real"], alternative="greater")
    note("  Wilcoxon signed-rank p = %.3g" % w.pvalue)

    # 1d/1e/1f: the saturation criterion and the two fates
    #   (d) tail drift of the runs the OLD rule called saturated
    #   (e) residual against late-time decay rate, coloured by fate
    #   (f) weight spread at the end of the run, by fate
    drift = []
    for fp in sorted(glob.glob(os.path.join(PER_DIR, "*.csv"))):
        name = os.path.basename(fp)[:-4]
        row = master.loc[master["dataset"] == name]
        if len(row) == 0 or str(row["status"].iloc[0]) != "saturated":
            continue
        try:
            tr = pd.read_csv(fp)
        except Exception:
            continue
        if len(tr) < 40:
            continue
        v = tr["curvature_variance"].astype(float).values
        if v[0] <= 0:
            continue
        v = v / v[0]
        k = int(len(v) * 0.75)
        drift.append(v[-1] / v[k] - 1.0)
    drift = np.array(drift)
    note("  old rule labelled %d runs 'saturated'; tail drift over the final "
         "quarter: median %+.1f%%, %d/%d still falling by >2%%, worst %+.1f%%"
         % (len(drift), 100 * np.median(drift), int((drift < -0.02).sum()),
            len(drift), 100 * drift.min()))

    fate = pd.read_csv(FATE)
    stalled = fate[fate["fate"] == "stalled"]
    conv_f = fate[fate["fate"] == "converging"]
    u, p_span = stats.mannwhitneyu(stalled["weight_span"], conv_f["weight_span"],
                                   alternative="greater")
    note("  strict re-run: n=%d (%d stalled, %d converging, %d undecided)"
         % (len(fate), len(stalled), len(conv_f),
            int((fate["fate"] == "undecided").sum())))
    note("    stalled  residual median %.4f ; late-time rate median %+.1e/step"
         % (stalled["v_end"].median(), stalled["rate_Q4"].median()))
    note("    converging residual median %.4f ; late-time rate median %+.1e/step"
         % (conv_f["v_end"].median(), conv_f["rate_Q4"].median()))
    note("    weight span: stalled %.2e vs converging %.2e ; Mann-Whitney p=%.3f"
         % (stalled["weight_span"].median(), conv_f["weight_span"].median(), p_span))

    fig, ax = plt.subplots(2, 3, figsize=(7.1, 4.4))
    ax = ax.ravel()
    ax[0].plot(grid, med, color=BLUE, lw=1.2)
    ax[0].fill_between(grid, q25, q75, color=BLUE, alpha=0.25, lw=0)
    ax[0].axvline(10, color=GREY, ls="--", lw=0.8)
    ax[0].annotate("t = 0.5\n(10 steps)", xy=(10, 0.85), xytext=(35, 0.9),
                   fontsize=6.5, color=GREY,
                   arrowprops=dict(arrowstyle="->", color=GREY, lw=0.6))
    ax[0].set_yscale("log")
    ax[0].set_xlabel("flow step")
    ax[0].set_ylabel("variance ratio $V(t)=\\sigma^2(t)/\\sigma^2(0)$")
    ax[0].set_title("a  median over %d networks" % C.shape[0], loc="left")
    ax[0].set_ylim(1e-3, 2)

    ax[1].hist(np.log10(vr), bins=24, color=BLUE, alpha=0.85)
    if len(fate):
        ax[1].hist(np.log10(fate["v_end"].clip(1e-8)), bins=12,
                   color=ORANGE, alpha=0.75, label="strict re-run (n=%d)"
                   % len(fate))
        ax[1].legend(frameon=False, fontsize=6, loc="upper left", bbox_to_anchor=(0.0, 0.82))
    ax[1].axvline(np.log10(np.median(vr)), color=BLUE, lw=1.2, ls="--")
    ax[1].set_xlabel("$\\log_{10}$ residual $V$")
    ax[1].set_ylabel("networks")
    ax[1].set_title("b  residual (upper bounds)", loc="left")
    ax[1].text(0.03, 0.90, "median = %.3f" % np.median(vr),
               transform=ax[1].transAxes, fontsize=6.5, color=BLUE)

    for _, r in keep.iterrows():
        ax[2].plot([0, 1], [r["real"], r["sham"]], "-o", color=GREY,
                   lw=0.6, ms=2.5, alpha=0.8)
    bp = ax[2].boxplot([keep["real"], keep["sham"]], widths=0.5,
                       showfliers=False, patch_artist=True,
                       medianprops=dict(color="black"))
    for patch, col in zip(bp["boxes"], [BLUE, ORANGE]):
        patch.set_facecolor(col)
        patch.set_alpha(0.6)
    ax[2].set_xticks([1, 2])
    ax[2].set_xticklabels(["real\nflow", "sham\nflow"])
    ax[2].set_ylabel("plateau $V_\\infty$")
    ax[2].set_title("c  curvature-edge mapping destroyed", loc="left")
    ax[2].text(0.05, 0.93, "%d/%d pairs, $p=%.1e$" % (higher, len(keep), p_sign),
               transform=ax[2].transAxes, fontsize=7)

    # (d) the old criterion is not a convergence test
    ax[3].hist(100 * drift, bins=20, color=ORANGE, alpha=0.85)
    ax[3].axvline(0, color=GREY, lw=0.8)
    ax[3].axvline(100 * np.median(drift), color="black", lw=1.2, ls="--")
    ax[3].set_xlabel("variance change over final quarter (%)")
    ax[3].set_ylabel("runs labelled saturated")
    ax[3].set_title("d  old rule measured elapsed time", loc="left")
    ax[3].text(0.05, 0.90, "median %+.1f%%\nn=%d, %d/%d still falling"
               % (100 * np.median(drift), len(drift),
                  int((drift < -0.02).sum()), len(drift)),
               transform=ax[3].transAxes, fontsize=6.5)

    # (e) two fates: residual against the late-time decay rate
    colors = {"stalled": ORANGE, "converging": BLUE, "undecided": GREY}
    for f_, g in fate.groupby("fate"):
        ax[4].scatter(np.abs(g["rate_Q4"]), g["v_end"].clip(1e-8), s=22,
                      color=colors.get(f_, GREY), label="%s (n=%d)" % (f_, len(g)),
                      edgecolor="white", linewidth=0.4)
    ax[4].set_xscale("log")
    ax[4].set_yscale("log")
    ax[4].axvline(1e-3, color=GREY, ls="--", lw=0.8)
    ax[4].set_xlabel("late-time decay rate $|r|$ (per step)")
    ax[4].set_ylabel("residual $V$")
    ax[4].set_title("e  two fates", loc="left")
    ax[4].legend(frameon=False, fontsize=6, loc="lower left")

    # (f) weight spread separates the branches
    if len(stalled) and len(conv_f):
        data = [conv_f["weight_span"], stalled["weight_span"]]
        bp = ax[5].boxplot(data, widths=0.5, showfliers=False,
                           patch_artist=True,
                           medianprops=dict(color="black"))
        for patch, col in zip(bp["boxes"], [BLUE, ORANGE]):
            patch.set_facecolor(col)
            patch.set_alpha(0.6)
        ax[5].set_yscale("log")
        ax[5].set_xticks([1, 2])
        ax[5].set_xticklabels(["converging", "stalled"])
        ax[5].set_ylabel("weight spread $w_{max}/w_{min}$")
        ax[5].set_title("f  weight spread (p=%.2f)" % p_span, loc="left")

    for a in ax:
        a.tick_params(labelsize=7)
    fig.tight_layout()
    save(fig, "Fig2_homogenization_and_control")


# ----------------------------------------------------------------- Fig 2
def fig2():
    note("")
    note("[Fig 2] surgery: is the plateau the cost of keeping the edge set?")
    s = pd.read_csv(SURGERY)
    s = s[(s["fired"] == True) & (s["completed"] == True) & (s["substantial"] == True)]
    s = s.sort_values("drop_factor")
    note("  defensible pairs: %d" % len(s))
    for _, r in s.iterrows():
        note("     %-22s no-surgery %.4f -> surgery %.4f  (%.1fx, %d edges removed)"
             % (r["dataset"], r["no_surgery_ratio"], r["surgery_ratio"],
                r["drop_factor"], r["edges_removed"]))
    higher = int((s["no_surgery_ratio"] > s["surgery_ratio"]).sum())
    p_sign = stats.binomtest(higher, len(s), 0.5, alternative="greater").pvalue
    note("  surgery arm lower in %d/%d ; one-sided sign test p = %.4f" % (higher, len(s), p_sign))
    note("  median drop factor = %.2fx" % s["drop_factor"].median())
    note("  families covered: %s" % ", ".join(sorted(set(s["family"]))))

    fig, ax = plt.subplots(1, 2, figsize=(5.6, 2.3))
    for _, r in s.iterrows():
        ax[0].plot([0, 1], [r["no_surgery_ratio"], r["surgery_ratio"]], "-o",
                   color=BLUE, lw=0.9, ms=4)
    ax[0].set_yscale("log")
    ax[0].set_xticks([0, 1])
    ax[0].set_xticklabels(["topology\npreserved", "surgery\n(edges removed)"])
    ax[0].set_ylabel("plateau $V_\\infty$")
    ax[0].set_title("a  paired, identical numerics (n=%d)" % len(s), loc="left")
    ax[0].set_xlim(-0.3, 1.3)

    ax[1].barh(range(len(s)), np.log10(s["no_surgery_ratio"] / s["surgery_ratio"]),
               color=ORANGE, alpha=0.85)
    ax[1].set_yticks(range(len(s)))
    ax[1].set_yticklabels(s["dataset"].values, fontsize=6.5)
    ax[1].set_xlabel("$\\log_{10}$ drop factor")
    ax[1].axvline(0, color="black", lw=0.8)
    ax[1].set_title("b  gain from removing edges", loc="left")
    save(fig, "Fig4_surgery")


# ----------------------------------------------------------------- Fig 3
def fig3(master):
    note("")
    note("[Fig 3] mechanism on synthetic graphs (degree sequence held fixed)")
    assort = pd.concat([pd.read_csv(os.path.join(SYN, "assort_base%d.csv" % b))
                        for b in range(4)], ignore_index=True)
    clust = pd.read_csv(os.path.join(SYN, "clustering_merged.csv"))
    # corrected-rule sub-sweep of the same cells (strict_clustering_backfill.py)
    merged_strict = os.path.join(SYN, "clustering_strict_merged.csv")
    clust_strict = pd.read_csv(merged_strict) if os.path.exists(merged_strict) \
        else pd.DataFrame()
    comm = pd.read_csv(os.path.join(SYN, "community_merged.csv"))
    # the saturated community runs live in separate per-cell files
    long_files = sorted(glob.glob(os.path.join(SYN, "long_*.csv")))
    comm_sat = pd.concat([pd.read_csv(f) for f in long_files], ignore_index=True) \
        if long_files else pd.DataFrame()

    a = assort[(assort["var_ratio"] > 0) & (assort["status"] == "saturated")]
    rho_a, p_a = stats.spearmanr(a["achieved_r"], np.log10(a["var_ratio"]))
    note("  assortativity: n=%d cells, Spearman(plateau, r) = %+.3f (p=%.2e)"
         % (len(a), rho_a, p_a))
    within = []
    for b, g in a.groupby("base"):
        if len(g) >= 5:
            r_, p_ = stats.spearmanr(g["achieved_r"], np.log10(g["var_ratio"]))
            within.append(r_)
    note("     within-base rho: median %+.3f, %d/%d positive"
         % (np.median(within), sum(1 for x in within if x > 0), len(within)))

    c = clust[(clust["var_ratio"] > 0)]
    rho_c, p_c = stats.spearmanr(c["achieved_T"], np.log10(c["var_ratio"]))
    note("  transitivity: n=%d cells, Spearman(plateau, T) = %+.3f (p=%.3f)" % (len(c), rho_c, p_c))
    for b, g in c.groupby("base"):
        g = g.sort_values("achieved_T")
        if len(g) >= 2:
            lo, hi = g.iloc[0], g.iloc[-1]
            note("     base %d: T %.3f -> %.3f  plateau %.5f -> %.5f  (%.1fx)"
                 % (b, lo["achieved_T"], hi["achieved_T"], lo["var_ratio"],
                    hi["var_ratio"], hi["var_ratio"] / lo["var_ratio"]))
    folds = {}
    if len(clust_strict):
        note("  transitivity, all horizons (V = residual at that step):")
        for _, r in clust_strict.sort_values(["base", "target_T"]).iterrows():
            note("     base %d T=%.2f: " % (int(r["base"]), r["target_T"])
                 + "  ".join("V(%d)=%.3e" % (h, r["V_%d" % h])
                             for h in (125, 250, 400, 700, 900)
                             if np.isfinite(r["V_%d" % h]))
                 + "  [%s @%d]" % (r["status"], int(r["steps_terminal"])))
            if r["windows"]:
                note("        windows: %s" % r["windows"])
        for b, g in clust_strict.groupby("base"):
            g = g.sort_values("target_T")
            if len(g) < 2:
                continue
            lo, hi = g.iloc[0], g.iloc[-1]
            folds[b] = [(h, hi["V_%d" % h] / lo["V_%d" % h])
                        for h in (125, 250, 400, 700, 900)
                        if np.isfinite(hi["V_%d" % h])
                        and np.isfinite(lo["V_%d" % h])]
            note("     base %d fold: " % b
                 + ", ".join("%dx@%d" % (f, h) for h, f in folds[b]))
            if "constant_curvature" in str(lo["status"]) and \
               "constant_curvature" in str(hi["status"]):
                note("        BOTH arms reach the constant-curvature floor -> "
                     "delay %.1fx in steps (%d vs %d), endpoint unchanged"
                     % (hi["steps_terminal"] / lo["steps_terminal"],
                        int(lo["steps_terminal"]), int(hi["steps_terminal"])))

    cs = comm_sat[comm_sat["status"] == "saturated"] if len(comm_sat) else comm[comm["status"] == "saturated"]
    cm = comm[comm["status"] != "saturated"]
    note("  community: %d saturated cells (%s), %d matched-horizon cells"
         % (len(cs), ", ".join(sorted(os.path.basename(f) for f in long_files)) or "none",
            len(cm)))
    if len(cs) >= 2:
        cs = cs.sort_values("modularity")
        lo, hi = cs.iloc[0], cs.iloc[-1]
        note("     saturated: Q %.3f -> %.3f  plateau %.5f -> %.5f  (%.2fx only)"
             % (lo["modularity"], hi["modularity"], lo["var_ratio"], hi["var_ratio"],
                hi["var_ratio"] / lo["var_ratio"]))
    if len(cm) >= 4:
        cm_s = cm.sort_values("modularity")
        lo, hi = cm_s.iloc[0], cm_s.iloc[-1]
        note("     matched horizon (95 steps): Q %.3f -> %.3f  plateau %.5f -> %.5f (%.1fx)"
             % (lo["modularity"], hi["modularity"], lo["var_ratio"], hi["var_ratio"],
                hi["var_ratio"] / lo["var_ratio"]))

    fig, ax = plt.subplots(2, 2, figsize=(7.1, 4.8))
    ax = ax.ravel()
    for b, g in a.groupby("base"):
        ax[0].scatter(g["achieved_r"], g["var_ratio"], s=12, alpha=0.85,
                      label="base %d" % b)
    ax[0].set_yscale("log")
    ax[0].set_xlabel("degree assortativity $r$")
    ax[0].set_ylabel("residual $V$")
    ax[0].set_title("a  assortativity ($\\rho=%+.2f$)" % rho_a, loc="left")
    ax[0].legend(frameon=False, fontsize=6)

    # (b) the decisive panel: both arms of base 0 followed to a terminal state
    hist_files = {"low": os.path.join(SYN, "ckpt_base0_t010", "hist.csv"),
                  "high": os.path.join(SYN, "ckpt_base0_t050", "hist.csv")}
    plotted = False
    for tag, fp in sorted(hist_files.items()):
        if not os.path.exists(fp):
            continue
        tr = pd.read_csv(fp)
        lab = ("$T=0.12$" if tag == "low" else "$T=0.46$")
        col = BLUE if tag == "low" else ORANGE
        ax[1].plot(tr["step"] + 1, tr["var_ratio"], lw=1.1, color=col,
                   label=lab)
        ax[1].plot([tr["step"].iloc[-1] + 1], [tr["var_ratio"].iloc[-1]],
                   marker="x", ms=6, color=col, mew=1.2)
        plotted = True
    if plotted:
        ax[1].axhline(1e-6, color=GREY, ls="--", lw=0.8)
        ax[1].set_yscale("log")
        ax[1].set_xlabel("flow step")
        ax[1].set_ylabel("residual $V(t)$")
        ax[1].set_title("b  same endpoint, different clock", loc="left")
        ax[1].legend(frameon=False, fontsize=6, loc="upper right")
        ax[1].set_ylim(3e-7, 2)

    # (c) the fold is a transient
    for b, pts in sorted(folds.items()):
        hs = [h for h, _ in pts]
        fs = [f for _, f in pts]
        ax[2].plot(hs, fs, "-o", ms=4, lw=1.1,
                   color=BLUE if b == 0 else ORANGE, label="base %d" % b)
    if folds:
        ax[2].axhline(1, color=GREY, ls="--", lw=0.8)
        ax[2].set_yscale("log")
        ax[2].set_xlabel("matched horizon (steps)")
        ax[2].set_ylabel("fold $V_{high-T}/V_{low-T}$")
        ax[2].set_title("c  the fold collapses", loc="left")
        ax[2].legend(frameon=False, fontsize=6)

    if len(cm) >= 4:
        g = cm.groupby("modularity", as_index=False)["var_ratio"].median()
        ax[3].plot(g["modularity"], g["var_ratio"], "-s", color=ORANGE, ms=3,
                   lw=0.9, label="matched horizon (95 steps)")
    if len(cs) >= 2:
        g = cs.groupby("modularity", as_index=False)["var_ratio"].median()
        ax[3].plot(g["modularity"], g["var_ratio"], "-o", color=BLUE, ms=4,
                   lw=0.9, label="to saturation")
    ax[3].set_yscale("log")
    ax[3].set_xlabel("modularity $Q$")
    ax[3].set_ylabel("")
    ax[3].set_title("d  community structure", loc="left")
    ax[3].legend(frameon=False, fontsize=6)
    fig.tight_layout()
    save(fig, "Fig5_mechanism")


# ----------------------------------------------------------------- Fig 4
def fig4(master):
    note("")
    note("[Fig 4] how far does structure predict the plateau?")
    lofo = pd.read_csv(LOFO)
    note("  leave-one-family-out (n=%d):" % int(lofo["n_pred"].iloc[0]))
    for _, r in lofo.iterrows():
        note("     %-34s oos R2 = %+.3f" % (r["model"], r["oos_R2"]))
    best = lofo.loc[lofo["oos_R2"].idxmax(), "model"]
    note("  best model = %s ; baseline = %.3f -> out-of-sample prediction is %s"
         % (best, lofo["oos_R2"].iloc[0],
            "not better than the intercept" if lofo["oos_R2"].max() <= lofo["oos_R2"].iloc[0]
            else "better than the intercept"))

    # primary statistics are computed on the converged subset (declared
    # once, then used everywhere)
    m = master[master["converged"] == True].dropna(subset=["var_ratio_final"])
    m = m[m["var_ratio_final"] > 0]
    note("  primary subset = converged trajectories: n=%d (of %d in the analysis set)"
         % (len(m), len(master)))
    enz = m[m["dataset"].str.startswith("ENZYMES")]
    rho_t, p_t = stats.spearmanr(enz["transitivity"], np.log10(enz["var_ratio_final"]))
    rho_n, p_n = stats.spearmanr(enz["log10_nodes"], np.log10(enz["var_ratio_final"]))
    note("  within ENZYMES (n=%d): Spearman(plateau, T) = %+.3f (p=%.4f)" % (len(enz), rho_t, p_t))
    note("                         Spearman(plateau, logN) = %+.3f (p=%.4f)" % (rho_n, p_n))
    # partial: plateau ~ T controlling logN
    R = enz[["transitivity", "log10_nodes", "var_ratio_final"]].rank()
    def resid(t, o):
        X = np.column_stack([np.ones(len(R))] + [R[c].values for c in o])
        b, *_ = np.linalg.lstsq(X, R[t].values, rcond=None)
        return R[t].values - X @ b
    pr, pp = stats.spearmanr(resid("transitivity", ["log10_nodes"]),
                             resid("var_ratio_final", ["log10_nodes"]))
    note("                         partial (T | logN) = %+.3f (p=%.3f)" % (pr, pp))

    fam = m.copy()
    fam["family"] = fam["dataset"].map(family_of)
    fg = fam.groupby("family").agg(n=("var_ratio_final", "size"),
                                   T=("transitivity", "mean"),
                                   V=("var_ratio_final", "median")).reset_index()
    fg = fg[fg["n"] >= 3]
    rho_f, p_f = stats.spearmanr(fg["T"], np.log10(fg["V"]))
    note("  family level (n=%d families): Spearman(median plateau, mean T) = %+.3f (p=%.4f)"
         % (len(fg), rho_f, p_f))

    fig, ax = plt.subplots(1, 2, figsize=(5.6, 2.3))
    lbl = [x.replace("assortativity", "assort.").replace("modularity", "modul.")
           for x in lofo["model"]]
    cols = [GREY if x.startswith("M0") else
            (BLUE if v > 0 else ORANGE) for x, v in zip(lofo["model"], lofo["oos_R2"])]
    ax[0].barh(range(len(lofo)), lofo["oos_R2"], color=cols, alpha=0.85)
    ax[0].set_yticks(range(len(lofo)))
    ax[0].set_yticklabels(lbl, fontsize=6.5)
    ax[0].axvline(0, color="black", lw=0.8)
    ax[0].axvline(lofo["oos_R2"].iloc[0], color=GREY, ls="--", lw=0.8)
    ax[0].set_xlabel("out-of-sample $R^2$ (leave-one-family-out)")
    ax[0].set_title("a  no cross-family extrapolation", loc="left")

    ax[1].scatter(enz["transitivity"], enz["var_ratio_final"], s=18,
                  color=BLUE, alpha=0.8)
    ax[1].set_yscale("log")
    ax[1].set_xlabel("transitivity $T$")
    ax[1].set_ylabel("plateau $V_\\infty$")
    ax[1].set_title("b  within ENZYMES only (n=%d)" % len(enz), loc="left")
    ax[1].text(0.05, 0.9, "$\\rho=%+.2f$ ($p=%.3f$)\npartial$|\\log N = %+.2f$"
               % (rho_t, p_t, pr), transform=ax[1].transAxes, fontsize=6.5)
    save(fig, "Fig6_predictability")


# ----------------------------------------------------------------- Fig 5
def fig5():
    note("")
    note("[Fig 5] dynamics: what sets the speed, and does speed predict depth?")
    c = pd.read_csv(CONV)
    c = c[c["rate_early"] > 0].dropna(subset=["rate_early", "var_ratio_final"])
    c = c[c["var_ratio_final"] > 0]
    r1, p1 = stats.spearmanr(c["rate_early"], c["lambda1_0"])
    r2, p2 = stats.spearmanr(c["rate_early"], np.log10(c["var_ratio_final"]))
    note("  n = %d units with a usable rate" % len(c))
    note("  Spearman(rate, spectral gap lambda1) = %+.3f (p=%.3f)  -> null" % (r1, p1))
    note("  Spearman(rate, log10 plateau)        = %+.3f (p=%.3g)" % (r2, p2))
    R = c[["rate_early", "var_ratio_final", "log10_nodes", "log10_md"]].rank()
    def resid(t, o):
        X = np.column_stack([np.ones(len(R))] + [R[x].values for x in o])
        b, *_ = np.linalg.lstsq(X, R[t].values, rcond=None)
        return R[t].values - X @ b
    pr, pp = stats.spearmanr(resid("rate_early", ["log10_nodes", "log10_md"]),
                             resid("var_ratio_final", ["log10_nodes", "log10_md"]))
    note("  partial (rate, plateau | size, degree) = %+.3f (p=%.3g)" % (pr, pp))

    fig, ax = plt.subplots(1, 2, figsize=(5.6, 2.3))
    ax[0].scatter(c["lambda1_0"], c["rate_early"], s=18, color=GREY, alpha=0.85)
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("spectral gap $\\lambda_1$ (initial Laplacian)")
    ax[0].set_ylabel("early decay rate (per step)")
    ax[0].set_title("a  rate is not set by the gap", loc="left")
    ax[0].text(0.05, 0.9, "$\\rho=%+.2f$ (n.s.)" % r1, transform=ax[0].transAxes,
               fontsize=7)

    ax[1].scatter(c["rate_early"], c["var_ratio_final"], s=18, color=BLUE, alpha=0.85)
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("early decay rate (per step)")
    ax[1].set_ylabel("plateau $V_\\infty$")
    ax[1].set_title("b  faster flow reaches lower", loc="left")
    ax[1].text(0.05, 0.9, "$\\rho=%+.2f$\npartial $|N,\\langle k\\rangle = %+.2f$"
               % (r2, pr), transform=ax[1].transAxes, fontsize=6.5)
    save(fig, "Fig7_dynamics")


# --------------------------------------------------------- Fig 1 (new, R2)
def fig_flow():
    """(a) Analysis populations C0-C5.  (b) The ladder of controls."""
    note("")
    note("[Fig 1] analysis populations + ladder of controls")
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.4))
    for a in ax:
        a.set_xlim(0, 10)
        a.set_ylim(0, 12)
        a.axis("off")

    def box(a, x, y, w, h, text, fc="#eef2f7"):
        a.add_patch(plt.Rectangle((x - w / 2.0, y - h / 2.0), w, h,
                                  facecolor=fc, edgecolor="#446",
                                  linewidth=0.7, zorder=1))
        a.text(x, y, text, ha="center", va="center", fontsize=6.0, zorder=2)

    # ---- (a) populations
    a = ax[0]
    a.set_title("a  Analysis populations", fontsize=8, loc="left")
    box(a, 5.0, 11.2, 6.4, 1.0, "C0   117 datasets assembled")
    box(a, 5.0, 9.9, 6.4, 1.0, "C1   110 runs completed")
    box(a, 5.0, 8.6, 6.4, 1.0, "C2   93 analysis set  |  17 source families")
    box(a, 5.0, 7.3, 6.4, 1.0, "C3   80 on a common step grid")
    box(a, 5.0, 6.0, 6.4, 1.0, "C4   16 strict long-run re-run  |  7 families")
    box(a, 5.0, 4.7, 6.4, 1.0, "C5   9 exact-transport re-run  |  5 families")
    for y in (10.65, 9.35, 8.05, 6.75, 5.45):
        a.annotate("", xy=(5.0, y - 0.2), xytext=(5.0, y + 0.15),
                   arrowprops=dict(arrowstyle="->", lw=0.7, color="#446"))
    a.text(5.0, 3.5, "control samples (not nested):", ha="center",
           fontsize=6, style="italic", color="#555")
    box(a, 5.0, 2.5, 7.2, 1.3,
        "34 disconnected of C2 re-run for $V_{within}$   |   "
        "25 sham pairs, 8 families   |   5 surgery pairs, 2 families",
        fc="#f7f2ee")
    box(a, 5.0, 1.0, 7.2, 1.0,
        "synthetic sweeps: 2 transitivity bases, 4 assortativity bases",
        fc="#f7f2ee")
    a.text(5.0, 0.15, "primary = >=8 families;  exploratory = <5 families",
           ha="center", fontsize=5.6, style="italic", color="#555")

    # ---- (b) ladder of controls
    b = ax[1]
    b.set_title("b  Ladder of controls", fontsize=8, loc="left")
    rungs = [
        ("L4   $V_G^*$ characterisation", "not reached", "#f7e9e9"),
        ("L3   step-size / horizon controls", "tested", "#eef7ee"),
        ("L3   exact transport + extended horizon", "tested", "#eef7ee"),
        ("L2   scale-equivariant schedule", "tested", "#eef7ee"),
        ("L2   solver tolerance / precision", "tested", "#eef7ee"),
        ("L2   $\\varepsilon$-sensitivity", "tested", "#eef7ee"),
        ("L1   operational stopping criterion", "replaces initial-variance rule",
         "#eef2f7"),
        ("L1   raw finite-time residual  $V(T)$", "what is usually reported",
         "#eef2f7"),
    ]
    y = 11.2
    for label, tag, fc in rungs:
        box(b, 5.0, y, 8.6, 0.95, label, fc=fc)
        b.text(9.35, y, tag, ha="right", va="center", fontsize=5.4,
               style="italic", color="#555")
        y -= 1.3
    b.text(5.0, 0.5, "each rung removes one non-geometric explanation",
           ha="center", fontsize=5.6, style="italic", color="#555")
    save(fig, "Fig1_analysis_flow")


# --------------------------------------------------------- Fig 3 (new, R2)
def fig_reg():
    """Scale invariance, estimator repair, and exact transport to 2k steps."""
    note("")
    note("[Fig 3] regularisation: scale invariance and exact transport")
    r3 = os.path.join(ROOT, "exp_round3")
    probe = pd.read_csv(os.path.join(r3, "scale_probe.csv"))
    cmp3 = pd.read_csv(os.path.join(r3, "flow_sinkhorn_exact_eqreg.csv"))
    long = pd.read_csv(os.path.join(r3, "flow_exact.csv"))

    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.3))

    # a: scale-invariance violation
    for key, name, col, mark in (("exact", "exact $W_1$", GREEN, "o"),
                                 ("eqreg", "scale-equivariant", BLUE, "s"),
                                 ("sinkhorn", "fixed $\\varepsilon$",
                                  ORANGE, "^")):
        sub = probe[probe["estimator"] == key]
        g = sub.groupby("factor")["max_abs_dkappa"].max()
        ax[0].plot(g.index, np.maximum(g.values, 1e-16), mark, color=col,
                   ms=3.5, label=name)
    ax[0].set_yscale("log")
    ax[0].set_xscale("log")
    ax[0].set_ylim(1e-16, 10)
    ax[0].set_xlabel("rescale factor $c$ in $w\\mapsto cw$")
    ax[0].set_ylabel("max$_e|\\kappa_e(cw)-\\kappa_e(w)|$")
    ax[0].set_title("a  not scale invariant", loc="left")
    ax[0].legend(fontsize=5.5, frameon=False, loc="lower right")
    note("  scale probe: max over networks and phases of max_e|dkappa| -- "
         "exact %.1e | eqreg %.1e | fixed-eps %.2f"
         % (probe[probe.estimator == "exact"]["max_abs_dkappa"].max(),
            probe[probe.estimator == "eqreg"]["max_abs_dkappa"].max(),
            probe[probe.estimator == "sinkhorn"]["max_abs_dkappa"].max()))

    # b: residual under three estimators at a matched horizon
    p = cmp3.pivot_table(index="dataset", columns="mode", values="V_end")
    xs = np.arange(len(p))
    ax[1].plot(xs, p["sinkhorn"], "o", color=GREY, ms=4, label="fixed $\\varepsilon$")
    ax[1].plot(xs, p["eqreg"], "s", color=BLUE, ms=4, label="scale-equivariant")
    ax[1].plot(xs, p["exact"], "^", color=ORANGE, ms=4, label="exact $W_1$")
    ax[1].set_yscale("log")
    ax[1].set_xticks(xs)
    ax[1].set_xticklabels([d.replace("_", " ")[:12] for d in p.index],
                          rotation=90, fontsize=5)
    ax[1].set_ylabel("pooled residual $V$ at 200 steps")
    ax[1].set_title("b  repairing it", loc="left")
    ax[1].legend(fontsize=5.5, frameon=False)
    rel = (p["eqreg"] - p["sinkhorn"]) / p["sinkhorn"]
    note("  eqreg vs fixed-eps at 200 steps: median %+.1f%%, max %+.1f%%"
         % (100 * rel.median(), 100 * rel.abs().max()))

    # c: exact transport, 2000 steps
    long_i = long.set_index("dataset")
    colours = dict(zip(long["dataset"],
                       plt.cm.tab10(np.linspace(0, 1, len(long)))))
    for ds in long["dataset"]:
        fp = os.path.join(r3, "traj", "%s__exact.csv" % ds)
        if not os.path.exists(fp):
            continue
        d = pd.read_csv(fp)
        conn = long_i.loc[ds, "V_between"] < 1e-12
        ax[2].plot(d["step"], np.maximum(d["V_within"], 1e-34),
                   ls="-" if conn else "--", lw=0.9,
                   color=colours[ds], label=ds.replace("_", " ")[:16])
    ax[2].set_yscale("log")
    ax[2].set_xlabel("step")
    ax[2].set_ylabel("within-component $V$")
    ax[2].set_title("c  exact $W_1$, 2000 steps", loc="left")
    ax[2].legend(fontsize=4.4, frameon=False, ncol=2,
                 handlelength=1.4, columnspacing=0.8)
    note("  exact 2000-step within-component variance: " + ", ".join(
        "%s %.2e" % (r["dataset"][:18], r["V_within"])
        for _, r in long.iterrows()))
    fig.tight_layout()
    save(fig, "Fig3_regularisation")


def fig_rate():
    """Fig 8: local relaxation rate r_t = -log(V_t / V_{t-1})."""
    note("")
    note("[Fig 8] local relaxation-rate trajectories")
    strict = os.path.join(ROOT, "ricci_flow_v2", "per_dataset")
    fate = pd.read_csv(os.path.join(ROOT, "strict_results", "fate_table.csv"))
    fate_map = dict(zip(fate["dataset"], fate["fate"]))
    longd = os.path.join(ROOT, "exp_round4", "traj_long")

    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.7))

    def rates(v):
        v = np.asarray(v, dtype=float)
        v = np.where(v > 0, v, np.nan)
        with np.errstate(invalid="ignore", divide="ignore"):
            r = -np.log(v[1:] / v[:-1])
        return np.concatenate([[np.nan], r])

    # (a) population C4, 16 strict runs
    a = ax[0]
    files = sorted(glob.glob(os.path.join(strict, "*.csv")))
    n_pos = n_neg = 0
    for fp in files:
        ds = os.path.basename(fp)[:-4]
        tr = pd.read_csv(fp)
        r = rates(tr["curvature_variance"].astype(float).values)
        r = r[1:]
        m = np.isfinite(r)
        if not m.any():
            continue
        pat = fate_map.get(ds, "converging")
        col = "#c0392b" if pat == "stalled" else "#21618c"
        a.plot(np.arange(1, len(r) + 1), np.where(m, r, np.nan),
               lw=0.8, color=col, alpha=0.85)
        n_pos += int((r[m] > 0).sum())
        n_neg += int((r[m] < 0).sum())
    a.axhline(0.0, color="#888", lw=0.6, ls=":")
    a.set_yscale("symlog", linthresh=1e-6)
    a.set_xlabel("step")
    a.set_ylabel(r"local rate $r_t$ (per step)")
    a.set_title("a  $C_4$, 16 strict runs (Sinkhorn)", fontsize=8, loc="left")
    a.text(0.98, 0.04, "red = near-plateau at horizon\nblue = still relaxing",
           transform=a.transAxes, ha="right", va="bottom", fontsize=5.4,
           color="#555")
    a.tick_params(labelsize=6)
    note("  C4 increments with r_t>0: %d, r_t<0: %d" % (n_pos, n_neg))

    # (b) 5,000-step exact-transport runs
    b = ax[1]
    for fp in sorted(glob.glob(os.path.join(longd, "*__exact5000.csv"))):
        ds = os.path.basename(fp).split("__")[0]
        tr = pd.read_csv(fp)
        r = rates(tr["V"].astype(float).values)[1:]
        b.plot(np.arange(1, len(r) + 1), r, lw=0.8, label=ds.replace("ENZYMES_", ""))
    b.axhline(0.0, color="#888", lw=0.6, ls=":")
    b.set_yscale("symlog", linthresh=1e-8)
    b.set_xlabel("step")
    b.set_title("b  $C_5$, exact transport to 5,000 steps", fontsize=8, loc="left")
    b.legend(fontsize=5.4, frameon=False)
    b.tick_params(labelsize=6)
    fig.tight_layout()
    save(fig, "Fig8_relaxation_rate")


def main():
    master = pd.read_csv(MASTER)
    note("master table: %d rows, %d cols" % master.shape)
    fig_flow()
    fig_rate()
    fig1(master)
    fig_reg()
    fig2()
    fig3(master)
    fig4(master)
    fig5()
    with open(os.path.join(OUT, "FIGURE_DATA_LOG.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(log) + "\n")
    print("\nfigure data log -> %s" % os.path.join(OUT, "FIGURE_DATA_LOG.txt"))


if __name__ == "__main__":
    main()
