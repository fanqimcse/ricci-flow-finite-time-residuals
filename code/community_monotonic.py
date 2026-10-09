# -*- coding: utf-8 -*-
"""Monotonicity addendum for experiment C.

With only 3 community levels per base graph, a Spearman rho inside a base is
barely interpretable, but ORDERING is: if the plateau rises monotonically as
community structure strengthens, that is 1 of 6 possible orderings of three
values, i.e. 1/6 under the null (1/3 if either direction counts). Across
independent base graphs the probability multiplies.
"""

import glob
import os

import numpy as np
import pandas as pd

SRC = r"D:\2026codes\2026821100data\synthetic_results"


def main():
    frames = [pd.read_csv(c)
              for c in sorted(glob.glob(os.path.join(SRC,
                                                     "community_b*.csv")))]
    d = pd.concat(frames, ignore_index=True)
    d = d.drop_duplicates(subset=["base", "target_mu", "rep"], keep="last")
    d = d[np.isfinite(d["var_ratio"].astype(float))]

    lines = []
    lines.append("EXPERIMENT C -- MONOTONICITY ADDENDUM")
    lines.append("=" * 56)
    lines.append("")
    lines.append("Rule: within a base graph the degree sequence is FIXED; the")
    lines.append("only thing that changes is the wiring. If the plateau rises")
    lines.append("monotonically as community structure strengthens (mu falls),")
    lines.append("that ordering has probability 1/6 under the null for three")
    lines.append("levels (1/3 if either direction is allowed).")
    lines.append("")

    n_incr = 0
    n_decr = 0
    bases = []
    for base, sub in d.groupby("base"):
        sub = sub.sort_values("target_mu", ascending=False)
        vals = sub["var_ratio"].astype(float).tolist()
        mus = sub["target_mu"].tolist()
        incr = all(vals[i] < vals[i + 1] for i in range(len(vals) - 1))
        decr = all(vals[i] > vals[i + 1] for i in range(len(vals) - 1))
        n_incr += int(incr)
        n_decr += int(decr)
        bases.append(len(vals))
        lines.append("  base {}  mu={}  plateau={}   {}".format(
            base,
            "[{}]".format(", ".join("{:.2f}".format(m) for m in mus)),
            "[{}]".format(", ".join("{:.5f}".format(v) for v in vals)),
            "INCREASING" if incr else ("decreasing" if decr else "non-monotone")))

    k = len(bases)
    lines.append("")
    lines.append("  monotone increasing in {}/{} base graphs".format(n_incr, k))
    lines.append("  monotone decreasing in {}/{} base graphs".format(n_decr, k))
    if k:
        p_one = (1.0 / 6.0) ** n_incr
        p_two = (1.0 / 3.0) ** n_incr
        lines.append("")
        lines.append("  under H0: P(all increasing by chance) = (1/6)^{} = "
                     "{:.4g}".format(n_incr, p_one))
        lines.append("            P(all monotone by chance)   = (1/3)^{} = "
                     "{:.4g}".format(n_incr, p_two))
    lines.append("")
    lines.append("CAVEAT that must travel with every number in experiment C:")
    lines.append("  all runs stopped at a MATCHED HORIZON of 95 steps")
    lines.append("  (status = horizon_reached, NOT saturated). The comparison")
    lines.append("  therefore shows that community structure slows and raises")
    lines.append("  the residual at equal flow time; it does not prove the")
    lines.append("  asymptotic plateau is higher. Synthetic values wouldkeep")
    lines.append("  falling if run longer, so the remaining gap to real")
    lines.append("  networks (which ARE saturated) is a LOWER bound on the")
    lines.append("  true gap.")

    with open(os.path.join(SRC, "COMMUNITY_MONOTONICITY.txt"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
