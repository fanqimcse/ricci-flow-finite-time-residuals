"""Report the weight spread reached by the real (no-surgery) long run.

Needed to choose a surgery threshold that actually triggers: if the flow
renormalises by the mean weight every step, the weights may never collapse,
in which case "surgery" is vacuous for this flow and the referee objection
has a different answer.
"""

import os
import sys

import numpy as np
import pandas as pd

SOURCE_DIR = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(SOURCE_DIR, "ricci_flow_longrun_results", "per_dataset")
OUT = os.path.join(SOURCE_DIR, "_weight_spread.txt")
NAMES = ["ENZYMES_g1", "ENZYMES_g103", "bn-mouse_visual-cortex_1",
         "reptilia-tortoise-network-pv", "mammalia-raccoon-proximity",
         "scc_rt_damascus", "eco-stmarks"]


def main():
    lines = []
    lines.append("{:<32s} {:>6s} {:>9s} {:>9s} {:>9s} {:>9s}".format(
        "dataset", "steps", "min_w_end", "min_w_ever", "max_w_ever",
        "min/mean_end"))
    for name in NAMES:
        path = os.path.join(REF, name + ".csv")
        if not os.path.exists(path):
            lines.append("{:<32s} (missing)".format(name))
            continue
        frame = pd.read_csv(path)
        mn = frame["min_weight"].astype(float).values
        mx = frame["max_weight"].astype(float).values
        med = frame["median_weight"].astype(float).values
        lines.append("{:<32s} {:>6d} {:>9.4f} {:>9.4f} {:>9.4f} {:>9.4f}".format(
            name[:32], len(frame), mn[-1], float(np.min(mn)),
            float(np.max(mx)), mn[-1] / med[-1]))
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
