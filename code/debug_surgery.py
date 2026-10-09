"""Debug helper: run one dataset with surgery and print a full traceback."""

import os
import sys
import traceback

SOURCE_DIR = os.path.dirname(os.path.abspath(__file__))
if SOURCE_DIR not in sys.path:
    sys.path.insert(0, SOURCE_DIR)

import ricci_flow_program_100data as R                    # noqa: E402
from ricci_flow_longrun import compute_curvature, build_base_args  # noqa: E402
from surgery_check import surgery_pass, CFG               # noqa: E402

DATA_DIR = r"D:\100data"
OUT = os.path.join(SOURCE_DIR, "_surg_debug.txt")


def main():
    dataset = sys.argv[1] if len(sys.argv) > 1 else "ENZYMES_g1"
    tol = float(sys.argv[2]) if len(sys.argv) > 2 else 0.30
    steps = int(sys.argv[3]) if len(sys.argv) > 3 else 150

    lines = []
    base = build_base_args(["--data_dir", DATA_DIR])
    files = R.discover_graph_files(DATA_DIR)
    path = [f for f in files
            if R.safe_dataset_id(f, DATA_DIR) == dataset][0]
    graph = R.load_graph(path, base)
    workspace = R.CurvatureWorkspace(
        graph, alpha=CFG["alpha"], backend=base.distance_backend,
        max_distance_mb=base.max_distance_mb)

    var0 = None
    try:
        for step in range(steps):
            mean_k, var_k, edge_k, _, _ = compute_curvature(
                workspace, CFG, CFG["solver"])
            if step == 0:
                var0 = float(var_k)
            R.ricci_flow_step(graph, edge_k, dt=CFG["dt"])
            removed = surgery_pass(graph, tol)
            if removed:
                lines.append("step {} removed {} ratio={:.4f} edges_left={}"
                             .format(step, removed, float(var_k) / var0,
                                     graph.number_of_edges()))
                workspace = R.CurvatureWorkspace(
                    graph, alpha=CFG["alpha"],
                    backend=base.distance_backend,
                    max_distance_mb=base.max_distance_mb)
        lines.append("completed {} steps; edges left {}; final ratio {:.5f}"
                     .format(steps, graph.number_of_edges(),
                             float(var_k) / var0))
    except Exception:
        lines.append("EXCEPTION:")
        lines.append(traceback.format_exc())

    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
