"""Command line: python -m cvqec.experiments <preset|all> [--quick] [--out results]

Writes results/<preset>.csv (+ .json) and results/<preset>.png.
"""

import argparse
from pathlib import Path

from cvqec.experiments.presets import PRESETS
from cvqec.experiments.runner import run_all, save_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("preset", choices=[*PRESETS, "all"])
    ap.add_argument("--quick", action="store_true", help="coarser sweeps, for a fast first look")
    ap.add_argument("--out", default="results")
    args = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    from cvqec.viz.sweeps import plot_sweep

    for name in PRESETS if args.preset == "all" else [args.preset]:
        for title, x, exps in PRESETS[name](quick=args.quick):
            rows = run_all(exps)
            path = save_rows(rows, Path(args.out) / f"{name}.csv")
            fig = plot_sweep(rows, x, title=title)
            fig.savefig(path.with_suffix(".png"), dpi=150, bbox_inches="tight")
            print(f"saved {path} and {path.with_suffix('.png')}")


if __name__ == "__main__":
    main()
