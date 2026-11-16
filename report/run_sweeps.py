"""Stage 1 of the report: run every preset at full resolution and store rows in results/.

    python report/run_sweeps.py            # all presets (~10 min)
    python report/run_sweeps.py loss       # one preset
"""

import sys
from pathlib import Path

from cvqec.experiments.presets import PRESETS
from cvqec.experiments.runner import run_all, save_rows

RESULTS = Path(__file__).resolve().parent.parent / "results"

if __name__ == "__main__":
    for name in sys.argv[1:] or PRESETS:
        for _, _, exps in PRESETS[name]():
            save_rows(run_all(exps), RESULTS / f"{name}.csv")
            print(f"== saved {name}", flush=True)
