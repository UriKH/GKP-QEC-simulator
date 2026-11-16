"""Run experiments and save rows. Codes and decoders are cached, so sweeping the noise is cheap."""

import csv
import json
import time
from pathlib import Path

import numpy as np

from cvqec.circuits.teleportation import effective_channel
from cvqec.codes import make_code
from cvqec.decoders import make_decoder
from cvqec.engines.effective import run_logical
from cvqec.experiments.config import Experiment
from cvqec.metrics import compute

_CACHE: dict = {}


def _code_and_decoder(exp: Experiment):
    key = (exp.code, json.dumps(exp.code_params, sort_keys=True), exp.decoder,
           json.dumps(exp.decoder_params, sort_keys=True))
    if key not in _CACHE:
        code = make_code(exp.code, **exp.code_params)
        _CACHE[key] = (code, make_decoder(exp.decoder, code, **exp.decoder_params))
    return _CACHE[key]


def run_experiment(exp: Experiment) -> dict:
    t0 = time.time()
    code, decoder = _code_and_decoder(exp)
    result = run_logical(code, effective_channel(exp.noise), decoder, exp.N)
    row = exp.flat()
    row["decoder"] = decoder.name
    row.update(compute(exp.metrics, result))
    row["seconds"] = round(time.time() - t0, 3)
    return row


def run_all(exps: list[Experiment], progress: bool = True) -> list[dict]:
    rows = []
    for i, exp in enumerate(exps, 1):
        rows.append(run_experiment(exp))
        if progress:
            print(f"[{i}/{len(exps)}] {rows[-1]['label']}  F_avg = {rows[-1].get('average_fidelity', np.nan):.5f}")
    return rows


def save_rows(rows: list[dict], path: str | Path) -> Path:
    """CSV (union of all keys) next to a JSON copy, so the exact configuration is always kept."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    path.with_suffix(".json").write_text(json.dumps(rows, indent=1, default=float))
    return path


def load_rows(path: str | Path) -> list[dict]:
    return json.loads(Path(path).with_suffix(".json").read_text())


def match_mean_photon(make, target: float, lo: float, hi: float, N: int, tol: float = 1e-3) -> float:
    """Bisect a code parameter so its n-bar equals `target` (fair comparisons at equal energy)."""
    for _ in range(60):
        mid = (lo + hi) / 2
        if make(mid).mean_photon(N) < target:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return (lo + hi) / 2
