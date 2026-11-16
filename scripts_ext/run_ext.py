"""Extension experiments E1-E3 (stored in results_ext/, separate from the report data in results/).

    conda run -n qcsim python scripts_ext/run_ext.py            # E1, E2, E3 and the convergence check
    conda run -n qcsim python scripts_ext/run_ext.py e1 conv    # a subset

E1  dephasing complementarity: pure dephasing e^{-i phi n}, phi ~ N(0, sigma^2), no displacement noise.
E2  the same codes under the effective teleportation channel vs squeezing (eta = eta_d = 1, sigma_phi = 0).
E3  feed-forward gain mismatch g != sqrt(eta) for square GKP at 10 dB, eta = eta_d = 1 (gain_channel.py).
conv  the same points recomputed at a larger Fock cutoff.

Codes are at equal mean photon number: GKP Delta = 0.3 (n-bar ~ 5.07 for every lattice), 4-leg cat with
alpha bisected to the square code's n-bar (as cvqec.experiments.presets.code_lineup(include_cat=True) does),
bare Fock qubit {|0>, |1>} as the break-even reference. Decoders: each code's standard decoder
("auto": gkp_ideal / cat_parity / projection) and Petz (noise-aware benchmark) for every code.
Each row stores its full configuration (code, parameters, decoder, channel, N).
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

from cvqec.circuits.teleportation import TeleportationParams, effective_channel
from cvqec.codes import FourLegCat, GKPCode, make_code
from cvqec.conventions import r_from_db
from cvqec.decoders import make_decoder
from cvqec.engines.effective import run_logical
from cvqec.experiments.runner import match_mean_photon, save_rows
from cvqec.metrics import compute
from cvqec.noise import Dephasing
from gain_channel import general_gain_channel

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results_ext"
N = 140
N_CONV = 180
DELTA = 0.3
METRICS = ("average_fidelity", "entanglement_fidelity", "pauli_errors", "mean_photon")


def lineup():
    nbar = GKPCode("square", DELTA).mean_photon(N)
    alpha = round(match_mean_photon(lambda a: FourLegCat(a), nbar, 1.0, 4.0, N), 4)
    return [
        ("GKP square", "gkp", {"lattice": "square", "delta": DELTA}),
        ("GKP hexagonal", "gkp", {"lattice": "hexagonal", "delta": DELTA}),
        (f"4-leg cat (α={alpha:.2f})", "cat4", {"alpha": alpha}),
        ("Fock qubit (break-even)", "fock", {}),
    ]


_CODES = {}


def code_of(key, params):
    k = (key, json.dumps(params, sort_keys=True))
    if k not in _CODES:
        _CODES[k] = make_code(key, **params)
    return _CODES[k]


def run_point(label, key, params, decoder, channel, channel_cfg, n=N, extra=None):
    t0 = time.time()
    code = code_of(key, params)
    dec = make_decoder(decoder, code)
    res = run_logical(code, channel, dec, n)
    row = {"label": label, "code": key, "decoder": dec.name, "N": n}
    row.update({f"code.{k}": v for k, v in params.items()})
    row.update({f"channel.{k}": v for k, v in channel_cfg.items()})
    row["channel.repr"] = channel.name
    row.update(compute(METRICS, res))
    row["code_nbar"] = code.mean_photon(n)
    row.update(extra or {})
    row["seconds"] = round(time.time() - t0, 3)
    return row


def e1(n=N, sigmas=None):
    sigmas = np.round(np.arange(0.02, 0.401, 0.02), 4).tolist() if sigmas is None else sigmas
    rows = []
    for label, key, params in lineup():
        for dec in ("auto", "petz"):
            for s in sigmas:
                rows.append(run_point(label, key, params, dec, Dephasing(s), {"type": "dephasing", "sigma_phi": s}, n))
            print(f"E1 {label} {dec} done", flush=True)
    return rows


def e2(n=N, dbs=None):
    dbs = list(range(2, 17)) if dbs is None else dbs
    rows = []
    for label, key, params in lineup():
        for dec in ("auto", "petz"):
            for db in dbs:
                p = TeleportationParams.from_db(db)
                cfg = {"type": "teleportation", "db": db, **p.as_dict()}
                rows.append(run_point(label, key, params, dec, effective_channel(p), cfg, n))
            print(f"E2 {label} {dec} done", flush=True)
    return rows


def e3(n=N, gains=None, db=10.0):
    gains = np.round(np.arange(0.85, 1.151, 0.01), 4).tolist() if gains is None else gains
    rows = []
    for label, key, params in [lineup()[0], lineup()[3]]:
        for dec in ("auto", "petz"):
            for g in gains:
                p = TeleportationParams.from_db(db, gain=g)
                cfg = {"type": "teleportation_general_gain", "db": db, **p.as_dict()}
                rows.append(run_point(label, key, params, dec, general_gain_channel(p), cfg, n))
            print(f"E3 {label} {dec} done", flush=True)
    return rows


def conv():
    """Recompute the hardest points of E1-E3 at N = 140 and N = N_CONV and store both."""
    rows = []
    for n in (N, N_CONV):
        _CODES.clear()
        for r in e1(n, sigmas=[0.2, 0.4]):
            rows.append(dict(r, experiment="E1"))
        for r in e2(n, dbs=[2, 10, 16]):
            rows.append(dict(r, experiment="E2"))
        for r in e3(n, gains=[0.85, 1.0, 1.15]):
            rows.append(dict(r, experiment="E3"))
    return rows


if __name__ == "__main__":
    todo = sys.argv[1:] or ["e1", "e2", "e3", "conv"]
    for name in todo:
        t = time.time()
        rows = {"e1": e1, "e2": e2, "e3": e3, "conv": conv}[name]()
        save_rows(rows, OUT / f"{name}.csv")
        print(f"== saved {name} ({len(rows)} rows, {time.time() - t:.0f} s)", flush=True)
