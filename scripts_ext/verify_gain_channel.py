"""Check gain_channel.general_gain_channel against the full three-mode Fock simulation of the circuit.

    conda run -n qcsim python scripts_ext/verify_gain_channel.py
Same method as tests/test_engines.py::test_full_circuit_equals_effective_channel (N = 25 per mode), but with
feed-forward gain g != sqrt(eta), including g < 1 (loss branch) and g > 1 (amplifier branch).
Writes results_ext/verify_gain_channel.json.
"""
import json
from pathlib import Path

import numpy as np
import qutip as qt

from cvqec.circuits import TeleportationParams, effective_channel, teleportation_circuit
from cvqec.engines.fock import FockEngine
from gain_channel import general_gain_channel

N = 25
OUT = Path(__file__).resolve().parent.parent / "results_ext"


def td(a, b):
    return 0.5 * np.abs(np.linalg.eigvalsh((a - b).full())).sum()


cases = [
    dict(r1=0.6, r2=0.6, gain=0.9),
    dict(r1=0.6, r2=0.6, gain=1.1),
    dict(r1=0.6, r2=0.4, eta=0.9, eta_d=0.95, sigma_phi=0.15, gain=0.8),
    dict(r1=0.5, r2=0.7, eta=0.95, eta_d=0.9, gain=1.15),
    dict(r1=0.6, r2=0.4, eta=0.9, eta_d=0.95, sigma_phi=0.15, gain=None),   # g = sqrt(eta): vs library too
]
psis = [(qt.basis(N, 0) + 0.7j * qt.basis(N, 1)).unit(), qt.coherent(N, 0.5 - 0.3j, method="analytic")]
rows = []
for c in cases:
    p = TeleportationParams(**c)
    for k, psi in enumerate(psis):
        eng = FockEngine(N)
        full = eng.run(teleportation_circuit(p), {"in": psi})
        mine = general_gain_channel(p).apply(qt.ket2dm(psi))
        row = dict(c, g=p.g, input=["fock01", "coherent"][k], trace_distance=td(full, mine),
                   engine_trace=eng.last_trace, out_trace=float(mine.tr().real))
        if c["gain"] is None:
            row["td_vs_library"] = td(mine, effective_channel(p).apply(qt.ket2dm(psi)))
        rows.append(row)
        print(row)
OUT.mkdir(exist_ok=True)
(OUT / "verify_gain_channel.json").write_text(json.dumps(rows, indent=1, default=float))
worst = max(r["trace_distance"] for r in rows)
print("worst trace distance:", worst)
assert worst < 1e-5, worst
