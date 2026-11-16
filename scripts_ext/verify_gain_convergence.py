"""Supplement to verify_gain_channel.py: residuals shrink with the engine cutoff, and the check is sensitive."""
import json
from pathlib import Path
import numpy as np
import qutip as qt
from cvqec.circuits import TeleportationParams, teleportation_circuit
from cvqec.engines.fock import FockEngine
from cvqec.noise import Compose, GaussianDisplacement, PhotonLoss
from gain_channel import general_gain_channel, Amplifier

def td(a, b):
    return 0.5 * np.abs(np.linalg.eigvalsh((a - b).full())).sum()

out = []
for c in (dict(r1=0.5, r2=0.7, eta=0.95, eta_d=0.9, gain=1.15), dict(r1=0.6, r2=0.6, gain=0.9)):
    p = TeleportationParams(**c)
    for N in (20, 25, 30):
        psi = (qt.basis(N, 0) + 0.7j * qt.basis(N, 1)).unit()
        full = FockEngine(N).run(teleportation_circuit(p), {"in": psi})
        d = td(full, general_gain_channel(p).apply(qt.ket2dm(psi)))
        out.append(dict(c, N=N, trace_distance=d)); print(out[-1], flush=True)
# sensitivity: pretend the gain were ideal (unit gain, g=1) for the g=0.9 case
N = 25
p = TeleportationParams(0.6, 0.6, gain=0.9)
psi = (qt.basis(N, 0) + 0.7j * qt.basis(N, 1)).unit()
full = FockEngine(N).run(teleportation_circuit(p), {"in": psi})
wrong = GaussianDisplacement(np.exp(-1.2)).apply(qt.ket2dm(psi))
out.append({"sensitivity": "g=0.9 circuit vs unit-gain channel", "trace_distance": td(full, wrong)}); print(out[-1])
Path("../results_ext/verify_gain_convergence.json").write_text(json.dumps(out, indent=1, default=float))
