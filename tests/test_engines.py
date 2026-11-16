"""The full three-mode simulation vs analytic results and vs the derived single-mode channel."""

import numpy as np
import pytest
import qutip as qt

from cvqec.circuits import ChannelOp, FeedForward, Measure, TeleportationParams, effective_channel, teleportation_circuit
from cvqec.codes import FockQubit
from cvqec.decoders import make_decoder
from cvqec.engines.effective import run_logical
from cvqec.engines.fock import FockEngine, defer_channels
from cvqec.noise import Identity

N = 25


def trace_distance(a, b):
    return 0.5 * np.abs(np.linalg.eigvalsh((a - b).full())).sum()


def test_channels_are_deferred_past_other_modes_only():
    ops = defer_channels(teleportation_circuit(TeleportationParams(0.5, 0.5, eta=0.9, sigma_phi=0.1)).ops)
    kinds = [type(op).__name__ for op in ops]
    loss = kinds.index("ChannelOp")
    # Bob's loss moves behind both of Alice's measurements but stays before Bob's feed-forward.
    assert loss > max(i for i, k in enumerate(kinds) if k == "Measure")
    assert loss < kinds.index("FeedForward")


@pytest.mark.parametrize("r", [0.3, 0.5])
def test_coherent_state_teleportation_fidelity(r):
    # Primer section G: unit-gain teleportation of |alpha> has F = 1 / (1 + e^{-2r}); F -> 1/2 without entanglement.
    alpha = 0.4 + 0.3j
    eng = FockEngine(N)
    rho = eng.run(teleportation_circuit(TeleportationParams(r, r)), {"in": qt.coherent(N, alpha, method="analytic")})
    assert eng.last_trace == pytest.approx(1.0, abs=1e-6)
    assert qt.expect(rho, qt.coherent(N, alpha, method="analytic")) == pytest.approx(1 / (1 + np.exp(-2 * r)), abs=1e-6)


def test_full_circuit_equals_effective_channel():
    # The central consistency check of the project: every imperfection at once, simulated mode by mode,
    # equals Dephasing o Loss o Displacement(e^{-2r} + (1 - eta_d)/eta_d) on a single mode.
    p = TeleportationParams(0.6, 0.4, eta=0.9, eta_d=0.95, sigma_phi=0.15)
    psi = (qt.basis(N, 0) + 0.7j * qt.basis(N, 1)).unit()
    full = FockEngine(N).run(teleportation_circuit(p), {"in": psi})
    eff = effective_channel(p).apply(qt.ket2dm(psi))
    assert trace_distance(full, eff) < 1e-5


def test_effective_channel_parameters():
    p = TeleportationParams.from_db(10, 7, eta_d=0.8)
    var_q, var_p = p.noise_variances()
    assert var_q == pytest.approx(0.1 + 0.25)
    assert var_p == pytest.approx(10 ** -0.7 + 0.25)
    with pytest.raises(ValueError):
        effective_channel(TeleportationParams(0.5, 0.5, eta=0.9, gain=1.0))


def test_circuit_component_table():
    c = teleportation_circuit(TeleportationParams(0.5, 0.5, eta=0.9, eta_d=0.9, sigma_phi=0.1))
    assert len(c.ops) == 11   # + the encoder (input) and the Steane round (decoder) = 13-14 components
    assert sum(isinstance(op, Measure) for op in c.ops) == 2
    assert sum(isinstance(op, FeedForward) for op in c.ops) == 2
    assert "B(pi/4) EPR source" in c.describe()


def test_noiseless_logical_channel_is_identity():
    code = FockQubit()
    res = run_logical(code, Identity(), make_decoder("projection", code), 10)
    assert res.logical.average_fidelity() == pytest.approx(1.0)
