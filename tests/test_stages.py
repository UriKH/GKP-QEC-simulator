"""The stage-by-stage figure (report Section 2): every panel must be a state the simulated circuit really produces."""

import numpy as np
import pytest
import qutip as qt

from cvqec.circuits import TeleportationParams, bob_before_feedforward, effective_channel, stage_params, teleportation_circuit
from cvqec.codes import GKPCode
from cvqec.conventions import quad_p, quad_q
from cvqec.decoders import GKPIdealDecoder
from cvqec.engines.fock import FockEngine
from cvqec.noise import PhotonLoss
from cvqec.ops.gates import beam_splitter, squeezer
from cvqec.viz.wigner import wigner

N = 25
P = TeleportationParams(0.6, 0.4, eta=0.85, eta_d=0.95, sigma_phi=0.15)


def trace_distance(a, b):
    return 0.5 * np.abs(np.linalg.eigvalsh((a - b).full())).sum()


def test_stages_switch_imperfections_on_in_order():
    labels = [label for label, _ in stage_params(P)]
    assert labels == ["shift noise", "+ loss", "+ dephasing"]
    first, _, last = (s for _, s in stage_params(P))
    assert (first.eta, first.sigma_phi, first.eta_d) == (1.0, 0.0, P.eta_d)   # detectors belong to the shift noise
    assert last == P


@pytest.mark.parametrize("k", [0, 1, 2])
def test_each_stage_equals_the_full_circuit(k):
    # The panel "after stage k" is computed with the single-mode channel; here the same stage is run through
    # the full three-mode simulation (squeezers, beam splitters, homodynes, feed-forward) on a small input.
    _, params = stage_params(P)[k]
    psi = (qt.basis(N, 0) + 0.6 * qt.basis(N, 1) - 0.4j * qt.basis(N, 2)).unit()
    full = FockEngine(N).run(teleportation_circuit(params), {"in": psi})
    assert trace_distance(full, effective_channel(params).apply(qt.ket2dm(psi))) < 1e-5


def test_bob_before_feedforward_matches_partial_trace():
    # Two-mode EPR source, loss on b, trace out a: Bob's state must equal the closed form used in the figure.
    n = 30
    vac = qt.tensor(qt.basis(n, 0), qt.basis(n, 0))
    epr = qt.ket2dm(beam_splitter(np.pi / 4, n) * qt.tensor(squeezer(P.r1, n), squeezer(-P.r2, n)) * vac)
    lossy = sum((qt.tensor(qt.qeye(n), K) * epr * qt.tensor(qt.qeye(n), K).dag() for K in PhotonLoss(P.eta).kraus(n)),
                0 * epr)
    assert trace_distance(lossy.ptrace(1), bob_before_feedforward(P, n)) < 1e-6


def test_bob_alone_is_thermal_for_symmetric_squeezing():
    # r1 = r2 = r: Bob holds a thermal state with nbar = eta sinh^2 r, whatever the input (no information).
    p, n = TeleportationParams(0.5, 0.5, eta=0.9), 40
    assert trace_distance(bob_before_feedforward(p, n), qt.thermal_dm(n, 0.9 * np.sinh(0.5) ** 2)) < 1e-10


def test_shift_noise_stage_adds_exactly_the_derived_variances():
    # On the figure's GKP input the first stage must add sigma_q^2 and sigma_p^2 to Var q and Var p.
    n = 140
    code = GKPCode("square", 0.3)
    plus = qt.ket2dm(code.encode(qt.Qobj(np.ones(2) / np.sqrt(2)), n))
    _, params = stage_params(TeleportationParams.from_db(10, eta=0.9, sigma_phi=0.1))[0]
    out = effective_channel(params).apply(plus)
    var = lambda op, r: qt.expect(op * op, r) - qt.expect(op, r) ** 2
    vq, vp = params.noise_variances()
    assert var(quad_q(n), out) - var(quad_q(n), plus) == pytest.approx(vq, rel=1e-6)
    assert var(quad_p(n), out) - var(quad_p(n), plus) == pytest.approx(vp, rel=1e-6)


def test_figure_window_holds_the_whole_wigner_function():
    # Every panel uses extent 8: the Wigner function must integrate to 1 inside it, or the plot hides part of the
    # state. (At extent 6.5 a 10.5 dB GKP state loses 2% of its weight outside the window; this test caught that.)
    n = 140
    code = GKPCode("square", 0.3)
    plus = qt.ket2dm(code.encode(qt.Qobj(np.ones(2) / np.sqrt(2)), n))
    full = effective_channel(TeleportationParams.from_db(10, eta=0.9, sigma_phi=0.1)).apply(plus)
    for rho in (plus, full, bob_before_feedforward(TeleportationParams.from_db(10, eta=0.9), n)):
        xs, W = wigner(rho, extent=8.0, points=201)
        assert W.sum() * (xs[1] - xs[0]) ** 2 == pytest.approx(1.0, abs=1e-3)


def test_reencoded_panel_is_the_decoded_qubit():
    # The last panel shows the decoded logical state put back into the code (V rho_L V^dag) for display only.
    # For a noiseless codeword this round trip must return the codeword itself.
    n = 140
    code = GKPCode("square", 0.3)
    V = code.encoder(n)
    plus = code.encode(qt.Qobj(np.ones(2) / np.sqrt(2)), n)
    rho_L = qt.Qobj(GKPIdealDecoder(code).decode(qt.ket2dm(plus)))
    assert qt.expect(V * rho_L * V.dag(), plus) > 0.99
