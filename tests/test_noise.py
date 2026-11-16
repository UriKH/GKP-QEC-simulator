"""Noise channels: CPTP checks, closed-form actions, and Kraus/exact maps vs qutip.mesolve (primer section I)."""

import numpy as np
import pytest
import qutip as qt

from cvqec.conventions import quad_p, quad_q
from cvqec.noise import Compose, Dephasing, GaussianDisplacement, Identity, PhotonLoss, make_channel

N = 30


def random_state(n=N, n_max=8, seed=1):
    rng = np.random.default_rng(seed)
    v = np.zeros(n, complex)
    v[:n_max] = rng.normal(size=n_max) + 1j * rng.normal(size=n_max)
    return qt.ket2dm(qt.Qobj(v).unit())


def mesolve_final(channel, rho):
    c_ops, t = channel.lindblad(rho.shape[0])
    H = 0 * qt.qeye(rho.shape[0])
    return qt.mesolve(H, rho, [0, t], c_ops, options={"atol": 1e-12, "rtol": 1e-10}).states[-1]


CHANNELS = [PhotonLoss(0.8), Dephasing(0.3), GaussianDisplacement(0.1, 0.05)]


@pytest.mark.parametrize("ch", CHANNELS, ids=str)
def test_trace_preserving_and_positive(ch):
    out = ch.apply(random_state())
    assert out.tr() == pytest.approx(1.0, abs=1e-8)
    assert np.linalg.eigvalsh(out.full()).min() > -1e-9


@pytest.mark.parametrize("ch", CHANNELS, ids=str)
def test_matches_lindblad_evolution(ch):
    # The same physics written two ways: discrete map vs continuous-time master equation.
    rho = random_state()
    assert (ch.apply(rho) - mesolve_final(ch, rho)).norm() < 1e-6


@pytest.mark.parametrize("ch", CHANNELS, ids=str)
def test_adjoint_definition(ch):
    # Tr[X E(rho)] = Tr[E^dag(X) rho] defines the Heisenberg-picture map used by the Petz decoder.
    rho, X = random_state(seed=2), qt.Qobj(np.diag(np.arange(N, dtype=float) % 3))
    assert (X * ch.apply(rho)).tr() == pytest.approx((ch.adjoint(X) * rho).tr(), abs=1e-8)


def test_loss_shrinks_coherent_amplitude():
    # Beam-splitter picture: |alpha> -> |sqrt(eta) alpha>, still pure.
    eta, alpha = 0.7, 1.2 + 0.4j
    out = PhotonLoss(eta).apply(qt.ket2dm(qt.coherent(N, alpha, method="analytic")))
    target = qt.coherent(N, np.sqrt(eta) * alpha, method="analytic")
    assert qt.expect(out, target) == pytest.approx(1.0, abs=1e-10)


def test_loss_composes_multiplicatively():
    rho = random_state()
    two = PhotonLoss(0.9).apply(PhotonLoss(0.8).apply(rho))
    assert (two - PhotonLoss(0.72).apply(rho)).norm() < 1e-10


def test_dephasing_kraus_matches_exact_mask():
    ch, rho = Dephasing(0.2), random_state()
    kraus_out = sum((K * rho * K.dag() for K in ch.kraus(N)), 0 * rho)
    assert (kraus_out - ch.apply(rho)).norm() < 1e-8


def test_displacement_noise_adds_variance():
    # The defining property: Var q grows by var_q, Var p by var_p, means unchanged.
    ch = GaussianDisplacement(0.15, 0.05)
    rho = qt.ket2dm(qt.coherent(N, 0.5 + 0.2j))
    out = ch.apply(rho)
    var = lambda op, r: qt.expect(op * op, r) - qt.expect(op, r) ** 2
    assert var(quad_q(N), out) == pytest.approx(0.5 + 0.15, abs=1e-6)
    assert var(quad_p(N), out) == pytest.approx(0.5 + 0.05, abs=1e-6)
    assert qt.expect(quad_q(N), out) == pytest.approx(qt.expect(quad_q(N), rho), abs=1e-8)


def test_compose_order_and_registry():
    rho = random_state()
    a, b = PhotonLoss(0.8), Dephasing(0.3)
    assert (Compose([a, Identity(), b]).apply(rho) - b.apply(a.apply(rho))).norm() < 1e-12
    assert make_channel("loss", eta=0.5).eta == 0.5
