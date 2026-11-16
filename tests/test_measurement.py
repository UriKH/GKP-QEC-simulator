"""Homodyne amplitudes, exact displacement matrix elements and Gaussian unitaries."""

import numpy as np
import pytest
import qutip as qt

from cvqec.conventions import quad_p, quad_q
from cvqec.ops.gates import displacement_exact, gaussian_unitary
from cvqec.ops.measurement import homodyne_distribution, quadrature_amplitudes

N = 40


@pytest.mark.parametrize("theta", [0.0, np.pi / 2, 0.7])
def test_homodyne_of_coherent_state(theta):
    # Course p.54: measuring q_theta on |alpha> gives a Gaussian at cos(theta) q0 + sin(theta) p0, variance 1/2.
    alpha = 0.8 - 0.5j
    q0, p0 = np.sqrt(2) * alpha.real, np.sqrt(2) * alpha.imag
    x = np.linspace(-8, 8, 1601)
    P = homodyne_distribution(qt.coherent(N, alpha, method="analytic"), x, theta)
    dx = x[1] - x[0]
    mean = np.sum(x * P) * dx
    assert np.sum(P) * dx == pytest.approx(1.0, abs=1e-10)
    assert mean == pytest.approx(np.cos(theta) * q0 + np.sin(theta) * p0, abs=1e-8)
    assert np.sum((x - mean) ** 2 * P) * dx == pytest.approx(0.5, abs=1e-8)


def test_momentum_amplitudes_are_fourier_transform():
    # <p|n> = (-i)^n psi_n(p): the theta = pi/2 homodyne basis is the Fourier transform of the q basis.
    x = np.linspace(-12, 12, 1201)
    dx = x[1] - x[0]
    Aq, Ap = quadrature_amplitudes(6, x), quadrature_amplitudes(6, x, np.pi / 2)
    F = np.exp(-1j * np.outer(x, x)) * dx / np.sqrt(2 * np.pi)
    np.testing.assert_allclose(F @ Aq, Ap, atol=1e-8)


@pytest.mark.parametrize("alpha", [0.3 + 0.2j, -1.1 + 0.7j, 2.5j, 0])
def test_exact_displacement_matches_padded_expm(alpha):
    padded = qt.displace(300, alpha).full()[:30, :30]
    np.testing.assert_allclose(displacement_exact(alpha, 30), padded, atol=1e-12)


def test_exact_displacement_large_amplitude_against_mpmath():
    mp = pytest.importorskip("mpmath")
    mp.mp.dps = 50
    a = 23 * np.exp(0.4j)                       # |alpha|^2 = 529: hopeless for a truncated expm
    D = displacement_exact(a, 140)
    am = mp.mpc(a.real, a.imag)
    for m, n in [(139, 139), (120, 130), (130, 120), (100, 50), (0, 139)]:
        lo, hi = min(m, n), max(m, n)
        ref = (mp.sqrt(mp.factorial(lo) / mp.factorial(hi)) * mp.e ** (-abs(am) ** 2 / 2)
               * mp.laguerre(lo, hi - lo, abs(am) ** 2) * (am if m >= n else -mp.conj(am)) ** (hi - lo))
        assert abs(complex(ref) - D[m, n]) < 1e-12 * max(1.0, abs(complex(ref)))


def test_exact_displacement_group_law_with_headroom():
    # D(-a) D(a) = 1 on low levels, provided the intermediate sum has room (column n spreads to ~|a|^2 + ...).
    a = 15 + 9j
    D1, D2 = displacement_exact(a, 600), displacement_exact(-a, 600)
    np.testing.assert_allclose((D2 @ D1)[:20, :20], np.eye(20), atol=1e-10)


def test_gaussian_unitary_implements_symplectic_matrix():
    M = np.array([[1.5, 0.3], [0.2, (1 + 0.3 * 0.2) / 1.5]])
    U = gaussian_unitary(M, 80)
    c = qt.coherent(80, 0.4 + 0.3j, method="analytic")
    x = np.array([qt.expect(quad_q(80), c), qt.expect(quad_p(80), c)])
    y = np.array([qt.expect(quad_q(80), U * c), qt.expect(quad_p(80), U * c)])
    np.testing.assert_allclose(y, M @ x, atol=1e-8)
    with pytest.raises(ValueError):
        gaussian_unitary(np.diag([2.0, 1.0]), 10)
