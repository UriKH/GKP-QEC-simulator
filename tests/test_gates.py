"""Physics checks of each gate against the Heisenberg actions in the course book (p.50-53)."""

import numpy as np
import pytest
import qutip as qt

from cvqec.conventions import quad_p, quad_q
from cvqec.ops.gates import (
    beam_splitter, displacement, on_modes, phase_shifter, squeezer, sum_gate, x_shift, z_shift,
)
from cvqec.truncation import TruncationWarning, check_truncation, tail_weight

N = 30
N2 = 20  # per-mode truncation for two-mode tests (400-dim space)


def two_mode(op1, op2, n=N2):
    return qt.tensor(op1, op2)


def q_(k, n=N2):
    return two_mode(quad_q(n), qt.qeye(n)) if k == 1 else two_mode(qt.qeye(n), quad_q(n))


def p_(k, n=N2):
    return two_mode(quad_p(n), qt.qeye(n)) if k == 1 else two_mode(qt.qeye(n), quad_p(n))


def var(op, state):
    return qt.expect(op * op, state) - qt.expect(op, state) ** 2


@pytest.mark.parametrize("U", [
    phase_shifter(0.7, N), squeezer(0.5, N), displacement(0.8 + 0.3j, N),
    beam_splitter(np.pi / 4, N2), sum_gate(1.0, N2),
])
def test_gates_are_unitary(U):
    dim = U.shape[0]
    np.testing.assert_allclose((U.dag() * U).full(), np.eye(dim), atol=1e-10)


# --- single mode ---------------------------------------------------------------------------------


def test_phase_shifter_rotates_quadratures():
    # Heisenberg: q -> cos(phi) q + sin(phi) p. So after R(pi/2), <q> equals the old <p>.
    state = qt.coherent(N, 1.0 + 0.5j)
    out = phase_shifter(np.pi / 2, N) * state
    assert qt.expect(quad_q(N), out) == pytest.approx(qt.expect(quad_p(N), state), abs=1e-8)
    # a -> e^{-i phi} a
    phi = 0.4
    out = phase_shifter(phi, N) * state
    assert qt.expect(qt.destroy(N), out) == pytest.approx(np.exp(-1j * phi) * (1.0 + 0.5j), abs=1e-8)


def test_squeezer_sign_and_variances():
    # Course p.39: V = diag(e^{-2r}, e^{2r}) / 2 for S(r)|0>; r > 0 squeezes q.
    r = 0.6
    sq = squeezer(r, N) * qt.basis(N, 0)
    assert var(quad_q(N), sq) == pytest.approx(0.5 * np.exp(-2 * r), rel=1e-6)
    assert var(quad_p(N), sq) == pytest.approx(0.5 * np.exp(2 * r), rel=1e-6)
    # Course p.31 Fock expansion: only even photon numbers.
    assert np.allclose(sq.full()[1::2], 0)


def test_quadrature_shifts():
    vac = qt.basis(N, 0)
    assert qt.expect(quad_q(N), x_shift(1.3, N) * vac) == pytest.approx(1.3, abs=1e-8)
    assert qt.expect(quad_p(N), x_shift(1.3, N) * vac) == pytest.approx(0.0, abs=1e-8)
    assert qt.expect(quad_p(N), z_shift(-0.9, N) * vac) == pytest.approx(-0.9, abs=1e-8)
    assert qt.expect(quad_q(N), z_shift(-0.9, N) * vac) == pytest.approx(0.0, abs=1e-8)


def test_displacement_commutation_phase():
    # Course p.29: D(b) D(a) = e^{b a* - b* a} D(a) D(b). This identity is the basis of GKP (p.78).
    a, b = 0.6 + 0.2j, -0.3 + 0.5j
    lhs = displacement(b, N) * displacement(a, N) * qt.basis(N, 0)
    rhs = np.exp(b * np.conj(a) - np.conj(b) * a) * displacement(a, N) * displacement(b, N) * qt.basis(N, 0)
    assert (lhs - rhs).norm() < 1e-8


# --- two modes -----------------------------------------------------------------------------------


def test_beam_splitter_mixes_coherent_amplitudes():
    # Heisenberg (p.52): a1 -> cos a1 + sin a2, a2 -> -sin a1 + cos a2.
    theta, a1, a2 = 0.3, 1.0, 0.5j
    state = beam_splitter(theta, N2) * qt.tensor(qt.coherent(N2, a1), qt.coherent(N2, a2))
    A1 = qt.tensor(qt.destroy(N2), qt.qeye(N2))
    A2 = qt.tensor(qt.qeye(N2), qt.destroy(N2))
    assert qt.expect(A1, state) == pytest.approx(np.cos(theta) * a1 + np.sin(theta) * a2, abs=1e-8)
    assert qt.expect(A2, state) == pytest.approx(-np.sin(theta) * a1 + np.cos(theta) * a2, abs=1e-8)


def test_hong_ou_mandel_dip():
    # Two single photons on a 50:50 splitter always leave together: P(|1,1>) = 0.
    out = beam_splitter(np.pi / 4, N2) * qt.tensor(qt.basis(N2, 1), qt.basis(N2, 1))
    assert abs(qt.tensor(qt.basis(N2, 1), qt.basis(N2, 1)).overlap(out)) ** 2 == pytest.approx(0, abs=1e-12)
    # ... and end up in (|2,0> - |0,2>)/sqrt2 for this sign convention.
    noon = (qt.tensor(qt.basis(N2, 2), qt.basis(N2, 0)) - qt.tensor(qt.basis(N2, 0), qt.basis(N2, 2))).unit()
    assert abs(noon.overlap(out)) == pytest.approx(1.0, abs=1e-10)


def test_sum_gate_heisenberg_action():
    # q2 -> q2 + g q1, p1 -> p1 - g p2, on a low-energy displaced input.
    g = 1.0
    state = qt.tensor(qt.coherent(N2, 0.4 + 0.2j), qt.coherent(N2, -0.3 + 0.35j))
    out = sum_gate(g, N2) * state
    assert qt.expect(q_(2), out) == pytest.approx(qt.expect(q_(2), state) + g * qt.expect(q_(1), state), abs=1e-6)
    assert qt.expect(p_(1), out) == pytest.approx(qt.expect(p_(1), state) - g * qt.expect(p_(2), state), abs=1e-6)
    assert qt.expect(q_(1), out) == pytest.approx(qt.expect(q_(1), state), abs=1e-6)


def epr_state(r1, r2, n):
    vac = qt.tensor(qt.basis(n, 0), qt.basis(n, 0))
    return beam_splitter(np.pi / 4, n) * qt.tensor(squeezer(r1, n), squeezer(-r2, n)) * vac


@pytest.mark.parametrize("r1, r2", [(0.5, 0.5), (0.7, 0.3)])
def test_epr_pair_correlations(r1, r2):
    # Components 2-4 of the circuit: S(r1) x S(-r2), then B(pi/4).
    # Expected (primer section E): Var(q1 - q2) = e^{-2 r1}, Var(p1 + p2) = e^{-2 r2}.
    n = 40
    state = epr_state(r1, r2, n)
    assert check_truncation(state, tol=1e-6) < 1e-6
    assert var(q_(1, n) - q_(2, n), state) == pytest.approx(np.exp(-2 * r1), rel=1e-6)
    assert var(p_(1, n) + p_(2, n), state) == pytest.approx(np.exp(-2 * r2), rel=1e-6)
    # Each half alone looks noisier than vacuum: entanglement hides the information in correlations.
    assert var(q_(1, n), state) > 0.5


def test_squeezing_error_converges_with_truncation():
    # A tiny tail in the *output* is not enough: exponentiating a truncated generator makes the
    # gate itself slightly wrong. At N=30 the tail is ~1e-12 yet Var(q1-q2) is off by ~1e-5.
    # The real check is convergence: the error must shrink rapidly as N grows.
    errors = []
    for n in (20, 30, 40):
        s = epr_state(0.7, 0.3, n)
        errors.append(abs(var(q_(1, n) - q_(2, n), s) / np.exp(-1.4) - 1))
    assert errors[0] > 10 * errors[1] > 100 * errors[2]
    assert errors[2] < 1e-6


def test_on_modes_embeds_into_register():
    dims = [4, 5, 6]
    U = on_modes(qt.tensor(phase_shifter(0.3, 5), phase_shifter(0.1, 6)), [1, 2], dims)
    assert U.dims == [dims, dims]
    expected = qt.tensor(qt.qeye(4), phase_shifter(0.3, 5), phase_shifter(0.1, 6))
    np.testing.assert_allclose(U.full(), expected.full(), atol=1e-12)


def test_truncation_guard_flags_overflowing_state():
    big = qt.coherent(10, 3.0)  # n-bar = 9 in a 10-level space: clearly truncated
    assert tail_weight(big) > 1e-2
    with pytest.warns(TruncationWarning):
        check_truncation(big)


# --- report claims about the EPR pair and the Bell measurement (Section 2, added 2026-11-14) -----------


def loss_on_second_mode(rho, eta, n):
    from cvqec.noise import PhotonLoss
    return sum((qt.tensor(qt.qeye(n), K) * rho * qt.tensor(qt.qeye(n), K).dag() for K in PhotonLoss(eta).kraus(n)),
               0 * rho)


def log_negativity(rho):
    pt = qt.partial_transpose(rho, [0, 1])
    return float(np.log2(np.abs(np.linalg.eigvalsh(pt.full())).sum()))


@pytest.mark.parametrize("eta", [1.0, 0.6])
def test_epr_pair_is_entangled_even_after_loss(eta):
    # Report: for any separable state Var(q_a - q_b) + Var(p_a + p_b) >= 2 (vacuum variance 1/2);
    # the EPR pair goes below 2, so it is entangled. With loss the state is mixed, and the bound still holds
    # (variance of a mixture >= mean of the variances). Cross-check with the log-negativity (> 0 iff entangled
    # for two-mode Gaussian states).
    n, r = 25, 0.5
    rho = loss_on_second_mode(qt.ket2dm(epr_state(r, r, n)), eta, n)
    duan = var(q_(1, n) - q_(2, n), rho) + var(p_(1, n) + p_(2, n), rho)
    # Heisenberg picture with loss on b: q_b -> sqrt(eta) q_b + sqrt(1-eta) q_v.
    c2, s2 = np.cosh(2 * r), np.sinh(2 * r)
    # Each sum then has variance (1+eta) cosh(2r)/2 - sqrt(eta) sinh(2r) + (1-eta)/2 (= e^{-2r} at eta = 1).
    expected = (1 + eta) * c2 - 2 * np.sqrt(eta) * s2 + (1 - eta)
    assert duan == pytest.approx(expected, rel=1e-5)
    assert duan < 2
    assert log_negativity(rho) > 0.1


def test_product_of_squeezed_states_respects_separable_bound():
    # Without the beam splitter the two squeezed modes are a product state: the sum is >= 2, negativity 0.
    n = 25
    vac = qt.tensor(qt.basis(n, 0), qt.basis(n, 0))
    prod = qt.ket2dm(qt.tensor(squeezer(0.5, n), squeezer(-0.5, n)) * vac)
    assert var(q_(1, n) - q_(2, n), prod) + var(p_(1, n) + p_(2, n), prod) >= 2 - 1e-9
    assert log_negativity(prod) == pytest.approx(0, abs=1e-9)


@pytest.mark.parametrize("r1, r2", [(0.4, 0.4), (0.5, 0.25)])
def test_bell_measurement_outcome_variances(r1, r2):
    # Components 6-9. m_q = (q_a - q_in)/sqrt2 is q of mode a after B(pi/4) on (in, a); m_p = (p_in + p_a)/sqrt2
    # is q of mode in after the extra R(pi/2). Report: Var m_q = (Var q_in + (e^{-2r1} + e^{2r2})/4) / 2,
    # Var m_p = (Var p_in + (e^{2r1} + e^{-2r2})/4) / 2: dominated by the anti-squeezed noise.
    n = 25
    rho_a = qt.ket2dm(epr_state(r1, r2, n)).ptrace(0)          # Alice's half of the pair
    rho_in = qt.ket2dm(qt.coherent(n, 0.3 + 0.2j, method="analytic"))
    out = beam_splitter(np.pi / 4, n) * qt.tensor(rho_in, rho_a) * beam_splitter(np.pi / 4, n).dag()
    var_mq = var(q_(2, n), out)
    R = qt.tensor(phase_shifter(np.pi / 2, n), qt.qeye(n))
    var_mp = var(q_(1, n), R * out * R.dag())
    assert var_mq == pytest.approx(0.5 * (0.5 + (np.exp(-2 * r1) + np.exp(2 * r2)) / 4), rel=1e-4)
    assert var_mp == pytest.approx(0.5 * (0.5 + (np.exp(2 * r1) + np.exp(-2 * r2)) / 4), rel=1e-4)
