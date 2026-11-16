import numpy as np
import pytest
import qutip as qt

from cvqec.conventions import (
    alpha_to_qp, gkp_delta_to_db, qp_to_alpha, quad_p, quad_q, quad_theta, r_from_db, squeezing_db,
)

N = 30


def test_canonical_commutator_below_cutoff():
    # [q, p] = i holds exactly except in the last Fock level, where truncation breaks it.
    comm = (quad_q(N) * quad_p(N) - quad_p(N) * quad_q(N)).full()
    np.testing.assert_allclose(comm[:-1, :-1], 1j * np.eye(N - 1), atol=1e-12)


def test_vacuum_variance_is_half():
    vac = qt.basis(N, 0)
    assert qt.expect(quad_q(N) ** 2, vac) == pytest.approx(0.5)
    assert qt.expect(quad_p(N) ** 2, vac) == pytest.approx(0.5)


def test_coherent_state_position_in_phase_space():
    # Course p.38: W_|alpha> is the vacuum Gaussian centered at sqrt2 (Re alpha, Im alpha).
    alpha = 1.2 - 0.7j
    state = qt.coherent(N, alpha)
    q0, p0 = alpha_to_qp(alpha)
    assert qt.expect(quad_q(N), state) == pytest.approx(q0, abs=1e-8)
    assert qt.expect(quad_p(N), state) == pytest.approx(p0, abs=1e-8)
    assert qp_to_alpha(q0, p0) == pytest.approx(alpha)


def test_quad_theta_interpolates():
    np.testing.assert_allclose(quad_theta(N, 0).full(), quad_q(N).full())
    np.testing.assert_allclose(quad_theta(N, np.pi / 2).full(), quad_p(N).full(), atol=1e-15)


def test_db_conversions():
    assert squeezing_db(r_from_db(10.0)) == pytest.approx(10.0)
    assert squeezing_db(1.0) == pytest.approx(8.6859, abs=1e-4)
    assert gkp_delta_to_db(np.sqrt(0.1)) == pytest.approx(10.0)
