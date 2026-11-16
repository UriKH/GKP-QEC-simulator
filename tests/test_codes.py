"""Physics checks of the bosonic codes against the course book (p.76-82)."""

import numpy as np
import pytest
import qutip as qt

from cvqec.codes import CODES, FockQubit, FourLegCat, GKPCode, Lattice, make_code
from cvqec.truncation import tail_weight

N = 140
DELTA = 0.3
ALL_CODES = [GKPCode("square", DELTA), GKPCode("rectangular", DELTA, lam=1.5), GKPCode("hexagonal", DELTA),
             FourLegCat(2.0), FockQubit()]


def course_square_gkp(mu, delta, n):
    """The course's q-space form (p.82): peaks of width Delta at (2s + mu) sqrt(pi), envelope kappa = Delta."""
    q = np.linspace(-14, 14, 6001)
    centers = (2 * np.arange(-12, 13) + mu) * np.sqrt(np.pi)
    psi = sum(np.exp(-0.5 * delta ** 2 * c ** 2) * np.exp(-(q - c) ** 2 / (2 * delta ** 2)) for c in centers)
    # Project onto Fock wavefunctions psi_n(q) (p.35) with the stable three-term recurrence.
    basis = np.zeros((n, q.size))
    basis[0] = np.pi ** -0.25 * np.exp(-q ** 2 / 2)
    basis[1] = np.sqrt(2) * q * basis[0]
    for k in range(2, n):
        basis[k] = np.sqrt(2 / k) * q * basis[k - 1] - np.sqrt((k - 1) / k) * basis[k - 2]
    return qt.Qobj(basis @ psi * (q[1] - q[0])).unit()


# --- generic code properties ---------------------------------------------------------------------


@pytest.mark.parametrize("code", ALL_CODES, ids=lambda c: c.name)
def test_encoder_is_isometry_and_paulis_anticommute(code):
    V = code.encoder(N)
    np.testing.assert_allclose((V.dag() * V).full(), np.eye(2), atol=1e-10)
    P = code.logical_paulis(N)
    assert (P["X"] * P["Z"] + P["Z"] * P["X"]).norm() < 1e-10
    assert (P["X"] * P["X"] - code.projector(N)).norm() < 1e-10


@pytest.mark.parametrize("code", ALL_CODES, ids=lambda c: c.name)
def test_codewords_fit_in_truncation(code):
    for k in code.codewords(N):
        assert tail_weight(k, n_top=20) < 1e-8


def test_registry():
    assert {"gkp", "cat4", "fock"} <= set(CODES)
    assert make_code("gkp", lattice="hexagonal", delta=0.25).name == "gkp_hexagonal"
    with pytest.raises(KeyError):
        make_code("toric")


def test_mean_photon_numbers():
    assert FockQubit().mean_photon(N) == pytest.approx(0.5)
    # Delta = 0.3 (~10.5 dB) costs about 5 photons; |1_L> is slightly more energetic than |0_L>.
    assert GKPCode("square", DELTA).mean_photon(N) == pytest.approx(5.07, abs=0.01)
    # Even cats have n-bar = |a|^2 tanh|a|^2, odd-like ones coth; the 4-leg average is close to |a|^2.
    assert FourLegCat(2.0).mean_photon(N) == pytest.approx(4.0, abs=0.01)


# --- GKP -----------------------------------------------------------------------------------------


@pytest.mark.parametrize("lattice", [Lattice.square(), Lattice.rectangular(1.5), Lattice.hexagonal()],
                         ids=lambda l: l.name)
def test_lattices_satisfy_course_condition(lattice):
    # beta alpha* - beta* alpha = i pi (p.80), i.e. the unit cell has area pi in (q, p).
    assert abs(np.linalg.det(lattice.generator_matrix)) == pytest.approx(np.pi)


def test_invalid_lattice_rejected():
    with pytest.raises(ValueError):
        Lattice(1.0 + 0j, 1.0j)


def test_square_gkp_matches_course_wavefunction():
    # Our e^{-Delta^2 n} construction equals the course's peaks-and-envelope form to O(Delta^4).
    code = GKPCode("square", DELTA)
    for mu in (0, 1):
        assert abs(course_square_gkp(mu, DELTA, N).overlap(code.codeword(mu, N))) > 0.999


@pytest.mark.parametrize("lattice", ["square", "hexagonal"])
def test_stabilizer_expectation_is_gaussian_in_length(lattice):
    # Finite energy smears each peak by Delta, so <D(g)> on the code ~ exp(-|g|^2 Delta^2 / 2).
    # For the square code |2 beta|^2 = 2 pi gives the familiar e^{-pi Delta^2} (p.82).
    code = GKPCode(lattice, DELTA)
    k0 = code.codeword(0, N)
    for g, S in zip(code.lattice.stabilizers, code.stabilizer_ops(N)):
        assert qt.expect(S, k0).real == pytest.approx(np.exp(-abs(g) ** 2 * DELTA ** 2 / 2), rel=2e-3)


def test_displacement_paulis_act_logically():
    # Course p.80: X = D(alpha) maps |0_L> to |1_L>, Z = D(beta) leaves |0_L> invariant (ideally).
    # With finite Delta both hold up to the same Gaussian factor exp(-|alpha|^2 Delta^2 / 2).
    code = GKPCode("square", DELTA)
    k0, k1 = code.codewords(N)
    P = code.displacement_paulis(N)
    expected = np.exp(-np.pi * DELTA ** 2 / 4)
    assert qt.expect(P["Z"], k0).real == pytest.approx(expected, rel=2e-3)
    assert abs(k1.overlap(P["X"] * k0)) == pytest.approx(expected, rel=2e-3)


def test_gkp_approaches_ideal_as_delta_shrinks():
    errs = [1 - qt.expect(GKPCode("square", d).stabilizer_ops(N)[1], GKPCode("square", d).codeword(0, N)).real
            for d in (0.4, 0.3, 0.2)]
    assert errs[0] > errs[1] > errs[2]


def test_rectangular_trades_x_for_z_protection():
    # Larger lam: longer X shift (|alpha| = lam sqrt(pi/2)) means rarer bit flips, shorter Z shift more phase flips.
    lat = Lattice.rectangular(1.5)
    assert abs(lat.alpha) > abs(Lattice.square().alpha) > abs(lat.beta)


# --- 4-leg cat -----------------------------------------------------------------------------------


def mod4_weights(ket):
    p = np.abs(ket.full().ravel()) ** 2
    return np.array([p[r::4].sum() for r in range(4)])


def test_cat_course_basis_is_mod4_superposition():
    # Course |0_L> ~ |a> + |-a> lives on even n; |0_L> +/- |1_L> live on n = 0 and n = 2 (mod 4).
    k0, k1 = FourLegCat(2.0).logical_basis(N)
    assert mod4_weights(k0)[[1, 3]].sum() < 1e-12
    np.testing.assert_allclose(mod4_weights((k0 + k1).unit()), [1, 0, 0, 0], atol=1e-12)
    np.testing.assert_allclose(mod4_weights((k0 - k1).unit()), [0, 0, 1, 0], atol=1e-12)


def test_single_loss_is_detected_by_parity():
    # p.77: a lowers n by one, so parity flips from +1 to -1 and the damaged state leaves the code space.
    code = FourLegCat(2.0)
    Pi, parity = code.projector(N), code.parity_op(N)
    psi = code.encode((qt.basis(2, 0) + 0.6j * qt.basis(2, 1)).unit(), N)
    damaged = (qt.destroy(N) * psi).unit()
    assert qt.expect(parity, psi) == pytest.approx(1.0)
    assert qt.expect(parity, damaged) == pytest.approx(-1.0)
    assert qt.expect(Pi, damaged) == pytest.approx(0.0, abs=1e-12)
