"""Decoders: independent cross-checks and the physics each one is supposed to capture (primer section J)."""

import numpy as np
import pytest
import qutip as qt

from cvqec.codes import FockQubit, FourLegCat, GKPCode
from cvqec.decoders import (GKPIdealDecoder, GKPSteaneDecoder, SquareGrid, cell_polygon, make_decoder,
                            polygon_fourier)
from cvqec.engines.effective import run_logical
from cvqec.noise import GaussianDisplacement, Identity, PhotonLoss
from cvqec.ops.measurement import homodyne_distribution

N = 140


def fidelity(psi2, rho2):
    v = psi2.full().ravel()
    return float(np.real(v.conj() @ rho2 @ v))


@pytest.fixture(scope="module")
def noisy_square():
    code = GKPCode("square", 0.3)
    k0, _ = code.logical_basis(N)
    return code, GaussianDisplacement(0.05).apply(qt.ket2dm(k0))


# --- geometry --------------------------------------------------------------------------------------


def test_cell_fourier_transform_of_square_is_sinc():
    verts = cell_polygon(GKPCode("square").lattice)
    s = np.sqrt(np.pi)
    for k in ([0.7, 0.0], [1.3, -2.1]):
        expected = s * s * np.sinc(k[0] * s / (2 * np.pi)) * np.sinc(k[1] * s / (2 * np.pi))
        assert polygon_fourier(verts, np.array(k)).real == pytest.approx(expected, abs=1e-12)


def test_hexagonal_voronoi_cell():
    verts = cell_polygon(GKPCode("hexagonal").lattice)
    assert len(verts) == 6
    assert polygon_fourier(verts, np.zeros(2)).real == pytest.approx(np.pi)   # cell area = pi for every lattice


# --- ideal decoder ---------------------------------------------------------------------------------


def test_ideal_decoder_equals_direct_q_binning(noisy_square):
    # For the square code <Z_L> = <sgn cos(sqrt(pi) q)>: integrate the homodyne distribution directly.
    code, rho = noisy_square
    x = np.linspace(-25, 25, 200001)
    P = homodyne_distribution(rho, x)
    even = np.floor(x / np.sqrt(np.pi) + 0.5) % 2 == 0
    direct = np.sum(P * even) * (x[1] - x[0])
    assert GKPIdealDecoder(code).decode(rho)[0, 0].real == pytest.approx(direct, abs=1e-5)


def test_grid_decoder_agrees_with_characteristic_function_decoder(noisy_square):
    # Two independent derivations of the same readout; the grid converges as O(h^2) (2e-4 at M = 16).
    code, rho = noisy_square
    char = GKPIdealDecoder(code).decode(rho)
    g = SquareGrid(32)
    np.testing.assert_allclose(g.subsystem_trace(g.from_fock(rho.full())), char, atol=1e-4)


def test_hexagonal_frame_change_matches_parallelogram_cells():
    # Gaussian unitary to the square frame + square cells == parallelogram cells of the hexagonal lattice.
    code = GKPCode("hexagonal", 0.3)
    k0, k1 = code.logical_basis(N)
    rho = GaussianDisplacement(0.05).apply(qt.ket2dm((k0 + 1j * k1).unit()))
    par = GKPIdealDecoder(code, cell="parallelogram").decode(rho)
    np.testing.assert_allclose(GKPSteaneDecoder(code).decode_without_round(rho), par, atol=1e-3)


@pytest.mark.parametrize("lattice", ["square", "rectangular", "hexagonal"])
def test_ideal_decoder_recovers_codewords(lattice):
    code = GKPCode(lattice, 0.3)
    dec = GKPIdealDecoder(code)
    for amps in ([1, 0], [0, 1], [1, 1], [1, 1j]):
        psi = qt.Qobj(np.array(amps, complex)).unit()
        assert fidelity(psi, dec.decode(qt.ket2dm(code.encode(psi, N)))) > 0.99


def test_decoded_state_is_physical():
    code = GKPCode("hexagonal", 0.3)
    rng = np.random.default_rng(3)
    v = rng.normal(size=60) + 1j * rng.normal(size=60)
    rho = qt.ket2dm(qt.Qobj(np.concatenate([v, np.zeros(N - 60)])).unit())
    out = GKPIdealDecoder(code).decode(rho)
    assert np.trace(out).real == pytest.approx(1.0, abs=1e-8)
    assert np.linalg.eigvalsh((out + out.conj().T) / 2).min() > -1e-8


def test_voronoi_cells_make_the_hexagonal_code_symmetric_and_better():
    # Minimum-distance (Voronoi) decoding is what gives the hexagonal lattice its edge (primer section H):
    # better than parallelogram cells, and equally good against X, Y and Z errors.
    code = GKPCode("hexagonal", 0.3)
    ch = GaussianDisplacement(0.1)
    vor = run_logical(code, ch, GKPIdealDecoder(code), N).logical
    par = run_logical(code, ch, GKPIdealDecoder(code, cell="parallelogram"), N).logical
    assert vor.average_fidelity() > par.average_fidelity() + 0.005
    rates = vor.pauli_error_rates()
    assert max(rates.values()) - min(rates.values()) < 1e-3
    sq = GKPCode("square", 0.3)
    a = run_logical(sq, ch, GKPIdealDecoder(sq), N).logical.average_fidelity()
    b = run_logical(sq, ch, GKPIdealDecoder(sq, cell="parallelogram"), N).logical.average_fidelity()
    assert a == pytest.approx(b, abs=1e-12)    # for the square lattice the two cells coincide


# --- Steane ----------------------------------------------------------------------------------------


def test_steane_approaches_ideal_as_ancillas_improve(noisy_square):
    code, _ = noisy_square
    k0, k1 = code.logical_basis(N)
    rho = GaussianDisplacement(0.05).apply(qt.ket2dm((k0 + k1).unit()))
    ideal = GKPSteaneDecoder(code).decode_without_round(rho)[0, 1].real
    coherences = [GKPSteaneDecoder(code, delta_ancilla=d).decode(rho)[0, 1].real for d in (0.15, 0.25, 0.4)]
    assert coherences[0] > coherences[1] > coherences[2]
    assert ideal - coherences[0] < 0.01


# --- cat, projection, Petz -------------------------------------------------------------------------


def test_cat_parity_decoder_undoes_a_single_loss():
    cat = FourLegCat(2.25)
    psi = qt.Qobj(np.array([0.6, 0.8j]))
    damaged = (qt.destroy(N) * cat.encode(psi, N)).unit()
    assert fidelity(psi, make_decoder("cat_parity", cat).decode(qt.ket2dm(damaged))) > 0.9999


def test_petz_is_within_barnum_knill_bound():
    # F_e(Petz) >= F_e(optimal)^2 >= F_e(any decoder)^2   (Barnum & Knill 2002).
    for code, other in ((FourLegCat(2.25), "cat_parity"), (FockQubit(), "projection")):
        ch = PhotonLoss(0.9)
        petz = run_logical(code, ch, make_decoder("petz", code), 60).logical.entanglement_fidelity()
        ref = run_logical(code, ch, make_decoder(other, code), 60).logical.entanglement_fidelity()
        assert petz >= ref ** 2 - 1e-9


def test_auto_decoder_choice():
    assert make_decoder("auto", GKPCode()).name == "gkp_ideal"
    assert make_decoder("auto", FourLegCat()).name == "cat_parity"
    assert make_decoder("auto", FockQubit()).name == "projection"


# --- report claims about the Steane round and the lattices (Sections 1-2, added 2026-11-14) -----------


@pytest.mark.parametrize("lattice", ["square", "rectangular", "hexagonal"])
def test_steane_round_does_not_read_the_logical_state(lattice):
    # The |+> / |0> ancillas have peaks at every multiple of sqrt(pi) in the measured quadrature, so the
    # syndrome reveals only the shift, never the logical value: with no channel noise and good ancillas the
    # round must leave all six Pauli eigenstates (the 2-design behind average_fidelity) intact.
    code = GKPCode(lattice, 0.3, lam=1.5) if lattice == "rectangular" else GKPCode(lattice, 0.3)
    res = run_logical(code, Identity(), GKPSteaneDecoder(code, delta_ancilla=0.15), N).logical
    assert res.average_fidelity() > 0.99


def test_finite_ancillas_add_their_own_error():
    # Report: the ancillas are finite-Delta GKP states, so the round itself causes errors, even without a channel.
    code = GKPCode("square", 0.3)
    f = [run_logical(code, Identity(), GKPSteaneDecoder(code, delta_ancilla=d), N).logical.average_fidelity()
         for d in (0.15, 0.3)]
    assert f[0] > f[1] + 0.01


def test_steane_hexagon_pays_for_separate_decoding():
    # With joint=False the Steane round corrects Q and P separately (parallelogram cell), so even with good ancillas
    # it cannot beat the separate-decoding ideal readout, which in turn is clearly below the joint (hexagon) one.
    code, ch = GKPCode("hexagonal", 0.3), GaussianDisplacement(0.1)
    joint = run_logical(code, ch, GKPIdealDecoder(code), N).logical.average_fidelity()
    separate = run_logical(code, ch, GKPIdealDecoder(code, cell="parallelogram"), N).logical.average_fidelity()
    steane = run_logical(code, ch, GKPSteaneDecoder(code, delta_ancilla=0.15, joint=False), N).logical.average_fidelity()
    assert steane < separate < joint - 0.01


def test_rectangular_code_wins_only_under_matching_bias():
    # Rectangular lambda = 1.5 has a longer cell along q: it loses to the square code under symmetric shift
    # noise and wins when the q-noise is the larger one (the asymmetric-squeezing result of Section 4).
    sq, rect = GKPCode("square", 0.3), GKPCode("rectangular", 0.3, lam=1.5)

    def F(code, vq, vp):
        return run_logical(code, GaussianDisplacement(vq, vp), GKPIdealDecoder(code), N).logical.average_fidelity()

    assert F(sq, 0.1, 0.1) > F(rect, 0.1, 0.1)
    assert F(rect, 0.25, 0.05) > F(sq, 0.25, 0.05)
    assert F(rect, 0.05, 0.25) < F(sq, 0.05, 0.25)


# --- joint decoding and inefficient syndrome detectors (added 2026-11-14) ---------------------------------


@pytest.mark.parametrize("lattice", ["square", "rectangular"])
def test_joint_equals_separate_for_rectangular_cells(lattice):
    # For square and rectangular lattices the nearest-point cell IS the rectangle: the lookup table must be empty.
    code = GKPCode(lattice, 0.3, lam=1.5) if lattice == "rectangular" else GKPCode(lattice, 0.3)
    assert not GKPSteaneDecoder(code).joint_table().any()


def test_joint_table_is_inversion_symmetric():
    # The hexagonal lattice is symmetric under e -> -e, so syndrome -s needs the opposite extra shift
    # (except on the cell boundary s = -M/2, which has no partner in the grid).
    dec = GKPSteaneDecoder(GKPCode("hexagonal", 0.3))
    t, M = dec.joint_table(), dec.grid.M
    assert t.any()
    for i in range(1, M):
        for j in range(1, M):
            np.testing.assert_array_equal(t[i, j], -t[M - i, M - j])


def test_joint_steane_recovers_the_hexagonal_advantage():
    # Same two syndromes, decoded together: better than separate decoding, and with good ancillas even better than
    # separate decoding with PERFECT syndromes. So the hexagon's loss under Steane came from the decoding rule.
    code, ch = GKPCode("hexagonal", 0.3), GaussianDisplacement(0.1)
    F = lambda dec: run_logical(code, ch, dec, N).logical.average_fidelity()
    joint = F(GKPSteaneDecoder(code, delta_ancilla=0.15))
    separate = F(GKPSteaneDecoder(code, delta_ancilla=0.15, joint=False))
    ideal_separate = F(GKPIdealDecoder(code, cell="parallelogram"))
    ideal_joint = F(GKPIdealDecoder(code))
    assert separate + 0.01 < joint
    assert ideal_separate < joint < ideal_joint


def test_inefficient_syndrome_detector():
    # eta_syndrome = 1 is the old round exactly; lower efficiency only adds noise; the round stays trace preserving.
    code = GKPCode("square", 0.3)
    k0, k1 = code.logical_basis(N)
    rho = GaussianDisplacement(0.05).apply(qt.ket2dm((k0 + k1).unit()))
    perfect = GKPSteaneDecoder(code).decode(rho)
    np.testing.assert_allclose(GKPSteaneDecoder(code, eta_syndrome=1.0).decode(rho), perfect, atol=1e-12)
    coh = [GKPSteaneDecoder(code, eta_syndrome=e).decode(rho) for e in (0.98, 0.9)]
    assert perfect[0, 1].real > coh[0][0, 1].real > coh[1][0, 1].real
    for out in coh:
        assert np.trace(out).real == pytest.approx(1.0, abs=1e-6)
    with pytest.raises(ValueError):
        GKPSteaneDecoder(code, eta_syndrome=0.0)
