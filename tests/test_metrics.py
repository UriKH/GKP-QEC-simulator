"""Logical-channel metrics against hand-computable channels."""

import numpy as np
import pytest

from cvqec.metrics import PAULIS, LogicalChannel


def pauli_channel(px, py, pz):
    out = np.zeros((2, 2, 2, 2), complex)
    for i in range(2):
        for j in range(2):
            E = np.outer(np.eye(2)[i], np.eye(2)[j])
            out[i, j] = sum(p * P @ E @ P.conj().T for p, P in
                            ((1 - px - py - pz, PAULIS["I"]), (px, PAULIS["X"]), (py, PAULIS["Y"]), (pz, PAULIS["Z"])))
    return LogicalChannel(out)


def test_identity_channel():
    ch = LogicalChannel.identity()
    assert ch.entanglement_fidelity() == pytest.approx(1)
    np.testing.assert_allclose(ch.ptm(), np.eye(4), atol=1e-12)


def test_pauli_channel_metrics():
    ch = pauli_channel(0.02, 0.01, 0.05)
    assert ch.entanglement_fidelity() == pytest.approx(0.92)          # F_e = probability of no error
    rates = ch.pauli_error_rates()
    assert (rates["p_x"], rates["p_y"], rates["p_z"]) == pytest.approx((0.02, 0.01, 0.05))
    np.testing.assert_allclose(np.diag(ch.ptm()), [1, 1 - 2 * 0.06, 1 - 2 * 0.07, 1 - 2 * 0.03], atol=1e-12)


def test_average_fidelity_formula_matches_2_design_average():
    # The six Pauli eigenstates form a 2-design, so their mean fidelity is the exact Haar average.
    ch = pauli_channel(0.03, 0.02, 0.07)
    states = [np.array(v, complex) / np.linalg.norm(v) for v in
              ([1, 0], [0, 1], [1, 1], [1, -1], [1, 1j], [1, -1j])]
    mean = np.mean([np.real(s.conj() @ ch.apply(np.outer(s, s.conj())) @ s) for s in states])
    assert ch.average_fidelity() == pytest.approx(mean, abs=1e-12)
