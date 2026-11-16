"""The logical channel: what the whole encode -> noise -> decode chain does to one qubit.

It is stored as the four outputs Lambda(|i><j|), i, j in {0, 1}. Linearity makes these
enough to know the channel on every input, and every metric below is computed from them.
"""

import numpy as np

PAULIS = {"I": np.eye(2), "X": np.array([[0, 1], [1, 0]]), "Y": np.array([[0, -1j], [1j, 0]]),
          "Z": np.diag([1.0, -1.0])}


class LogicalChannel:
    def __init__(self, outputs: np.ndarray):
        self.outputs = np.asarray(outputs, dtype=complex)  # shape (2, 2, 2, 2): outputs[i, j] = Lambda(|i><j|)

    def apply(self, rho: np.ndarray) -> np.ndarray:
        return np.einsum("ij,ijab->ab", rho, self.outputs)

    def choi(self) -> np.ndarray:
        """J = (1/2) sum_ij |i><j| (x) Lambda(|i><j|), the state the channel makes from half a Bell pair."""
        J = np.zeros((4, 4), complex)
        for i in range(2):
            for j in range(2):
                J[2 * i:2 * i + 2, 2 * j:2 * j + 2] = self.outputs[i, j] / 2
        return J

    def entanglement_fidelity(self) -> float:
        """F_e = <Phi+| J |Phi+>: how well the channel preserves entanglement with a reference qubit."""
        return float(np.real(sum(self.outputs[i, j][i, j] for i in range(2) for j in range(2))) / 4)

    def average_fidelity(self) -> float:
        """Mean fidelity over all pure inputs (Horodecki/Nielsen): F_avg = (2 F_e + 1) / 3 for a qubit."""
        return (2 * self.entanglement_fidelity() + 1) / 3

    def ptm(self) -> np.ndarray:
        """Pauli transfer matrix R_ab = Tr[P_a Lambda(P_b)] / 2 (rows/cols ordered I, X, Y, Z). Identity = no error."""
        P = list(PAULIS.values())
        return np.real(np.array([[np.trace(Pa @ self.apply(Pb)) / 2 for Pb in P] for Pa in P]))

    def pauli_error_rates(self) -> dict[str, float]:
        """Error probabilities of the Pauli-twirled channel, from the PTM diagonal."""
        _, rx, ry, rz = np.diag(self.ptm())
        return {"p_x": (1 + rx - ry - rz) / 4, "p_y": (1 - rx + ry - rz) / 4, "p_z": (1 - rx - ry + rz) / 4}

    @classmethod
    def identity(cls) -> "LogicalChannel":
        out = np.zeros((2, 2, 2, 2), complex)
        for i in range(2):
            for j in range(2):
                out[i, j][i, j] = 1
        return cls(out)
