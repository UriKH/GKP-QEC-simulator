"""The interface every bosonic code implements, plus a registry so experiments can say "gkp_hex".

A code is fully described by its two logical codewords |0_L>, |1_L> in a truncated Fock space.
Everything else (encoder, projector, logical Paulis, n-bar) is derived from them here, so a new
code only has to implement `codewords`.
"""

from abc import ABC, abstractmethod

import numpy as np
import qutip as qt


class BosonicCode(ABC):
    """A qubit encoded in one bosonic mode."""

    name: str = "code"

    @abstractmethod
    def codewords(self, N: int) -> tuple[qt.Qobj, qt.Qobj]:
        """The code's |0_L>, |1_L> as normalized kets on N Fock levels (possibly not orthogonal)."""

    def logical_basis(self, N: int) -> tuple[qt.Qobj, qt.Qobj]:
        """Orthonormal logical basis: the codewords after symmetric (Loewdin) orthonormalization.

        Finite-energy codewords overlap slightly (e.g. ~3e-4 for square GKP at Delta = 0.3).
        Symmetric orthonormalization removes the overlap while moving each codeword as little
        as possible, so the encoder below is an exact isometry.
        """
        k0, k1 = self.codewords(N)
        A = np.column_stack([k0.full().ravel(), k1.full().ravel()])
        eigval, eigvec = np.linalg.eigh(A.conj().T @ A)
        A = A @ eigvec @ np.diag(eigval ** -0.5) @ eigvec.conj().T
        return qt.Qobj(A[:, 0]), qt.Qobj(A[:, 1])

    def encoder(self, N: int) -> qt.Qobj:
        """Isometry V = |0_L><0| + |1_L><1| from a qubit into the mode (dims [[N], [2]])."""
        k0, k1 = self.logical_basis(N)
        return qt.Qobj(np.column_stack([k0.full().ravel(), k1.full().ravel()]), dims=[[N], [2]])

    def encode(self, qubit_state: qt.Qobj, N: int) -> qt.Qobj:
        """Map a qubit ket or density matrix into the code space."""
        V = self.encoder(N)
        return V * qubit_state if qubit_state.isket else V * qubit_state * V.dag()

    def projector(self, N: int) -> qt.Qobj:
        """Projector onto the code space."""
        V = self.encoder(N)
        return V * V.dag()

    def logical_paulis(self, N: int) -> dict[str, qt.Qobj]:
        """Logical X, Y, Z as operators on the mode, acting inside the code space (zero outside)."""
        V = self.encoder(N)
        return {k: V * P * V.dag() for k, P in (("X", qt.sigmax()), ("Y", qt.sigmay()), ("Z", qt.sigmaz()))}

    def mean_photon(self, N: int) -> float:
        """n-bar of the maximally mixed code state, the fair 'energy cost' for comparing codes."""
        return float(qt.expect(qt.num(N), self.projector(N) / 2))


CODES: dict[str, type[BosonicCode]] = {}


def register_code(key: str):
    """Class decorator: `@register_code("gkp")` makes the code constructible by name."""
    def wrap(cls):
        CODES[key] = cls
        return cls
    return wrap


def make_code(key: str, **params) -> BosonicCode:
    """Build a registered code, e.g. make_code("gkp", lattice="hexagonal", delta=0.3)."""
    try:
        return CODES[key](**params)
    except KeyError:
        raise KeyError(f"unknown code {key!r}; registered: {sorted(CODES)}") from None
