"""The channel interface. Every noise source is a CPTP map on a single mode's density matrix.

A channel must implement `apply`. The other methods have sensible defaults:
  - `adjoint` (the Heisenberg-picture map E^dag, needed by the Petz decoder) defaults to Kraus form;
  - `kraus` is needed only by the full multi-mode engine for channels it cannot defer;
  - `lindblad` gives the continuous-time form (collapse operators, duration) used to cross-check
    `apply` against qutip.mesolve in the tests (primer section I).
"""

from abc import ABC, abstractmethod

import numpy as np
import qutip as qt


class Channel(ABC):
    name: str = "channel"

    @abstractmethod
    def apply(self, rho: qt.Qobj) -> qt.Qobj:
        """E(rho). Must be linear, so it also accepts non-Hermitian operators such as |0_L><1_L|."""

    def adjoint(self, X: qt.Qobj) -> qt.Qobj:
        """E^dag(X) = sum_k K_k^dag X K_k."""
        return sum((K.dag() * X * K for K in self.kraus(X.shape[0])), 0 * X)

    def kraus(self, N: int) -> list[qt.Qobj]:
        raise NotImplementedError(f"{type(self).__name__} has no Kraus form")

    def lindblad(self, N: int) -> tuple[list[qt.Qobj], float]:
        raise NotImplementedError(f"{type(self).__name__} has no Lindblad form")

    def apply_many(self, rhos: np.ndarray) -> np.ndarray:
        """E applied to a batch of matrices with shape (B, N, N). Used by the Fock engine."""
        out = np.zeros_like(rhos)
        for K in self.kraus(rhos.shape[-1]):
            K = K.full()
            out += K @ rhos @ K.conj().T  # batched matmul; much faster than a three-operand einsum
        return out

    def __repr__(self):
        return self.name


class Compose(Channel):
    """Channels applied in list order: Compose([A, B])(rho) = B(A(rho))."""

    def __init__(self, channels: list[Channel]):
        self.channels = [c for c in channels if not isinstance(c, Identity)]
        self.name = " -> ".join(c.name for c in self.channels) or "identity"

    def apply(self, rho):
        for c in self.channels:
            rho = c.apply(rho)
        return rho

    def adjoint(self, X):
        for c in reversed(self.channels):
            X = c.adjoint(X)
        return X


class Identity(Channel):
    name = "identity"

    def apply(self, rho):
        return rho

    def adjoint(self, X):
        return X

    def kraus(self, N):
        return [qt.qeye(N)]


CHANNELS: dict[str, type[Channel]] = {}


def register_channel(key: str):
    def wrap(cls):
        CHANNELS[key] = cls
        return cls
    return wrap


def make_channel(key: str, **params) -> Channel:
    try:
        return CHANNELS[key](**params)
    except KeyError:
        raise KeyError(f"unknown channel {key!r}; registered: {sorted(CHANNELS)}") from None
