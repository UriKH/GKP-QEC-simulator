"""Homodyne detection in the Fock basis (course p.53-55).

A homodyne detector at angle theta measures q_theta = cos(theta) q + sin(theta) p. Its eigenstates are
rotated position eigenstates |x>_theta = R(theta)^dag |x>, so in the Fock basis
    <x_theta | n> = e^{-i theta n} psi_n(x),
where psi_n are the Hermite functions (course p.35). Everything here is built on that one formula.
"""

import numpy as np
import qutip as qt


def fock_wavefunctions(n_max: int, x: np.ndarray) -> np.ndarray:
    """psi_n(x) = <x|n> for n < n_max, shape (n_max, len(x)), via the stable three-term recurrence."""
    x = np.asarray(x, dtype=float)
    psi = np.zeros((n_max, x.size))
    psi[0] = np.pi ** -0.25 * np.exp(-x ** 2 / 2)
    if n_max > 1:
        psi[1] = np.sqrt(2) * x * psi[0]
    for n in range(2, n_max):
        psi[n] = np.sqrt(2 / n) * x * psi[n - 1] - np.sqrt((n - 1) / n) * psi[n - 2]
    return psi


def quadrature_amplitudes(N: int, x: np.ndarray, theta: float = 0.0) -> np.ndarray:
    """Matrix A[i, n] = <x_i, theta | n>. A @ ket gives the quadrature wavefunction on the grid x."""
    return (fock_wavefunctions(N, x) * np.exp(-1j * theta * np.arange(N))[:, None]).T


def homodyne_distribution(state: qt.Qobj, x: np.ndarray, theta: float = 0.0) -> np.ndarray:
    """Outcome density P(x) of an ideal homodyne at angle theta."""
    rho = qt.ket2dm(state) if state.isket else state
    A = quadrature_amplitudes(rho.shape[0], x, theta)
    return np.real(np.einsum("in,nm,im->i", A, rho.full(), A.conj()))


def position_matrix(rho: np.ndarray, x: np.ndarray) -> np.ndarray:
    """<x_i| rho |x_j> on a grid, for a Fock-basis matrix rho (any operator, not only states)."""
    A = quadrature_amplitudes(rho.shape[0], x)
    return A @ rho @ A.conj().T
