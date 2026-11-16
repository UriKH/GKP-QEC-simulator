"""Gaussian gates in the course book's conventions, as QuTiP unitaries on a truncated Fock space.

Each gate is written as U = exp(generator) with the course's sign convention and documented by
its Heisenberg action U^dag A U (the form the course uses on p.50-53).

Truncation note: exponentiating a *truncated* anti-Hermitian generator gives an exactly unitary
matrix, but it agrees with the true infinite-dimensional gate only on states whose photon-number
support stays well below N. Number-conserving gates (phase shifter, beam splitter) are exact on
every state with total photon number < N; the others (squeezer, displacement, SUM) need headroom.
Use `cvqec.truncation` to check the headroom.
"""

import numpy as np
import qutip as qt
from scipy.special import gammaln

from cvqec.conventions import quad_p, quad_q

# --- single-mode gates -------------------------------------------------------------------------


def phase_shifter(phi: float, N: int) -> qt.Qobj:
    """R(phi) = exp(-i phi a^dag a)  (p.52).

    Heisenberg: a -> e^{-i phi} a,  q -> cos(phi) q + sin(phi) p.
    So R(pi/2) followed by a q-homodyne measures p.
    """
    return (-1j * phi * qt.num(N)).expm()


def squeezer(r: float, N: int) -> qt.Qobj:
    """S(r) = exp(r/2 (a^2 - a^dag^2))  (p.31, p.53).

    Heisenberg: a -> cosh(r) a - sinh(r) a^dag,  q -> e^{-r} q,  p -> e^{+r} p.
    r > 0 squeezes q. QuTiP's squeeze(N, z) = exp((z* a^2 - z a^dag^2)/2) is the same operator.
    """
    return qt.squeeze(N, r)


def displacement(alpha: complex, N: int) -> qt.Qobj:
    """D(alpha) = exp(alpha a^dag - alpha* a)  (p.29).

    Heisenberg: a -> a + alpha,  (q, p) -> (q + sqrt2 Re alpha, p + sqrt2 Im alpha).
    """
    return qt.displace(N, alpha)


def displacement_exact(alpha: complex, N: int) -> np.ndarray:
    """Exact matrix elements <m|D(alpha)|n>, m, n < N, of the *infinite-dimensional* displacement.

    `displacement` exponentiates a truncated generator: exactly unitary, but wrong near the cutoff and
    useless for |alpha| comparable to sqrt(N). This block is the opposite: every entry is exact, but the
    block is not unitary because amplitude pushed above level N is (correctly) missing. Used where large
    displacements matter (GKP decoding, feed-forward). Formula (Cahill-Glauber), for m = n + k >= n:
        <m|D|n> = sqrt(n!/m!) alpha^k e^{-|alpha|^2/2} L_n^{(k)}(|alpha|^2),   <n|D|m> = (-alpha*)^k / alpha^k * <m|D|n>,
    with the Laguerre polynomials from their forward recurrence in n (stable, vectorized over k).
    """
    x = abs(alpha) ** 2
    if x > 1300:
        raise ValueError("|alpha|^2 > 1300 would overflow the Laguerre recurrence")
    k = np.arange(N)
    L = np.zeros((N, N))  # L[n, k] = L_n^{(k)}(x)
    L[0] = 1.0
    if N > 1:
        L[1] = 1.0 + k - x
    for n in range(1, N - 1):
        L[n + 1] = ((2 * n + 1 + k - x) * L[n] - (n + k) * L[n - 1]) / (n + 1)
    n_idx, k_idx = np.meshgrid(np.arange(N), k, indexing="ij")
    valid = n_idx + k_idx < N
    lg = gammaln(np.arange(2 * N) + 1)
    log_pref = -x / 2 + 0.5 * (lg[n_idx] - lg[n_idx + k_idx])
    if alpha != 0:
        log_pref = log_pref + k_idx * np.log(abs(alpha))
    else:
        log_pref[k_idx > 0] = -np.inf  # D(0) = identity: only the diagonal survives
    lower = np.where(valid, np.exp(log_pref) * L, 0.0) * np.exp(1j * k_idx * np.angle(alpha))
    D = np.zeros((N, N), dtype=complex)
    rows, cols = n_idx[valid] + k_idx[valid], n_idx[valid]
    D[rows, cols] = lower[valid]                                        # m >= n
    D[cols, rows] = lower[valid] * (-1.0) ** k_idx[valid] * np.exp(-2j * k_idx[valid] * np.angle(alpha))
    return D


def gaussian_unitary(M: np.ndarray, N: int) -> qt.Qobj:
    """Single-mode Gaussian unitary U whose Heisenberg action is U^dag (q, p) U = M (q, p).

    Euler decomposition M = R(phi1) diag(e^{-r}, e^{r}) R(phi2) from the SVD, so U = R(phi1) S(r) R(phi2)
    (Heisenberg matrices compose in the same order as the operators). The squeezer is the only
    truncation-sensitive factor: give it headroom by choosing N well above the state's support.
    """
    if not np.isclose(np.linalg.det(M), 1.0):
        raise ValueError("M must be symplectic (det = 1)")
    U, sv, Vt = np.linalg.svd(M)
    if np.linalg.det(U) < 0:
        F = np.diag([-1.0, 1.0])
        U, Vt = U @ F, F @ Vt
    phi1, phi2 = np.arctan2(U[0, 1], U[0, 0]), np.arctan2(Vt[0, 1], Vt[0, 0])
    return phase_shifter(phi1, N) * squeezer(-np.log(sv[0]), N) * phase_shifter(phi2, N)


def x_shift(s: float, N: int) -> qt.Qobj:
    """X(s) = exp(-i s p): shifts q by s  (p.40). Equal to D(s/sqrt2) (p.81)."""
    return displacement(s / np.sqrt(2), N)


def z_shift(s: float, N: int) -> qt.Qobj:
    """Z(s) = exp(i s q): shifts p by s  (p.40). Equal to D(i s/sqrt2)."""
    return displacement(1j * s / np.sqrt(2), N)


# --- two-mode gates (mode order: (first, second)) ----------------------------------------------


def beam_splitter(theta: float, N: int) -> qt.Qobj:
    """B(theta) = exp(-theta (a1 a2^dag - a1^dag a2))  (p.52).

    Heisenberg: a1 -> cos(theta) a1 + sin(theta) a2,
                a2 -> -sin(theta) a1 + cos(theta) a2.
    theta = pi/4 is the balanced (50:50) splitter used for the EPR pair and the Bell measurement.
    """
    a = qt.destroy(N)
    a1 = qt.tensor(a, qt.qeye(N))
    a2 = qt.tensor(qt.qeye(N), a)
    return (-theta * (a1 * a2.dag() - a1.dag() * a2)).expm()


def sum_gate(g: float, N: int) -> qt.Qobj:
    """SUM(g) = exp(-i g q1 p2), the CV analogue of CNOT (bridge material, not in the course).

    Heisenberg: q2 -> q2 + g q1,  p1 -> p1 - g p2  (q1, p2 unchanged).
    Used in Steane-type GKP error correction to copy a data-mode shift onto an ancilla.
    """
    return (-1j * g * qt.tensor(quad_q(N), quad_p(N))).expm()


# --- placing gates inside a multi-mode register ------------------------------------------------


def on_modes(gate: qt.Qobj, targets: int | list[int], dims: list[int]) -> qt.Qobj:
    """Embed a 1- or 2-mode gate into a register with Fock dimensions `dims`, acting on `targets`."""
    return qt.expand_operator(gate, dims=dims, targets=targets)
