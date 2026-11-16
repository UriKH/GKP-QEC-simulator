"""Single source of truth for units, quadratures and phase-space conventions.

Everything here follows the course book (page numbers in comments):
    hbar = omega = 1                       (p.28)
    q = (a + a^dag) / sqrt(2)              (p.28)
    p = -i (a - a^dag) / sqrt(2)           (p.28)
    [q, p] = i,  vacuum variance = 1/2     (p.28, p.31)
    D(alpha) shifts (q, p) by sqrt(2) * (Re alpha, Im alpha)   (p.37-38)

If a convention ever needs to change, it changes here and nowhere else.
"""

import numpy as np
import qutip as qt

HBAR = 1.0
VACUUM_VARIANCE = 0.5

# qutip.wigner(..., g=WIGNER_G) uses x = q and y = p in the course's convention.
WIGNER_G = np.sqrt(2)


def quad_q(N: int) -> qt.Qobj:
    """Position quadrature q = (a + a^dag)/sqrt(2) on a Fock space truncated at N levels."""
    a = qt.destroy(N)
    return (a + a.dag()) / np.sqrt(2)


def quad_p(N: int) -> qt.Qobj:
    """Momentum quadrature p = -i(a - a^dag)/sqrt(2)."""
    a = qt.destroy(N)
    return -1j * (a - a.dag()) / np.sqrt(2)


def quad_theta(N: int, theta: float) -> qt.Qobj:
    """Generalized quadrature q_theta = cos(theta) q + sin(theta) p, what a homodyne measures (p.54)."""
    return np.cos(theta) * quad_q(N) + np.sin(theta) * quad_p(N)


def alpha_to_qp(alpha: complex) -> tuple[float, float]:
    """Phase-space shift (dq, dp) produced by D(alpha)."""
    return np.sqrt(2) * alpha.real, np.sqrt(2) * alpha.imag


def qp_to_alpha(dq: float, dp: float) -> complex:
    """The alpha with D(alpha) shifting q by dq and p by dp."""
    return (dq + 1j * dp) / np.sqrt(2)


def squeezing_db(r: float) -> float:
    """Squeezing of S(r) in dB: the squeezed variance is e^{-2r} times vacuum."""
    return -10 * np.log10(np.exp(-2 * r))


def r_from_db(db: float) -> float:
    return db * np.log(10) / 20


def gkp_delta_to_db(delta: float) -> float:
    """Finite-energy GKP squeezing: each spike has variance Delta^2/2, i.e. Delta^2 times vacuum."""
    return -10 * np.log10(delta ** 2)
