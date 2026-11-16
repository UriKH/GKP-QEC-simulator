"""Wigner functions in the course's (q, p) coordinates (course p.38: vacuum variance 1/2)."""

import numpy as np
import qutip as qt

from cvqec.conventions import WIGNER_G


def wigner(state: qt.Qobj, extent: float = 6.0, points: int = 301) -> tuple[np.ndarray, np.ndarray]:
    """Return (xs, W) with W[i, j] = W(q = xs[j], p = xs[i])."""
    xs = np.linspace(-extent, extent, points)
    return xs, qt.wigner(state, xs, xs, g=WIGNER_G)


def negativity(W: np.ndarray, xs: np.ndarray) -> float:
    """Wigner negativity  delta = integral |W| - 1  (zero for every Gaussian state)."""
    d = xs[1] - xs[0]
    return float(np.sum(np.abs(W)) * d * d - 1)


def plot_wigner(ax, state: qt.Qobj, extent: float = 6.0, points: int = 301, title: str | None = None):
    """Draw W on a symmetric red/blue scale (blue = negative), so interference fringes stand out."""
    xs, W = wigner(state, extent, points)
    vmax = np.abs(W).max()
    im = ax.imshow(W, origin="lower", extent=[-extent, extent, -extent, extent],
                   cmap="RdBu_r", vmin=-vmax, vmax=vmax, interpolation="bilinear")
    ax.set_xlabel("q")
    ax.set_ylabel("p")
    ax.set_aspect("equal")
    if title:
        ax.set_title(title)
    return im


def plot_code(code, N: int, extent: float = 6.0, points: int = 201):
    """Wigner functions of |0_L>, |1_L>, |+_L> for a code, side by side."""
    import matplotlib.pyplot as plt

    k0, k1 = code.logical_basis(N)
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.4))
    for ax, (label, st) in zip(axes, [("|0_L>", k0), ("|1_L>", k1), ("|+_L>", (k0 + k1).unit())]):
        plot_wigner(ax, st, extent, points, title=f"{code.name}  {label}")
    fig.tight_layout()
    return fig
