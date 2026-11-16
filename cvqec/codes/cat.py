"""Cat codes (course p.76-78) and the unencoded Fock qubit used as the break-even reference."""

import numpy as np
import qutip as qt

from cvqec.codes.base import BosonicCode, register_code


def _coherent(alpha: complex, N: int) -> qt.Qobj:
    # Analytic Fock coefficients: exact for any N, unlike displacing the truncated vacuum.
    return qt.coherent(N, alpha, method="analytic")


@register_code("cat4")
class FourLegCat(BosonicCode):
    """4-leg cat in the course's basis:  |0_L> ~ |a> + |-a>,   |1_L> ~ |ia> + |-ia>  (p.77).

    In photon-number terms |0_L> = (|C0> + |C2>)/sqrt2 and |1_L> = (|C0> - |C2>)/sqrt2, where |C0>, |C2>
    are the four-coherent-state superpositions supported on n = 0 and n = 2 (mod 4). A photon loss
    maps n mod 4 from {0, 2} to {3, 1}: it changes parity without leaving a trace on the logical
    information, so measuring parity detects it and one loss is correctable (N/2 - 1 = 1 for N = 4 legs).
    """

    def __init__(self, alpha: float = 2.0):
        self.alpha = alpha
        self.name = f"cat4(alpha={alpha})"

    def codewords(self, N: int) -> tuple[qt.Qobj, qt.Qobj]:
        a = self.alpha
        k0 = (_coherent(a, N) + _coherent(-a, N)).unit()
        k1 = (_coherent(1j * a, N) + _coherent(-1j * a, N)).unit()
        return k0, k1

    def parity_op(self, N: int) -> qt.Qobj:
        """Photon-number parity (-1)^n, the error syndrome for a single loss."""
        return qt.Qobj(np.diag((-1.0) ** np.arange(N)))


@register_code("fock")
class FockQubit(BosonicCode):
    """Unencoded qubit |0>, |1> in the lowest Fock states. A code 'beats break-even' if it outlives this."""

    name = "fock"

    def codewords(self, N: int) -> tuple[qt.Qobj, qt.Qobj]:
        return qt.basis(N, 0), qt.basis(N, 1)
