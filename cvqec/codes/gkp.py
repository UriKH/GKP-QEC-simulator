"""GKP codes on an arbitrary lattice (course p.78-82).

The course defines a GKP code by two complex numbers (alpha, beta):
    logical  X = D(alpha),   Z = D(beta),
    stabilizers   S_X = D(2 alpha),   S_Z = D(2 beta),
with  beta alpha* - beta* alpha = i pi  so that X and Z anticommute (via D(b)D(a) = e^{b a* - b* a} D(a)D(b)).
Square, rectangular and hexagonal codes differ only in (alpha, beta), so they share one class.

Construction. The ideal |mu_L> (mu = 0, 1) is the code-space projection of any seed state with the
right Z-eigenvalue: summing the abelian group {D(2 m alpha) D(n beta)} over all m, n projects onto
Z = +1, and one extra D(alpha) flips it to |1_L>. With the vacuum as seed,
    |mu_L>  ~  sum_{m,n}  D((2m + mu) alpha) D(n beta) |0>,
a lattice of coherent states. Finite energy (p.82) is the Gaussian envelope e^{-Delta^2 n},
which acts on a coherent state exactly:  e^{-Delta^2 n}|g> = e^{-|g|^2 (1 - e^{-2 Delta^2})/2} |g e^{-Delta^2}>.
So every codeword is built analytically in the Fock basis, without exponentiating truncated
operators. For the square code this equals the course's q-space form (peaks of width Delta,
envelope kappa = Delta) up to O(Delta^4): overlap 0.9997 at Delta = 0.3 (tests/test_codes.py).
"""

from dataclasses import dataclass

import numpy as np
import qutip as qt
from scipy.special import gammaln

from cvqec.codes.base import BosonicCode, register_code


@dataclass(frozen=True)
class Lattice:
    """A GKP lattice given by the logical displacements X = D(alpha), Z = D(beta)."""

    alpha: complex
    beta: complex
    name: str = "custom"

    def __post_init__(self):
        area = (self.beta * np.conj(self.alpha) - np.conj(self.beta) * self.alpha) / 1j
        if not np.isclose(area, np.pi):
            raise ValueError(f"need beta alpha* - beta* alpha = i pi, got i*{area:.6f}")

    @classmethod
    def square(cls) -> "Lattice":
        s = np.sqrt(np.pi / 2)
        return cls(s + 0j, 1j * s, "square")

    @classmethod
    def rectangular(cls, lam: float) -> "Lattice":
        """lam > 1 makes X longer (bit flips rarer) and Z shorter (phase flips more frequent)."""
        s = np.sqrt(np.pi / 2)
        return cls(lam * s + 0j, 1j * s / lam, f"rectangular(lam={lam})")

    @classmethod
    def hexagonal(cls) -> "Lattice":
        s = np.sqrt(np.pi / np.sqrt(3))
        return cls(s + 0j, s * np.exp(2j * np.pi / 3), "hexagonal")

    @property
    def generator_matrix(self) -> np.ndarray:
        """Columns are the logical shifts (dq, dp) of X and Z. det = pi (course units)."""
        return np.sqrt(2) * np.array([[self.alpha.real, self.beta.real],
                                      [self.alpha.imag, self.beta.imag]])

    @property
    def stabilizers(self) -> tuple[complex, complex]:
        return 2 * self.alpha, 2 * self.beta


def _damped_coherent(g: complex, delta: float, N: int) -> np.ndarray:
    """Fock coefficients of e^{-Delta^2 n} |g>, computed in log space so |g| ~ 20 is safe."""
    out = np.zeros(N, dtype=complex)
    gp = g * np.exp(-delta ** 2)
    if gp == 0:
        out[0] = 1.0
        return out
    k = np.arange(N)
    return np.exp(-abs(g) ** 2 / 2 + k * np.log(abs(gp)) - gammaln(k + 1) / 2 + 1j * k * np.angle(gp))


@register_code("gkp")
class GKPCode(BosonicCode):
    """Finite-energy GKP code. `lattice` is a Lattice or one of "square", "rectangular", "hexagonal"."""

    def __init__(self, lattice: Lattice | str = "square", delta: float = 0.3, lam: float = 1.5,
                 weight_cutoff: float = 1e-16):
        if isinstance(lattice, str):
            lattice = Lattice.rectangular(lam) if lattice == "rectangular" else getattr(Lattice, lattice)()
        self.lattice = lattice
        self.delta = delta
        self.weight_cutoff = weight_cutoff
        self.name = f"gkp_{lattice.name}"

    def codeword(self, mu: int, N: int) -> qt.Qobj:
        a, b = self.lattice.alpha, self.lattice.beta
        # Drop lattice points whose envelope weight e^{-|g|^2 (1 - e^{-2 Delta^2})/2} is negligible.
        g_max = np.sqrt(-2 * np.log(self.weight_cutoff) / (1 - np.exp(-2 * self.delta ** 2)))
        n_max = int(np.ceil(2 * g_max / min(abs(a), abs(b)))) + 2
        psi = np.zeros(N, dtype=complex)
        for m in range(-n_max, n_max + 1):
            x = (2 * m + mu) * a
            for n in range(-n_max, n_max + 1):
                z = n * b
                if abs(x + z) > g_max:
                    continue
                phase = np.exp((x * np.conj(z) - np.conj(x) * z) / 2)  # D(x) D(z) = phase * D(x + z)
                psi += phase * _damped_coherent(x + z, self.delta, N)
        return qt.Qobj(psi).unit()

    def codewords(self, N: int) -> tuple[qt.Qobj, qt.Qobj]:
        return self.codeword(0, N), self.codeword(1, N)

    def displacement_paulis(self, N: int) -> dict[str, qt.Qobj]:
        """The course's logical operators as displacements: X = D(alpha), Z = D(beta)."""
        return {"X": qt.displace(N, self.lattice.alpha), "Z": qt.displace(N, self.lattice.beta)}

    def stabilizer_ops(self, N: int) -> tuple[qt.Qobj, qt.Qobj]:
        sx, sz = self.lattice.stabilizers
        return qt.displace(N, sx), qt.displace(N, sz)
