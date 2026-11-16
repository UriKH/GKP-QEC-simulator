"""The physical noise channels of the project (primer section I).

| channel              | Kraus / action                                         | Lindblad (duration t)        |
|----------------------|--------------------------------------------------------|------------------------------|
| PhotonLoss(eta)      | K_l = sqrt((1-eta)^l / l!) eta^{n/2} a^l               | L = a,       eta = e^{-t}    |
| Dephasing(sigma)     | rho_mn -> rho_mn e^{-(m-n)^2 sigma^2 / 2}              | L = n,       sigma^2 = t     |
| GaussianDisplacement | rho -> int G(u, v) X(u) Z(v) rho Z(v)^dag X(u)^dag     | L = p (q-noise), L = q (p-noise), t = 1 |

Loss is a beam splitter of transmissivity eta to a vacuum environment (Stinespring); the
displacement channel is the classical-noise channel that finite squeezing and inefficient detectors
produce in teleportation, with variances given in the course's units (vacuum variance = 1/2).
"""

import numpy as np
import qutip as qt
from scipy.sparse.linalg import expm_multiply
from scipy.special import gammaln

from cvqec.conventions import quad_p, quad_q
from cvqec.noise.base import Channel, register_channel


@register_channel("loss")
class PhotonLoss(Channel):
    """Pure loss with transmissivity eta (eta = 1: no loss). Number-lowering, so exact under truncation."""

    def __init__(self, eta: float, tol: float = 1e-16):
        if not 0 <= eta <= 1:
            raise ValueError("eta must lie in [0, 1]")
        self.eta, self.tol = eta, tol
        self.name = f"loss(eta={eta:g})"

    def kraus(self, N):
        # <n-l| K_l |n> = sqrt(C(n, l)) eta^{(n-l)/2} (1-eta)^{l/2}: the binomial photon-survival amplitude.
        n = np.arange(N)
        Ks = []
        for l in range(N):
            m = n[l:]
            if self.eta == 1:
                amp = np.where(l == 0, 1.0, 0.0) * np.ones(m.size)
            else:
                log_amp = 0.5 * (gammaln(m + 1) - gammaln(l + 1) - gammaln(m - l + 1)) + 0.5 * l * np.log1p(-self.eta)
                log_amp = log_amp + (0.5 * (m - l) * np.log(self.eta) if self.eta > 0 else np.where(m == l, 0, -np.inf))
                amp = np.exp(log_amp)
            if amp.max() ** 2 < self.tol and l > 0:
                break
            K = np.zeros((N, N))
            K[m - l, m] = amp
            Ks.append(qt.Qobj(K))
        return Ks

    def apply(self, rho):
        return sum((K * rho * K.dag() for K in self.kraus(rho.shape[0])), 0 * rho)

    def lindblad(self, N):
        return [qt.destroy(N)], -np.log(self.eta)


@register_channel("dephasing")
class Dephasing(Channel):
    """Random rotation e^{-i phi n} with phi ~ Normal(0, sigma^2). sigma^2 = kappa_phi * t (primer section I)."""

    def __init__(self, sigma: float, nodes: int = 40):
        self.sigma, self.nodes = sigma, nodes
        self.name = f"dephasing(sigma={sigma:g})"

    def _mask(self, N):
        d = np.subtract.outer(np.arange(N), np.arange(N))
        return np.exp(-d ** 2 * self.sigma ** 2 / 2)

    def apply(self, rho):
        return qt.Qobj(rho.full() * self._mask(rho.shape[0]), dims=rho.dims)

    adjoint = apply  # the noise distribution is symmetric, so E^dag = E

    def apply_many(self, rhos):
        return rhos * self._mask(rhos.shape[-1])

    def kraus(self, N):
        # Gauss-Hermite quadrature of the phase average: accurate while (N * sigma) is modest.
        x, w = np.polynomial.hermite.hermgauss(self.nodes)
        n = np.arange(N)
        return [qt.Qobj(np.diag(np.sqrt(wk / np.sqrt(np.pi)) * np.exp(-1j * np.sqrt(2) * self.sigma * xk * n)))
                for xk, wk in zip(x, w)]

    def lindblad(self, N):
        return [qt.num(N)], self.sigma ** 2


@register_channel("displacement")
class GaussianDisplacement(Channel):
    """Random shifts q -> q + u, p -> p + v with Var u = var_q, Var v = var_p.

    Applied exactly (up to truncation) as the Lindblad evolution dρ/dt = var_q D[p]ρ + var_p D[q]ρ for
    unit time: D[p] generates random X(u) kicks and D[q] random Z(v) kicks (primer section I). The two
    generators commute, so asymmetric noise (r1 != r2 in teleportation) needs no special treatment.
    """

    def __init__(self, var_q: float, var_p: float | None = None):
        self.var_q = var_q
        self.var_p = var_q if var_p is None else var_p
        self.name = f"displacement(var_q={self.var_q:.3g}, var_p={self.var_p:.3g})"
        self._liouvillians = {}

    def lindblad(self, N):
        return [np.sqrt(self.var_q) * quad_p(N), np.sqrt(self.var_p) * quad_q(N)], 1.0

    def _liouvillian(self, N):
        if N not in self._liouvillians:
            c_ops, _ = self.lindblad(N)
            L = qt.liouvillian(None, c_ops)
            self._liouvillians[N] = L.to("CSR").data.as_scipy()
        return self._liouvillians[N]

    def apply(self, rho):
        if self.var_q == 0 and self.var_p == 0:
            return rho
        N = rho.shape[0]
        vec = rho.full().ravel(order="F")  # qutip's column-stacking convention
        out = expm_multiply(self._liouvillian(N), vec)
        return qt.Qobj(out.reshape((N, N), order="F"), dims=rho.dims)

    adjoint = apply  # D[L] with Hermitian L is self-adjoint

    def apply_many(self, rhos):
        return np.array([self.apply(qt.Qobj(r)).full() for r in rhos])
