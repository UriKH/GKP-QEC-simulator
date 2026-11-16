"""Effective teleportation channel for a general feed-forward gain g (the library's effective_channel needs g = sqrt(eta)).

Derivation (Heisenberg picture; course units, vacuum variance 1/2). The EPR source makes
    q_a = (q_A + q_B)/sqrt2,  q_b = (q_B - q_A)/sqrt2   (and the same for p),
from a q-squeezed mode A (Var q_A = e^{-2 r1}/2, Var p_A = e^{2 r1}/2) and a p-squeezed mode B
(Var q_B = e^{2 r2}/2, Var p_B = e^{-2 r2}/2), so Var(q_a - q_b) = e^{-2 r1}, Var(p_a + p_b) = e^{-2 r2}.
Bob's half passes loss eta; the homodynes add (1 - eta_d)/eta_d of variance to each outcome (after rescaling).
After the feed-forward X(-sqrt2 g m_q), Z(sqrt2 g m_p):
    q_out = g q_in + nu_q,   nu_q = sqrt(eta) q_b - g q_a + sqrt(1-eta) q_env - sqrt2 g (detector noise)
    p_out = g p_in + nu_p,   nu_p = sqrt(eta) p_b + g p_a + sqrt(1-eta) p_env + sqrt2 g (detector noise)
    Var nu_q = [(sqrt(eta) - g)^2 e^{2 r2} + (sqrt(eta) + g)^2 e^{-2 r1}] / 4 + (1 - eta)/2 + g^2 (1 - eta_d)/eta_d
    Var nu_p = [(sqrt(eta) - g)^2 e^{2 r1} + (sqrt(eta) + g)^2 e^{-2 r2}] / 4 + (1 - eta)/2 + g^2 (1 - eta_d)/eta_d
(nu_q and nu_p are uncorrelated and independent of the input.) A Gaussian channel x -> g x + nu is
    g <= 1:  Loss(g^2)      followed by classical displacement noise Var nu - (1 - g^2)/2,
    g >  1:  Amplifier(g^2) (quantum-limited) followed by classical noise  Var nu - (g^2 - 1)/2,
and the added classical noise is >= 0 (for eta = eta_d = 1, r1 = r2 = r it is (g cosh r - sinh r)^2 resp.
(g sinh r - cosh r)^2). At g = sqrt(eta) this is exactly cvqec's effective_channel. The dephasing
component (12) acts after the feed-forward, so it is applied last.

Verified against the full three-mode Fock engine in verify_gain_channel.py.
"""
import numpy as np
import qutip as qt
from scipy.special import gammaln

from cvqec.noise import Compose, Dephasing, GaussianDisplacement, PhotonLoss
from cvqec.noise.base import Channel


class Amplifier(Channel):
    """Quantum-limited phase-insensitive amplifier, power gain G >= 1:  a -> sqrt(G) a + sqrt(G-1) e^dag.

    Kraus: <n+k| A_k |n> = sqrt(C(n+k, k)) (1 - 1/G)^{k/2} G^{-(n+1)/2}. It raises photon number, so under a
    cutoff N the trace of high-n inputs leaks; callers must check convergence (the GKP state here has
    negligible weight near the cutoff).
    """

    def __init__(self, G: float, tol: float = 1e-16):
        if G < 1:
            raise ValueError("G must be >= 1")
        self.G, self.tol = G, tol
        self.name = f"amplifier(G={G:g})"
        self._kraus = {}

    def kraus(self, N):
        if N in self._kraus:
            return self._kraus[N]
        n = np.arange(N)
        Ks = [qt.Qobj(np.diag(self.G ** (-(n + 1) / 2)))]
        if self.G > 1:
            x = 1 - 1 / self.G
            for k in range(1, N):
                m = n[: N - k]
                log_amp = 0.5 * (gammaln(m + k + 1) - gammaln(k + 1) - gammaln(m + 1)) + 0.5 * k * np.log(x) \
                    - 0.5 * (m + 1) * np.log(self.G)
                amp = np.exp(log_amp)
                if amp.max() ** 2 < self.tol:
                    break
                K = np.zeros((N, N))
                K[m + k, m] = amp
                Ks.append(qt.Qobj(K))
        self._kraus[N] = Ks
        return Ks

    def apply(self, rho):
        return sum((K * rho * K.dag() for K in self.kraus(rho.shape[0])), 0 * rho)


def noise_variances_general(p) -> tuple[float, float]:
    g, eta = p.g, p.eta
    det = g ** 2 * (1 - p.eta_d) / p.eta_d
    se = np.sqrt(eta)
    vq = ((se - g) ** 2 * np.exp(2 * p.r2) + (se + g) ** 2 * np.exp(-2 * p.r1)) / 4 + (1 - eta) / 2 + det
    vp = ((se - g) ** 2 * np.exp(2 * p.r1) + (se + g) ** 2 * np.exp(-2 * p.r2)) / 4 + (1 - eta) / 2 + det
    return vq, vp


def general_gain_channel(p) -> Compose:
    """The single-mode channel of the teleportation circuit for any gain g (p.gain), any eta, eta_d, r1, r2."""
    g = p.g
    vq, vp = noise_variances_general(p)
    quantum = abs(1 - g ** 2) / 2
    aq, ap = vq - quantum, vp - quantum
    if min(aq, ap) < -1e-12:
        raise ValueError(f"negative classical noise ({aq}, {ap}): check the derivation")
    parts = []
    if g < 1:
        parts.append(PhotonLoss(g ** 2))
    elif g > 1:
        parts.append(Amplifier(g ** 2))
    parts.append(GaussianDisplacement(max(aq, 0.0), max(ap, 0.0)))
    if p.sigma_phi > 0:
        parts.append(Dephasing(p.sigma_phi))
    return Compose(parts)
