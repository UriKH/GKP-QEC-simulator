"""The project's main circuit: CV teleportation with every imperfection we study (primer sections A, F, G).

Modes: "in" (Alice's input, carries the encoded qubit), "a" (Alice's EPR half), "b" (Bob's EPR half).

     #  component                          physics
     1  encoder (the input state)          |psi_L> of the chosen bosonic code
     2  S(r1) on a                         squeeze q
     3  S(-r2) on b                        squeeze p
     4  B(pi/4) on (a, b)                  EPR pair: Var(q_a - q_b) = e^{-2 r1}, Var(p_a + p_b) = e^{-2 r2}
     5  loss eta on b                      Bob's half travels through a lossy channel
     6  B(pi/4) on (in, a)                 Bell-measurement beam splitter
     7  R(pi/2) on in                      so that a q-homodyne reads p
     8  homodyne q on a  (eff. eta_d)      m_q = (q_a - q_in)/sqrt2
     9  homodyne q on in (eff. eta_d)      m_p = (p_in + p_a)/sqrt2
    10  X(-sqrt2 g m_q) on b               feed-forward
    11  Z(+sqrt2 g m_p) on b               feed-forward
    12  dephasing sigma_phi on b           e.g. storage or a phase-noisy channel
   13-14 SUM + syndrome homodyne + shift   Steane-type GKP correction: cvqec.decoders.gkp.SteaneDecoder

Heisenberg result (primer section G): with gain g = sqrt(eta) the output is the input after the
single-mode channel   E = Dephasing(sigma_phi) o Loss(eta) o Displacement(var_q, var_p),   where
    var_q = e^{-2 r1} + (1 - eta_d)/eta_d,     var_p = e^{-2 r2} + (1 - eta_d)/eta_d.
`effective_channel` returns that E; tests check it against the full simulation of this circuit.
"""

from dataclasses import asdict, dataclass, replace

import numpy as np
import qutip as qt

from cvqec.circuits.base import ChannelOp, Circuit, FeedForward, Gate, Measure
from cvqec.conventions import r_from_db
from cvqec.noise import Compose, Dephasing, GaussianDisplacement, PhotonLoss
from cvqec.ops.gates import beam_splitter, phase_shifter, squeezer


@dataclass
class TeleportationParams:
    r1: float                  # squeezing of the q-correlations (r1 = r2 is the symmetric case)
    r2: float
    eta: float = 1.0           # transmissivity of Bob's channel
    eta_d: float = 1.0         # homodyne detector efficiency
    sigma_phi: float = 0.0     # dephasing strength (std of the random phase)
    gain: float | None = None  # feed-forward gain; None means sqrt(eta), the value that keeps E simple

    @classmethod
    def from_db(cls, db: float, db_p: float | None = None, **kw) -> "TeleportationParams":
        """Squeezing in dB (course p.39): r = db ln10 / 20. Pass db_p for asymmetric squeezing."""
        return cls(r_from_db(db), r_from_db(db if db_p is None else db_p), **kw)

    @property
    def g(self) -> float:
        return np.sqrt(self.eta) if self.gain is None else self.gain

    def noise_variances(self) -> tuple[float, float]:
        detector = (1 - self.eta_d) / self.eta_d
        return np.exp(-2 * self.r1) + detector, np.exp(-2 * self.r2) + detector

    def as_dict(self) -> dict:
        return asdict(self)


def teleportation_circuit(p: TeleportationParams) -> Circuit:
    g = p.g
    c = Circuit(["in", "a", "b"], input_mode="in", output_mode="b", name="CV teleportation")
    c.add(Gate("S(r1) squeeze q", ("a",), lambda N: squeezer(p.r1, N)))
    c.add(Gate("S(-r2) squeeze p", ("b",), lambda N: squeezer(-p.r2, N)))
    c.add(Gate("B(pi/4) EPR source", ("a", "b"), lambda N: beam_splitter(np.pi / 4, N)))
    if p.eta < 1:
        c.add(ChannelOp(f"loss eta={p.eta:g}", "b", PhotonLoss(p.eta)))
    c.add(Gate("B(pi/4) Bell measurement", ("in", "a"), lambda N: beam_splitter(np.pi / 4, N)))
    c.add(Gate("R(pi/2) rotate p into q", ("in",), lambda N: phase_shifter(np.pi / 2, N)))
    c.add(Measure("m_q", "a", efficiency=p.eta_d))
    c.add(Measure("m_p", "in", efficiency=p.eta_d))
    # X(s) = D(s/sqrt2) and Z(s) = D(i s/sqrt2) (p.81), so the amplitudes are -g m_q and i g m_p.
    c.add(FeedForward("X(-sqrt2 g m_q)", "b", lambda o: -g * o["m_q"] + 0j))
    c.add(FeedForward("Z(+sqrt2 g m_p)", "b", lambda o: 1j * g * o["m_p"]))
    if p.sigma_phi > 0:
        c.add(ChannelOp(f"dephasing sigma={p.sigma_phi:g}", "b", Dephasing(p.sigma_phi)))
    return c


def effective_channel(p: TeleportationParams):
    """The single-mode channel the whole circuit implements (valid for gain = sqrt(eta))."""
    if p.gain is not None and not np.isclose(p.gain, np.sqrt(p.eta)):
        raise ValueError("the effective channel is derived for gain = sqrt(eta); use the Fock engine otherwise")
    var_q, var_p = p.noise_variances()
    parts = [GaussianDisplacement(var_q, var_p)]
    if p.eta < 1:
        parts.append(PhotonLoss(p.eta))
    if p.sigma_phi > 0:
        parts.append(Dephasing(p.sigma_phi))
    return Compose(parts)


def stage_params(p: TeleportationParams) -> list[tuple[str, TeleportationParams]]:
    """The circuit with its imperfections switched on one at a time (for the stage-by-stage figure).

    Stage "shift noise" is the circuit with perfect channel and storage (only finite squeezing and the
    detectors); "+ loss" adds Bob's channel loss; "+ dephasing" is the full circuit p. Each stage is a real
    circuit, so tests/test_stages.py can check effective_channel(stage) against the full three-mode simulation.
    """
    if p.gain is not None:
        raise ValueError("stages are defined for the default gain sqrt(eta)")
    stages = [("shift noise", replace(p, eta=1.0, sigma_phi=0.0))]
    if p.eta < 1:
        stages.append(("+ loss", replace(p, sigma_phi=0.0)))
    if p.sigma_phi > 0:
        stages.append(("+ dephasing", p))
    return stages


def bob_before_feedforward(p: TeleportationParams, N: int) -> qt.Qobj:
    """Bob's mode on its own after the EPR source and the channel loss, before any feed-forward.

    Heisenberg picture: q_b = (-e^{-r1} q_a0 + e^{r2} q_b0)/sqrt2 and p_b = (-e^{r1} p_a0 + e^{-r2} p_b0)/sqrt2, so
    the state is Gaussian with zero mean and Var q = (e^{-2r1} + e^{2r2})/4, Var p = (e^{2r1} + e^{-2r2})/4; loss maps
    V -> eta V + (1 - eta)/2. It does not depend on the input at all: before the feed-forward Bob holds no information.
    Built as S(s) rho_thermal S(s)^dag with matching variances (tests/test_stages.py checks it against ptrace).
    """
    vq = p.eta * (np.exp(-2 * p.r1) + np.exp(2 * p.r2)) / 4 + (1 - p.eta) / 2
    vp = p.eta * (np.exp(2 * p.r1) + np.exp(-2 * p.r2)) / 4 + (1 - p.eta) / 2
    nbar = np.sqrt(vq * vp) - 0.5                       # (2 nbar + 1)/2 = sqrt(Var q Var p)
    s = 0.25 * np.log(vp / vq)                          # S(s) maps Var q -> e^{-2s} Var q
    S = squeezer(s, N)
    return S * qt.thermal_dm(N, nbar) * S.dag()

