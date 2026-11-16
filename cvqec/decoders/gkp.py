"""GKP decoders (primer section J).

GKPIdealDecoder: ideal error correction followed by a perfect logical readout.
    Shift errors are only known modulo the logical lattice L = {a alpha + b beta}: the syndrome is
    the shift's position inside one cell V of L. Ideal correction assumes the *smallest* shift, so
    V is the Voronoi cell (the points closer to 0 than to any other lattice point). For square and
    rectangular codes V is a rectangle; for the hexagonal code it is a hexagon, which is where its
    advantage comes from. Tracing out the syndrome (the stabilizer subsystem decomposition) leaves
    a qubit whose Pauli expectations are sums of the characteristic function chi(l) = Tr[rho D(l)]
    over lattice points:
        <P_L> = sum_{l in coset(P)}  w_V(l) s_P(l) chi(l),
        w_V(l) = (1/|V|) int_V e^{i Omega(l, u)} d^2u        (Fourier transform of the cell),
    where coset(P) holds the points a alpha + b beta with (a mod 2, b mod 2) = (1,0) for X, (0,1) for Z,
    (1,1) for Y, and s_P = +-1 is the sign with D(l) = s_P P_L on the ideal code. For the square code
    this reduces to binning q into cells of width sqrt(pi): <Z_L> = < sgn cos(sqrt(pi) q) >.

GKPSteaneDecoder: one realistic Steane round with finite-energy ancillas (components 13-14).
    Simulated exactly on a position grid (any lattice is first mapped to the square frame by a
    Gaussian unitary). The q-round couples a |+~> ancilla with SUM (q_anc -> q_anc + q_data),
    homodynes q_anc and shifts the data back. The p-round does the same in p with a |0~>
    ancilla. The ancilla's own finite width enters the syndrome and kicks back into the data, as
    in the lab. After the round, the qubit is read out ideally (square cells in the square frame).

    joint=True (default): the two syndromes are measured first and the correction is chosen from both together,
    as the shortest shift (in the code's own phase space) consistent with them. That is the nearest-lattice-point
    rule of the ideal decoder, so the hexagonal code keeps its hexagonal cell. In the square frame it equals the
    separate correction followed by a syndrome-dependent logical X, Z or Y shift (a lattice vector), so it costs
    only a lookup table. joint=False corrects q right after the q-round and p after the p-round: the parallelogram.
    For square and rectangular lattices the two are identical.

    eta_syndrome < 1: inefficient syndrome homodynes, modelled exactly like Alice's detectors in the circuit (loss
    before an ideal detector, outcome rescaled by 1/sqrt(eta)): the read syndrome gets extra Gaussian noise of
    variance (1 - eta)/(2 eta). The ancilla is discarded, so the data's conditional state is the ideal one averaged
    over that noise: each outcome's weight becomes a Gaussian-smoothed bin.
"""

import numpy as np
import qutip as qt
from scipy.spatial import Voronoi

from cvqec.codes.gkp import GKPCode, Lattice
from cvqec.decoders.base import Decoder, register_decoder
from cvqec.ops.gates import displacement_exact, gaussian_unitary
from cvqec.ops.measurement import fock_wavefunctions, quadrature_amplitudes

PAULI = {"X": np.array([[0, 1], [1, 0]]), "Y": np.array([[0, -1j], [1j, 0]]), "Z": np.diag([1.0, -1.0])}


def phase_space(gamma: complex) -> np.ndarray:
    """(dq, dp) of D(gamma) (course p.38)."""
    return np.sqrt(2) * np.array([gamma.real, gamma.imag])


def cell_polygon(lattice: Lattice, kind: str = "voronoi") -> np.ndarray:
    """Vertices (counter-clockwise) of a fundamental cell of the logical lattice in (q, p)."""
    gx, gz = phase_space(lattice.alpha), phase_space(lattice.beta)
    if kind == "parallelogram":
        verts = np.array([(-gx - gz) / 2, (gx - gz) / 2, (gx + gz) / 2, (gz - gx) / 2])
    elif kind == "voronoi":
        pts = np.array([i * gx + j * gz for i in range(-3, 4) for j in range(-3, 4)])
        vor = Voronoi(pts)
        origin = int(np.argmin(np.linalg.norm(pts, axis=1)))
        verts = vor.vertices[vor.regions[vor.point_region[origin]]]
    else:
        raise ValueError(kind)
    order = np.argsort(np.arctan2(verts[:, 1], verts[:, 0]))
    return verts[order]


def polygon_fourier(verts: np.ndarray, k: np.ndarray) -> complex:
    """int_P e^{i k.u} d^2u for a CCW polygon, exactly, via the divergence theorem (one term per edge)."""
    k = np.asarray(k, float)
    kk = k @ k
    if kk == 0:
        x, y = verts[:, 0], verts[:, 1]
        return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    total = 0j
    for A, B in zip(verts, np.roll(verts, -1, axis=0)):
        t = B - A
        z = k @ t
        edge = 1.0 if abs(z) < 1e-12 else (np.exp(1j * z) - 1) / (1j * z)
        total += -1j * (k @ np.array([t[1], -t[0]])) / kk * np.exp(1j * (k @ A)) * edge
    return total


@register_decoder("gkp_ideal")
class GKPIdealDecoder(Decoder):
    name = "gkp_ideal"

    def __init__(self, code: GKPCode, cell: str = "voronoi"):
        super().__init__(code)
        self.cell = cell
        self._paulis = {}

    def logical_paulis(self, N: int) -> dict[str, np.ndarray]:
        """The N x N matrices P_L with <P_L> = Tr[rho P_L] (built once per cutoff)."""
        if N not in self._paulis:
            lat = self.code.lattice
            verts = cell_polygon(lat, self.cell)
            area = polygon_fourier(verts, np.zeros(2)).real
            G = np.column_stack([phase_space(lat.alpha), phase_space(lat.beta)])
            R = 2 * np.sqrt(2 * N + 1) + 2   # beyond this, D(l) has no matrix elements inside N levels
            n_max = int(np.ceil(R * np.abs(np.linalg.inv(G)).sum(axis=1).max())) + 1
            ops = {P: np.zeros((N, N), complex) for P in "XYZ"}
            for a in range(-n_max, n_max + 1):
                for b in range(-n_max, n_max + 1):
                    if a % 2 == 0 and b % 2 == 0:
                        continue                         # stabilizers and identity: no logical content
                    l = a * G[:, 0] + b * G[:, 1]
                    if np.linalg.norm(l) > R:
                        continue
                    P = "X" if b % 2 == 0 else "Z" if a % 2 == 0 else "Y"
                    # D(a alpha + b beta) = e^{i pi ab/2} X^a Z^b on the code, and X Z = -i Y.
                    s = np.exp(1j * np.pi * a * b / 2) * (-1j if P == "Y" else 1)
                    w = polygon_fourier(verts, np.array([-l[1], l[0]])).real / area
                    gamma = a * lat.alpha + b * lat.beta
                    ops[P] += w * np.real(s) * displacement_exact(gamma, N)
            self._paulis[N] = ops
        return self._paulis[N]

    def decode(self, rho):
        r = rho.full()
        out = np.trace(r) * np.eye(2) / 2
        for P, op in self.logical_paulis(r.shape[0]).items():
            out = out + np.trace(r @ op) * PAULI[P] / 2
        return out


class SquareGrid:
    """Position/momentum grid aligned with the square GKP cells, sqrt(pi) = M * h.

    G = 2 M^2 points make the discrete Fourier transform between the q and p grids exactly unitary
    (G h^2 = 2 pi). Points sit at half-integer multiples of h, so no point lies on a cell boundary.
    """

    def __init__(self, M: int = 16):
        if M % 2:
            raise ValueError("M must be even")
        self.M, self.h, self.G = M, np.sqrt(np.pi) / M, 2 * M * M
        self.x0 = (-M - 0.5) * np.sqrt(np.pi)
        self.x = self.x0 + (np.arange(self.G) + 0.5) * self.h
        self.F = np.exp(-1j * np.outer(self.x, self.x)) / np.sqrt(self.G)   # <p_k|x_j> h

    def from_fock(self, rho: np.ndarray) -> np.ndarray:
        A = quadrature_amplitudes(rho.shape[0], self.x) * np.sqrt(self.h)
        return A @ rho @ A.conj().T

    def to_p(self, R):
        return self.F @ R @ self.F.conj().T

    def to_q(self, R):
        return self.F.conj().T @ R @ self.F

    def fold(self, l: np.ndarray) -> np.ndarray:
        """Shift index r in [-M/2, M/2) equivalent to l modulo M: the syndrome in grid units."""
        return (l + self.M // 2) % self.M - self.M // 2

    def subsystem_trace(self, R: np.ndarray) -> np.ndarray:
        """Ideal square-cell readout: rho_L[mu, nu] = sum_{c = mu, c' = nu mod 2} sinc((c'-c)/2) B[c, c']."""
        M = self.M
        B = np.einsum("aibi->ab", R.reshape(2 * M, M, 2 * M, M))
        c = np.arange(-M, M)
        W = np.sinc(np.subtract.outer(c, c) / 2)
        out = np.zeros((2, 2), complex)
        for mu in range(2):
            for nu in range(2):
                sel = np.ix_(c % 2 == mu, c % 2 == nu)
                out[mu, nu] = np.sum(W[sel] * B[sel])
        return out


@register_decoder("gkp_steane")
class GKPSteaneDecoder(Decoder):
    name = "gkp_steane"

    def __init__(self, code: GKPCode, delta_ancilla: float | None = None, M: int = 16, pad: int = 60,
                 N_ancilla: int = 200, joint: bool = True, eta_syndrome: float = 1.0):
        super().__init__(code)
        if not 0 < eta_syndrome <= 1:
            raise ValueError("eta_syndrome must lie in (0, 1]")
        self.joint, self.eta_s = joint, eta_syndrome
        self.delta_a = code.delta if delta_ancilla is None else delta_ancilla
        self.grid = SquareGrid(M)
        self.pad = pad
        self.N_anc = N_ancilla
        self._frames = {}
        self._rounds = None

    # -- frame: map the code's lattice onto the square lattice ---------------------------------------
    def _to_square_frame(self, r: np.ndarray) -> np.ndarray:
        N = r.shape[0]
        M_lat = np.column_stack([phase_space(self.code.lattice.alpha),
                                 phase_space(self.code.lattice.beta)]) / np.sqrt(np.pi)
        if np.allclose(M_lat, np.eye(2)):
            return r
        Np = N + self.pad
        if Np not in self._frames:
            self._frames[Np] = gaussian_unitary(np.linalg.inv(M_lat), Np).full()
        U = self._frames[Np]
        big = np.zeros((Np, Np), complex)
        big[:N, :N] = r
        return U @ big @ U.conj().T

    # -- the two Steane rounds as precomputed conditional-state weights ------------------------------
    def _prepare_rounds(self):
        g = self.grid
        anc = GKPCode("square", self.delta_a)
        k0, k1 = anc.logical_basis(self.N_anc)
        plus, zero = (k0 + k1).unit().full().ravel(), k0.full().ravel()
        # Outcomes m_l = l h. Ancilla wavefunctions are needed at m_l - x_j (q-round) and m_l + p_k
        # (p-round); both are half-integer multiples of h, so tabulate them once in 1D and index.
        extent = np.sqrt(2 * self.N_anc + 1) + 2               # the ancilla vanishes beyond this
        L = int(np.ceil((np.abs(g.x).max() + extent) / g.h))
        l = np.arange(-L, L + 1)
        j = np.arange(g.G)

        def tabulate(ket, theta, idx, offset):
            vals = np.arange(idx.min(), idx.max() + 1)
            wf = quadrature_amplitudes(self.N_anc, (vals + offset) * g.h, theta) @ ket
            return wf[idx - vals[0]]

        a = tabulate(plus, 0.0, l[:, None] - j[None, :], -g.x0 / g.h - 0.5)
        b = tabulate(zero, np.pi / 2, l[:, None] + j[None, :], g.x0 / g.h + 0.5)
        r = g.fold(l)
        w = self._syndrome_bins(r)                          # w[s][l]: probability that true bin l reads as s
        q_weights = {s: g.h * (a.T * w[s]) @ a.conj() for s in np.unique(r)}
        p_weights = {s: g.h * (b.T * w[s]) @ b.conj() for s in np.unique(r)}
        self._rounds = (q_weights, p_weights)

    def _syndrome_bins(self, r: np.ndarray) -> dict:
        """Perfect detector: bin l is read as fold(l). Inefficient: Gaussian noise moves it to fold(l + k)."""
        g = self.grid
        if self.eta_s == 1:
            return {s: (r == s).astype(float) for s in np.unique(r)}
        var = (1 - self.eta_s) / (2 * self.eta_s)
        k = np.arange(-int(np.ceil(8 * np.sqrt(var) / g.h)) - 1, int(np.ceil(8 * np.sqrt(var) / g.h)) + 2)
        kernel = np.exp(-(k * g.h) ** 2 / (2 * var))
        kernel /= kernel.sum()
        return {s: sum(kk * (g.fold(r + kj) == s) for kj, kk in zip(k, kernel)) for s in np.unique(r)}

    def steane_round(self, R: np.ndarray) -> np.ndarray:
        """One q-round then one p-round on a q-grid density matrix."""
        if self._rounds is None:
            self._prepare_rounds()
        qw, pw = self._rounds
        g = self.grid
        # q-round: conditional state R(x, x') a(m - x) a*(m - x'), then shift the data by -fold(m).
        R = sum(np.roll(R * W, -s, axis=(0, 1)) for s, W in qw.items())
        # p-round in momentum space: R(p, p') b(m + p) b*(m + p'), then shift p by +fold(m).
        P = g.to_p(R)
        P = sum(np.roll(P * W, s, axis=(0, 1)) for s, W in pw.items())
        return g.to_q(P)

    def joint_table(self) -> np.ndarray:
        """Extra lattice shift (a, b), in units of sqrt(pi) in the square frame, for each syndrome pair.

        Syndromes s1 (q-round) and s2 (p-round) are grid indices in [-M/2, M/2). The separate correction shifts q by
        -s1 h and p by +s2 h; the joint rule adds (a sqrt(pi), b sqrt(pi)) with a, b in {-1, 0, 1} so that the implied
        error, mapped back to the code's own phase space, is as short as possible (nearest lattice point).
        """
        M, h = self.grid.M, self.grid.h
        M_lat = np.column_stack([phase_space(self.code.lattice.alpha),
                                 phase_space(self.code.lattice.beta)]) / np.sqrt(np.pi)
        s = np.arange(-M // 2, M // 2)
        table = np.zeros((M, M, 2), int)
        for i, s1 in enumerate(s):
            for j, s2 in enumerate(s):
                best = None
                for a in (0, -1, 1):          # 0 first: ties keep the separate correction
                    for b in (0, -1, 1):
                        err = np.array([(s1 - a * M) * h, -(s2 + b * M) * h])   # minus the applied shift
                        d = np.linalg.norm(M_lat @ err)
                        if best is None or d < best[0] - 1e-12:
                            best = (d, a, b)
                table[i, j] = best[1:]
        return table

    def steane_round_joint(self, R: np.ndarray) -> np.ndarray:
        """Both syndromes first, then one correction chosen from the pair (see joint_table)."""
        if self._rounds is None:
            self._prepare_rounds()
        qw, pw = self._rounds
        g, M = self.grid, self.grid.M
        table = self.joint_table()
        out = np.zeros_like(R)
        for s1, Wq in qw.items():
            P = g.to_p(R * Wq)                      # q-round outcome s1, correction postponed
            groups = {}
            for s2, Wp in pw.items():
                a, b = table[s1 + M // 2, s2 + M // 2]
                term = np.roll(P * Wp, s2 + b * M, axis=(0, 1))       # p-shift incl. the extra lattice step
                groups[a] = groups.get(a, 0) + term
            for a, Pa in groups.items():            # q-shift depends on (s1, a) only: at most 3 transforms
                out += np.roll(g.to_q(Pa), -s1 + a * M, axis=(0, 1))
        return out

    def decode(self, rho):
        r = self._to_square_frame(rho.full())
        R = self.grid.from_fock(r)
        if self.joint and self.joint_table().any():
            return self.grid.subsystem_trace(self.steane_round_joint(R))
        return self.grid.subsystem_trace(self.steane_round(R))

    def decode_without_round(self, rho):
        """Perfect readout in the square frame with no Steane round (parallelogram cells)."""
        r = self._to_square_frame(rho.full())
        return self.grid.subsystem_trace(self.grid.from_fock(r))
