"""Exact multi-mode simulation in a truncated Fock space: the reference every shortcut is checked against.

Strategy (memory is the constraint: three modes at N = 25 is already a 15625-dim ket):
  1. Channels are moved later past operations on other modes (exact: they commute). In teleportation
     this moves Bob's loss and dephasing behind Alice's measurements.
  2. Until the last measurement the state is an ensemble of unnormalized kets ("branches"). A channel
     multiplies the branches by its Kraus operators. A homodyne contracts the measured mode with
     <x_i|n> on an outcome grid, adding an outcome axis (weights sqrt(dx), so sums are integrals).
     Detector inefficiency eta_d is loss before an ideal detector, with outcomes rescaled by 1/sqrt(eta_d).
  3. After the last measurement the branches are summed into one density matrix per outcome
     record. Feed-forward displacements (exact matrix elements, see ops.gates.displacement_exact)
     and remaining channels act on those, and the outcome axes are finally summed.
The trace of the result (`last_trace`) should be 1: a deficit flags too small a grid or cutoff.
"""

import numpy as np
import qutip as qt

from cvqec.circuits.base import ChannelOp, Circuit, FeedForward, Gate, Measure
from cvqec.noise import PhotonLoss
from cvqec.ops.gates import displacement_exact, on_modes
from cvqec.ops.measurement import quadrature_amplitudes


def defer_channels(ops: list) -> list:
    """Bubble every ChannelOp past later operations that act on other modes."""
    ops = list(ops)
    moved = True
    while moved:
        moved = False
        for j in range(len(ops) - 1):
            a, b = ops[j], ops[j + 1]
            if isinstance(a, ChannelOp) and not isinstance(b, ChannelOp) and not set(a.modes) & set(b.modes):
                ops[j], ops[j + 1] = b, a
                moved = True
    return ops


class FockEngine:
    def __init__(self, N: int):
        self.N = N
        self.last_trace = None

    def run(self, circuit: Circuit, inputs: dict[str, qt.Qobj] | None = None) -> qt.Qobj:
        N, inputs = self.N, inputs or {}
        ops = defer_channels(circuit.ops)
        last_measure = max((i for i, op in enumerate(ops) if isinstance(op, Measure)), default=-1)

        modes = list(circuit.modes)
        ket = qt.tensor([inputs.get(m, qt.basis(N, 0)) for m in modes]).full().reshape([N] * len(modes))
        T = ket[None]                         # axes: branch, outcome axes..., mode axes...
        labels, values = [], []               # measurement records

        def mode_axis(m):
            return 1 + len(labels) + modes.index(m)

        for op in ops[: last_measure + 1]:
            if isinstance(op, Gate):
                T = self._gate_on_ket(T, op, mode_axis, N)
            elif isinstance(op, ChannelOp):
                T = self._kraus_on_ket(T, op.channel.kraus(N), mode_axis(op.mode))
            elif isinstance(op, Measure):
                if op.efficiency < 1:
                    T = self._kraus_on_ket(T, PhotonLoss(op.efficiency).kraus(N), mode_axis(op.mode))
                x = np.linspace(*op.grid)
                A = quadrature_amplitudes(N, x, op.theta) * np.sqrt(x[1] - x[0])
                T = np.tensordot(T, A, axes=([mode_axis(op.mode)], [1]))
                T = np.moveaxis(T, -1, 1 + len(labels))
                modes.remove(op.mode)
                labels.append(op.label)
                values.append(x / np.sqrt(op.efficiency))
            else:
                raise NotImplementedError("feed-forward before the last measurement is not supported yet")

        # Density phase: one matrix per outcome record, branches summed.
        dim = N ** len(modes)
        T = T.reshape(T.shape[: 1 + len(labels)] + (dim,))
        T = np.moveaxis(T, 0, -2)                                 # (...outcomes, branch, dim)
        R = np.conj(np.swapaxes(T, -1, -2)) @ T                   # sum over branches of |v><v|, transposed
        R = np.swapaxes(R, -1, -2)
        out_shape = R.shape[:-2]
        R = R.reshape((-1, dim, dim))
        dims = [N] * len(modes)
        grids = np.meshgrid(*values, indexing="ij") if values else []
        outcomes = {lab: g.ravel() for lab, g in zip(labels, grids)}

        for op in ops[last_measure + 1:]:
            if isinstance(op, Gate):
                U = on_modes(op.build(N), [modes.index(m) for m in op.modes], dims).full()
                R = U @ R @ U.conj().T
            elif isinstance(op, ChannelOp):
                if len(modes) == 1:
                    R = op.channel.apply_many(R)
                else:
                    Ks = [on_modes(K, modes.index(op.mode), dims).full() for K in op.channel.kraus(N)]
                    R = sum(K @ R @ K.conj().T for K in Ks)
            elif isinstance(op, FeedForward):
                alphas = np.broadcast_to(np.asarray(op.alpha(outcomes)), (R.shape[0],))
                uniq, inv = np.unique(np.round(alphas, 12), return_inverse=True)
                Ds = np.array([on_modes(qt.Qobj(displacement_exact(a, N)), modes.index(op.mode), dims).full()
                               for a in uniq])
                D = Ds[inv]
                R = D @ R @ np.conj(np.swapaxes(D, 1, 2))
            else:
                raise NotImplementedError(type(op).__name__)

        rho = R.sum(axis=0)
        self.last_trace = float(np.real(np.trace(rho)))
        self.outcome_shape = out_shape
        return qt.Qobj(rho, dims=[dims, dims])

    @staticmethod
    def _gate_on_ket(T, op, mode_axis, N):
        k = len(op.modes)
        U = op.build(N).full().reshape([N] * (2 * k))
        axes = [mode_axis(m) for m in op.modes]
        T = np.tensordot(U, T, axes=(list(range(k, 2 * k)), axes))
        return np.moveaxis(T, list(range(k)), axes)

    @staticmethod
    def _kraus_on_ket(T, kraus, ax):
        Ks = np.array([K.full() for K in kraus])
        T = np.tensordot(Ks, T, axes=([2], [ax]))      # (K, d', branch, ...)
        T = np.moveaxis(T, 1, ax + 1)
        return T.reshape((-1,) + T.shape[2:])
