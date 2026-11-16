"""Guards for the Fock-space truncation, the main numerical approximation of the project.

A truncated simulation is trustworthy only if every mode keeps negligible population in its
top Fock levels. These helpers measure that population so engines and tests can refuse to
report numbers from a state that has hit the ceiling.

A small tail is necessary but not sufficient: gates built by exponentiating truncated
generators (squeezers, displacements) are themselves slightly wrong, even when the output
tail is ~1e-12. The decisive check is convergence, meaning results unchanged when N grows
(see tests/test_gates.py::test_squeezing_error_converges_with_truncation).
"""

import warnings

import numpy as np
import qutip as qt


class TruncationWarning(UserWarning):
    pass


def tail_weight(state: qt.Qobj, n_top: int = 5, mode: int = 0) -> float:
    """Population in the top `n_top` Fock levels of one mode (ket or density matrix)."""
    rho = state.ptrace(mode) if len(state.dims[0]) > 1 else state
    if rho.isket:
        rho = qt.ket2dm(rho)
    populations = np.real(rho.diag())
    return float(populations[-n_top:].sum())


def check_truncation(state: qt.Qobj, tol: float = 1e-6, n_top: int = 5) -> float:
    """Warn if any mode has more than `tol` population in its top levels; return the worst tail."""
    n_modes = len(state.dims[0])
    worst = max(tail_weight(state, n_top, m) for m in range(n_modes))
    if worst > tol:
        warnings.warn(
            f"Fock truncation reached: {worst:.2e} population in the top {n_top} levels "
            f"(tol {tol:.0e}). Increase N.",
            TruncationWarning,
            stacklevel=2,
        )
    return worst
