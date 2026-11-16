"""Metric registry. A metric is a function of a RunResult that returns a float or a dict of floats.

Add one with:
    @metric("my_metric")
    def my_metric(result): ...
Experiments list metrics by name, so no other file changes.
"""

from dataclasses import dataclass
from typing import Callable

import numpy as np
import qutip as qt

from cvqec.metrics.logical import LogicalChannel


@dataclass
class RunResult:
    logical: LogicalChannel
    noisy_plus: qt.Qobj        # E(|+_L><+_L|): the physical state the decoder receives
    N: int


METRICS: dict[str, Callable[[RunResult], float | dict]] = {}


def metric(name: str):
    def wrap(fn):
        METRICS[name] = fn
        return fn
    return wrap


def compute(names, result: RunResult) -> dict[str, float]:
    row = {}
    for name in names:
        value = METRICS[name](result)
        row.update(value if isinstance(value, dict) else {name: value})
    return row


@metric("average_fidelity")
def _average_fidelity(r):
    return r.logical.average_fidelity()


@metric("entanglement_fidelity")
def _entanglement_fidelity(r):
    return r.logical.entanglement_fidelity()


@metric("logical_infidelity")
def _infidelity(r):
    return 1 - r.logical.average_fidelity()


@metric("pauli_errors")
def _pauli(r):
    return r.logical.pauli_error_rates()


@metric("mean_photon")
def _nbar(r):
    return float(qt.expect(qt.num(r.N), r.noisy_plus))


@metric("wigner_negativity")
def _negativity(r):
    from cvqec.viz.wigner import negativity, wigner  # slow: only computed when requested
    xs, W = wigner(r.noisy_plus, extent=8.0, points=241)
    return negativity(W, xs)


@metric("truncation_tail")
def _tail(r):
    p = np.real(r.noisy_plus.diag())
    return float(p[-10:].sum())
