"""A circuit is a list of operations on named modes. It says *what* happens, never *how* to simulate it.

Engines (cvqec/engines) decide how to simulate. Four operation types are enough for every circuit in this project:
    Gate         a unitary on one or two modes, built lazily for a given Fock cutoff N;
    ChannelOp    a noise channel (cvqec.noise) on one mode;
    Measure      a homodyne detection of q_theta with efficiency eta_d; the mode is consumed;
    FeedForward  a displacement whose amplitude is a function of earlier outcomes.
"""

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
import qutip as qt

from cvqec.noise.base import Channel


@dataclass
class Gate:
    label: str
    modes: tuple[str, ...]
    build: Callable[[int], qt.Qobj]


@dataclass
class ChannelOp:
    label: str
    mode: str
    channel: Channel

    @property
    def modes(self):
        return (self.mode,)


@dataclass
class Measure:
    label: str
    mode: str
    theta: float = 0.0
    efficiency: float = 1.0
    grid: tuple[float, float, int] = (-7.0, 7.0, 81)  # outcome discretization (min, max, points)

    @property
    def modes(self):
        return (self.mode,)


@dataclass
class FeedForward:
    """Displace `mode` by alpha(outcomes), where `outcomes` maps measurement labels to arrays."""

    label: str
    mode: str
    alpha: Callable[[dict[str, np.ndarray]], np.ndarray]

    @property
    def modes(self):
        return (self.mode,)


Op = Gate | ChannelOp | Measure | FeedForward


@dataclass
class Circuit:
    modes: list[str]
    ops: list[Op] = field(default_factory=list)
    input_mode: str | None = None
    output_mode: str | None = None
    name: str = "circuit"

    def add(self, op: Op) -> "Circuit":
        unknown = set(op.modes) - set(self.modes)
        if unknown:
            raise ValueError(f"{op.label}: unknown modes {unknown}")
        self.ops.append(op)
        return self

    def describe(self) -> str:
        """Numbered component table, the form used in the report."""
        rows = [f"{i:>2}. {type(op).__name__:<11} {op.label:<34} on {', '.join(op.modes)}"
                for i, op in enumerate(self.ops, 1)]
        return f"{self.name}  (input: {self.input_mode}, output: {self.output_mode})\n" + "\n".join(rows)
