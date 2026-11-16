"""Single-mode engine: encode, apply the circuit's effective channel, decode.

Using the derived channel instead of the three-mode circuit is what makes GKP feasible
(N ~ 140 per mode). tests/test_engines.py checks it against the full Fock engine.
"""

import numpy as np
import qutip as qt

from cvqec.codes.base import BosonicCode
from cvqec.decoders.base import Decoder
from cvqec.metrics.logical import LogicalChannel
from cvqec.metrics.registry import RunResult
from cvqec.noise.base import Channel


def run_logical(code: BosonicCode, channel: Channel, decoder: Decoder, N: int) -> RunResult:
    """Push the four operators V|i><j|V^dag through channel and decoder; return the logical channel."""
    V = code.encoder(N)
    decoder.bind(channel, N)
    outputs = np.zeros((2, 2, 2, 2), complex)
    for i in range(2):
        for j in range(2):
            ij = qt.Qobj(np.outer(np.eye(2)[i], np.eye(2)[j]))
            outputs[i, j] = decoder.decode(channel.apply(V * ij * V.dag()))
    plus = code.encode(qt.Qobj(np.ones(2) / np.sqrt(2)), N)
    return RunResult(LogicalChannel(outputs), channel.apply(plus * plus.dag()), N)
