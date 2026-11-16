"""Petz (transpose) recovery: a channel-adapted decoder that is close to optimal for any code.

    R(sigma) = P E^dag( E(P)^{-1/2} sigma E(P)^{-1/2} ) P,    P = code projector.

It knows the noise exactly, so it is not a realistic decoder. It is a benchmark: it shows how much
protection the *code* could offer, separating that from the quality of a specific decoder
(Barnum & Knill 2002 prove it is within a factor 2 of optimal in entanglement infidelity).
"""

import numpy as np
import qutip as qt

from cvqec.decoders.base import Decoder, register_decoder


@register_decoder("petz")
class PetzDecoder(Decoder):
    name = "petz"

    def __init__(self, code, rcond: float = 1e-12):
        super().__init__(code)
        self.rcond = rcond
        self.channel = None

    def bind(self, channel, N):
        self.channel = channel
        self.V = self.code.encoder(N)
        EP = channel.apply(self.V * self.V.dag()).full()
        w, U = np.linalg.eigh((EP + EP.conj().T) / 2)
        keep = w > self.rcond * w.max()
        self.EP_inv_sqrt = qt.Qobj((U[:, keep] / np.sqrt(w[keep])) @ U[:, keep].conj().T)

    def decode(self, rho):
        if self.channel is None:
            raise RuntimeError("PetzDecoder needs bind(channel, N) first")
        X = self.channel.adjoint(self.EP_inv_sqrt * rho * self.EP_inv_sqrt)
        return (self.V.dag() * X * self.V).full()
