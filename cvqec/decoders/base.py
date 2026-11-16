"""Decoder interface: a decoder turns the noisy mode back into a qubit (a linear map to 2x2 matrices).

The engine calls `bind(channel, N)` once before decoding, so a decoder that needs to know the
noise (Petz) can prepare. Decoders that don't need it ignore the call.
"""

from abc import ABC, abstractmethod

import numpy as np
import qutip as qt

from cvqec.codes.base import BosonicCode


class Decoder(ABC):
    name = "decoder"

    def __init__(self, code: BosonicCode):
        self.code = code

    def bind(self, channel, N: int) -> None:
        pass

    @abstractmethod
    def decode(self, rho: qt.Qobj) -> np.ndarray:
        """Recover the qubit from the noisy mode. Must be linear in rho."""


DECODERS: dict[str, type[Decoder]] = {}


def register_decoder(key: str):
    def wrap(cls):
        DECODERS[key] = cls
        cls.key = key
        return cls
    return wrap


def make_decoder(key: str, code: BosonicCode, **params) -> Decoder:
    if key == "auto":
        key = default_decoder_key(code)
    try:
        return DECODERS[key](code, **params)
    except KeyError:
        raise KeyError(f"unknown decoder {key!r}; registered: {sorted(DECODERS)}") from None


def default_decoder_key(code: BosonicCode) -> str:
    from cvqec.codes import FockQubit, FourLegCat, GKPCode
    if isinstance(code, GKPCode):
        return "gkp_ideal"
    if isinstance(code, FourLegCat):
        return "cat_parity"
    if isinstance(code, FockQubit):
        return "projection"
    return "petz"
