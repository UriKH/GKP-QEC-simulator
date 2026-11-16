"""Bosonic codes. Add a code by subclassing BosonicCode and decorating it with @register_code."""

from cvqec.codes.base import CODES, BosonicCode, make_code, register_code
from cvqec.codes.cat import FockQubit, FourLegCat
from cvqec.codes.gkp import GKPCode, Lattice

__all__ = ["CODES", "BosonicCode", "FockQubit", "FourLegCat", "GKPCode", "Lattice", "make_code", "register_code"]
