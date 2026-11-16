"""Truncation convergence of the headline numbers (plan: rerun at ~1.2 N and require agreement)."""

import pytest

from cvqec.circuits import TeleportationParams, effective_channel
from cvqec.codes import FourLegCat, GKPCode
from cvqec.decoders import make_decoder
from cvqec.engines.effective import run_logical


@pytest.mark.parametrize("code", [GKPCode("hexagonal", 0.3), GKPCode("square", 0.3), FourLegCat(2.25)],
                         ids=lambda c: c.name)
def test_logical_fidelity_converged_in_cutoff(code):
    ch = effective_channel(TeleportationParams.from_db(12, eta=0.95, sigma_phi=0.05))
    f = [run_logical(code, ch, make_decoder("auto", code), n).logical.average_fidelity() for n in (140, 170)]
    assert f[0] == pytest.approx(f[1], abs=1e-6)


def test_steane_grid_converged():
    """E2 (2026-11-15): the Steane round's position grid (M points per syndrome cell) is converged. At the default
    point M = 16 (used everywhere) and M = 24 agree to 1.4 % in 1-F (square; hexagonal 0.4 %, too slow for the suite)."""
    from cvqec.decoders.gkp import GKPSteaneDecoder
    ch = effective_channel(TeleportationParams.from_db(12))
    code = GKPCode("square", 0.3)
    inf = [1 - run_logical(code, ch, GKPSteaneDecoder(code, M=M), 140).logical.average_fidelity() for M in (16, 24)]
    assert inf[0] == pytest.approx(inf[1], rel=0.03)
