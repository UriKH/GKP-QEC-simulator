"""Experiment plumbing: configs, sweeps, saving, energy matching."""

import pytest

from cvqec.circuits import TeleportationParams
from cvqec.codes import FourLegCat
from cvqec.experiments import Experiment, grid, load_rows, match_mean_photon, run_experiment, save_rows
from cvqec.experiments.presets import PRESETS


def test_with_and_grid_do_not_mutate_base():
    base = Experiment("fock", {}, "projection", N=12)
    exps = grid(base, **{"noise.eta": [0.9, 1.0], "decoder": ["projection", "petz"]})
    assert len(exps) == 4
    assert base.noise.eta == 1.0 and base.decoder == "projection"
    assert {e.noise.eta for e in exps} == {0.9, 1.0}


def test_run_and_save_roundtrip(tmp_path):
    exp = Experiment("fock", {}, "projection", noise=TeleportationParams.from_db(10, eta=0.9), N=20)
    row = run_experiment(exp)
    assert 0.5 < row["average_fidelity"] < 1
    assert {"p_x", "p_y", "p_z", "noise.eta", "mean_photon"} <= set(row)
    path = save_rows([row], tmp_path / "r.csv")
    assert load_rows(path)[0]["average_fidelity"] == pytest.approx(row["average_fidelity"])


def test_energy_matching():
    alpha = match_mean_photon(lambda a: FourLegCat(a), 5.07, 1.0, 4.0, 100)
    assert FourLegCat(alpha).mean_photon(100) == pytest.approx(5.07, abs=0.01)


def test_presets_build():
    for name, make in PRESETS.items():
        for title, x, exps in make(quick=True):
            assert exps and all(isinstance(e, Experiment) for e in exps), name


def test_budget_adds_imperfections_cumulatively():
    # Error budget (report Section 4): each stage keeps every imperfection of the previous one and adds one more;
    # the last stage swaps the ideal decoder for the Steane round (ideal syndrome homodynes: eta_s is not reported).
    from cvqec.experiments.presets import BUDGET_STAGES, REALISTIC
    (_, _, exps), = PRESETS["budget"]()
    sq = {e.label.split(" :: ")[1]: e for e in exps if e.label.startswith("GKP square ::")}
    assert list(sq) == BUDGET_STAGES
    n = [sq[s].noise for s in BUDGET_STAGES]
    assert (n[0].eta, n[0].eta_d, n[0].sigma_phi) == (1.0, 1.0, 0.0)
    assert (n[1].eta, n[1].eta_d) == (REALISTIC["eta"], 1.0)
    assert (n[2].eta, n[2].eta_d, n[2].sigma_phi) == (REALISTIC["eta"], REALISTIC["eta_d"], 0.0)
    assert n[3] == n[4] and n[3].sigma_phi == REALISTIC["sigma_phi"]
    assert [sq[s].decoder for s in BUDGET_STAGES[3:]] == ["auto", "gkp_steane"]
    assert sq[BUDGET_STAGES[4]].decoder_params == {}
    # the bare qubit has no Steane round
    assert sum(e.label.startswith("bare qubit") for e in exps) == 4


def test_reported_runs_use_joint_steane_with_ideal_syndrome_detectors():
    # Every Steane run in the report is the JOINT round (default) and never sets eta_syndrome (not modelled in the
    # report); the separate round appears only in the decoder-split table, explicitly labelled "separate".
    for name, make in PRESETS.items():
        for _, _, exps in make():
            for e in exps:
                if e.decoder == "gkp_steane":
                    assert "eta_syndrome" not in e.decoder_params, (name, e.label)
                    assert e.decoder_params.get("joint", True) or "separate" in e.label, (name, e.label)


def test_default_operating_point():
    # The two squeezings are different parameters: code Delta = 0.3 (10.5 dB) and teleportation 12 dB.
    from cvqec.conventions import squeezing_db
    from cvqec.experiments.presets import DELTA, TELE_DB
    assert DELTA == 0.3 and TELE_DB == 12.0
    for name in ("loss", "detector", "dephasing", "steane_ancilla", "decoder_split"):
        for _, _, exps in PRESETS[name]():
            for e in exps:
                assert squeezing_db(e.noise.r1) == pytest.approx(TELE_DB) == squeezing_db(e.noise.r2), name
                assert e.code != "gkp" or e.code_params["delta"] == DELTA, name
