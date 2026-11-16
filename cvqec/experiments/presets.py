"""The sweeps behind the report figures. Each preset returns a list of (title, x-axis, experiments).

The report studies one code family (GKP) on three lattices at equal energy, with the unencoded
Fock qubit as the break-even reference. Pass include_cat=True to code_lineup to add the 4-leg cat.
"""

import numpy as np

from cvqec.circuits.teleportation import TeleportationParams
from cvqec.codes import FourLegCat, GKPCode
from cvqec.experiments.config import Experiment, grid
from cvqec.experiments.runner import match_mean_photon

N = 140
# Two DIFFERENT squeezing parameters appear in every run (do not confuse them):
#   DELTA   = the GKP code's own width (code / state-preparation quality): Delta = 0.3 = 10.5 dB, nbar = 5.07.
#   TELE_DB = the teleportation (EPR resource) squeezing r1 = r2, in dB: 12 dB, r = 1.38.
# Every sweep that does not vary a parameter keeps it at these defaults (eta = eta_d = 1, sigma_phi = 0, ideal decoder).
DELTA = 0.3
TELE_DB = 12.0
LAM = 1.5          # rectangular lattice stretch used in the code line-up


def code_lineup(delta: float = DELTA, include_cat: bool = False) -> list[Experiment]:
    lineup = [
        Experiment("gkp", {"lattice": "square", "delta": delta}, label="GKP square"),
        Experiment("gkp", {"lattice": "rectangular", "delta": delta, "lam": LAM}, label=f"GKP rectangular (λ={LAM:g})"),
        Experiment("gkp", {"lattice": "hexagonal", "delta": delta}, label="GKP hexagonal"),
        Experiment("fock", {}, label="bare qubit (break-even)"),
    ]
    if include_cat:
        nbar = GKPCode("square", delta).mean_photon(N)
        alpha = match_mean_photon(lambda a: FourLegCat(a), nbar, 1.0, 4.0, N)
        lineup.insert(3, Experiment("cat4", {"alpha": round(alpha, 4)}, label=f"4-leg cat (α={alpha:.2f})"))
    return lineup


def _sweep(path, values, base_noise=None, lineup=None):
    exps = []
    for e in lineup or code_lineup():
        if base_noise is not None:
            e = _set_noise(e, base_noise)
        exps += [x.with_("label", e.label) for x in grid(e, **{path: values})]
    return exps


def _set_noise(e, noise):
    e = e.with_("label", e.label)  # a copy
    e.noise = noise
    return e


def squeezing(quick=False):
    dbs = np.arange(2, 17, 2 if quick else 1).tolist()
    exps = []
    for e in code_lineup():
        for db in dbs:
            x = _set_noise(e, TeleportationParams.from_db(db))
            x.label = e.label
            exps.append(x)
    return [("Resource squeezing", "noise.r1", exps)]


def loss(quick=False):
    etas = np.linspace(0.8, 1.0, 5 if quick else 11).tolist()
    return [("Channel loss at 12 dB", "noise.eta", _sweep("noise.eta", etas, TeleportationParams.from_db(TELE_DB)))]


def detector(quick=False):
    etas = np.linspace(0.8, 1.0, 5 if quick else 11).tolist()
    return [("Detector efficiency at 12 dB", "noise.eta_d", _sweep("noise.eta_d", etas, TeleportationParams.from_db(TELE_DB)))]


def dephasing(quick=False):
    sig = np.linspace(0.0, 0.2, 5 if quick else 11).tolist()
    return [("Dephasing at 12 dB", "noise.sigma_phi", _sweep("noise.sigma_phi", sig, TeleportationParams.from_db(TELE_DB)))]


def asymmetry(quick=False):
    """Same total squeezing budget (24 dB), split unevenly between q and p correlations.

    The rectangular code (lam = 1.5) is matched to var_q / var_p = lam^4, i.e. q-noise ~7 dB stronger.
    """
    splits = np.arange(-12, 13, 4 if quick else 1)
    exps = []
    for e in code_lineup():
        for s in splits:
            x = _set_noise(e, TeleportationParams.from_db(TELE_DB + s / 2, TELE_DB - s / 2))
            x.label = e.label
            exps.append(x)
    return [("Asymmetric squeezing (dB_q - dB_p)", "noise.r1", exps)]


def decoders(quick=False):
    dbs = np.arange(8, 17, 2 if quick else 1).tolist()
    exps = []
    names = {"gkp_ideal": "ideal", "gkp_steane": "Steane (joint)", "petz": "Petz"}
    for lat in ("square", "hexagonal"):
        for dec, extra in (("gkp_ideal", {}), ("gkp_steane", {}), ("petz", {})):
            for db in dbs:
                exps.append(Experiment("gkp", {"lattice": lat, "delta": DELTA}, dec, extra,
                                       TeleportationParams.from_db(db), label=f"GKP {lat}, {names[dec]}"))
    return [("Decoders compared", "noise.r1", exps)]


def steane_ancilla(quick=False):
    deltas = np.linspace(0.15, 0.4, 4 if quick else 11).tolist()
    exps = [Experiment("gkp", {"lattice": lat, "delta": DELTA}, "gkp_steane", {"delta_ancilla": d},
                       TeleportationParams.from_db(TELE_DB), label=f"GKP {lat}, Steane (joint)")
            for lat in ("square", "hexagonal") for d in deltas]
    return [("Steane EC vs ancilla quality at 12 dB", "decoder.delta_ancilla", exps)]


def decoder_split(quick=False):
    """Separates the two costs of a realistic correction round at 12 dB (the decoder table of the report).

    ideal joint = perfect syndromes, nearest lattice point; ideal separate = perfect syndromes, Q and P corrected
    separately (parallelogram); Steane separate / joint = finite-Delta ancillas, decoded separately / jointly.
    The syndrome homodynes are ideal in every reported run (the report lists their inefficiency as not modelled;
    GKPSteaneDecoder(eta_syndrome=...) exists in the library and is tested, but is not used here).
    """
    rows = [("gkp_ideal", {}, "ideal, joint"), ("gkp_ideal", {"cell": "parallelogram"}, "ideal, separate"),
            ("gkp_steane", {"joint": False}, "Steane, separate"), ("gkp_steane", {}, "Steane, joint")]
    exps = [Experiment("gkp", {"lattice": lat, "delta": DELTA}, dec, extra, TeleportationParams.from_db(TELE_DB),
                       label=f"GKP {lat}, {name}")
            for lat in ("square", "hexagonal") for dec, extra, name in rows]
    return [("Decoder split at 12 dB", "noise.r1", exps)]


# One illustrative "realistic" operating point with every imperfection on at once (not taken from a specific experiment):
REALISTIC = dict(eta=0.95, eta_d=0.95, sigma_phi=0.05)
BETTER_ANCILLA = 0.15   # Delta_a = 0.15 = 16.5 dB, the best ancilla of the steane_ancilla sweep


def combined(quick=False):
    """All imperfections at once (REALISTIC), vs teleportation squeezing: ideal decoder for every code, plus a
    realistic Steane round (joint, ancillas Delta_a = DELTA and BETTER_ANCILLA, ideal syndrome homodynes) for the
    square and hexagonal codes."""
    dbs = np.arange(4, 17, 2 if quick else 1).tolist()
    exps = []
    for e in code_lineup():
        for db in dbs:
            x = _set_noise(e, TeleportationParams.from_db(db, **REALISTIC))
            x.label = e.label
            exps.append(x)
    for lat in ("square", "hexagonal"):
        for db in dbs:
            exps.append(Experiment("gkp", {"lattice": lat, "delta": DELTA}, "gkp_steane", {},
                                   TeleportationParams.from_db(db, **REALISTIC), label=f"GKP {lat}, Steane (joint)"))
            exps.append(Experiment("gkp", {"lattice": lat, "delta": DELTA}, "gkp_steane",
                                   {"delta_ancilla": BETTER_ANCILLA}, TeleportationParams.from_db(db, **REALISTIC),
                                   label=f"GKP {lat}, Steane, better ancillas"))
    return [("All imperfections at once", "noise.r1", exps)]


SENSITIVITY_ETAS = [0.90, 0.97]   # eta = eta_d around REALISTIC (0.95); sigma_phi kept at REALISTIC


def sensitivity(quick=False):
    """E3 (2026-11-15): is the all-at-once conclusion robust? The `combined` runs repeated with eta = eta_d = 0.90 and
    0.97 (dephasing as in REALISTIC). Labels end with ' @ eta=<value>'; the 0.95 point is the `combined` preset."""
    dbs = np.arange(4, 17, 2 if quick else 1).tolist()
    exps = []
    for eta in SENSITIVITY_ETAS:
        noise = dict(REALISTIC, eta=eta, eta_d=eta)
        tag = f" @ eta={eta:.2f}"
        for e in code_lineup():
            for db in dbs:
                x = _set_noise(e, TeleportationParams.from_db(db, **noise))
                x.label = e.label + tag
                exps.append(x)
        for lat in ("square", "hexagonal"):
            for db in dbs:
                exps.append(Experiment("gkp", {"lattice": lat, "delta": DELTA}, "gkp_steane", {},
                                       TeleportationParams.from_db(db, **noise), label=f"GKP {lat}, Steane (joint)" + tag))
                exps.append(Experiment("gkp", {"lattice": lat, "delta": DELTA}, "gkp_steane",
                                       {"delta_ancilla": BETTER_ANCILLA}, TeleportationParams.from_db(db, **noise),
                                       label=f"GKP {lat}, Steane, better ancillas" + tag))
    return [("All imperfections at once, other eta = eta_d", "noise.r1", exps)]


BUDGET_STAGES = ["squeezing only", "+ channel loss", "+ detectors", "+ dephasing", "+ Steane round"]


def budget(quick=False):
    """Error budget at TELE_DB: the imperfections of REALISTIC switched on one at a time (label = "code :: stage"),
    then the ideal decoder replaced by the realistic Steane round (Delta_a = DELTA) for the GKP codes."""
    steps = [{}, {"eta": REALISTIC["eta"]}, {"eta_d": REALISTIC["eta_d"]}, {"sigma_phi": REALISTIC["sigma_phi"]}]
    exps = []
    for e in code_lineup():
        noise = {}
        for k, step in enumerate(steps):
            noise = {**noise, **step}
            x = _set_noise(e, TeleportationParams.from_db(TELE_DB, **noise))
            x.label = f"{e.label} :: {BUDGET_STAGES[k]}"
            exps.append(x)
        if e.code == "gkp":
            exps.append(Experiment("gkp", e.code_params, "gkp_steane", {}, TeleportationParams.from_db(TELE_DB, **noise),
                                   label=f"{e.label} :: {BUDGET_STAGES[4]}"))
    return [("Error budget at 12 dB", "noise.r1", exps)]


def code_quality(quick=False):
    """Logical error vs the GKP code's own squeezing (Delta) at a fixed 12 dB teleportation channel."""
    code_dbs = np.arange(8, 15, 2 if quick else 0.5).tolist()
    exps = [Experiment("gkp", {"lattice": lat, "delta": float(10 ** (-db / 20))}, noise=TeleportationParams.from_db(TELE_DB), N=200,
                       metrics=("average_fidelity", "pauli_errors", "mean_photon", "truncation_tail"),
                       label=f"GKP {lat}")
            for lat in ("square", "hexagonal") for db in code_dbs]
    exps += [Experiment("gkp", {"lattice": "rectangular", "delta": float(10 ** (-db / 20)), "lam": LAM},
                        noise=TeleportationParams.from_db(TELE_DB), N=200,
                        metrics=("average_fidelity", "pauli_errors", "mean_photon", "truncation_tail"),
                        label=f"GKP rectangular (λ={LAM:g})") for db in code_dbs]
    return [("Code quality at 12 dB teleportation", "code.delta", exps)]


DELTA_DEPH_SIGMAS = [0.0, 0.1]   # E5: dephasing (rad) at which the code width is swept


def delta_dephasing(quick=False):
    """E5: is there a best code width under dephasing? Logical error vs the code squeezing (Delta) at a fixed 12 dB
    teleportation channel, without and with sigma_phi = 0.1 rad. Expected trade-off: a smaller Delta (more photons)
    narrows the peaks (fewer shift errors) but a rotation moves a peak at distance x by about phi*x, so far peaks
    leave their cell; hence an optimum Delta under dephasing. Bare qubit added at each sigma as the reference."""
    code_dbs = np.arange(6, 16.5, 1 if quick else 0.5).tolist()
    metrics = ("average_fidelity", "pauli_errors", "mean_photon", "truncation_tail")
    exps = []
    for sig in DELTA_DEPH_SIGMAS:
        noise = TeleportationParams.from_db(TELE_DB, sigma_phi=sig)
        for lat, extra in (("square", {}), ("rectangular", {"lam": LAM}), ("hexagonal", {})):
            name = f"GKP rectangular (λ={LAM:g})" if lat == "rectangular" else f"GKP {lat}"
            exps += [Experiment("gkp", {"lattice": lat, "delta": float(10 ** (-db / 20)), **extra}, noise=noise, N=200,
                                metrics=metrics, label=f"{name} @ sigma={sig:g}") for db in code_dbs]
        exps.append(Experiment("fock", {}, noise=noise, metrics=metrics, label=f"bare qubit @ sigma={sig:g}"))
    return [("Code width under dephasing", "code.delta", exps)]


RECT_LAMBDAS = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
RECT_BIASES = [0, -6, -12]   # dB_q - dB_p around TELE_DB (total 2 TELE_DB): symmetric, medium, strong q-noise


def rect_lambda(quick=False):
    """Why lambda = 1.5? The rectangular stretch lambda swept at three noise biases (same split as the asymmetry
    sweep, negative = noisier q). Expected: the best lambda grows with the bias; lambda = 1 (the square code) is best
    for symmetric noise. Same Delta for every lambda (nbar grows slowly with lambda, recorded as mean_photon)."""
    lams = RECT_LAMBDAS[::2] if quick else RECT_LAMBDAS
    return [("Rectangular stretch lambda", "code.lam",
             [Experiment("gkp", {"lattice": "rectangular", "delta": DELTA, "lam": lam},
                         noise=TeleportationParams.from_db(TELE_DB + s / 2, TELE_DB - s / 2),
                         metrics=("average_fidelity", "pauli_errors", "mean_photon"), label=f"bias {s:+d} dB")
              for s in RECT_BIASES for lam in lams])]


def cat_reference(quick=False):
    """Backs one sentence of Section 2.2: a 4-leg cat at the same nbar as GKP (Delta = 0.3) sent through the same
    circuit, vs the bare qubit and square GKP. Not plotted; macros only."""
    exps = []
    for e in code_lineup(include_cat=True):
        if e.code in ("cat4", "fock") or e.code_params.get("lattice") == "square":
            for db in (6, 12, 16):
                x = _set_noise(e, TeleportationParams.from_db(db))
                x.label = e.label
                exps.append(x)
    return [("4-leg cat reference", "noise.r1", exps)]


PRESETS = {f.__name__: f for f in (squeezing, loss, detector, dephasing, asymmetry, decoders, steane_ancilla,
                                   decoder_split, code_quality, combined, budget, rect_lambda, cat_reference,
                                   sensitivity, delta_dephasing)}
