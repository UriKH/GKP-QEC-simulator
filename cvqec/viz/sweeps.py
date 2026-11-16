"""Sweep plots for the report: logical infidelity vs one noise parameter, one line per code/decoder.

Series colors follow a fixed categorical order (colorblind-validated palette) and every series also
gets its own marker, so identity never depends on color alone. Infidelity is on a log axis because
the interesting differences between codes span orders of magnitude.
"""

import numpy as np

from cvqec.conventions import squeezing_db

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
MARKERS = ["o", "s", "^", "D", "v", "P"]
# Color follows the entity in every figure: each code keeps its color and marker; the decoder sets the line style.
ENTITY_STYLE = [("square", 0), ("rect", 1), ("hexagonal", 2), ("cat", 3), ("Fock", 4), ("bare", 4)]
DECODER_LINESTYLE = [("better ancillas", "-."), ("Steane", "--"), ("Petz", ":")]


def style_for(label: str, index: int) -> dict:
    slot = next((k for key, k in ENTITY_STYLE if key in label), index)
    ls = next((v for key, v in DECODER_LINESTYLE if key in label), "-")
    return {"color": SERIES[slot % len(SERIES)], "marker": MARKERS[slot % len(MARKERS)], "ls": ls}
AXIS_LABELS = {"noise.r1": "teleportation squeezing (dB)", "noise.eta": r"channel transmissivity $\eta$",
               "noise.eta_d": r"detector efficiency $\eta_d$", "noise.sigma_phi": r"dephasing $\sigma_\phi$ (rad)",
               "decoder.delta_ancilla": r"Steane ancilla $\Delta_a$", "code.delta": "GKP code squeezing (dB)"}


def _x_values(rows, x):
    if x == "noise.r1":  # show squeezing in dB; for asymmetric runs show the q-p split instead
        if any(abs(r["noise.r1"] - r["noise.r2"]) > 1e-9 for r in rows):
            return [squeezing_db(r["noise.r1"]) - squeezing_db(r["noise.r2"]) for r in rows], r"$\mathrm{dB}_q-\mathrm{dB}_p$"
        return [squeezing_db(r["noise.r1"]) for r in rows], AXIS_LABELS[x]
    if x == "code.delta":
        return [-20 * np.log10(r[x]) for r in rows], AXIS_LABELS[x]
    return [r[x] for r in rows], AXIS_LABELS.get(x, x)


def plot_sweep(rows, x, y="average_fidelity", title="", ax=None, legend=True):
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
    else:
        fig = ax.figure
    xs, xlabel = _x_values(rows, x)
    labels = list(dict.fromkeys(r["label"] for r in rows))  # first-seen order = fixed color order
    for i, lab in enumerate(labels):
        pts = sorted((xv, r[y]) for xv, r in zip(xs, rows) if r["label"] == lab)
        px, py = np.array(pts).T
        if y in ("average_fidelity", "entanglement_fidelity"):
            py = 1 - py
        ax.plot(px, py, lw=1.8, ms=4.5, label=lab, **style_for(lab, i))
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(f"1 - {y.replace('_', ' ')}" if "fidelity" in y else y)
    ax.set_title(title, loc="left", fontsize=11)
    ax.grid(True, which="major", color="#d9d9d4", lw=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    if legend:
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    return fig
