"""Stage 2 of the report: figures, tables and number macros from results/*.json (see docs/RUNBOOK.md).

    python report/run_sweeps.py                          # stage 1, ~10 min -> results/
    python report/make_figures.py                        # stage 2, ~2 min  -> report/figures, report/generated
    python report/make_figures.py stages                 # one job only

Every panel is its own file (figures/<name>.pdf, vector; figures/png/<name>.png, 600 dpi), so LaTeX places them.
One visual language in every sweep figure:
    colour + marker = the code (square blue o, rectangular orange s, hexagonal green ^, bare qubit pink v),
    line style      = the decoder (solid = ideal, dashed = Steane round, dash-dot = Steane with better ancillas,
                      dotted = Petz), thin black solid = the analytic Gaussian model,
    dotted grey vertical line = the default operating point used by the other sweeps.
Numbers quoted in the text are macros in generated/numbers.tex (each with a comment saying what it is); the same
list, with full precision and the source file, is written to docs/results_numbers.md.
"""

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import qutip as qt
from matplotlib.lines import Line2D
from scipy.optimize import brentq
from scipy.special import erfc

from cvqec.codes import GKPCode
from cvqec.conventions import r_from_db, squeezing_db
from cvqec.decoders.gkp import cell_polygon, phase_space
from cvqec.experiments import load_rows
from cvqec.experiments import presets as P
from cvqec.viz.sweeps import MARKERS, SERIES, style_for
from cvqec.viz.wigner import wigner

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
RESULTS = REPO / "results"
FIG = GEN = None          # set in main (--out)
WIDTH = 6.5               # inches, = \textwidth of main.tex
PANEL = (3.15, 2.35)      # one sweep panel (two side by side fill \textwidth)
YLABEL = r"logical error $1-\bar F$"
GREY = "#6b6b6b"
MODEL = dict(color="#18212c", lw=0.8, zorder=5, label="Gaussian model (square)")   # thin black, drawn on top

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["STIX Two Text", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 8.5, "axes.titlesize": 9, "legend.fontsize": 7,
    "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03, "savefig.dpi": 600, "lines.linewidth": 1.5, "lines.markersize": 4,
})

CODE_LABELS = {"GKP square": "GKP square", "GKP rectangular (λ=1.5)": r"GKP rectangular ($\lambda=1.5$)",
               "GKP hexagonal": "GKP hexagonal", "bare qubit (break-even)": "bare qubit (break-even)"}
DECODER_LS = {"ideal": "-", "Steane": "--", "better": "-.", "Petz": ":"}


def rows(name):
    return load_rows(RESULTS / f"{name}.csv")


# ---------------------------------------------------------------------------------------------- figure helpers
# Paired panels share one fixed axes box (inches from the figure edges), so side-by-side plots line up exactly; their
# legend is a separate file (figures/legend_<figure>.pdf) placed once under the pair in LaTeX.
MARGINS = dict(left=0.62, right=0.24, bottom=0.48, top=0.27)


def panel(title, size=PANEL, top=None, fixed=True):
    fig = plt.figure(figsize=size)
    if fixed:
        m = {**MARGINS, **({"top": top} if top is not None else {})}
        W, H = size
        ax = fig.add_axes([m["left"] / W, m["bottom"] / H, 1 - (m["left"] + m["right"]) / W,
                           1 - (m["bottom"] + m["top"]) / H])
        # every title at the same height from the top edge, starting at the axes' left edge
        fig.text(m["left"] / W, 1 - 0.04 / H, title, ha="left", va="top", fontsize=plt.rcParams["axes.titlesize"])
    else:
        ax = fig.add_subplot()
        ax.set_title(title, loc="left", pad=4)
    return fig, ax


def save_legend(handles, name, ncol):
    """A legend on its own, shared by the panels of one figure."""
    fig = plt.figure(figsize=(WIDTH, 0.4))
    fig.legend(handles, [h.get_label() for h in handles], loc="center", ncol=ncol, frameon=False,
               handlelength=2.8, columnspacing=1.6)
    fig.savefig(FIG / f"legend_{name}.pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(FIG / "png" / f"legend_{name}.png", dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print(f"  legend {name}")


def rows_of_codes_and_decoders(codes, decoders):
    """Order handles so a 2-row legend shows codes in row 1 and decoders in row 2 (matplotlib fills by column)."""
    blank = Line2D([], [], ls="none", label=" ")
    n = max(len(codes), len(decoders))
    codes, decoders = codes + [blank] * (n - len(codes)), decoders + [blank] * (n - len(decoders))
    return [h for pair in zip(codes, decoders) for h in pair], n


def finish(ax, xlabel, ylabel=YLABEL, log=True):
    if log:
        ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, which="major", color="#dcdcd6", lw=0.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def legend_below(ax, handles=None, labels=None, ncol=2, y=-0.27):
    if handles is None:
        handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, y), ncol=ncol, frameon=False,
              handlelength=2.8, columnspacing=1.2, borderaxespad=0)


def legend_codes_decoders(ax, codes, decoders, y=-0.27):
    """Codes (colour + marker) in the left column, decoders (line style, black) in the right column."""
    blank = Line2D([], [], ls="none", label=" ")
    n = max(len(codes), len(decoders))
    hs = codes + [blank] * (n - len(codes)) + decoders + [blank] * (n - len(decoders))
    legend_below(ax, hs, [h.get_label() for h in hs], ncol=2, y=y)


def code_handle(code):
    st = style_for(code, 0)
    return Line2D([], [], color=st["color"], marker=st["marker"], ls="-", label=CODE_LABELS.get(code, code))


def decoder_handle(ls, text):
    return Line2D([], [], color="#18212c", ls=ls, lw=1.3, label=text)


def default_line(ax, x, text, ypos=0.98, ha="left"):
    """The default operating point the other sweeps use."""
    ax.axvline(x, color=GREY, ls=(0, (1, 2)), lw=0.9, zorder=0)
    ax.text(x, ypos, " " + text + " ", transform=ax.get_xaxis_transform(), fontsize=6.5, color=GREY, ha=ha,
            va="top", bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.85), zorder=7)


def series(ax, rs, xfun, labels=None, display=None):
    """One line per label (rows are already 1-F as y); colour/marker from the code, line style from the decoder."""
    labs = labels or list(dict.fromkeys(r["label"] for r in rs))
    for i, lab in enumerate(labs):
        pts = sorted((xfun(r), 1 - r["average_fidelity"]) for r in rs if r["label"] == lab)
        x, y = np.array(pts).T
        ax.plot(x, y, label=(display or CODE_LABELS).get(lab, lab), zorder=4 if lab.startswith("bare") else 3,
                **style_for(lab, i))


def save(fig, name, tight=False):
    # tight: crop to the content (stand-alone figures); otherwise keep the full fixed-size canvas so paired panels align
    with plt.rc_context({"savefig.bbox": "tight" if tight else "standard"}):
        fig.savefig(FIG / f"{name}.pdf")
        fig.savefig(FIG / "png" / f"{name}.png", dpi=600)
    plt.close(fig)
    print(f"  figure {name}")


def db_of(r):
    return squeezing_db(r["noise.r1"])


def bias_of(r):
    return squeezing_db(r["noise.r1"]) - squeezing_db(r["noise.r2"])


# ---------------------------------------------------------------------------------------------- models
def gaussian_model(var_noise, delta=P.DELTA):
    """Ideal square GKP + Gaussian shifts of variance var_noise + Delta^2/2 per quadrature.

    A quadrature fails when its shift passes sqrt(pi)/2: P = erfc(sqrt(pi) / (2 sqrt(2 s2))).
    q- and p-failures are independent: p_x = p_z = P(1-P), p_y = P^2, so 1 - F = (2/3)(2P - P^2).
    """
    s2 = var_noise + delta ** 2 / 2
    Pf = erfc(np.sqrt(np.pi) / (2 * np.sqrt(2 * s2)))
    return (2 / 3) * (2 * Pf - Pf ** 2)


def rect_matched_lambda(bias_db, tele_db=P.TELE_DB, delta=P.DELTA):
    """Gaussian-model prediction: the rectangular cell (half-widths lam sqrt(pi)/2 in q, sqrt(pi)/(2 lam) in p)
    is best when both quadratures fail equally often, lam^4 = (s_q^2 + Delta^2/2) / (s_p^2 + Delta^2/2)."""
    vq = 10 ** (-(tele_db + bias_db / 2) / 10) + delta ** 2 / 2
    vp = 10 ** (-(tele_db - bias_db / 2) / 10) + delta ** 2 / 2
    return (vq / vp) ** 0.25


# ---------------------------------------------------------------------------------------------- Wigner figures
def fig_wigner():
    N, ext = 140, 5.2
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 2.35))
    for ax, lat, title in zip(axes, ["square", "rectangular", "hexagonal"],
                              ["(a) square", r"(b) rectangular, $\lambda=1.5$", "(c) hexagonal"]):
        code = GKPCode(lat, P.DELTA, **({"lam": P.LAM} if lat == "rectangular" else {}))
        xs, W = wigner(code.logical_basis(N)[0], extent=ext, points=261)
        v = np.abs(W).max()
        ax.imshow(W, origin="lower", extent=[-ext, ext, -ext, ext], cmap="RdBu_r", vmin=-v, vmax=v,
                  interpolation="bilinear", rasterized=True)
        gx, gz = phase_space(code.lattice.alpha), phase_space(code.lattice.beta)
        pts = np.array([i * gx + j * gz for i in range(-8, 9) for j in range(-8, 9)])
        pts = pts[(np.abs(pts) < ext).all(axis=1)]
        ax.plot(pts[:, 0], pts[:, 1], ".", color="#18212c", ms=1.6)
        cell = cell_polygon(code.lattice)
        ax.fill(cell[:, 0], cell[:, 1], fill=False, ec="#18212c", lw=1.0)
        ax.set_title(title, loc="left")
        ax.set_xlabel("$q$")
        ax.set_xticks([-4, 0, 4])
        ax.set_yticks([-4, 0, 4])
    axes[0].set_ylabel("$p$")
    fig.tight_layout(w_pad=0.6)
    save(fig, "wigner_codes")


# The stage figure uses the default 12 dB teleportation squeezing with 10% channel loss and sigma_phi = 0.1 rad: the
# same values as the loss and dephasing macros (\lossNinety..., \dephTenth...), large enough to be seen by eye.
STAGE_PARAMS = dict(db=P.TELE_DB, eta=0.9, sigma_phi=0.1)
STAGE_EXTENT = 8.0   # tests/test_stages.py checks that every panel's Wigner function integrates to 1 in this window


def stage_states():
    """The panels of the stage figure: (title, state) and the fidelities quoted in the text."""
    from cvqec.circuits import TeleportationParams, bob_before_feedforward, effective_channel, stage_params
    from cvqec.decoders import GKPIdealDecoder
    from cvqec.engines.effective import run_logical
    from cvqec.noise import Identity

    N = 140
    code = GKPCode("square", P.DELTA)
    dec = GKPIdealDecoder(code)
    p = TeleportationParams.from_db(STAGE_PARAMS["db"], eta=STAGE_PARAMS["eta"], sigma_phi=STAGE_PARAMS["sigma_phi"])
    plus_L = qt.Qobj(np.ones(2) / np.sqrt(2))
    plus = qt.ket2dm(code.encode(plus_L, N))

    def F(channel):
        return run_logical(code, channel, dec, N).logical.average_fidelity()

    fids = {"A": F(Identity())}
    panels = [(f"(a) input $|+_L\\rangle$\n$\\bar F={fids['A']:.3f}$", plus)]
    nbar = p.eta * np.sinh(p.r1) ** 2
    fids["Nbar"] = nbar
    panels.append((f"(b) Bob before feed-forward\n$\\bar n={nbar:.2f}$", bob_before_feedforward(p, N)))
    names = {"shift noise": ("C", "(c) after feed-forward"), "+ loss": ("D", f"(d) + loss $\\eta={p.eta:g}$"),
             "+ dephasing": ("E", f"(e) + dephasing $\\sigma_\\phi={p.sigma_phi:g}$")}
    for label, params in stage_params(p):
        ch = effective_channel(params)
        key, text = names[label]
        fids[key] = F(ch)
        panels.append((f"{text}\n$\\bar F={fids[key]:.3f}$", ch.apply(plus)))
    final = effective_channel(p).apply(plus)
    rho_L = qt.Qobj(dec.decode(final))
    V = code.encoder(N)
    fids["Plus"] = float(np.real(plus_L.dag() * rho_L * plus_L))
    panels.append((f"(f) after ideal correction\n(re-encoded), $\\langle +|\\rho_L|+\\rangle={fids['Plus']:.3f}$",
                   V * rho_L * V.dag()))
    return panels, fids


def fig_stages():
    """The teleported square-GKP |+_L> as the circuit's imperfections are switched on one at a time.

    Every panel is a state computed by the simulator (tests/test_stages.py checks each stage against the full
    three-mode simulation). Panel titles give the logical average fidelity after ideal GKP correction.
    """
    ext = STAGE_EXTENT
    panels, _ = stage_states()
    fig, axes = plt.subplots(2, 3, figsize=(WIDTH, 5.0))
    for ax, (title, rho) in zip(axes.ravel(), panels):
        xs, W = wigner(rho, extent=ext, points=241)
        v = np.abs(W).max()
        ax.imshow(W, origin="lower", extent=[-ext, ext, -ext, ext], cmap="RdBu_r", vmin=-v, vmax=v,
                  interpolation="bilinear", rasterized=True)
        ax.set_title(title, loc="left", fontsize=8)
        ax.set_xticks([-6, 0, 6])
        ax.set_yticks([-6, 0, 6])
    for ax in axes[1]:
        ax.set_xlabel("$q$")
    for ax in axes[:, 0]:
        ax.set_ylabel("$p$")
    fig.tight_layout(h_pad=1.8, w_pad=0.4)
    save(fig, "stages")


# ---------------------------------------------------------------------------------------------- sweep figures
def fig_squeezing():
    rs = rows("squeezing")
    fig, ax = panel("(a) teleportation squeezing")
    series(ax, rs, db_of)
    db = np.linspace(2, 16, 200)
    ax.plot(db, gaussian_model(10 ** (-db / 10)), **MODEL)
    default_line(ax, P.TELE_DB, "default 12 dB", ha="right")
    finish(ax, "teleportation squeezing $r_1=r_2$ (dB)")
    ax.set_xticks(range(2, 17, 2))
    save(fig, "squeezing_tele")

    cq = rows("code_quality")
    fig, ax = panel("(b) code squeezing, teleportation at 12 dB")
    code_db = lambda r: -20 * np.log10(r["code.delta"])
    series(ax, cq, code_db, labels=["GKP square", "GKP rectangular (λ=1.5)", "GKP hexagonal"])
    xs = sorted({code_db(r) for r in cq})
    cdb = np.linspace(xs[0], xs[-1], 200)
    ax.plot(cdb, gaussian_model(10 ** (-P.TELE_DB / 10), 10 ** (-cdb / 20)), **MODEL)
    fock = [1 - r["average_fidelity"] for r in rs
            if r["label"].startswith("bare") and abs(db_of(r) - P.TELE_DB) < 1e-9][0]
    st = style_for("Fock", 0)
    ax.axhline(fock, color=st["color"], lw=1.2, label="bare qubit (break-even)")
    default_line(ax, -20 * np.log10(P.DELTA), r"default $\Delta=0.3$", ypos=0.78, ha="right")
    finish(ax, r"GKP code squeezing $-20\log_{10}\Delta$ (dB)")
    save(fig, "squeezing_code")
    save_legend([code_handle(c) for c in CODE_LABELS] + [Line2D([], [], **MODEL)], "squeezing", ncol=3)


def fig_noise():
    specs = [("loss", "noise_loss", lambda r: r["noise.eta"], "(a) channel loss, teleportation 12 dB",
              r"channel transmissivity $\eta$"),
             ("detector", "noise_detector", lambda r: r["noise.eta_d"], "(b) detector inefficiency, teleportation 12 dB",
              r"homodyne detector efficiency $\eta_d$"),
             ("dephasing", "noise_dephasing", lambda r: r["noise.sigma_phi"], "(c) dephasing, teleportation 12 dB",
              r"dephasing $\sigma_\phi$ (rad)"),
             ("asymmetry", "noise_asymmetry", bias_of, r"(d) asymmetric squeezing, $\mathrm{dB}_q+\mathrm{dB}_p=24$",
              r"$\mathrm{dB}_q-\mathrm{dB}_p$   (negative: noisier $q$)")]
    for name, out, xfun, title, xlabel in specs:
        fig, ax = panel(title)
        series(ax, rows(name), xfun)
        if name == "asymmetry":
            default_line(ax, 0, "symmetric", ha="left")
        finish(ax, xlabel)
        save(fig, out)
    save_legend([code_handle(c) for c in CODE_LABELS], "noise", ncol=4)


def fig_rect_lambda():
    rs = rows("rect_lambda")
    fig, ax = panel(r"rectangular stretch $\lambda$, $\mathrm{dB}_q+\mathrm{dB}_p=24$", fixed=False)
    shades = ["#f3b48f", "#eb6834", "#9c3b10"]
    for k, bias in enumerate(P.RECT_BIASES):
        lab = f"bias {bias:+d} dB"
        pts = sorted((r["code.lam"], 1 - r["average_fidelity"]) for r in rs if r["label"] == lab)
        x, y = np.array(pts).T
        text = "symmetric noise" if bias == 0 else rf"$\mathrm{{dB}}_q-\mathrm{{dB}}_p={bias}$"
        ax.plot(x, y, color=shades[k], marker=MARKERS[k], label=text)
        lam_star = rect_matched_lambda(bias)
        if x[0] <= lam_star <= x[-1]:
            ax.plot([lam_star], [np.interp(lam_star, x, y)], marker="*", ms=8, mfc="white", mec=shades[k], ls="none",
                    zorder=6)
    ax.plot([], [], marker="*", ms=8, mfc="white", mec="#18212c", ls="none", label="Gaussian-model optimum")
    default_line(ax, P.LAM, r"$\lambda=1.5$ (report)", ha="left")
    ax.text(1.0, 0.98, r" $\lambda=1$: square", transform=ax.get_xaxis_transform(), fontsize=6.5, color=GREY, va="top")
    finish(ax, r"rectangular stretch $\lambda$")
    legend_below(ax)
    save(fig, "rect_lambda", tight=True)


def delta_dephasing_curves():
    """E5 rows -> {sigma: {code label: (code dB array, 1-F array)}} plus {sigma: bare-qubit 1-F}."""
    curves, bare = {}, {}
    code_db = lambda r: -20 * np.log10(r["code.delta"])
    for r in rows("delta_dephasing"):
        base, sig = r["label"].split(" @ sigma=")
        if base.startswith("bare"):
            bare[float(sig)] = 1 - r["average_fidelity"]
        else:
            curves.setdefault(float(sig), {}).setdefault(base, []).append((code_db(r), 1 - r["average_fidelity"]))
    curves = {s: {lab: tuple(np.array(sorted(p)).T) for lab, p in c.items()} for s, c in curves.items()}
    return curves, bare


def best_point(x, y):
    """Minimum of a sampled curve, refined by a parabola through the three points around it (in log y)."""
    i = int(np.argmin(y))
    if 0 < i < len(y) - 1:
        a, b, c = np.polyfit(x[i - 1:i + 2], np.log(y[i - 1:i + 2]), 2)
        if a > 0:
            xm = -b / (2 * a)
            return xm, float(np.exp(np.polyval([a, b, c], xm))), True
    return x[i], y[i], False   # minimum at the edge of the sweep: no interior optimum


def fig_delta_dephasing():
    """E5: logical error vs the code width at sigma_phi = 0.1 rad (solid) and without dephasing (faint)."""
    curves, bare = delta_dephasing_curves()
    sig = max(curves)
    fig, ax = panel(rf"code width at dephasing $\sigma_\phi={sig:g}$ rad", fixed=False)
    for lab in curves[sig]:
        st = style_for(lab, 0)
        x0, y0 = curves[0.0][lab]
        ax.plot(x0, y0, color=st["color"], lw=0.9, alpha=0.45, ls="-", zorder=2)
        x, y = curves[sig][lab]
        ax.plot(x, y, color=st["color"], marker=st["marker"], ls="-", label=CODE_LABELS[lab], zorder=3)
        xm, ym, interior = best_point(x, y)
        if interior:
            ax.plot([xm], [ym], marker="*", ms=11, mfc=st["color"], mec="#18212c", mew=0.8, ls="none", zorder=6)
    xs = np.linspace(*curves[sig]["GKP square"][0][[0, -1]], 200)
    dl = 10 ** (-xs / 20)   # Gaussian model, square code: rotation adds a shift variance sigma_phi^2 / (2 Delta^2)
    ax.plot(xs, gaussian_model(10 ** (-P.TELE_DB / 10) + sig ** 2 / (2 * dl ** 2), dl), **MODEL)
    st = style_for("bare", 0)
    ax.axhline(bare[sig], color=st["color"], lw=1.2, label=rf"bare qubit, $\sigma_\phi={sig:g}$")
    ax.plot([], [], color="#18212c", lw=0.9, alpha=0.45, label=r"same code, no dephasing")
    ax.plot([], [], marker="*", ms=11, mfc="white", mec="#18212c", mew=0.8, ls="none", label=r"best $\Delta$ (fitted minimum)")
    default_line(ax, -20 * np.log10(P.DELTA), r"default $\Delta=0.3$", ypos=0.98, ha="right")
    finish(ax, r"GKP code squeezing $-20\log_{10}\Delta$ (dB)")
    legend_below(ax, ncol=2, y=-0.2)
    save(fig, "delta_dephasing", tight=True)


def delta_dephasing_numbers(M):
    """E5 macros: best code width per lattice at sigma_phi = 0.1, its infidelity, and the default-Delta point."""
    if not (RESULTS / "delta_dephasing.csv").exists():
        print("  (no results/delta_dephasing.csv: run python -m cvqec.experiments delta_dephasing)")
        return
    curves, bare = delta_dephasing_curves()
    sig = max(curves)
    src = "results/delta_dephasing.json"
    M.add("eFiveSigma", f"{sig:g}", "E5: dephasing (rad) of the code-width sweep", src)
    M.inf("eFiveFock", bare[sig], f"bare qubit, sigma_phi = {sig:g}, 12 dB", src, 3)
    for lab, short in list(SHORT.items())[:3]:
        x, y = curves[sig][lab]
        xm, ym, interior = best_point(x, y)
        M.add(f"eFiveBestDB{short}", f"{xm:.1f}" if interior else "none",
              f"E5: code squeezing (dB) that minimises 1-F of {lab} at sigma_phi = {sig:g}", src, float(xm))
        M.add(f"eFiveBestDelta{short}", f"{10 ** (-xm / 20):.2f}", f"E5: the same optimum as Delta", src,
              float(10 ** (-xm / 20)))
        M.inf(f"eFiveBest{short}", ym, f"E5: smallest 1-F of {lab} at sigma_phi = {sig:g}", src, 3)
        d0 = -20 * np.log10(P.DELTA)
        M.inf(f"eFiveDefault{short}", float(np.exp(np.interp(d0, x, np.log(y)))),
              f"E5: 1-F of {lab} at the default Delta = 0.3, sigma_phi = {sig:g}", src, 3)
        x0, y0 = curves[0.0][lab]
        M.inf(f"eFiveTop{short}", float(y[-1]), f"E5: 1-F of {lab} at the largest code squeezing "
              f"({x[-1]:g} dB), sigma_phi = {sig:g}", src, 3)
        M.inf(f"eFiveTopNoDeph{short}", float(y0[-1]), f"E5: same, without dephasing", src, 3)
        up = x >= xm   # above the optimum: where the code becomes worse than the bare qubit again
        worse = crossing(x[up], y[up], np.full(up.sum(), bare[sig]))
        M.add(f"eFiveWorseDB{short}", "none" if np.isnan(worse) else f"{worse:.1f}",
              f"E5: code squeezing (dB) above which {lab} is worse than the bare qubit, sigma_phi = {sig:g}", src)
    M.add("eFiveTopDB", f"{max(curves[sig][next(iter(SHORT))][0]):g}", "E5: largest code squeezing swept (dB)", src)


def fig_decoders():
    rs = rows("decoders")
    print(rs)
    fig, ax = panel("(a) decoders vs teleportation squeezing", top=0.52)   # same box as (b), which has a top axis
    series(ax, rs, db_of)
    default_line(ax, P.TELE_DB, "12 dB", ha="right")
    finish(ax, "teleportation squeezing $r_1=r_2$ (dB)")
    save(fig, "decoders_squeezing")

    sa = rows("steane_ancilla")
    fig, ax = panel("(b) Steane round vs ancilla width, 12 dB", top=0.52)
    series(ax, sa, lambda r: r["decoder.delta_ancilla"], display={"GKP square, Steane (joint)": "GKP square",
                                                                   "GKP hexagonal, Steane (joint)": "GKP hexagonal"})
    ideal = {lab: 1 - r["average_fidelity"] for r in rows("decoder_split")
             for lab in ("GKP square", "GKP hexagonal") if r["label"] == f"{lab}, ideal, joint"}
    for lab, v in ideal.items():
        ax.axhline(v, color=style_for(lab, 0)["color"], lw=1.0, ls="-", alpha=0.7)
    default_line(ax, P.DELTA, r"$\Delta_a$ = code $\Delta$", ha="left")
    finish(ax, r"ancilla width $\Delta_a$   (top axis: ancilla squeezing, dB)")
    top = ax.secondary_xaxis("top", functions=(lambda d: -20 * np.log10(np.maximum(d, 1e-3)),
                                               lambda s: 10 ** (-s / 20)))
    top.set_xticks([8, 10, 12, 14, 16])
    top.tick_params(labelsize=7)
    save(fig, "decoders_ancilla")
    save_legend([code_handle("GKP square"), code_handle("GKP hexagonal"), decoder_handle("-", "ideal decoder"),
                 decoder_handle("--", "Steane round (joint)"), decoder_handle(":", "Petz (benchmark)")], "decoders", ncol=5)


def fig_combined():
    rs = rows("combined")
    fig, ax = panel(r"(a) all imperfections at once")
    series(ax, rs, db_of)
    default_line(ax, P.TELE_DB, "12 dB", ha="right")
    finish(ax, "teleportation squeezing $r_1=r_2$ (dB)")
    save(fig, "combined_squeezing")

    budget = rows("budget")
    fig, ax = panel("(b) error budget, teleportation 12 dB")
    for i, code in enumerate(CODE_LABELS):
        pts = sorted((P.BUDGET_STAGES.index(r["label"].split(" :: ")[1]), 1 - r["average_fidelity"])
                     for r in budget if r["label"].split(" :: ")[0] == code)
        x, y = np.array(pts).T
        st = style_for(code, i)
        k = min(len(x), 4)            # stages 0-3: ideal decoder (solid); stage 4: Steane round (dashed)
        ax.plot(x[:k], y[:k], color=st["color"], marker=st["marker"], ls="-", label=CODE_LABELS[code])
        if len(x) > 4:
            ax.plot(x[3:], y[3:], color=st["color"], marker=st["marker"], ls="--", markevery=[1])
    ax.set_xticks(range(len(P.BUDGET_STAGES)))
    ax.set_xticklabels(["squeezing\nonly", "+ channel\nloss", "+ detec-\ntors", "+ de-\nphasing", "+ Steane\nround"])
    finish(ax, "")
    save(fig, "combined_budget")
    hs, n = rows_of_codes_and_decoders(
        [code_handle(c) for c in CODE_LABELS],
        [decoder_handle("-", "ideal decoder"), decoder_handle("--", r"Steane round, $\Delta_a=0.3$"),
         decoder_handle("-.", rf"Steane round, $\Delta_a={P.BETTER_ANCILLA:g}$")])
    save_legend(hs, "combined", ncol=n)


# ---------------------------------------------------------------------------------------------- numbers
def powers_of_ten(text):
    """'7e-07' -> 7\\times10^{-7}; every power-of-ten number is wrapped in \\ensuremath so it works in text and math."""
    m = re.fullmatch(r"(-?\d+(?:\.\d+)?)e([+-]\d+)", text)
    if m:
        mant, exp = m.group(1), int(m.group(2))
        text = (f"{mant}\\times" if mant not in ("1", "1.0") else "") + f"10^{{{exp}}}"
    if "10^{" in text and not text.startswith("\\ensuremath"):
        text = f"\\ensuremath{{{text}}}"
    return text


class Macros:
    """name -> (LaTeX text, exact value, meaning, source). Written to numbers.tex and docs/results_numbers.md."""

    def __init__(self):
        self.items = {}

    def add(self, name, text, meaning, source="", exact=None):
        assert name.isalpha(), name
        text = powers_of_ten(str(text))
        meaning = meaning.replace("λ", "lambda").replace("α", "alpha")   # numbers.tex must stay ASCII (pdfLaTeX)
        assert name not in self.items, name
        self.items[name] = (str(text), exact if exact is not None else text, meaning, source)

    def inf(self, name, value, meaning, source, digits=2):
        """An infidelity 1-F (2 significant digits) and the matching fidelity F (4 decimals) as <name>Fid."""
        self.add(name, f"{value:.2g}" if value < 0.01 else f"{value:.3f}", "1-F: " + meaning, source, value)
        self.add(name + "Fid", f"{1 - value:.4f}", "F: " + meaning, source, 1 - value)

    def write(self, tex: Path, md: Path):
        tex.write_text("% Generated by report/make_figures.py -- do not edit by hand (docs/results_numbers.md).\n"
                       + "".join(f"\\newcommand{{\\{k}}}{{{t}}}  % {m}\n" for k, (t, _, m, _) in self.items.items()))
        lines = ["# Result numbers (generated)", "",
                 "Generated by `report/make_figures.py` from `results/*.json`. Each row is a LaTeX macro in "
                 "`generated/numbers.tex`. Write `\\name` in the report, never the number.", "",
                 "| macro | in the PDF | exact value | meaning | source |", "|---|---|---|---|---|"]
        for k, (t, v, m, s) in self.items.items():
            exact = f"{v:.6g}" if isinstance(v, float) else v
            lines.append(f"| `\\{k}` | {t} | {exact} | {m} | {s} |")
        md.write_text("\n".join(lines) + "\n")


def crossing(xs, f_code, f_ref):
    """First x where the code's infidelity crosses the reference's (log-linear interpolation)."""
    d = np.log(np.asarray(f_code)) - np.log(np.asarray(f_ref))
    for i in range(len(d) - 1):
        if d[i] * d[i + 1] < 0:
            return xs[i] + (xs[i + 1] - xs[i]) * d[i] / (d[i] - d[i + 1])
    return float("nan")


def break_even_text(xs, f_code, f_ref, lo):
    x = crossing(xs, f_code, f_ref)
    if np.isnan(x):   # no crossing in the range: always better ("below lo") or never better ("none")
        return f"below {lo:g}" if np.all(np.asarray(f_code) < np.asarray(f_ref)) else "none"
    return f"{x:.1f}"


def by_label(rs, x_of):
    out = {}
    for r in rs:
        out.setdefault(r["label"], []).append((x_of(r), 1 - r["average_fidelity"]))
    return {k: np.array(sorted(v)).T for k, v in out.items()}


SHORT = {"GKP square": "Square", "GKP rectangular (λ=1.5)": "Rect", "GKP hexagonal": "Hex",
         "bare qubit (break-even)": "Fock"}


def parameter_numbers(M):
    src = "cvqec/experiments/presets.py"
    M.add("parCodeDelta", f"{P.DELTA:g}", "code width Delta (all codes, all runs except the code-quality sweep)", src)
    M.add("parCodeDB", f"{-20 * np.log10(P.DELTA):.1f}", "code squeezing in dB, -20 log10(Delta)", src)
    for lat, short, kw in (("square", "Square", {}), ("hexagonal", "Hex", {}), ("rectangular", "Rect", {"lam": P.LAM})):
        M.add(f"parNbar{short}", f"{GKPCode(lat, P.DELTA, **kw).mean_photon(P.N):.2f}", f"mean photon number, {lat} code",
              "GKPCode.mean_photon(N=140)")
    M.add("parLam", f"{P.LAM:g}", "rectangular stretch lambda in the code line-up", src)
    M.add("parTeleDB", f"{P.TELE_DB:g}", "default teleportation squeezing r1 = r2 (dB)", src)
    r = r_from_db(P.TELE_DB)
    M.add("parTeleR", f"{r:.2f}", "default teleportation squeezing parameter r", src)
    M.add("parTeleVar", f"{np.exp(-2 * r):.3f}", "shift-noise variance e^{-2r} at the default squeezing", src)
    M.add("parBiasTotal", f"{2 * P.TELE_DB:g}", "total squeezing dB_q + dB_p in the asymmetry and lambda sweeps", src)
    M.add("parN", f"{P.N}", "Fock cutoff of every run", src)
    M.add("parNcodeQuality", "200", "Fock cutoff of the code-quality sweep (smaller Delta needs more photons)", src)
    M.add("parRealEta", f"{P.REALISTIC['eta']:g}", "all-at-once runs: channel transmissivity", src)
    M.add("parRealEtaD", f"{P.REALISTIC['eta_d']:g}", "all-at-once runs: detector efficiency", src)
    M.add("parRealSigmaPhi", f"{P.REALISTIC['sigma_phi']:g}", "all-at-once runs: dephasing (rad)", src)
    M.add("parAncillaGood", f"{P.BETTER_ANCILLA:g}", "better Steane ancilla width Delta_a", src)
    M.add("parAncillaGoodDB", f"{-20 * np.log10(P.BETTER_ANCILLA):.1f}", "better Steane ancilla squeezing (dB)", src)
    M.add("parStageDB", f"{STAGE_PARAMS['db']:g}", "stage figure: teleportation squeezing (dB)", "make_figures.py")
    M.add("parStageEta", f"{STAGE_PARAMS['eta']:g}", "stage figure: channel transmissivity", "make_figures.py")
    M.add("parStageSigmaPhi", f"{STAGE_PARAMS['sigma_phi']:g}", "stage figure: dephasing (rad)", "make_figures.py")
    M.add("parSteaneM", "16", "Steane grid: grid points per syndrome cell sqrt(pi)", "cvqec/decoders/gkp.py")


def squeezing_numbers(M):
    src = "results/squeezing.json"
    sq = by_label(rows("squeezing"), db_of)
    fock = sq["bare qubit (break-even)"]
    for lab, short in SHORT.items():
        for db, word in ((6, "Six"), (12, "Twelve"), (16, "Sixteen")):
            i = int(np.argmin(np.abs(sq[lab][0] - db)))
            M.inf(f"inf{word}{short}", sq[lab][1][i], f"{lab}, ideal decoder, teleportation {db} dB", src)
            if short != "Fock":
                M.add(f"gain{word}{short}", f"{fock[1][i] / sq[lab][1][i]:.1f}",
                      f"(1-F bare qubit)/(1-F {lab}) at {db} dB", src)
        if short != "Fock":
            M.add(f"breakEven{short}", break_even_text(sq[lab][0], sq[lab][1], fock[1], 2),
                  f"teleportation dB above which {lab} beats the bare qubit (2-16 dB sweep)", src)


def code_quality_numbers(M):
    src = "results/code_quality.json"
    cq = rows("code_quality")
    d = by_label(cq, lambda r: -20 * np.log10(r["code.delta"]))
    for lab, short in list(SHORT.items())[:3]:
        for db, word in ((8, "Eight"), (10.5, "Default"), (14, "Fourteen")):
            M.inf(f"codeQ{word}{short}", float(np.interp(db, *d[lab])), f"{lab} at code {db} dB, teleportation 12 dB",
                  src)
    narrow = cq + rows("delta_dephasing")   # the two N = 200 sweeps (codes narrower than the default)
    M.add("maxTail", f"{max(r.get('truncation_tail', 0) for r in narrow):.0e}",
          "largest population in the top Fock levels (code-quality and delta_dephasing sweeps, N = 200)", src)


def noise_numbers(M):
    loss = by_label(rows("loss"), lambda r: r["noise.eta"])
    det = by_label(rows("detector"), lambda r: r["noise.eta_d"])
    deph = by_label(rows("dephasing"), lambda r: r["noise.sigma_phi"])
    for lab, short in SHORT.items():
        lo, de = float(np.interp(0.9, *loss[lab])), float(np.interp(0.9, *det[lab]))
        M.inf(f"lossNinety{short}", lo, f"{lab}, channel transmissivity 0.9, 12 dB", "results/loss.json", 3)
        M.inf(f"detNinety{short}", de, f"{lab}, detector efficiency 0.9, 12 dB", "results/detector.json", 3)
        M.inf(f"dephTenth{short}", float(np.interp(0.1, *deph[lab])), f"{lab}, sigma_phi = 0.1 rad, 12 dB",
              "results/dephasing.json", 3)
        if short != "Fock":
            M.add(f"detOverLoss{short}", f"{de / lo:.1f}", f"(1-F detectors 0.9)/(1-F channel 0.9), {lab}",
                  "results/loss.json, results/detector.json")
            x = crossing(deph[lab][0], deph[lab][1], deph["bare qubit (break-even)"][1])
            M.add(f"dephasingCross{short}", f"{x:.3f}", f"sigma_phi (rad) where {lab} falls to the bare qubit",
                  "results/dephasing.json", x)
    fock = "bare qubit (break-even)"
    x = crossing(loss[fock][0], loss["GKP rectangular (λ=1.5)"][1], loss[fock][1])
    M.add("lossCrossRect", f"{x:.2f}", "channel transmissivity below which GKP rectangular is worse than the bare qubit",
          "results/loss.json", x)
    rise = float(np.interp(0.1, *deph[fock])) - float(deph[fock][1][0])
    M.add("dephFockRiseSim", f"{rise:.4f}", "bare qubit: increase of 1-F from sigma_phi = 0 to 0.1 rad (simulated)",
          "results/dephasing.json", rise)
    pz = (1 - np.exp(-0.1 ** 2 / 2)) / 2   # a random rotation phi multiplies |1> by e^{-i phi}: Z error probability
    M.add("dephFockPz", f"{pz:.4f}", "bare qubit: Z-error probability (1-e^{-sigma^2/2})/2 at sigma_phi = 0.1", "formula", pz)
    M.add("dephFockRiseModel", f"{2 * pz / 3:.4f}", "bare qubit: predicted rise of 1-F, (2/3) p_z at sigma_phi = 0.1",
          "formula", 2 * pz / 3)


def asymmetry_numbers(M):
    src = "results/asymmetry.json"
    asym = by_label(rows("asymmetry"), bias_of)
    for lab, short in SHORT.items():
        i = int(np.argmin(asym[lab][1]))
        M.add(f"asymBest{short}", f"{asym[lab][0][i]:+.0f}", f"bias dB_q - dB_p where {lab} is best (1 dB grid)", src)
        M.inf(f"asymBestInf{short}", asym[lab][1][i], f"{lab} at its best bias", src, 3)
    rect = asym["GKP rectangular (λ=1.5)"]
    best = float(rect[0][int(np.argmin(rect[1]))])
    for lab, short in SHORT.items():
        M.inf(f"asymAtRectBest{short}", float(np.interp(best, *asym[lab])), f"{lab} at the rectangular code's best bias",
              src, 3)
    hexa = asym["GKP hexagonal"]
    x = crossing(rect[0], rect[1], np.interp(rect[0], *hexa))
    M.add("asymRectBeatsHex", f"{x:+.1f}", "bias dB_q - dB_p below which GKP rectangular beats the hexagonal code", src, x)
    pred = brentq(lambda s: rect_matched_lambda(s) - P.LAM, -30, 0)
    M.add("asymPredRect", f"{pred:+.1f}", "Gaussian-model bias at which lambda = 1.5 is matched", "make_figures.py", pred)


def rect_lambda_numbers(M):
    src = "results/rect_lambda.json"
    rs = rows("rect_lambda")
    words = {0: "Zero", -6: "Six", -12: "Twelve"}
    for bias in P.RECT_BIASES:
        pts = sorted((r["code.lam"], 1 - r["average_fidelity"]) for r in rs if r["label"] == f"bias {bias:+d} dB")
        x, y = np.array(pts).T
        i = int(np.argmin(y))
        w = words[bias]
        M.add(f"rectBestLam{w}", f"{x[i]:.1f}", f"best lambda (0.1 grid) at bias {bias} dB", src, float(x[i]))
        M.inf(f"rectBestInf{w}", y[i], f"rectangular code at its best lambda, bias {bias} dB", src, 3)
        M.inf(f"rectLamOneInf{w}", float(y[0]), f"lambda = 1 (square) at bias {bias} dB", src, 3)
        M.inf(f"rectLamUsedInf{w}", float(np.interp(P.LAM, x, y)), f"lambda = 1.5 at bias {bias} dB", src, 3)
        M.add(f"rectPredLam{w}", f"{rect_matched_lambda(bias):.2f}", f"Gaussian-model best lambda at bias {bias} dB",
              "make_figures.py")
    lam_max = max(r["code.lam"] for r in rs)
    nb = [r["mean_photon"] for r in rs if r["code.lam"] == lam_max][0]
    M.add("rectNbarMax", f"{nb:.2f}", f"mean photon number of the rectangular code at lambda = {lam_max:g}", src)


def decoder_numbers(M):
    dec = by_label(rows("decoders"), db_of)
    names = (("GKP square, ideal", "SqIdeal"), ("GKP square, Steane (joint)", "SqSteane"), ("GKP square, Petz", "SqPetz"),
             ("GKP hexagonal, ideal", "HexIdeal"), ("GKP hexagonal, Steane (joint)", "HexSteane"),
             ("GKP hexagonal, Petz", "HexPetz"))
    for lab, short in names:
        for db, word in ((12, "Twelve"), (16, "")):
            M.inf(f"dec{word}{short}", float(np.interp(db, *dec[lab])), f"{lab}, teleportation {db} dB",
                  "results/decoders.json")
    for lab, short in (("GKP square, Steane (joint)", "Square"), ("GKP hexagonal, Steane (joint)", "Hex")):
        st = by_label(rows("steane_ancilla"), lambda r: r["decoder.delta_ancilla"])[lab]
        M.inf(f"steaneBest{short}", st[1][0], f"{lab}, ancilla Delta_a = {st[0][0]:g}, 12 dB",
              "results/steane_ancilla.json", 3)
        M.inf(f"steaneThirty{short}", float(np.interp(0.3, *st)), f"{lab}, ancilla Delta_a = 0.3, 12 dB",
              "results/steane_ancilla.json", 3)


def decoder_split_table(M):
    """Table + macros from results/decoder_split.csv: 1-F at 12 dB for each way of correcting (square, hexagonal)."""
    order = [("ideal, joint", "IdealJoint", "ideal syndromes, joint (nearest lattice point)"),
             ("ideal, separate", "IdealSep", "ideal syndromes, $Q$ and $P$ separately"),
             ("Steane, separate", "SteaneSep", r"Steane round ($\Delta_a=0.3$), $Q$ and $P$ separately"),
             ("Steane, joint", "SteaneJoint", r"Steane round ($\Delta_a=0.3$), joint")]
    inf = {r["label"]: 1 - r["average_fidelity"] for r in rows("decoder_split")}
    lines = [r"\begin{tabular}{|l|cc|cc|}", r"\toprule",
             r" & \multicolumn{2}{c|}{square} & \multicolumn{2}{c|}{hexagonal}\\",
             r"correction at 12\,dB & $1-\bar F$ & $\bar F$ & $1-\bar F$ & $\bar F$\\", r"\midrule"]
    for key, short, text in order:
        sq, hx = inf[f"GKP square, {key}"], inf[f"GKP hexagonal, {key}"]
        M.inf(f"split{short}Sq", sq, f"square, {key}, 12 dB", "results/decoder_split.json", 3)
        M.inf(f"split{short}Hex", hx, f"hexagonal, {key}, 12 dB", "results/decoder_split.json", 3)
        lines.append(f"{text} & \\num{{{sq:.3f}}} & \\num{{{1 - sq:.4f}}} & \\num{{{hx:.3f}}} & \\num{{{1 - hx:.4f}}} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    (GEN / "decoder_table.tex").write_text("\n".join(lines) + "\n")


def combined_numbers(M):
    """Macros for the all-at-once results (results/combined.json, results/budget.json)."""
    comb = by_label(rows("combined"), db_of)
    short = {**SHORT, "GKP square, Steane (joint)": "SqSteane", "GKP hexagonal, Steane (joint)": "HexSteane",
             "GKP square, Steane, better ancillas": "SqSteaneGood",
             "GKP hexagonal, Steane, better ancillas": "HexSteaneGood"}
    fock = comb["bare qubit (break-even)"]
    for lab, sh in short.items():
        M.inf(f"combTwelve{sh}", float(np.interp(12, *comb[lab])), f"{lab}, all imperfections, 12 dB",
              "results/combined.json", 3)
        if sh != "Fock":
            xs, f_code = comb[lab]
            M.add(f"combBreakEven{sh}", break_even_text(xs, f_code, np.interp(xs, *fock), 4),
                  f"teleportation dB above which {lab} beats the bare qubit, all imperfections (4-16 dB)",
                  "results/combined.json")
    words = ["Sqz", "Loss", "Det", "Deph", "Steane"]
    for r in rows("budget"):
        code, stage = r["label"].split(" :: ")
        M.inf(f"budget{SHORT[code]}{words[P.BUDGET_STAGES.index(stage)]}", 1 - r["average_fidelity"],
              f"{code}, error budget stage '{stage}', 12 dB", "results/budget.json", 3)


def validation_numbers(M):
    """The central check (tests/test_engines.py): full three-mode circuit vs the single-mode channel E."""
    from cvqec.circuits import TeleportationParams, effective_channel, teleportation_circuit
    from cvqec.engines.fock import FockEngine

    n = 25
    p = TeleportationParams(0.6, 0.4, eta=0.9, eta_d=0.95, sigma_phi=0.15)
    psi = (qt.basis(n, 0) + 0.7j * qt.basis(n, 1)).unit()
    diff = FockEngine(n).run(teleportation_circuit(p), {"in": psi}) - effective_channel(p).apply(qt.ket2dm(psi))
    td = 0.5 * np.abs(np.linalg.eigvalsh(diff.full())).sum()
    mant, exp = f"{td:.1e}".split("e")
    src = "make_figures.validation_numbers (as tests/test_engines.py)"
    M.add("valTraceDist", f"{mant}\\times10^{{{int(exp)}}}", "trace distance, full 3-mode circuit vs channel E", src, td)
    for k, v, m in (("valEta", p.eta, "eta"), ("valEtaD", p.eta_d, "eta_d"), ("valRone", p.r1, "r1"),
                    ("valRtwo", p.r2, "r2"), ("valSigmaPhi", p.sigma_phi, "sigma_phi")):
        M.add(k, f"{v:g}", f"validation run: {m}", src)
    M.add("valRoneDB", f"{squeezing_db(p.r1):.1f}", "validation run: r1 in dB", src)
    M.add("valRtwoDB", f"{squeezing_db(p.r2):.1f}", "validation run: r2 in dB", src)


def stage_numbers(M):
    _, f = stage_states()
    src = "make_figures.stage_states"
    for k, text in (("C", "(c) after feed-forward"), ("D", "(d) + loss"), ("E", "(e) + dephasing")):
        M.inf(f"stageInf{k}", 1 - f[k], f"stage figure panel {text}", src, 3)
    M.add("stageNbar", f"{f['Nbar']:.2f}", "stage figure (b): Bob's mean photon number before the feed-forward", src)
    M.add("stagePlus", f"{f['Plus']:.3f}", "stage figure (f): <+|rho_L|+> after ideal correction", src)


def pauli_table(M):
    target = r_from_db(P.TELE_DB)
    lines = []
    for r in rows("squeezing"):
        if abs(r["noise.r1"] - target) < 1e-9:
            lab = CODE_LABELS[r["label"]]
            vals = [r["p_x"], r["p_y"], r["p_z"], 1 - r["average_fidelity"]]
            lines.append(lab + "".join(f" & \\num{{{v:.2e}}}" for v in vals) + f" & \\num{{{r['average_fidelity']:.4f}}} \\\\")
            sh = SHORT[r["label"]]
            for k, v in zip(("px", "py", "pz"), vals):
                M.add(f"pauli{k.upper()}{sh}", f"{v:.2g}", f"{r['label']}: {k} at 12 dB", "results/squeezing.json", v)
    (GEN / "pauli_table.tex").write_text("\\begin{tabular}{|l|c|c|c|c|c|}\n\\toprule\n"
                                         "code & $p_x$ & $p_y$ & $p_z$ & $1-\\bar F$ & $\\bar F$\\\\\n\\midrule\n"
                                         + "\n".join(lines) + "\n\\bottomrule\n\\end{tabular}\n")


def cat_numbers(M):
    rs = rows("cat_reference")
    for r in rs:
        if r["label"].startswith("4-leg cat"):
            db = round(db_of(r))
            word = {6: "Six", 12: "Twelve", 16: "Sixteen"}[db]
            M.inf(f"cat{word}", 1 - r["average_fidelity"], f"{r['label']}, same nbar as GKP, teleportation {db} dB",
                  "results/cat_reference.json", 3)
    M.add("catAlpha",
          [x["label"] for x in rs if x["label"].startswith("4-leg")][0].split("=")[1].rstrip(")"),
          "4-leg cat amplitude alpha with the same nbar as GKP Delta = 0.3", "results/cat_reference.json")


def appendix_numbers(M):
    """Closed-form numbers quoted in the derivations appendix (the derivations appendix of the report)."""
    src = "analytic (appendix_derivations.tex)"
    d, r, eta = P.DELTA, r_from_db(P.TELE_DB), STAGE_PARAMS["eta"]
    M.add("appPeakVar", f"{d ** 2 / 2:.3f}", "variance Delta^2/2 of one GKP peak", src)
    M.add("appNbarThermal", f"{1 / np.expm1(2 * d ** 2):.2f}", "GKP nbar from the thermal argument, 1/(e^{2 Delta^2}-1)", src)
    M.add("appNbarApprox", f"{1 / (2 * d ** 2) - 0.5:.2f}", "leading terms 1/(2 Delta^2) - 1/2", src)
    M.add("appRadius", f"{1 / d:.1f}", "typical distance 1/Delta of a GKP peak from the origin", src)
    M.add("appRotThreshold", f"{np.sqrt(np.pi) * d / 2:.2f}", "rotation angle sqrt(pi) Delta/2 that moves a peak at 1/Delta by sqrt(pi)/2", src)
    M.add("appBobNbarIdeal", f"{np.sinh(r) ** 2:.2f}", "Bob's thermal nbar sinh^2 r at the default squeezing, no loss", src)
    duan = (1 + eta) * np.cosh(2 * r) - 2 * np.sqrt(eta) * np.sinh(2 * r) + (1 - eta)
    M.add("appDuanSum", f"{duan:.2f}", "Duan sum after loss (stage eta) at the default squeezing; separable needs >= 2", src)
    M.add("appCohTeleFid", f"{1 / (1 + np.exp(-2 * r)):.3f}", "coherent-state teleportation fidelity 1/(1+e^{-2r}), default squeezing", src)
    M.add("appDetVarNinety", f"{(1 - 0.9) / 0.9:.3f}", "shift variance (1-eta_d)/eta_d added by detectors with eta_d = 0.9", src)
    M.add("appLossVarNinety", f"{(1 - 0.9) / (2 * 0.9):.3f}", "shift variance (1-eta)/(2 eta) added by channel loss eta = 0.9 (rescaled units)", src)
    s2 = np.exp(-2 * r) + d ** 2 / 2
    pf = erfc(np.sqrt(np.pi) / (2 * np.sqrt(2 * s2)))
    M.add("appGaussVar", f"{s2:.3f}", "Gaussian model: total shift variance e^{-2r} + Delta^2/2 at the default squeezing", src)
    M.add("appGaussFail", f"{pf:.4f}", "Gaussian model: one-quadrature failure probability at the default squeezing", src)
    M.inf("appGaussInf", float(gaussian_model(np.exp(-2 * r))), "Gaussian model, square code, default squeezing", src)
    # 2026-11-16: the model with every imperfection as an extra shift variance, against the simulation (square code)
    tv = np.exp(-2 * r)
    M.inf("modelDetNinety", float(gaussian_model(tv + (1 - 0.9) / 0.9)), "Gaussian model, detectors eta_d = 0.9", src, 3)
    M.inf("modelLossNinety", float(gaussian_model(tv + (1 - 0.9) / (2 * 0.9))),
          "Gaussian model, channel loss eta = 0.9 (vacuum noise only, the shrinking is ignored)", src, 3)
    M.inf("modelDephTenth", float(gaussian_model(tv + 0.1 ** 2 / (2 * d ** 2))),
          "Gaussian model, dephasing 0.1 rad as an extra shift variance sigma_phi^2/(2 Delta^2)", src, 3)
    M.add("modelRotVarTenth", f"{0.1 ** 2 / (2 * d ** 2):.3f}", "extra shift variance sigma_phi^2/(2 Delta^2), 0.1 rad, Delta = 0.3", src)
    M.add("modelDeltaOpt", f"{np.sqrt(0.1):.2f}", "Gaussian-model best Delta = sqrt(sigma_phi) at 0.1 rad", src)
    M.add("modelDeltaOptDB", f"{-20 * np.log10(np.sqrt(0.1)):.1f}", "the same best Delta in dB, -20 log10 sqrt(sigma_phi)", src)


def duan_sum(r, eta):
    return (1 + eta) * np.cosh(2 * r) - 2 * np.sqrt(eta) * np.sinh(2 * r) + (1 - eta)


def duan_numbers(M):
    """Largest Duan sum over every (squeezing, eta) the report simulates with r1 = r2 (separable needs >= 2)."""
    pts = [(db, 1.0) for db in range(2, 17)] + [(P.TELE_DB, e) for e in np.arange(0.8, 1.0001, 0.02)]
    pts += [(db, e) for db in range(4, 17) for e in [P.REALISTIC["eta"]] + list(P.SENSITIVITY_ETAS)]
    vals = [(duan_sum(r_from_db(db), e), db, e) for db, e in pts]
    v, db, e = max(vals)
    M.add("appDuanMax", f"{v:.2f}", f"largest Duan sum over the simulated points (at {db:g} dB, eta = {e:.2f})",
          "analytic (appendix)")


SENS_SERIES = [("bare qubit (break-even)", "bare qubit", "Fock"), ("GKP square", "square, ideal", "Square"),
               ("GKP hexagonal", "hexagonal, ideal", "Hex"),
               ("GKP square, Steane (joint)", r"square, Steane ($\Delta_a=0.3$)", "SqSteane"),
               ("GKP hexagonal, Steane (joint)", r"hexagonal, Steane ($\Delta_a=0.3$)", "HexSteane"),
               ("GKP square, Steane, better ancillas", r"square, Steane ($\Delta_a=0.15$)", "SqSteaneGood"),
               ("GKP hexagonal, Steane, better ancillas", r"hexagonal, Steane ($\Delta_a=0.15$)", "HexSteaneGood")]


def sensitivity_numbers(M):
    """E3: the all-at-once runs at eta = eta_d = 0.90 / 0.95 (combined) / 0.97 -> macros sens* and sensitivity_table.tex."""
    if not (RESULTS / "sensitivity.csv").exists():
        print("  (no results/sensitivity.csv: run report/run_sweeps.py sensitivity)")
        return
    data = {P.REALISTIC["eta"]: by_label(rows("combined"), db_of)}
    for lab, curve in by_label(rows("sensitivity"), db_of).items():
        base, eta = lab.split(" @ eta=")
        data.setdefault(float(eta), {})[base] = curve
    etas = sorted(data)
    word = {0.90: "Ninety", 0.95: "NinetyFive", 0.97: "NinetySeven"}
    lines = [r"\begin{tabular}{|l|" + "c" * len(etas) + "|}", r"\toprule",
             r"$1-\bar F$ at \parTeleDB\,dB & " + " & ".join(f"$\\eta=\\eta_d={e:.2f}$" for e in etas) + r"\\",
             r"\midrule"]
    for lab, text, sh in SENS_SERIES:
        vals = [float(np.interp(P.TELE_DB, *data[e][lab])) for e in etas]
        for e, v in zip(etas, vals):
            M.inf(f"sens{sh}{word[round(e, 2)]}", v, f"{lab}, all imperfections with eta = eta_d = {e:.2f}, 12 dB",
                  "results/sensitivity.json + combined.json", 3)
        lines.append(text + "".join(f" & \\num{{{v:.3f}}}" for v in vals) + r" \\")
    lines.append(r"\midrule")
    for lab, text, sh in SENS_SERIES[3:5]:
        cells = []
        for e in etas:
            xs, f = data[e][lab]
            be = break_even_text(xs, f, np.interp(xs, *data[e]["bare qubit (break-even)"]), 4)
            M.add(f"sensBreakEven{sh}{word[round(e, 2)]}", be,
                  f"dB above which {lab} beats the bare qubit, eta = eta_d = {e:.2f}", "results/sensitivity.json")
            cells.append(be)
        lines.append(f"break-even (dB), {text}" + "".join(f" & {c}" for c in cells) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    (GEN / "sensitivity_table.tex").write_text("\n".join(lines) + "\n")


def summary_table(M):
    """Summary of the single-error sweeps at one strength each (ideal decoder, 12 dB) + the all-at-once run.

    Writes generated/summary_table.tex (best GKP code per column in bold) and the sumNoisyQ*/sumNoisyP* macros.
    """
    src = "results/{squeezing,loss,detector,dephasing,asymmetry,combined}.json"
    sq = by_label(rows("squeezing"), db_of)
    loss = by_label(rows("loss"), lambda r: r["noise.eta"])
    det = by_label(rows("detector"), lambda r: r["noise.eta_d"])
    deph = by_label(rows("dephasing"), lambda r: r["noise.sigma_phi"])
    asym = by_label(rows("asymmetry"), bias_of)
    comb = by_label(rows("combined"), db_of)
    bias = P.TELE_DB   # dB_q - dB_p = -12 (6 dB in q, 18 dB in p: noisy q) and +12 (noisy p)
    cols = [(sq, P.TELE_DB), (loss, 0.9), (det, 0.9), (deph, 0.1), (asym, -bias), (asym, bias), (comb, P.TELE_DB)]
    order = [("bare qubit (break-even)", "bare qubit"), ("GKP square", "square"),
             ("GKP rectangular (λ=1.5)", r"rectangular ($\lambda=1.5$)"), ("GKP hexagonal", "hexagonal")]
    vals = {lab: [float(np.interp(x, *d[lab])) for d, x in cols] for lab, _ in order}
    for lab, short in SHORT.items():
        M.inf(f"sumNoisyQ{short}", vals[lab][4], f"{lab}, bias dB_q - dB_p = -{bias:g} (noisy q), 24 dB total", src, 3)
        M.inf(f"sumNoisyP{short}", vals[lab][5], f"{lab}, bias dB_q - dB_p = +{bias:g} (noisy p), 24 dB total", src, 3)
    best = [min(vals[lab][j] for lab, _ in order[1:]) for j in range(len(cols))]
    lines = [r"\begin{tabular}{|l|ccccccc|}", r"\toprule",
             r"$1-\bar F$ & squeezing & channel loss & detectors & dephasing & noisy $q$ & noisy $p$ & all at once\\",
             r"(ideal decoder) & only & $\eta=0.9$ & $\eta_d=0.9$ & $\sigma_\phi=0.1$ & ($-12$\,dB bias) & ($+12$\,dB bias)"
             r" & (Fig.~\ref{fig:combined}(a))\\", r"\midrule"]
    for lab, text in order:
        cells = []
        for j, v in enumerate(vals[lab]):
            c = f"\\num{{{v:.3f}}}"
            cells.append(f"\\textbf{{{v:.3f}}}" if lab != order[0][0] and v == best[j] else c)   # \\num drops bold
        lines.append(text + "".join(f" & {c}" for c in cells) + r" \\")
        if lab == order[0][0]:
            lines.append(r"\midrule")
    lines += [r"\bottomrule", r"\end{tabular}"]
    (GEN / "summary_table.tex").write_text("\n".join(lines) + "\n")


def numbers_and_tables():
    M = Macros()
    for f in (parameter_numbers, squeezing_numbers, code_quality_numbers, noise_numbers, asymmetry_numbers,
              rect_lambda_numbers, decoder_numbers, decoder_split_table, combined_numbers, validation_numbers,
              stage_numbers, pauli_table, cat_numbers, appendix_numbers, duan_numbers,
              sensitivity_numbers, delta_dephasing_numbers, summary_table):
        f(M)
    M.write(GEN / "numbers.tex", REPO / "docs" / "results_numbers.md")
    print(f"  {len(M.items)} macros -> {GEN / 'numbers.tex'} and docs/results_numbers.md")


JOBS = {"squeezing": fig_squeezing, "noise": fig_noise, "rect_lambda": fig_rect_lambda, "decoders": fig_decoders,
        "combined": fig_combined, "delta_dephasing": fig_delta_dephasing, "wigner": fig_wigner, "stages": fig_stages, "numbers": numbers_and_tables}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("jobs", nargs="*", help=f"any of {', '.join(JOBS)} (default: all)")
    ap.add_argument("--out", default="report", help="folder that gets figures/ and generated/ (LaTeX macros, tables)")
    a = ap.parse_args()
    if unknown := set(a.jobs) - set(JOBS):
        ap.error(f"unknown job(s) {sorted(unknown)}; choose from {list(JOBS)}")
    out = Path(a.out) if Path(a.out).is_absolute() else REPO / a.out
    FIG, GEN = out / "figures", out / "generated"
    (FIG / "png").mkdir(parents=True, exist_ok=True)
    GEN.mkdir(exist_ok=True)
    for name in a.jobs or JOBS:
        JOBS[name]()
