"""Figures and key numbers for the extension experiments (results_ext/*.json -> results_ext/figures/*.pdf).

    conda run -n qcsim python scripts_ext/make_ext_figures.py

Styling follows report/make_figures.py: same rcParams and width, and cvqec.viz.sweeps.style_for, so every code
keeps its colour and marker (square blue circles, hexagonal green triangles, cat yellow diamonds, Fock pink
down-triangles); the standard decoder is a solid line and Petz is dotted. (Series labels must not contain
"rect" unless they are the rectangular code: style_for matches substrings, so "no correction" would turn orange.)
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cvqec.experiments import load_rows
from cvqec.viz.sweeps import style_for

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results_ext"
FIG = RES / "figures"
WIDTH = 6.5
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["STIX Two Text", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 9, "axes.titlesize": 9.5, "legend.fontsize": 7.5,
    "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8, "savefig.bbox": "tight",
    "savefig.dpi": 600,   # raster parts (Wigner heat maps); lines and text stay vector in the PDF
})
STD_NAME = {"gkp_ideal": "ideal", "cat_parity": "parity", "projection": "projection", "petz": "Petz"}


def series(rows, xkey):
    out = {}
    for r in rows:
        lab = f"{r['label'].split(' (')[0] if 'Fock' in r['label'] else r['label']}, {STD_NAME[r['decoder']]}"
        out.setdefault(lab, []).append((r[xkey], 1 - r["average_fidelity"]))
    return {k: np.array(sorted(v)).T for k, v in out.items()}


def draw(ax, data, title, xlabel):
    for i, (lab, (x, y)) in enumerate(data.items()):
        ax.plot(x, y, lw=1.6, ms=3.8, label=lab, **style_for(lab, i))
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$1-\bar F$")
    ax.set_title(title, loc="left")
    ax.grid(True, which="major", color="#d9d9d4", lw=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def crossing(x, a, b):
    """First x where curve a crosses curve b (log-linear interpolation); nan if none."""
    d = np.log(a) - np.log(b)
    for i in range(len(d) - 1):
        if d[i] == 0:
            return float(x[i])
        if d[i] * d[i + 1] < 0:
            return float(x[i] + (x[i + 1] - x[i]) * d[i] / (d[i] - d[i + 1]))
    return float("nan")


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    e1, e2 = load_rows(RES / "e1.csv"), load_rows(RES / "e2.csv")
    d1, d2 = series(e1, "channel.sigma_phi"), series(e2, "channel.db")
    opt = load_rows(RES / "optimal.csv") if (RES / "optimal.json").exists() else []

    fig, (a, b) = plt.subplots(1, 2, figsize=(WIDTH, 2.9), sharey=True)
    draw(a, d2, "(a) shifts: teleportation channel", "teleportation squeezing (dB)")
    draw(b, d1, "(b) rotations: pure dephasing", r"dephasing $\sigma_\phi$ (rad)")
    b.set_ylabel("")
    for ax, d, kind, xk in ((a, d2, "E2", "db"), (b, d1, "E1", "sigma_phi")):
        # optimal recovery (SDP) for cat and Fock
        for lab in dict.fromkeys(o["label"] for o in opt):
            pts = sorted((o[xk], 1 - o["F_avg_opt"]) for o in opt if o["experiment"] == kind and o["label"] == lab)
            x, y = np.array(pts).T
            ax.plot(x, y, marker="*", ls="none", ms=5.5, mfc="none", mew=0.7, color=style_for(lab, 0)["color"])
        # GKP: Barnum-Knill bracket for the optimum,  F_e(opt) <= sqrt(F_e(Petz))
        for lab, (x, y) in d.items():
            if lab.startswith("GKP") and lab.endswith("Petz"):
                fe = 1 - 1.5 * y
                lower = (2 / 3) * (1 - np.sqrt(fe))
                ax.fill_between(x, lower, y, color=style_for(lab, 0)["color"], alpha=0.13, lw=0)
    a.plot([], [], marker="*", ls="none", ms=5.5, mfc="none", mew=0.7, color="#18212c", label="optimal recovery (SDP)")
    a.fill_between([], [], [], color="#18212c", alpha=0.13, lw=0, label="GKP optimum lies in band")
    a.set_ylim(1e-5, 1)
    h, l = a.get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout()
    fig.savefig(FIG / "complementarity.pdf")
    plt.close(fig)

    nums = {}
    fock1 = d1["Fock qubit, projection"]
    for lab, (x, y) in d1.items():
        nums[f"E1 sigma where '{lab}' crosses Fock (projection)"] = crossing(x, y, fock1[1])
        for s in (0.1, 0.2, 0.3, 0.4):
            nums[f"E1 1-F at sigma={s} '{lab}'"] = float(np.interp(s, x, y))
    fock2 = d2["Fock qubit, projection"]
    for lab, (x, y) in d2.items():
        nums[f"E2 dB where '{lab}' crosses Fock (projection)"] = crossing(x, y, fock2[1])
        for db in (6, 10, 16):
            nums[f"E2 1-F at {db} dB '{lab}'"] = float(np.interp(db, x, y))

    for kind, xk in (("E1", "sigma_phi"), ("E2", "db")):
        for lab in dict.fromkeys(o["label"] for o in opt):
            pts = sorted((o[xk], 1 - o["F_avg_opt"], o["F_avg_opt"] - o["F_avg_opt_upper"])
                         for o in opt if o["experiment"] == kind and o["label"] == lab)
            x, y, gap = np.array(pts).T
            nums[f"{kind} max certified gap (F_avg) '{lab}, optimal'"] = float(np.abs(gap).max())
            d = d1 if kind == "E1" else d2
            if "cat" in lab:
                fock_opt = np.array(sorted((o[xk], 1 - o["F_avg_opt"]) for o in opt
                                           if o["experiment"] == kind and "Fock" in o["label"])).T
                fy = np.interp(x, fock_opt[0], fock_opt[1])     # common grid (E2: cat only at 8-16 dB)
                nums[f"{kind} x where cat optimal crosses Fock optimal"] = crossing(x, y, fy)
                nums[f"{kind} min ratio cat_opt/Fock_opt"] = float(np.min(y / fy))
            for xv in ((0.1, 0.2, 0.3, 0.4) if kind == "E1" else (6, 10, 16)):
                if x.min() <= xv <= x.max():             # np.interp would silently clamp outside the grid
                    nums[f"{kind} 1-F at {xk}={xv} '{lab}, optimal'"] = float(np.interp(xv, x, y))
    for d, kind, xs in ((d1, "E1", (0.1, 0.2, 0.3, 0.4)), (d2, "E2", (6, 10, 16))):
        for lab, (x, y) in d.items():
            if lab.startswith("GKP") and lab.endswith("Petz"):
                low = (2 / 3) * (1 - np.sqrt(1 - 1.5 * y))
                for xv in xs:
                    nums[f"{kind} GKP optimum lower bound (Barnum-Knill) at {xv} '{lab}'"] = float(np.interp(xv, x, low))

    e3p = RES / "e3.csv"
    if e3p.with_suffix(".json").exists():
        e3 = load_rows(e3p)
        d3 = series(e3, "channel.gain")
        fig, ax = plt.subplots(figsize=(WIDTH * 0.6, 2.6))
        draw(ax, d3, r"feed-forward gain error (10 dB, $\eta=1$)", r"feed-forward gain $g$")
        ax.axvline(1.0, color="#18212c", lw=0.6, ls=(0, (2, 2)))
        ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.0, 1.0))
        fig.tight_layout()
        fig.savefig(FIG / "gain_mismatch.pdf")
        plt.close(fig)
        for lab, (x, y) in d3.items():
            i = int(np.argmin(y))
            nums[f"E3 '{lab}' best g"] = float(x[i])
            nums[f"E3 '{lab}' 1-F at best g"] = float(y[i])
            for g in (0.9, 0.95, 1.0, 1.05, 1.1):
                nums[f"E3 '{lab}' 1-F at g={g}"] = float(np.interp(g, x, y))
    (RES / "key_numbers.json").write_text(json.dumps(nums, indent=1))
    for k, v in nums.items():
        print(f"{k}: {v:.4g}")


if __name__ == "__main__":
    main()
