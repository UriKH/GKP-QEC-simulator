# report/figures/: the figures of the submitted report

Every figure used in the report, as vector PDF (one file per panel, as included in the LaTeX) and as a 600 dpi PNG
copy in `png/`. Figure numbers are those of the submitted report (last local draft, 2026-11-17). The PDFs are made
by `report/make_figures.py` from the data in `results/` (see `results/README.md`); the two drawings are TikZ.

| report figure | file(s) | content | data |
|---|---|---|---|
| Fig. 1 | `wigner_codes.pdf` | Wigner functions of the logical 0 of the square, rectangular and hexagonal codes (Delta = 0.3) | computed directly (job `wigner`) |
| Fig. 2 | `hexcell_tikz.tex` (source), `hexcell.pdf` | hexagonal lattice with the joint (hexagon) and separate (parallelogram) correction cells | TikZ drawing |
| Fig. 3 | `circuit_tikz.tex` (source), `circuit.pdf` | the teleportation circuit, components 1-14 | TikZ drawing |
| Fig. 4 | `stages.pdf` | square GKP through the circuit, stage by stage (12 dB, eta = 0.9, sigma_phi = 0.1) | computed directly (job `stages`) |
| Fig. 5(a), (b) | `squeezing_tele.pdf`, `squeezing_code.pdf`, legend `legend_squeezing.pdf` | 1 - F vs teleportation squeezing; vs the code's own squeezing | `squeezing`, `code_quality` |
| Fig. 6(a)-(d) | `noise_loss.pdf`, `noise_detector.pdf`, `noise_dephasing.pdf`, `noise_asymmetry.pdf`, legend `legend_noise.pdf` | channel loss, detector efficiency, dephasing, q/p squeezing split | `loss`, `detector`, `dephasing`, `asymmetry` |
| Fig. 7 | `delta_dephasing.pdf` | 1 - F vs code width at sigma_phi = 0.1 rad (best Delta) | `delta_dephasing` |
| Fig. 8 | `rect_lambda.pdf` | rectangular stretch lambda at three noise biases, with the Gaussian-model optimum | `rect_lambda` |
| Fig. 9(a), (b) | `decoders_squeezing.pdf`, `decoders_ancilla.pdf`, legend `legend_decoders.pdf` | ideal / Steane / Petz vs squeezing; Steane vs ancilla width | `decoders`, `steane_ancilla` |
| Fig. 10(a), (b) | `combined_squeezing.pdf`, `combined_budget.pdf`, legend `legend_combined.pdf` | every imperfection at once; error budget at 12 dB | `combined`, `budget` |

The tables (Pauli rates, decoders, sensitivity, summary) and the numbers quoted in the text are written by
`make_figures.py` into `report/generated/` (not stored in git; recreate with `python report/make_figures.py numbers`).

To compile a TikZ file on its own, wrap it in a standalone document, e.g.
`\documentclass[tikz]{standalone} \usepackage{bm} \usetikzlibrary{calc} \begin{document} \input{hexcell_tikz} \end{document}`
(the circuit also needs the libraries and colours listed in its first line). `circuit.pdf` and `hexcell.pdf` (and
their PNGs) were rendered this way for the repository README.
