# Runbook: how every number and figure in the report is produced

Last full rerun: 2026-11-14 (all presets). The report LaTeX is no longer in the repository (submitted 2026-11-17);
the figures it used are in `report/figures/`.
All commands run from the repository root. `PY` is the `qcsim` environment's Python:

```bash
PY=$(conda run -n qcsim which python)   # or: conda activate qcsim; PY=python
```

## 0. The two squeezings (read this first)

Every run has **two different squeezing parameters**. They are not related to each other:

| what | symbol | default | in dB | where it enters |
|---|---|---|---|---|
| the GKP code's own width (how well the code state is made) | $\Delta$ | 0.3 | 10.5 dB ($-20\log_{10}\Delta$) | component 1 (encoder); also the Steane ancillas ($\Delta_a=\Delta$ unless stated) |
| the teleportation resource (EPR pair) | $r_1=r_2=r$ | $r=1.38$ | **12 dB** ($20\,r/\ln 10$) | components 2-3; gives shift noise of variance $e^{-2r}=0.063$ |

So "12 dB" always means the teleportation squeezing, and "$\Delta=0.3$" always means the code. Figures say
"teleportation 12 dB" in the title. Macros: `\parTeleDB`, `\parTeleR`, `\parTeleVar`, `\parCodeDelta`, `\parCodeDB`.

**Default operating point** (every sweep keeps all of these except the one it varies):
teleportation 12 dB ($r_1=r_2$), code $\Delta=0.3$, $\eta=\eta_d=1$, $\sigma_\phi=0$, ideal decoder, Fock cutoff $N=140$.
Defined once in `cvqec/experiments/presets.py` (`DELTA`, `TELE_DB`, `LAM`, `N`); `tests/test_experiments.py::test_default_operating_point` checks it.

## 1. Commands

```bash
# tests (131 tests, ~65 s)
$PY -m pytest -q

# stage 1: all sweeps at full resolution -> results/<preset>.csv + .json  (~12 min; overwrites results/)
$PY report/run_sweeps.py
$PY report/run_sweeps.py rect_lambda code_quality        # only some presets

# quick smoke run that does NOT touch results/ (coarse grids, ~40 s)
$PY -m cvqec.experiments all --quick --out /tmp/cvqec_smoke

# stage 2: figures + tables + number macros -> report/figures, report/generated, docs/results_numbers.md (~2 min)
$PY report/make_figures.py
$PY report/make_figures.py squeezing numbers   # only some jobs
#   jobs: squeezing noise rect_lambda decoders combined delta_dephasing wigner stages numbers
#   (--out <folder> writes figures/ and generated/ somewhere else)
```

Before a rerun, keep a copy of the old data: `cp -r results results_backup_<date>` (the 2026-11-14 pre-rerun data is in
`results/archive_2026-11-14_before_rerun/`).

## 2. The experiments (presets)

| preset | varied | fixed | codes / decoder | figure file(s) | macros |
|---|---|---|---|---|---|
| `squeezing` | teleportation $r_1=r_2$: 2-16 dB (1 dB) | defaults | 3 lattices + bare qubit; ideal | `squeezing_tele` | `\inf{Six,Twelve,Sixteen}*`, `\gain*`, `\breakEven*`, `\pauli*`, `pauli_table.tex` |
| `code_quality` | code squeezing 8-14 dB (0.5 dB), $N=200$ | teleportation 12 dB | 3 lattices; ideal (bare qubit = line) | `squeezing_code` | `\codeQ*`, `\maxTail` |
| `loss` | $\eta$: 0.8-1 (0.02) | defaults | 3 lattices + bare; ideal | `noise_loss` | `\lossNinety*`, `\detOverLoss*` |
| `detector` | $\eta_d$: 0.8-1 (0.02) | defaults | 3 lattices + bare; ideal | `noise_detector` | `\detNinety*` |
| `dephasing` | $\sigma_\phi$: 0-0.2 rad (0.02) | defaults | 3 lattices + bare; ideal | `noise_dephasing` | `\dephTenth*`, `\dephasingCross*` |
| `asymmetry` | $\mathrm{dB}_q-\mathrm{dB}_p$: -12..12 (1 dB), sum 24 dB | defaults | 3 lattices + bare; ideal | `noise_asymmetry` | `\asymBest*`, `\asymAtRectBest*`, `\asymPredRect` |
| `rect_lambda` (new) | rectangular $\lambda$: 1.0-2.0 (0.1) at bias 0, -6, -12 dB | defaults | rectangular; ideal | `rect_lambda` | `\rectBestLam*`, `\rectPredLam*`, `\rectLam*Inf*`, `\rectNbarMax` |
| `decoders` | teleportation 8-16 dB | defaults | square, hex; ideal, Steane (joint, $\Delta_a=0.3$), Petz | `decoders_squeezing` | `\dec*`, `\decTwelve*` |
| `steane_ancilla` | ancilla $\Delta_a$: 0.15-0.4 | teleportation 12 dB | square, hex; Steane (joint) | `decoders_ancilla` | `\steaneBest*`, `\steaneThirty*` |
| `decoder_split` | 4 ways of correcting | teleportation 12 dB | square, hex | `decoder_table.tex` | `\split*` |
| `combined` | teleportation 4-16 dB, with $\eta=\eta_d=0.95$, $\sigma_\phi=0.05$ | REALISTIC | 3 lattices + bare (ideal); square, hex Steane (joint) with $\Delta_a=0.3$ and $0.15$ | `combined_squeezing` | `\combTwelve*`, `\combBreakEven*` |
| `sensitivity` (E3, 2026-11-15) | as `combined`, with $\eta=\eta_d\in\{0.90, 0.97\}$ | REALISTIC otherwise | as `combined` | `sensitivity_table.tex` | `\sens*` |
| `delta_dephasing` (E5, 2026-11-16) | code width $\Delta$: 6-16 dB code squeezing, step 0.5, N = 200 | $\sigma_\phi\in\{0, 0.1\}$, 12 dB, rest ideal | 3 lattices + bare qubit | `figures/delta_dephasing.pdf` (job `delta_dephasing`) | `\eFive*` |
| `budget` | imperfections added one at a time | teleportation 12 dB | 3 lattices + bare; ideal, last stage Steane | `combined_budget` | `\budget*` |

Paired panels (Figs. squeezing, noise, decoders, combined) share one fixed axes box and have no legend of their own;
each figure's legend is a separate file `figures/legend_<squeezing|noise|decoders|combined>.pdf`, placed once under
the panels in `results_section.tex`. `rect_lambda` is stand-alone and keeps its own legend.

Not from a preset: `wigner_codes` (Section 1), `stages` (Section 2, `\stage*`, `\parStage*`), the validation
macros `\val*` (full three-mode circuit vs the channel $\mathcal E$, same check as `tests/test_engines.py`).

Every macro also exists as `<name>Fid` = the fidelity $\bar F=1-(1-\bar F)$ where the macro is an infidelity.
The full list with exact values and meanings is generated into `docs/results_numbers.md`.

## 3. Choices and why (each must be stated in the report)

| choice | value | why |
|---|---|---|
| code width | $\Delta=0.3$ (10.5 dB, $\bar n\approx5.07$) | a realistic, near-term GKP quality (similar to recent experiments; cite in iteration 4); same $\Delta$ for every lattice gives (nearly) equal energy: $\bar n$ = 5.07 square/hex, 5.11 rectangular; fits the Fock cutoff $N=140$ with a negligible tail (`\maxTail`). |
| teleportation squeezing | 12 dB | (i) the channel noise $e^{-2r}=0.063$ is of the same size as the code's own peak variance $\Delta^2/2=0.045$, so neither hides the other; (ii) the logical error is ~1 %: error correction clearly works (7-10x better than the bare qubit), yet every added imperfection is still visible; (iii) inside what squeezed-light sources have shown (record ~15 dB; cite in iteration 4). The squeezing sweep shows every other value. |
| rectangular $\lambda$ | 1.5 | a clearly non-square cell (q-tolerance 2.25x the p-tolerance); the `rect_lambda` sweep shows that the best $\lambda$ follows the noise bias and that $\lambda=1$ (square) is best for symmetric noise. |
| rectangular code left out of the decoder / Steane sweeps | | for symmetric noise it is the worst code (close to the bare qubit); its role is the biased-noise test. Its Steane round is the same as the square one in a stretched frame, so it adds no new information. It is included in the all-at-once runs. |
| Steane round | joint, ideal syndrome homodynes | joint decoding keeps the hexagonal cell (Section 1, step 3); the separate version appears only in the decoder table. Syndrome-detector inefficiency, a noisy SUM, loss/dephasing during the round and repeated rounds are NOT in any reported run (listed as not modelled in Section 2.3.5). |
| Steane ancillas | $\Delta_a=\Delta=0.3$ (and 0.15 = 16.5 dB as "better ancillas") | the ancilla is made by the same process as the code state; 0.15 is the best ancilla of the `steane_ancilla` sweep. |
| all-at-once operating point | $\eta=\eta_d=0.95$, $\sigma_\phi=0.05$ rad | illustrative, inside the ranges of the single-imperfection sweeps; NOT taken from one experiment. |
| stage figure | 12 dB, $\eta=0.9$, $\sigma_\phi=0.1$, $\eta_d=1$ | the default squeezing, with the loss and dephasing values used by `\lossNinety*` and `\dephTenth*`, large enough to see by eye. |
| figure of merit | average fidelity $\bar F=(2F_e+1)/3$ over all input states | one number per run that does not depend on which logical state is sent; Pauli rates $p_x,p_y,p_z$ show the kind of error. |
| break-even | bare qubit $\{\lvert0\rangle,\lvert1\rangle\}$ through the same circuit | the standard "is encoding worth it" reference. |

## 4. Accuracy settings (numerical approximations)

- Fock cutoff $N=140$ ($N=200$ for the code-quality sweep): `tests/test_convergence.py` checks that results do not change up to $N=170$.
- Steane round: position grid with $M=16$ points per syndrome cell ($\sqrt\pi$), $G=2M^2=512$ points; ancillas built with $N=200$.
- Shift-noise channel: Lindblad equation solved exactly for unit time (`expm_multiply` on the Liouvillian).
- Validation: full three-mode circuit vs single-mode channel $\mathcal E$, trace distance `\valTraceDist`.

Macros in powers of ten are written as `\ensuremath{a\times10^{b}}` (usable in text and math). The generated tables use
the vertical-rule style of the user's report_v6. Closed-form appendix numbers: `\app*` (incl. `\appDuanMax`).
