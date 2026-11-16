# results/: the data behind the submitted report

Every file here is the full-resolution output of one experiment preset. The presets are defined in
`cvqec/experiments/presets.py`, all files are written by `report/run_sweeps.py` (about 12 min for all of them),
and every figure, table and quoted number is made from them by `report/make_figures.py`. Nothing in the report
was typed by hand.

Figure and table numbers below are those of the submitted report (the last local draft, 2026-11-17).

## Defaults

Every preset keeps these values unless it varies them (`presets.py`: `DELTA`, `TELE_DB`, `LAM`, `N`):
- teleportation squeezing r1 = r2 = 12 dB (r = 1.38, shift-noise variance e^{-2r} = 0.063);
- code width Delta = 0.3 (10.5 dB, mean photon number about 5.07, the same for every lattice);
- rectangular stretch lambda = 1.5;
- eta = eta_d = 1, sigma_phi = 0, ideal decoder, Fock cutoff N = 140.

The two squeezings are different things. "12 dB" always means the teleportation resource; "Delta = 0.3" always
means the code.

## Files

Each preset has a `.csv` (flat table) and a `.json` (the same rows). One row = one simulated point.

| preset (file) | what it varies | codes / decoders | used in the report |
|---|---|---|---|
| `squeezing` | teleportation squeezing 2-16 dB (1 dB steps) | square, rectangular, hexagonal GKP, bare qubit; ideal | Fig. 5(a); Table 1 (Pauli rates at 12 dB); Section 4.1-4.2 numbers; Table 4 |
| `code_quality` | the code's own squeezing 8-14.5 dB (0.5 dB steps), N = 200 | the three lattices; ideal | Fig. 5(b) |
| `loss` | channel transmissivity eta = 0.80-1.00 | the three lattices + bare qubit | Fig. 6(a); Table 4 |
| `detector` | detector efficiency eta_d = 0.80-1.00 | same | Fig. 6(b); Table 4 |
| `dephasing` | dephasing sigma_phi = 0-0.2 rad | same | Fig. 6(c); crossing points in Section 4.3; Table 4 |
| `asymmetry` | split of 24 dB total squeezing between q and p, dB_q - dB_p = -12 ... +12 | same | Fig. 6(d); Table 4 (noisy q / noisy p columns) |
| `delta_dephasing` | code squeezing 6-16 dB at sigma_phi = 0 and 0.1 rad, N = 200 | the three lattices + bare qubit | Fig. 7 (best code width under dephasing) |
| `rect_lambda` | rectangular stretch lambda = 1.0-2.0 at biases 0, -6, -12 dB | rectangular GKP | Fig. 8 |
| `decoders` | teleportation squeezing 8-16 dB | square and hexagonal; ideal, Steane (joint), Petz | Fig. 9(a); Petz fidelity in Section 4.1 |
| `steane_ancilla` | Steane ancilla width Delta_a = 0.15-0.40 at 12 dB | square and hexagonal; Steane (joint) | Fig. 9(b) |
| `decoder_split` | four ways of correcting at 12 dB (ideal or Steane, joint or separate) | square and hexagonal | Table 2 |
| `combined` | teleportation squeezing 4-16 dB with every imperfection on (eta = eta_d = 0.95, sigma_phi = 0.05) | all codes, ideal; square and hexagonal also Steane with Delta_a = 0.3 and 0.15 | Fig. 10(a); middle column of Table 3; last column of Table 4 |
| `budget` | at 12 dB, the imperfections switched on one at a time (label `code :: stage`) | all codes; last stage = Steane round | Fig. 10(b); the error-budget numbers in Section 4.6 |
| `sensitivity` | the `combined` run repeated at eta = eta_d = 0.90 and 0.97 (label ends `@ eta=...`) | as `combined` | Table 3 |
| `cat_reference` | 6, 12, 16 dB | 4-leg cat (same energy), square GKP, bare qubit | one number in Section 2.2 (why not the cat code); not plotted |

Not from this folder: Fig. 1 (Wigner functions of the codes), Fig. 4 (the circuit stage by stage) and the
effective-channel check (trace distance 8.5e-7) are computed directly by `make_figures.py`. Figs. 2 and 3 are TikZ
drawings (`report/figures/*_tikz.tex`).

## Columns

| column | meaning |
|---|---|
| `label` | series name shown in the legend (code, decoder, and for some presets a stage or `@ eta=` / `@ sigma=` tag) |
| `code`, `code.lattice`, `code.delta`, `code.lam` | code family (`gkp`, `fock` = bare qubit, `cat4`), lattice, width Delta, rectangular stretch |
| `decoder`, `decoder.delta_ancilla` | `gkp_ideal`, `gkp_steane`, `petz`, `projection` (bare qubit); Steane ancilla width |
| `N` | Fock cutoff |
| `noise.r1`, `noise.r2` | teleportation squeezing of the q and p correlations (natural units; dB = 20 r / ln 10) |
| `noise.eta`, `noise.eta_d`, `noise.sigma_phi`, `noise.gain` | channel transmissivity, detector efficiency, dephasing (rad), feed-forward gain (empty = sqrt(eta)) |
| `average_fidelity` | F = (2F_e + 1)/3; the plots show 1 - F |
| `entanglement_fidelity` | F_e |
| `p_x`, `p_y`, `p_z` | logical Pauli error rates |
| `mean_photon` | mean photon number of the encoded state |
| `truncation_tail` | population in the top Fock levels (N = 200 sweeps only; largest value 3e-5) |
| `seconds` | run time of the point |

## Archives (superseded, not used in the report)

- `archive_2026-11-14_before_rerun/`: the data before the full rerun of 2026-11-14. The rerun removed the
  syndrome-detector inefficiency from the reported runs (`combined`, `budget`, `decoder_split`), added the rectangular
  code to `code_quality`, and added the `rect_lambda` and `cat_reference` presets.
- `archive_2026-11-14_separate_steane/`: `decoders` and `steane_ancilla` with the old Steane round that corrects q and
  p separately. Kept because it shows why joint decoding was needed for the hexagonal code (now Table 2).

## Regenerate

```bash
conda run -n qcsim python report/run_sweeps.py              # all presets, overwrites this folder
conda run -n qcsim python report/run_sweeps.py loss detector # some presets
conda run -n qcsim python report/make_figures.py             # figures and numbers from this folder
```
Copy this folder first if you want to keep the current data. The meaning of every generated number is listed in
`docs/results_numbers.md`.
