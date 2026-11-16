# Extension experiments: summary (2026-11-10)

**Status (2026-11-17): not used in the submitted report.** The `report_v2/` template and the figures named below
were removed with the old report folders; `scripts_ext/make_ext_figures.py` recreates the figures. Data and file list:
`results_ext/README.md`.

What two helper agents produced while the report was being written, after Claude checked their output. Nothing
here is in the report yet; we decide at the end which parts to use.

Full details: `results_ext/README.md`. Every number: `results_ext/key_numbers.json`. Scripts: `scripts_ext/`.

---

## 1. Audit of the existing results: no problems
- All 53 numbers in the report reproduce exactly from `results/`. Saved results recompute identically with today's code.
- All 102 tests still pass. The library (`cvqec/`), `results/`, `report/` and the user's report files were not changed.
- Two useful details:
  - At 2 dB, GKP beats the bare qubit only barely (x1.05). The advantage grows with squeezing to x26-45 at 16 dB.
    Write "from 2 dB, growing to x26-45 at 16 dB", not "a large margin everywhere".
  - "The cat loses to the bare qubit under teleportation noise" holds even with the best possible correction for
    both: 0.850 vs 0.915 average fidelity at 10 dB. This strengthens the "why GKP" argument (Section 2.2).

## 2. Complementarity: GKP vs cat (figure `results_ext/figures/complementarity.pdf`)

All codes at equal energy (n = 5.07; the bare qubit has n = 0.5).

**(a) Shifts (teleportation noise), 10 dB, infidelity 1 - F:**

| code | decoder | 1 - F |
|---|---|---|
| square GKP | ideal | 0.026 |
| bare qubit | best possible | 0.085 |
| cat | best possible | 0.150 |

**(b) Rotations (pure dephasing), sigma = 0.2 rad:**

| code | decoder | 1 - F |
|---|---|---|
| cat | best possible | 0.003 |
| bare qubit | best possible | 0.007 |
| square GKP | Petz (noise-aware) | 0.025 (its true best is at least 0.013) |
| square GKP | ideal (standard) | 0.095 |
| cat | parity (standard) | 0.078 |

- From 0.04 rad upward the cat beats even GKP's best possible decoder.
- The decoder matters as much as the code. The cat's parity decoder catches photon loss only, so it is about 25x
  worse than the best correction under rotations.

**Correction to an earlier suggestion.** Petz is not a fair "best case" for every code. For the bare qubit under
dephasing it applies the noise a second time and doubles the error. So the agent computed the truly optimal
correction (an SDP with a certificate that nothing does better) for the cat and the bare qubit. GKP is too large for
that, so the figure shows a band where GKP's best decoder must lie (from the Barnum-Knill theorem).

**If used in the report:**
- Cite Grimsmo, Combes & Baragiola, PRX 10, 011058 (2020).
- Say the optimal and Petz decoders are benchmarks that know the noise exactly, not something an experiment runs.
- Say only one cat size was tested.

## 3. Feed-forward gain error (figure `results_ext/figures/gain_mismatch.pdf`)

The library only handled the exact gain g = sqrt(eta). The agent derived the general case and checked it against
the full three-mode simulation (agreement to 7.7e-6).

Square GKP, 10 dB, ideal decoder:

| gain g | 1 - F |
|---|---|
| 1.00 (correct) | 0.026 |
| 0.95 / 1.05 | 0.034 / 0.043 (1.3-1.7x worse) |
| 0.90 / 1.10 | 0.070 / 0.091 (about 3x worse) |

**Message:** the feed-forward electronics must be calibrated to within a few percent.

**Leave out:** the Petz curve in that figure keeps improving as g grows. The agent's explanation is its own
reasoning, not from a paper. Don't use it unless we check it together.

## 4. "Beyond this project" template: `report_v2/beyond_section.tex`

Same format as the section 2 template; not yet `\input` in `main.tex`. Four parts:
1. Complementarity, with the numbers above.
2. A 10-row table of real noise sources: where each enters the circuit and what error it causes.
3. Squeezed cat codes, described as in the paper. They are better against loss, not against displacement noise.
4. A short "what we would simulate next" list.

New references `grimsmo2020` and `schlegel2022` were verified and appended to `references.bib`.

**Check first:** the local-oscillator phase-noise hint in the table uses a formula the agent derived itself
(extra shift noise of about delta^2 e^{2r}/4). Verify it together before using it, or keep that row qualitative.

## 5. Practical notes
- The optimal-correction (SDP) run needs `cvxpy` and `clarabel`, which are not in `qcsim`. The README has the install line.
- To use a figure in the report, copy it from `results_ext/figures/` to `report_v2/figures/`.
