# Session handoff: start here

The one file a new session should read first. It says what the project is, what has been done, where
everything lives, and what comes next. Details are in [IMPLEMENTATION_LOG.md](IMPLEMENTATION_LOG.md)
(per-phase engineering log) and [TEST_PLAN.md](TEST_PLAN.md).

Last updated: 2026-11-17 (status below); sections 2-8 were last revised 2026-11-14.

**Final state 2026-11-17:** the report was submitted. The repository was cleaned: the report LaTeX versions
(`report_v2/`-`report_v15/`, the old `report/main.tex`) were removed; the final text is in the user's Overleaf project
and a full backup of the pre-cleanup tree is `../GKP-QEC-simulator_backup_2026-11-17.tar.gz`. What remains: the
library, tests, sweep data (`results/`, `results_ext/`), the scripts that make them (`report/run_sweeps.py`,
`report/make_figures.py`, `scripts_ext/`), the figures the report used (`report/figures/`), and these docs.
Paths to `report_vN/` further down are historical. Lessons for future projects: LESSONS_LEARNED.md. Also see `CLAUDE.md` (repo root) and the skills in `.claude/skills/`.

---

## 1. The project in one paragraph
This is the final project for *Modern Quantum Computing* (Technion 00460054), by Uri Kasher Hitin. We built a
continuous-variable (CV) quantum circuit simulator in QuTiP. It encodes a qubit in a bosonic mode with
**GKP codes** (square, rectangular, hexagonal lattices), sends it through **CV teleportation** with realistic
imperfections, corrects errors, and measures how well each lattice protects the qubit. The simulator is done
and tested. The current phase is **writing the report**: the user writes it themselves, and Claude explains
and reviews, so the user can defend every part in person.

## 2. How to work with this user
- **Mentor first.** Explain the physics and the reasoning before (or instead of) writing things for them. The
  user wants to understand everything well enough to explain it to the instructor in person. For the
  report, the user writes the text; Claude explains, suggests structure, and reviews drafts.
- **Course-based.** Ground explanations in the course book `context/course_materials/summary.pdf`
  (Hebrew, 128 pp.) with page numbers, and flag anything that goes beyond the course ("bridge material").
- **Background.** Besides this course, the user has taken a quantum-technologies course (Lindblad master
  equation, beam splitters and photon bunching/HOM) and a quantum-information course (channels, Kraus
  operators). Connect to those where useful.
- **Keep the logs current.** After each phase, update `docs/IMPLEMENTATION_LOG.md`, the HTML build log
  (artifact below, source copy `docs/build_log.html`), and this file.
- Nothing is committed to git yet. Commit only when the user asks.

## 3. The instructor's required report structure (given 2026-11-07)
1. Explanation of the codes: how they work, and what errors they can fix (using the lecture notes).
2. Why this circuit and this code? What does the circuit do?
3. What did we run (which experiments), and why?
4. Results (graphs, tables).
5. Conclusions from the results.

**Status (2026-11-09):** the new report lives in `report_v2/` (main.tex `\input`s one file per item).
Item 1 is done: `report_v2/code.tex`, archived as `report_v2/archive/code_2026-11-09_section1.tex`. Item 2 has a
template with hints and placeholders: `report_v2/circuit_section.tex`; the user is writing it. Edits are proposed as
new file versions with `\chg{}` colour and `% [C#]` tags (see memory "report-edit-workflow").
**Earlier status:** not started in this structure.
**Status (2026-11-14):** the user moved to `report_v3/` (main.tex now defines `\chg`, `\todo`, `\vt` and the math
macros). `circuit_section.tex` (item 2) is mostly written by the user; open `\todo{}`s remain (entanglement for mixed
states, Steane round, implementation note, validation and closing sentences). Claude fixed two compile errors only
(duplicate `\chg` in code.tex commented out; `\ad'` -> `\op{a}^{\dagger\prime}`). Notes for the open points and a
review list R1-R11: `docs/section2_notes.html`. Course check: the Steane CODE is in the course (p.74, 90), Steane-type
GKP correction and the SUM gate are not (bridge, gkp2001). Decoder split computed 2026-11-14 (not yet in results/):
1-F at 12 dB, square ideal/separate/Steane 0.0093/0.0093/0.0487, hexagonal 0.0070/0.0153/0.0585.
Agreed queue: test audit (DONE 2026-11-14: 122 tests, mutation check), stage-by-stage Wigner figure (DONE:
`report_v3/figures/stages.pdf`, tests/test_stages.py; not yet \input in the text), implementation subsection (user,
Section 3), choice of main graphs, final numerical check. Validation macros `\valTraceDist` etc. are in numbers.tex.
**Status (2026-11-14, later):** `report_v4/` = v3 + hint comments in Section 2 (`% NOTES s.N` -> docs/section2_notes.html,
`% R#` review items, `% CITE:` where to cite) + templates `runs_section.tex` (item 3, incl. Implementation),
`results_section.tex` (item 4, figures included, macro lists), `conclusions_section.tex` (item 5); compiles, 16 pp.
Figures regenerated with 600 dpi heat maps (vector otherwise); results_ext figures copied in. Open decision: Steane
extensions (joint Q/P decoding = separate decoding + syndrome-dependent logical Pauli, 16x16 branches on the 512-point
grid, est. half a day incl. tests; syndrome-detector inefficiency = Gaussian blur in the round's weights, ~1-2 h).
New stage agreed: citations + readability pass, section by section, with the user.
**Status (2026-11-14, end):** joint Steane decoding + syndrome-detector inefficiency implemented and tested (128 tests);
presets decoders, steane_ancilla rerun, new decoder_split (table `generated/decoder_table.tex`, macros `\split*`);
report_v4 updated (figure, numbers, table in results_section, hints in code.tex R-S1 + course pages, circuit_section (4)).
All-at-once runs added (presets combined, budget; figure combined.pdf; report_v4 Section 4.6 template):
realistic Steane (Delta_a = 0.3) never beats the bare qubit with everything on; Delta_a = 0.15 does.
Course check: ancilla stabilizer measurement IS in the course (p.70-72, 84, 86, 93-96); the physical GKP measurement
(Steane-type, SUM gate) is NOT. Next: the user writes; then the citations + readability pass, section by section.

**Status (2026-11-14, evening): full rerun for report_v5.** The user moved the report to `report_v5/` (v4 = previous).
- Presets: `DELTA` (code, 0.3 = 10.5 dB) and `TELE_DB` (teleportation, 12 dB) are now named constants: two different
  squeezings, never to be mixed. New presets `rect_lambda` (why lambda = 1.5: best lambda follows the bias; Gaussian
  model predicts it to 0.02) and `cat_reference` (backs the cat sentence of Section 2.2: cat worse than bare qubit).
  `code_quality` now includes the rectangular code. Syndrome-detector inefficiency removed from every reported run
  (the report says it is not modelled); `decoder_split` has 4 rows; budget has 5 stages.
- CHANGED RESULT: with all imperfections on, the realistic Steane round (Delta_a = 0.3) now beats the bare qubit above
  ~4 dB, narrowly (0.110 vs 0.121 at 12 dB); Delta_a = 0.15 gives 0.073. (Before, eta_s = 0.95 pushed it above.)
- `make_figures.py --out report_v5`: one PDF per panel (+ 600 dpi PNG), one legend scheme (colour = code, line style =
  decoder), 296 macros with meanings in `docs/results_numbers.md`. Stage figure moved to 12 dB (was 10 dB).
- Commands and the reasons for every choice: `docs/RUNBOOK.md`. Plan to finish: `docs/FINALIZATION_PLAN.md`.
- 130 tests (two new run-policy tests in test_experiments.py). Pre-rerun data: `results/archive_2026-11-14_before_rerun/`;
  the user's pre-rerun report files: `report_v5/archive/2026-11-14_before_rerun/`.

**Extensions (2026-11-10, not yet in the report):** `results_ext/` + `scripts_ext/` (README there). Audit: all report
numbers reproduce from `results/`. E1 pure dephasing and E2 teleportation noise for square/hex GKP, 4-leg cat and
Fock qubit (incl. SDP-optimal recovery for cat/Fock; GKP optimum bounded via Barnum-Knill) -> figure
`results_ext/figures/complementarity.pdf`: GKP wins under shifts, cat wins under rotations. E3 feed-forward gain
mismatch (general-gain channel in `scripts_ext/gain_channel.py`, checked vs the 3-mode engine) ->
`gain_mismatch.pdf`. Caveat: Petz is not a fair best case for the Fock qubit (it applies dephasing twice). SDP
needs cvxpy + clarabel (not in `qcsim`; reinstall line in the README). Template `report_v2/beyond_section.tex`
(not yet `\input`); bib entries grimsmo2020, schlegel2022 added (verified). The existing draft `report/main.tex` (9 pages) is organized
differently (intro, circuit, codes/decoders, results, summary, appendix). Its content covers all five items
and can be re-cut, but the plan is to go through the items **one at a time**, with the user writing.

## 4. Locked decisions
| Topic | Decision |
|---|---|
| Circuit | CV teleportation, 14 components: encoder, S(r₁), S(−r₂), BS(π/4) (EPR), loss η, Bell BS(π/4), R(π/2), 2 homodynes (efficiency η_d), feed-forward X(−√2 g m_q) and Z(√2 g m_p), dephasing, Steane EC (SUM + homodyne + shift) |
| Code family | **GKP only** in the report: square, rectangular (λ = 1.5), hexagonal. The 4-leg cat was dropped from the report (it stays in the library). Bare Fock qubit {\|0⟩,\|1⟩} is the break-even reference |
| Why GKP over cat | Teleportation noise is random displacements, which is GKP's native error model. GKP Paulis and Cliffords are Gaussian (same toolbox as the circuit). Our test: cat F̄ ≈ 0.77 vs bare qubit 0.90 at 10 dB. **Not** "dephasing": dephasing is *GKP's* weak spot too |
| Errors | Finite squeezing (intrinsic), channel loss η, detector inefficiency η_d, dephasing σ_φ, asymmetric squeezing r₁ ≠ r₂ |
| Code energy | Δ = 0.3 (10.5 dB), n̄ ≈ 5.07 for every lattice (n̄ ≈ 1/(2Δ²) − ½), so the comparison is at equal energy. A Δ sweep is included |
| Decoders | Ideal (nearest lattice point, Voronoi cell), Steane with finite-Δ ancillas (realistic), Petz (benchmark that knows the noise) |
| Metric | Average fidelity F̄ = (2F_e + 1)/3 of the logical qubit channel; plotted as 1 − F̄ = ⅔(p_x + p_y + p_z) |

## 5. Conventions (course book p.28–55, 78–82; single source: `cvqec/conventions.py`)
- ħ = 1, q = (a + a†)/√2, p = −i(a − a†)/√2, vacuum variance 1/2.
- S(r) = e^{r/2(a² − a†²)} squeezes q.
- B(θ) = e^{−θ(a₁a₂† − a₁†a₂)}.
- R(φ) = e^{−iφn}.
- X(s) = e^{−isp} = D(s/√2), Z(s) = e^{isq} = D(is/√2).
- GKP: X̄ = D(α), Z̄ = D(β), stabilizers D(2α), D(2β), with βα* − β*α = iπ.

## 6. Key physics results (what the report should say)
- **Teleportation = one simple channel.** With gain g = √η:
  E = Dephasing(σ_φ) ∘ Loss(η) ∘ Displacement(σ_q², σ_p²), where σ² = e^{−2r} + (1 − η_d)/η_d.
  The full three-mode simulation confirms it to a trace distance of 8.5×10⁻⁷. Coherent-state teleportation
  fidelity is 1/(1 + e^{−2r}), confirmed.
- **Squeezing.** Square and hexagonal GKP beat the bare qubit at every squeezing from 2 to 16 dB: by ×7–9 at
  12 dB and ×26–45 at 16 dB. (An early note claiming "~7 dB break-even" was wrong and has been corrected.)
- **Gaussian model for the square code**, accurate to within 0.4% over 2–16 dB:
  P = erfc(√π / (2√(2σ²))), with σ² = e^{−2r} + Δ²/2. Then p_x = p_z = P(1 − P), p_y = P², and
  1 − F̄ = ⅔(2P − P²).
- **Why hexagonal wins.** Every cell has area π, but the hexagon's inscribed radius is 0.952 against the
  square's 0.886 (densest packing). It also makes p_x = p_y = p_z. In the square code p_y ≪ p_x = p_z,
  because a Ȳ shift is √2 longer.
- **Rectangular (λ = 1.5)** has a long q-side and a short p-side, so p_x is tiny and p_z dominates. It loses
  under symmetric noise and wins under matched bias: at dB_q − dB_p = −10 it gives 0.018, against 0.049
  for square and 0.038 for hexagonal.
- **Loss vs detector.** At 12 dB and 90% efficiency, square GKP gives 0.038 (channel loss) vs 0.076
  (detector). The detector adds (1 − η_d)/η_d = 0.11, more than e^{−2r} = 0.063 itself.
- **Dephasing is GKP's weak spot.** A rotation by φ shifts a peak at radius R by about Rφ, and R ≈ 3.3. GKP
  loses to the bare qubit beyond σ_φ ≈ 0.14 rad (rectangular: 0.06). This is a trade-off: more energy means
  better shift protection but worse rotation sensitivity.
- **Realistic EC.** Steane with 10.5 dB ancillas levels off at a few percent (0.026 at 16 dB vs 0.0011
  ideal). Better ancillas help: 0.049 at Δ_a = 0.3 → 0.017 at Δ_a = 0.15. Hexagonal loses its advantage under
  Steane, because sequential q/p syndromes decode with a parallelogram cell.
- **Petz beats ideal** at 16 dB (2.5e-4 vs 1.1e-3), because the ideal decoder counts the code's own width Δ
  as noise. Petz is a benchmark, not a realistic decoder.

## 7. Where everything is
```
cvqec/                    the simulator (Python package, `pip install -e .`)
  conventions.py          course conventions
  ops/gates.py            gates; displacement_exact (Laguerre); gaussian_unitary
  ops/measurement.py      homodyne <x_theta|n>
  codes/                  BosonicCode ABC + registry; gkp.py (Lattice), cat.py (cat + Fock reference)
  noise/                  Channel ABC; loss, dephasing, Gaussian displacement (Kraus + Lindblad)
  circuits/               Circuit ops; teleportation.py (circuit + effective_channel + TeleportationParams)
  engines/                fock.py (exact 3-mode), effective.py (run_logical, single mode)
  decoders/               gkp_ideal, gkp_steane, cat_parity, projection, petz
  metrics/                LogicalChannel (F_e, F_avg, PTM, Pauli rates) + @metric registry
  experiments/            Experiment config, sweeps, presets, CLI (python -m cvqec.experiments)
  viz/                    wigner.py, sweeps.py
tests/                    130 tests, ~65 s (conda run -n qcsim pytest)
results/*.csv|json        full-resolution sweep data used by the report
report/                   main.tex, circuit.tex (TikZ), references.bib, figures/*.pdf,
                          generated/{numbers,pauli_table}.tex, run_sweeps.py, make_figures.py,
                          overleaf.zip, main.pdf, README.md
report_v2/                previous version (kept)
report_v3/                previous version (user's Section 2 draft, 2026-11-13)
report_v4/                CURRENT report: main.tex inputs code.tex (item 1), circuit_section.tex (item 2), ...;
                          archive/ holds accepted versions
context/course_materials/ course book summary.pdf
docs/                     HANDOFF.md (this), IMPLEMENTATION_LOG.md, TEST_PLAN.md, build_log.html
```
- **Environment:** conda env `qcsim` (`environment.yml`: Python 3.12, qutip 5.3.1, mpmath, tectonic, ...).
- **Rebuild the report:**
  1. `python report/run_sweeps.py` (~10 min)
  2. `python report/make_figures.py`
  3. `cd report && tectonic -X compile main.tex` (old draft). The current report is `report_v2/`: copy
     `report/figures/*.pdf` and `report/generated/*.tex` into `report_v2/` after step 2, then compile there.
  - Smoke runs: `python -m cvqec.experiments all --quick --out <scratch>`. Without `--out` it overwrites `results/`.

  The numbers in the text are macros from `generated/numbers.tex`.
- **Artifacts (private, claude.ai):**
  - Physics primer (circuit diagram, derivations, exercises): https://claude.ai/artifact/HTonK6KUPg4U1k3Mk3T3Dq
  - Build log: https://claude.ai/artifact/NuidExyFL7KAhb8kZGebqH

  To update one from a new session, read it first, then publish with its `url`.

## 8. Lessons and gotchas (do not relearn these)
- **Truncation.** A small Fock tail does not prove accuracy: exponentiating a truncated generator makes the
  gate wrong. Test convergence in N. Use `qt.coherent(..., method="analytic")` and `displacement_exact` for
  large shifts.
- **Speed.** Three-operand `np.einsum` was 20× slower than batched `@` in the Fock engine.
- **LaTeX.** `\input` inside a booktabs `tabular` fails, so the generator writes the whole table. Keep the
  sources ASCII-only for Overleaf's pdfLaTeX.
- **References.** `report_v2/references.bib` was verified on the web on 2026-11-09 (old unverified copy in
  `report_v2/archive/`). Verify every new entry before adding it. `report/references.bib` is the old, unverified one.

## 9. Next steps
1. Go through the instructor's five items one at a time. For each: Claude explains (with course page
   references) and suggests what the section should contain; the user writes; Claude reviews.
2. Re-cut `report/main.tex` into the five-part structure as the user's sections arrive.
3. Optional follow-ups noted earlier:
   - a rectangular λ matched to the noise bias;
   - a Steane decoder with joint q/p (Voronoi) decoding for the hexagonal code;
   - a small end-to-end GKP run through the three-mode engine.
