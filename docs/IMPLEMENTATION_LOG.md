# Implementation log: CV teleportation + bosonic QEC simulator

This is the running record of what was built, why, and how it was verified, one phase per section.
It is the markdown twin of the HTML tracker page. Update both when a phase changes.

- Plan: `~/.claude/plans/hi-claude-i-want-structured-grove.md` (milestones 1-9)
- Background primer (physics refresher): https://claude.ai/artifact/HTonK6KUPg4U1k3Mk3T3Dq
- Test plan: [TEST_PLAN.md](TEST_PLAN.md)
- Environment: conda env `qcsim` (`environment.yml`); `pip install -e .` once; run tests with `conda run -n qcsim pytest`

Status (2026-11-03): milestones 2-8 implemented; 102 tests pass in ~31 s.
Milestone 1 (the user's own teleportation derivation) and milestone 9 (final figures and report) are with the user.

---

## Conventions used everywhere (course book p.28-55, 78-82)

- ħ = 1, q = (a + a†)/√2, p = −i(a − a†)/√2, vacuum variance 1/2.
- R(φ) = e^{−iφn}: a → e^{−iφ}a.
- B(θ) = e^{−θ(a₁a₂† − a₁†a₂)}: a₁ → cos θ a₁ + sin θ a₂.
- S(r) = e^{r/2(a² − a†²)} squeezes q.
- X(s) = e^{−isp} = D(s/√2), Z(s) = e^{isq} = D(is/√2).
- GKP: X̄ = D(α), Z̄ = D(β), with βα* − β*α = iπ.
  - Square: α = √(π/2).
  - Rectangular: α = λ√(π/2), β = (i/λ)√(π/2).
  - Hexagonal: α = √(π/√3), β = α e^{2πi/3}.
- The only file that encodes these conventions is `cvqec/conventions.py`.

---

## Phase 1: environment, course reading, primer (pre-plan)
- Created the `qcsim` env and read the course book.
- Published the primer: circuit diagram, derivations, exercises, Wigner gallery.
- Locked decisions:
  - Teleportation circuit.
  - Codes: GKP square, rectangular and hexagonal; 4-leg cat; Fock break-even reference.
  - Errors: finite squeezing, loss, dephasing, η_d, r₁≠r₂.
  - Decoders, staged: ideal → Petz → Steane.

## Milestone 2: conventions and Gaussian gates
**Files:** `cvqec/conventions.py`, `cvqec/ops/gates.py`, `cvqec/truncation.py`, `tests/test_conventions.py`, `tests/test_gates.py`

**What:** quadratures, dB conversions, and the gates R, S, D, X, Z, B and SUM. Each gate's docstring states its Heisenberg action and course page. `on_modes` embeds a gate in a register.

**Verified:**
- unitarity;
- squeezer sign and covariance (p.39);
- beam-splitter action (p.52);
- Hong–Ou–Mandel dip;
- SUM Heisenberg action;
- EPR correlations Var(q_a − q_b) = e^{−2r₁} and Var(p_a + p_b) = e^{−2r₂}.

**Lesson (truncation):**
- At N = 30 the EPR tail weight was 3e-12, yet a variance was off by 1.4e-5. A truncated generator makes the *gate* wrong, not just the state.
- Convergence in N is the real check. Errors were 1.6e-3, 1.4e-5 and 1.1e-7 at N = 20, 30 and 40.
- This lesson recurs in milestones 4–5: exact displacement matrix elements, and analytic coherent states.

## Milestone 3: bosonic codes and Wigner plots
**Files:** `cvqec/codes/{base,gkp,cat}.py`, `cvqec/viz/wigner.py`, `tests/test_codes.py`

**Design:**
- `BosonicCode` only needs `codewords(N)`. The encoder (Löwdin-orthonormalized, so an exact isometry), projector, logical Paulis and n̄ are derived from it.
- Registry: `make_code("gkp", lattice="hexagonal", delta=0.3)`.

**GKP construction:**
- |μ_L⟩ ∝ e^{−Δ²n} Σ_{m,n} D((2m+μ)α) D(nβ)|0⟩.
- The group sum projects onto the code, and e^{−Δ²n} acts on coherent states analytically. So there's no truncated operator, and one class covers every lattice.

**Cat:** the course basis |0_L⟩ ∝ |α⟩+|−α⟩, |1_L⟩ ∝ |iα⟩+|−iα⟩. This equals (|C0⟩ ± |C2⟩)/√2.

**Verified:**
- overlap with the course's ψ(q) form is 0.9997;
- ⟨D(γ)⟩ ≈ e^{−|γ|²Δ²/2} (square: e^{−πΔ²});
- D(α) and D(β) act logically;
- the rectangular trade-off;
- the cat's mod-4 support;
- a single loss flips parity and leaves the code space.

**Numbers:** n̄ = 5.07 for all GKP lattices at Δ = 0.3. The cat needs α ≈ 2.25 for equal energy.

## Milestone 4: noise channels and homodyne measurement
**Files:** `cvqec/noise/{base,channels}.py`, `cvqec/ops/measurement.py`, additions to `cvqec/ops/gates.py`, `tests/test_noise.py`, `tests/test_measurement.py`

**Channel interface:**
- `apply` (linear, so it also acts on |0_L⟩⟨1_L|), `adjoint` (for Petz), `kraus` (for the Fock engine), `lindblad` (for cross-checks), `apply_many` (batched).
- `Compose`, plus a registry.

**Channels:**

| Channel | Implementation |
|---|---|
| `PhotonLoss(η)` | Exact binomial Kraus operators |
| `Dephasing(σ)` | Exact element-wise mask ρ_mn e^{−(m−n)²σ²/2}; Gauss–Hermite Kraus when needed |
| `GaussianDisplacement(var_q, var_p)` | Exact Lindblad evolution var_q·D[p] + var_p·D[q] for unit time (`expm_multiply` on a sparse Liouvillian); asymmetric noise is free |

**Measurement:** ⟨x_θ|n⟩ = e^{−iθn} ψ_n(x) is the single formula for all homodyne work (course p.54).

**New gate helpers:**
- `displacement_exact(α, N)`: exact ⟨m|D(α)|n⟩ via a Laguerre recurrence. It's valid at |α|² ≈ 500, where expm is hopeless; checked against mpmath to 1e-13 relative.
- `gaussian_unitary(M, N)`: any symplectic 2×2 matrix as R·S·R (Euler/SVD).

**Verified:**
- CPTP;
- each channel equals `qutip.mesolve` of its Lindblad form (< 1e-6);
- the adjoint identity;
- loss maps |α⟩ → |√η α⟩, and losses compose;
- displacement noise adds exactly var_q and var_p;
- the coherent-state homodyne distribution;
- ⟨p|n⟩ = (−i)ⁿψ_n(p).

## Milestone 5: circuits and the full Fock engine
**Files:** `cvqec/circuits/{base,teleportation}.py`, `cvqec/engines/fock.py`, `tests/test_engines.py`

**Circuit:** a list of `Gate | ChannelOp | Measure | FeedForward` on named modes. It describes *what* happens; engines decide how to simulate it.

`teleportation_circuit(params)` is components 2–12:
- S(r₁) and S(−r₂), then B(π/4): the EPR pair;
- loss η on Bob's mode;
- Bell B(π/4) and R(π/2);
- two homodynes with efficiency η_d;
- feed-forward X(−√2 g m_q) and Z(√2 g m_p);
- dephasing.

The encoder (component 1) is the input state; the Steane round (13–14) is a decoder.

**Fock engine:**
1. Channels are moved past operations on other modes, which is exact.
2. Up to the last measurement the state is an ensemble of kets; Kraus operators multiply the branches. A homodyne contracts a mode against ⟨x_i|n⟩·√dx on an outcome grid, and η_d is modeled as loss before an ideal detector, with outcomes rescaled by 1/√η_d.
3. After the last measurement there is one density matrix per outcome record. Feed-forward uses exact displacements, then the outcomes are summed.

**Verified:**
- Coherent-state teleportation F = 1/(1+e^{−2r}) to 1e-6 (N = 25, r = 0.3 and 0.5).
- **The central check:** with η = 0.9, η_d = 0.95, r₁ = 0.6, r₂ = 0.4 and σφ = 0.15, the full three-mode simulation equals the derived single-mode channel to a trace distance of 8.5e-7.

**Performance:** replacing three-operand `einsum` with batched matmul took the run from 63 s to 3.7 s.

## Milestone 6: effective engine
**Files:** `cvqec/circuits/teleportation.py::effective_channel`, `cvqec/engines/effective.py`

- E = Dephasing(σφ) ∘ Loss(η) ∘ Displacement(e^{−2r₁} + (1−η_d)/η_d, e^{−2r₂} + (1−η_d)/η_d), valid for gain g = √η. It raises an error otherwise; use the Fock engine for other gains.
- `run_logical(code, channel, decoder, N)` pushes V|i⟩⟨j|V† through E and the decoder, returning the 4 outputs as a `LogicalChannel` plus the noisy |+_L⟩.
- GKP at N = 140 takes about 1 s per point.

## Milestone 7: decoders
**Files:** `cvqec/decoders/{base,gkp,subspace,petz}.py`, `tests/test_decoders.py`

| Key | What it models | Notes |
|---|---|---|
| `gkp_ideal` | Ideal EC with minimum-distance (Voronoi) cells, then a perfect readout | Stabilizer subsystem decomposition: ⟨P_L⟩ = Σ_{ℓ∈coset(P)} w_V(ℓ) s_P(ℓ) Tr[ρD(ℓ)], with w_V the exact Fourier transform of the polygonal cell. Precomputed as three N×N operators (≈1 s). `cell="parallelogram"` is available for comparison. |
| `gkp_steane` | One Steane round with finite-Δ_a GKP ancillas (components 13–14), then a perfect readout | Exact on a position grid: SUM is an index shift, homodyne outcomes are summed, the correction is a roll. Non-square lattices are first mapped to the square frame by `gaussian_unitary` (so they get parallelogram cells). The p-round uses an exactly unitary DFT (G = 2M², grid aligned with √π cells). |
| `cat_parity` | Knill–Laflamme recovery for {1, a}: measure the parity sector, rotate back | Generic `ErrorSubspaceDecoder`; anything uncorrectable goes to I/2. |
| `projection` | No correction (keep the code-space part) | The baseline. |
| `petz` | Petz/transpose channel, which knows the noise | A benchmark only; within Barnum–Knill: F_e(Petz) ≥ F_e(opt)². |

**Verified:**
- The ideal decoder equals direct integration of the q-distribution (square, 1e-5).
- The grid decoder agrees with the characteristic-function decoder, converging as O(h²) (2e-4 at M = 16, below 1e-4 at M = 32).
- Hex frame change plus grid equals parallelogram cells.
- Codewords decode with F > 0.99 for every lattice.
- The decoded state is physical.
- Voronoi beats parallelogram for hex and makes p_x = p_y = p_z; for square the two cells coincide.
- Steane improves monotonically as Δ_a shrinks and approaches ideal.
- The cat undoes a single loss with F > 0.9999.
- The Petz Barnum–Knill bound holds.

## Milestone 8: metrics, experiments, sweep plots
**Files:** `cvqec/metrics/{logical,registry}.py`, `cvqec/experiments/{config,runner,presets,__main__}.py`, `cvqec/viz/sweeps.py`, `tests/test_metrics.py`, `tests/test_experiments.py`, `tests/test_convergence.py`

**Metrics** (`@metric` registry):
- average and entanglement fidelity, logical infidelity;
- Pauli-twirled error rates (from the PTM);
- n̄ of the noisy |+_L⟩;
- Wigner negativity (slow, opt-in);
- truncation tail.

**Experiments:**
- `Experiment` dataclass and `grid(base, **{"noise.eta": [...]})`.
- `run_all` caches codes and decoders.
- `save_rows` writes CSV plus a JSON copy. (The plan said parquet; CSV/JSON needs no extra dependency.)
- `match_mean_photon` gives equal-energy comparisons.

**Presets** (`python -m cvqec.experiments <name|all> [--quick]`): squeezing, loss, detector, dephasing, asymmetry, decoders, steane_ancilla. All of them in quick mode take about 36 s.

**Verified:**
- metrics on hand-computable Pauli channels;
- F_avg equals the exact 2-design average;
- config immutability and the save/load round trip;
- energy matching;
- all presets build;
- headline fidelities unchanged from N = 140 to 170 (1e-6).

## First results (quick sweeps, equal energy n̄ ≈ 5.07, Δ = 0.3; preliminary)
1. GKP (square, hexagonal) beats the Fock reference at every squeezing level in 2–16 dB. (An earlier note said "~7 dB break-even"; that was a misreading of the quick plot.) The advantage factor is about 7–9 at 12 dB and 26–45 at 16 dB. Hexagonal is best, then square; rectangular (λ = 1.5) is worse under symmetric noise. Square GKP follows the Gaussian model 1−F = (2/3)(2P−P²), P = erfc(√π/(2√(2σ²))), σ² = e^{−2r} + Δ²/2, to within 0.4%.
2. Hexagonal plus Voronoi gives p_x = p_y = p_z. The square code has a suppressed p_y.
3. The 4-leg cat (α = 2.25) is below break-even for teleportation noise. Gaussian displacement noise is not its error model. Under pure loss its parity decoder works (η = 0.95: 0.98 vs 0.88 without correction), but it still trails the Fock qubit. The decoder ignores the amplitude shrinking α → √η α. GKP is better under loss at equal energy (0.9992), consistent with Albert et al. 2018.
4. Dephasing is GKP's weak spot: rotations move far lattice points far. Beyond σφ ≈ 0.12 rad it falls below break-even.
5. Asymmetric squeezing: the rectangular code gains as q-noise dominates, but λ = 1.5 over-corrects for ±4 dB splits. The matched λ is a natural follow-up.
6. Steane with Δ_a = 0.3 (10.5 dB) ancillas is close to break-even. Ancilla quality dominates (Δ_a = 0.15 gives 1.7e-2). Hex under Steane is worse than square, because of the parallelogram decoding.
7. Petz is at or slightly above ideal at high squeezing, and can be slightly worse than simple decoders at low squeezing. It's a benchmark, not an upper bound.

## Report draft (2026-11-04)
- The cat was dropped from the report and figures at the user's request (one code family). `code_lineup(include_cat=True)` still adds it.
- New preset `code_quality` (Δ sweep, N = 200). Asymmetry is widened to ±12 dB; the rectangular optimum is at −10 dB (1−F = 0.018 vs square 0.049 and hex 0.038 there).
- `report/`: `run_sweeps.py` (stage 1), `make_figures.py` (stage 2: PDFs, `generated/numbers.tex` macros, Pauli table), `main.tex`, `circuit.tex` (TikZ), `references.bib`, `overleaf.zip`. Compiled locally with Tectonic (added to the env).
- Petz beats ideal at 16 dB (2.5e-4 vs 1.1e-3), because the ideal decoder counts the code's own width Δ as noise.

## Open questions for discussion with the user
- Choose the final report figure set and resolution (non-quick sweeps take a few minutes each).
- Rectangular λ: tie it to the noise asymmetry, λ⁴ = var_p/var_q?
- A cat decoder that tracks amplitude decay, or a cat with larger α?
- A Steane variant that decodes hex with Voronoi cells (joint q/p decision)?
- Should the Fock engine also run small GKP states (N ≈ 40) as an end-to-end check including EC?

## Test audit, stage figure, validation macros (2026-11-14)
- **Audit.** Every test file was read against what it claims; all 102 were meaningful. A mutation check (10
  deliberate bugs in a scratch copy) caught all ten, but one (Steane q-correction sign) only with a new test, so
  the old Steane test (monotonic in ancilla quality) was not enough. Details: `docs/TEST_PLAN.md`.
- **20 new tests (122 total):** EPR entanglement with and without loss (Duan sum vs closed form, log-negativity),
  separable bound for a product state, Bell-measurement outcome variances; Steane keeps all Pauli eigenstates (3
  lattices), finite ancillas add error, hex Steane < ideal separate < ideal joint, rectangular wins only under
  q-biased noise; `tests/test_stages.py` (9 tests) for the stage figure.
- **Library:** `stage_params(p)` (the circuit with imperfections switched on in order: shift noise, + loss,
  + dephasing) and `bob_before_feedforward(p, N)` (Bob's Gaussian state before any correction; thermal with
  n = eta sinh^2 r for r1 = r2) in `cvqec/circuits/teleportation.py`.
- **Figure `stages.pdf`** (`python report/make_figures.py stages`): square GKP |+_L> at 10 dB, eta = 0.9,
  sigma_phi = 0.1: input, Bob before feed-forward, after feed-forward, + loss, + dephasing, after ideal correction
  (decoded qubit re-encoded for display). Average fidelities 1.000, 0.974, 0.940, 0.904. Each stage is checked
  against the full three-mode simulation. The window must be extent 8: at 6.5 a Delta = 0.3 state loses 2%.
- **`make_figures.py`** now takes job names (`stages`, `numbers`, ...) and writes validation macros
  (`\valTraceDist` = 8.5e-07 and its parameters) by rerunning the central check. Existing macros unchanged.
- **Decoder split** (computed for the notes, not yet a preset): 1-F at 12 dB, square ideal/separate/Steane
  0.0093/0.0093/0.0487, hexagonal 0.0070/0.0153/0.0585; Steane with no channel noise: 1-F = 0.015 (Delta_a = 0.3).

## Joint Steane decoding and syndrome-detector inefficiency (2026-11-14)
- **Why:** the hexagon's loss under Steane came from correcting Q and P separately (parallelogram), not from the
  measurement. A Steane round measures both syndromes, so it can decode jointly (nearest lattice point).
- **How:** in the square frame, joint = separate correction + a syndrome-dependent logical shift (a, b) in {-1,0,1}
  units of sqrt(pi), chosen so the implied error is shortest in the code's own phase space (`joint_table`, 16 x 16 for
  M = 16). The q-shift is postponed until both outcomes are known (it commutes with the p-round), and terms are grouped
  by the q-shift, so at most 3 inverse transforms per q-outcome: ~1.2 s per point for hex. Square/rect: table empty,
  old fast path. Default `joint=True`.
- **Syndrome detectors:** `eta_syndrome` smooths each outcome bin with a Gaussian of variance (1-eta)/(2 eta)
  (POVM of an inefficient homodyne; ancilla discarded). eta = 1 reproduces the old round exactly.
- **Results (12 dB, 1-F, Delta_a = 0.3):** square ideal 0.009, Steane 0.049 (joint = separate); hexagonal ideal joint
  0.007, ideal separate 0.015, Steane separate 0.058, Steane joint 0.048; with eta_s = 0.9: 0.090 / 0.094. Ancilla sweep:
  hex joint 0.014 at Delta_a = 0.15 (square 0.017). Reran presets decoders, steane_ancilla, new decoder_split; old
  results archived in results/archive_2026-11-14_separate_steane/. 5 new tests (127), mutation M11 caught.

## All imperfections at once (2026-11-14)
- Until now every sweep varied ONE imperfection (others ideal). New presets: `combined` (4-16 dB, eta = eta_d = 0.95,
  sigma_phi = 0.05; ideal decoder for all codes; Steane joint with eta_s = 0.95 and Delta_a = 0.3 or 0.15 for square
  and hex) and `budget` (12 dB, imperfections added one at a time, then the Steane round and syndrome detectors).
  The operating point is illustrative (constants REALISTIC, REALISTIC_ETA_SYNDROME in presets.py), not from an experiment.
- Figure `combined.pdf` (make_figures job `combined`), macros `\comb*`, `\budget*`.
- Results at 12 dB (1-F): ideal square 0.059, hex 0.053, rect 0.122, bare qubit 0.121; realistic Steane 0.130 / 0.131
  (never better than the bare qubit in 4-16 dB); with Delta_a = 0.15 ancillas 0.091 / 0.088 (better everywhere).
  Budget (square): 0.009 squeezing, 0.020 + loss, 0.050 + detectors, 0.059 + dephasing, 0.110 + Steane, 0.130 + syndrome.
  Biggest steps: detectors and the Steane round. One new test (budget stages cumulative): 128 tests.


## 2026-11-14 (evening): full rerun for report_v5, separate figures, macros with meanings
- Presets: named DELTA / TELE_DB / LAM; eta_syndrome removed from combined, budget, decoder_split (not modelled in the
  report); code_quality + rectangular; new rect_lambda (lambda 1-2 at bias 0, -6, -12 dB) and cat_reference.
- Full rerun of all 13 presets (~12 min). Old data archived in results/archive_2026-11-14_before_rerun/.
- make_figures.py rewritten: --out folder, one file per panel (PDF + 600 dpi PNG), legends below the axes with codes
  (colour/marker) left and decoders (line style) right, dotted line at the default operating point, Gaussian model drawn
  on top; Macros class writes numbers.tex with a comment per macro and docs/results_numbers.md (exact values, sources).
  Every infidelity macro has an F twin (<name>Fid). Stage figure at 12 dB, eta = 0.9, sigma_phi = 0.1 (its panels c, d
  equal the 12 dB squeezing and eta = 0.9 loss values: a cross-check).
- Findings: rect_lambda optimum = Gaussian-model prediction (1.0, 1.2, 1.5 vs 1.00, 1.22, 1.49); cat (alpha 2.25) 1-F
  0.172 vs bare 0.066 at 12 dB; all-at-once Steane (Delta_a 0.3) 0.110 vs bare 0.121 at 12 dB (break-even ~4.2 dB).
- Tests: 130 (test_reported_runs_use_joint_steane_with_ideal_syndrome_detectors, test_default_operating_point).
- report_v5: results_section uses the panel files (+ new fig:rect_lambda); runs table rows/caption updated [C7-C14];
  stage-figure caption [C13]. Compiles (20 pages), ASCII-only.
- Follow-up: paired panels now use a fixed axes box (margins in inches, titles at a fixed height, full canvas saved)
  so side-by-side plots align exactly; one shared legend file per figure (legend_squeezing/noise/decoders/combined.pdf),
  included under the panels [C15].

## 2026-11-14 (late): derivations appendix

- `report_v5/appendix_derivations.tex` (new, Claude's draft, all `\color{chg}`): 18 items [A1]-[A18] deriving every
  formula the report uses without proof (displacements, phase shifter R(phi), squeezer/BS/dB, Wigner and characteristic
  functions, thermal state, finite GKP energy nbar = 1/(e^{2 Delta^2}-1), Kraus/Stinespring/Lindblad/Liouvillian, loss
  (BS and Lindblad routes), dephasing, shift noise, channel order, deferred measurement, Duan, validation numbers,
  fidelity/PTM/Pauli rates, Gaussian model and lambda^4 rule, ideal decoder via chi, Petz). Not yet in main.tex: add
  `\appendix` + `\input{appendix_derivations}` after the conclusions. Scratch build compiles (35 pages, no overfull boxes).
- `report/make_figures.py`: new `appendix_numbers` -> 14 closed-form macros `\app*` (existing macros unchanged).
- `report/appendix_to_html.py`: converts the appendix tex (its LaTeX subset) to `docs/appendix_derivations.html` with
  MathJax; macros from numbers.tex, section numbers from main.aux.
- Found: main.tex's `\bkt`/`\bkkt` do not compile (\middle inside a brace group); unused by the report, so the appendix
  defines `\abk`/`\abkk`. The course's rotation R(theta) (p.40) has the opposite sign to the phase shifter (p.50) we use.

## 2026-11-15: report_v6 round (bridge appendix, E2, E3, powers of ten)

- report_v6 (user's copy) backed up to `report_v6/archive/2026-11-15_before_round/`; previous \chg rounds accepted.
- `report_v6/appendix_bridge.tex` [B1]-[B10]: only non-course material the text uses (Kraus/Lindblad link, loss Kraus
  from BS and Lindblad, deferred measurement, Bob's thermal state, Duan sum after loss, finite GKP nbar, detector vs
  channel loss, fidelity/Pauli rates/trace distance, Gaussian model + lambda rule, Petz). Input from main.tex [C25].
- Main-text edits [C1]-[C25] (coloured): appendix pointers, fixed `\parCodeDB` -> `\parTeleDB` (12 dB choice), D1 bullet,
  Duan numbers (`\appDuanSum`, `\appDuanMax`), wrong dephasing numbers and the 19.1 "factor", bare-qubit dephasing
  explanation, anchor numbers, E3 table and text, conclusions nuance.
- E3: preset `sensitivity` (eta = eta_d = 0.90, 0.97; ~2 min) -> `generated/sensitivity_table.tex`, macros `\sens*`.
- E2: Steane grid M = 16 vs 24 converged (test added; 131 tests).
- make_figures.py: default --out report_v6; `powers_of_ten` (7e-07 -> \ensuremath{7\times10^{-7}}); tables in the user's
  vertical-rule style; `duan_numbers`, `sensitivity_numbers`.

### 2026-11-16: E5, code width under dephasing
- Preset `delta_dephasing` (~3 min): code squeezing 6-16 dB at sigma_phi = 0 and 0.1 rad, 12 dB teleportation, N = 200
  (tail <= 3e-5). At 0.1 rad square/hex have an optimum at ~9.5 dB (Delta ~ 0.33, 1-F = 0.040 / 0.036); above ~14 dB
  they are worse than the bare qubit (0.068); rectangular is worse than the bare qubit at every Delta. Default
  Delta = 0.3 is near the optimum. The 10.5 dB point matches the dephasing sweep.
- make_figures.py: job `delta_dephasing` -> figures/delta_dephasing.pdf; `delta_dephasing_numbers` -> 37 `\eFive*`
  macros (only additions to numbers.tex). Report text not edited: a paste-in snippet was handed to the user.

### 2026-11-16: report_v8 = review round on the user's report_v7 (proposals only)
- report_v8/ is a copy of report_v7 (the user's copy, untouched) with \chg proposals [C1]-[C25]; answers the user's todos.
- Gaussian model moved into the main text as Section 2.6 (sec:model): every imperfection is an extra shift variance
  (code Delta^2/2, detectors (1-eta_d)/eta_d, loss (1-eta)/(2 eta), rotation sigma_phi^2/(2 Delta^2)). Checked against the
  sweeps: exact for squeezing/detectors (4 digits), ~10-25% for loss/dephasing; predicts the E5 optimum Delta = sqrt(sigma_phi).
  New macros model* (make_figures appendix_numbers); thin black model line added to figures/delta_dephasing.pdf.
- Thermal state: in the course book's Appendix A (p.107-109); Gaussian states p.33. Bridge item rewritten as a 3-step argument.

### 2026-11-16: report_v10 (todo round on the user's report_v9)
- report_v10 = copy of report_v9; earlier \chg rounds uncoloured. New changes C30-C33 (blue):
  C30 Appendix A.4 (br:gkp) removed, kappa = Delta reason is one clause in code.tex; C31 component table removed,
  non-trivial components described in the Fig. 3 caption, runs_section points to the figure; C32 loss-vs-detector
  paragraph rewritten (shrinking, variance c^2/2, factor (sqrt2)^2); C33 Gaussian-model formulas and comparison moved
  to Appendix br:model, Section 2.6 keeps a short description plus the three uses. Builds clean (30 pp).

### 2026-11-16: validation iteration on report_v12 -> report_v13
- Tests: 131 passed (67 s). Ten headline macros re-run from the current code (square/hex/rect/bare at 12 dB, loss and
  detectors 0.9, dephasing 0.1 square/bare, Steane budget point, Petz): all reproduce to the printed digit.
- Every derivation in code.tex, circuit_section.tex and both appendices re-checked by hand; derived numbers recomputed in
  closed form (Gaussian model, Duan sums, Bob's nbar, lambda predictions, cell inradii). One error: the lambda derivation
  wrote P_p = 1 - erf(lambda c/s_p) instead of 1 - erf(c/(lambda s_p)) (result was right). Fixed as C38.
- Truncation: maxTail now covers both N = 200 sweeps (3e-5, was 7e-7 from code_quality only; make_figures.py changed).
  N = 200 vs 250 at the narrowest codes: <= 1.2 % relative change in 1-F (hex, 16 dB, no dephasing), else <= 1e-4.
- report_v13 = v12 with earlier rounds uncoloured and C38-C47 (accuracy fixes, typos). Builds clean, 28 pp.

### 2026-11-16: report_v13, second round (readability R1-R4, sim-vs-report, citations)
- C48-C51: rotation paragraph in code.tex shortened (details stay in the appendix); three-item list before E removed;
  the experiment list is one line per run; thermal-appendix meta sentence removed.
- C52-C54 (simulator details the report did not state): Steane result includes a perfect logical readout after the
  round; the bare qubit is not corrected (projection onto {|0>,|1>}, leakage -> maximally mixed); the 3-mode
  validation uses one test state (|0> + 0.7i|1>, N = 25).
- Citations, each claim located in the paper:
  gkp2001 Sec. X "Error Recovery", Eq. (105), Fig. 4 (|+> ancilla, SUM data->ancilla, q mod alpha; |0> ancilla for p);
  barnum2002 Theorem 3, Eq. (13) + corollary (error <= 2x optimal), map Eq. (12) with rho = P/d;
  nielsen2002 Eq. (3) F = (d F_e + 1)/(d + 1), due to Horodecki et al. PRA 60, 1888 (1999);
  duan2000 Theorem 1, Eq. (3), [x,p] = i, a = -1 gives the report's u, v and bound 2;
  albert2018 Eqs. (1.4)-(1.6): loss Kraus operators and beam splitter with vacuum (now cited for the loss model,
  replacing braunstein2005, which does not state it); braunstein2005 Sec. IV.A pp. 28-30: CV Bell detection with a
  beam splitter and two homodynes (restored sentence, C56); massar1995 bound (N+1)/(N+2) -> 2/3 (attributed in
  Derka-Buzek-Ekert quant-ph/9707028; also braunstein2005 p. 32); vahlbruch2016 15 dB; vaidman1994 CV generalization in
  the abstract; braunstein1998 finite squeezing + nonideal detection; furusawa1998 first unconditional CV teleportation
  (F = 0.58); menicucci2014 20.5 dB threshold for fault-tolerant CV-cluster MBQC (open \cnote: keep or drop).
- Citation round 2 (user rule: cite the paper that states a result; drop long derivations it already contains):
  removed menicucci2014, furusawa1998, massar1995 (2/3 sentence), qutip5 (QuTiP named, not cited) -> C58-C61.
  C62 fidelity appendix: Bloch-sphere derivation replaced by Nielsen Eq. (3); C63 Duan proof replaced by Theorem 1;
  C64 loss Kraus derivation replaced by Albert et al. Eqs. (1.4), (1.6); C65 Sec. 2.3 derivation credited to
  Braunstein-van Loock Sec. IV.A.1 Eqs. (161)-(164); C66 Var(m_q), Var(m_p) block removed. 9 references, 26 pp.

## 2026-11-16 (late): report_v14 round C67-C78 (results explanations, thermal appendix, instructor wording)
- Accepted C62-C66 (uncoloured); pre-round copy in `report_v14/archive/2026-11-16_before_C67/` (incl. make_figures.py.bak).
- C67 removed appendix "Bob's mode ... thermal state" (kept as comments); C68 Fig. 4(b) caption now gives the no-signalling
  reason (course p.115). Stage panel title no longer says "thermal" (`make_figures.py stages` rerun).
- C69 "see below" -> Section~\ref{sec:modeling_reality} (Energy item). C70 runs: E is the same circuit with modes a, b and
  the outcomes removed analytically; Steane round has a real finite-Delta_a ancilla. C71: F-bar averages over all inputs.
- C72-C76 results explanations checked against results/*.csv: Pauli caption (hex: three equivalent directions; square
  p_y ~ p_x p_z), rect crosses bare qubit at eta ~ \lossCrossRect, detectors: rect slightly better than bare everywhere,
  rect dephasing crossing, bare-qubit dephasing p_z = (1-e^{-s^2/2})/2 (predicted rise 0.0017, simulated 0.0015), rect beats
  hex only for bias < \asymRectBeatsHex dB.
- C77 new Section 4.7 "Summary: how each code corrects each error" + generated/summary_table.tex (new job function
  `summary_table` in make_figures.py; macros sumNoisyQ*/sumNoisyP*). C78 conclusions: rect only for noise along its wide
  side (q), pointer to the table. New macros: lossCrossRect, dephFockPz/RiseModel/RiseSim, asymRectBeatsHex (426 total).
- Compiles clean (27 pp), ASCII OK, no undefined refs.

## 2026-11-17: report_v15 Appendix B derivations (C79, C80)
- C79 B.1 rewritten step by step: one shift X(u) is a phase per p-cell, Gaussian average (completing the square) gives
  e^{-sigma_q^2 (p-p')^2/2}, meaning; Lindblad L1 = sigma_q p solved explicitly gives the same factor at t = 1 (why L1 is right).
- C80 B.3 uses the B.1 method with n in place of p (Fock cells), meaning (populations kept, far coherences lost), and the
  Lindblad solution with an explicit arrow to t = sigma_phi^2. Pre-round copy: report_v15/archive/2026-11-17_before_C79/.
- \xRightarrow needs mathtools (not loaded): used \overset{..}{\Longrightarrow}. Compiles clean, 28 pp, ASCII OK.
- (same day) C79 shortened at the user's request: no integral derivation; B.1 is now "why L1" + one arrow line (phase per p-cell,
  Gaussian average) + the Lindblad solution with an arrow. B.3 step 2 now says "the same Gaussian average as in B.1".
- (2026-11-17) C81 Appendix A/B cross-references, in place (no new lines): primes + Section 2.3.2 pointer in the Duan
  paragraph, loss substitution on Bob's mode q_b', variance sources (Eqs. var_quads, squeezed_quad_q), loss variance ->
  Section 2.4.2, B.2 beam splitter -> sec:loss_channel; new inline label sec:epr in circuit_section.tex. CLAUDE.md and
  HANDOFF.md now point at report_v15; CLAUDE.md corrected: the user did learn Lindblad (Quantum Technologies).

## 2026-11-17: submission and repository cleanup
- Report submitted (final text in the user's Overleaf). Removed `report_v2/`-`report_v15/`, the old report LaTeX/PDF/zip,
  unused figures (`results/*.png`, `results_ext/figures/`), notes (`docs/*.html`, FINALIZATION_PLAN.md), the audit script
  and its output (it compared against `report_v2/`), the original task prompt (now in the lessons zip). Full backup of the
  pre-cleanup tree: `../GKP-QEC-simulator_backup_2026-11-17.tar.gz`.
- Kept: library, tests (131 pass), all sweep data, `report/run_sweeps.py`, `report/make_figures.py` (default `--out report`),
  the 18 figures the report used + TikZ circuit in `report/figures/`, scripts_ext/results_ext data, docs.
- Added `docs/LESSONS_LEARNED.md`, `claude_harness_lessons.zip` (lessons, skills, memory, context files), a README;
  CLAUDE.md rewritten for the final layout.
- (same day) READMEs: `results/README.md` (each preset file -> report figure/table, columns, archives), `report/figures/README.md`
  (each file -> figure number), status header + file table in `results_ext/README.md` and `docs/EXTENSIONS_SUMMARY.md`
  (not used in the report). Recovered the Fig. 2 TikZ (hexagonal cell) from the backup into `report/figures/hexcell_tikz.tex`.
  `results_ext/figures/` ignored (recreated by make_ext_figures.py).
- (same day) README gallery: circuit and hexagonal cell rendered from TikZ (report/figures/circuit.pdf, hexcell.pdf + png/),
  README now shows circuit, code Wigner functions, hexagonal cell, stage figure, two result plots, key results, link to
  report/full report.pdf.
