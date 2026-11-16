# Test plan

## How to run
```bash
conda env create -f environment.yml      # once
conda activate qcsim && pip install -e . # once
pytest                                   # 130 tests, ~65 s
pytest tests/test_engines.py -v          # one area
pytest --durations=8                     # slowest tests
python -m cvqec.experiments all --quick --out /tmp/smoke  # smoke-run every sweep (~40 s); ALWAYS pass --out:
                                                          # the default (results/) overwrites the report's full-resolution data
```

## Principles
1. **Every test checks physics, not plumbing,** wherever possible: a statement from the course book or an analytic result, with its page or derivation in the test's comment.
2. **Independent cross-checks beat self-consistency.** The strongest tests compare two derivations that share no code:
   - full circuit vs effective channel;
   - characteristic-function decoder vs direct q-binning;
   - Kraus vs master equation;
   - Laguerre vs mpmath.
3. **Tolerances follow from the numerics, not from what passes.** Each tolerance is about 10× the measured discretization or truncation error, and that error is documented next to it.
4. **Truncation is tested by convergence** (results unchanged as N grows), never by tail weight alone.

## Coverage by milestone

| Area | File | What is verified | Kind of check |
|---|---|---|---|
| Conventions | `test_conventions.py` | [q,p] = i below the cutoff; vacuum variance 1/2; coherent-state position (p.38); q_θ; dB conversions | analytic |
| Gates | `test_gates.py` | unitarity; R, S, B, D, X, Z, SUM Heisenberg actions (p.50–53); squeezer variances (p.39); HOM dip; displacement commutation phase (p.29); EPR correlations; convergence in N; `on_modes`; truncation warning; **EPR entanglement** (Var(q_a-q_b)+Var(p_a+p_b) < 2 and log-negativity > 0, with and without loss, against the closed form); product state respects the separable bound; **Bell-measurement outcome variances** Var m_q, Var m_p (report Section 2) | analytic and convergence |
| Codes | `test_codes.py` | isometric encoder; anticommuting logical Paulis; codewords fit in N; registry; n̄; lattice condition βα*−β*α = iπ (p.80); match with the course's ψ(q) GKP (p.82); ⟨S⟩ = e^{−\|γ\|²Δ²/2}; D(α), D(β) act logically; Δ→0 trend; rectangular trade-off; cat mod-4 structure and parity syndrome (p.77) | analytic and independent construction |
| Measurement | `test_measurement.py` | homodyne of a coherent state at any θ (p.54); ⟨p\|n⟩ is the Fourier transform of ⟨q\|n⟩; exact displacement vs padded expm (small α) and vs mpmath (\|α\|² = 529); group law; Gaussian unitary from a symplectic matrix | analytic and high precision |
| Noise | `test_noise.py` | CPTP; each channel = `qutip.mesolve` of its Lindblad form; adjoint identity; loss on coherent states; loss composition; dephasing Kraus = exact mask; displacement adds var_q and var_p; `Compose` order | independent (master equation) |
| Engines | `test_engines.py` | channel deferral is legal; coherent teleportation F = 1/(1+e^{−2r}); **full three-mode circuit = effective channel** (all imperfections, trace distance < 1e-5); effective-channel parameters; component table; noiseless logical channel = identity | analytic and independent |
| Decoders | `test_decoders.py` | cell Fourier transform (square = sinc², hex area = π); ideal decoder = direct q-binning; grid decoder = characteristic-function decoder; hex frame change = parallelogram cells; codeword recovery for all lattices; physical output; Voronoi advantage and X/Y/Z symmetry for hex; Steane → ideal as Δ_a → 0; **Steane round keeps all Pauli eigenstates** (syndrome does not read the logical value; all lattices); **finite ancillas add their own error**; **hex Steane < ideal separate < ideal joint**; **rectangular wins only under q-biased noise**; **joint Steane**: empty lookup table for square/rectangular, inversion-symmetric table for hex, joint beats separate and even ideal-separate; **inefficient syndrome detector** (eta_s = 1 reproduces the old round, monotone, trace preserving); cat undoes one loss; Petz Barnum–Knill bound; auto decoder choice | independent and physics |
| Metrics | `test_metrics.py` | identity channel; Pauli channel (F_e, error rates, PTM); F_avg formula = exact 2-design average | analytic |
| Experiments | `test_experiments.py` | config copies don't mutate; sweep grid size; run, save and load round trip; energy matching; every preset builds; error-budget stages add imperfections cumulatively | plumbing |
| Stage figure | `test_stages.py` | stages switch imperfections on in order; **each stage = full three-mode circuit with only those imperfections** (trace distance < 1e-5); Bob's state before feed-forward = partial trace of the lossy EPR pair, thermal with n̄ = η sinh²r for r1 = r2; first stage adds exactly σ_q², σ_p² to a GKP state; the Wigner window (extent 8) holds the whole state; decode + re-encode round trip | independent and physics |
| Convergence | `test_convergence.py` | logical fidelity of hex GKP, square GKP and cat unchanged from N = 140 to 170 (1e-6) under a realistic channel | convergence |

## Tolerances and why

| Quantity | Tolerance | Measured error |
|---|---|---|
| EPR covariances (N = 40) | rel 1e-6 | 1.1e-7 |
| Channel vs mesolve | 1e-6 (norm) | ~1e-8 |
| Coherent teleportation fidelity (N = 25, 81-point outcome grids) | 1e-6 | 8e-8 |
| Full engine vs effective channel | trace distance 1e-5 | 8.5e-7 |
| Ideal decoder vs direct binning | 1e-5 | 2e-6 (limited by the direct quadrature) |
| Grid decoder (M = 32) vs characteristic-function decoder | 1e-4 | 5e-5 (O(h²); 2e-4 at the default M = 16) |
| Hex grid vs parallelogram (M = 16, padded squeezer) | 1e-3 | 2e-4 |
| Convergence N = 140 → 170 | 1e-6 | < 1e-7 |
| Bell-measurement variances (N = 25) | rel 1e-4 | passes |
| EPR Duan sum with loss vs closed form (N = 25) | rel 1e-5 | passes |
| Each figure stage vs full circuit (N = 25) | trace distance 1e-5 | passes |
| Wigner weight inside the figure window (extent 8) | 1e-3 | 6e-4 (extent 6.5 lost 2%) |

## Mutation check (2026-11-14): do the tests catch real bugs?
Ten deliberate bugs were put into a scratch copy of `cvqec/` (never into the repo) and the suite was run against each.
All ten were caught.

| Bug | Caught by |
|---|---|
| M1 squeezer sign flipped | squeezer variances, EPR correlations, coherent teleportation, full = effective |
| M2 detector noise (1-eta_d) instead of (1-eta_d)/eta_d | full = effective, effective-channel parameters |
| M3 loss applied before the shift noise in E | full = effective |
| M4 feed-forward sign of X flipped | coherent teleportation, full = effective |
| M5 dephasing exponent /4 instead of /2 | mesolve check, Kraus = mask |
| M6 q and p noise swapped | displacement adds variance, full = effective, every stage, rectangular bias |
| M7 Steane q-correction shifted the wrong way | **only the new** `test_steane_round_does_not_read_the_logical_state` (the original 102 tests missed it) |
| M8 default decoding cell = parallelogram | hex Voronoi advantage, hex Steane ordering |
| M9 F_avg = (F_e+1)/2 | 2-design average |
| M10 feed-forward gain 1 instead of sqrt(eta) | full = effective, stages with loss |
| M11 (2026-11-14) wrong p-sign in the joint lookup table | `test_joint_steane_recovers_the_hexagonal_advantage` |

## Manual and visual checks (not automated)
- Wigner gallery of every code (`viz.wigner.plot_code`): peak positions, lattice orientation, cat fringes.
- Sweep figures (`python -m cvqec.experiments all --quick --out <scratch>`):
  - curves are monotone in the noise parameter;
  - the Fock reference is smooth;
  - no curve sits above 0.5 infidelity.
- The `truncation_tail` metric column is below 1e-8 in every saved row.

## Known gaps and next tests to add
- Steane for non-square lattices uses parallelogram cells. There's no test yet for a Voronoi-aware Steane variant (none exists).
- The Fock engine is validated only on low-energy inputs (N = 25). An end-to-end GKP run through the three-mode engine is infeasible at N = 140 per mode.
- Report-resolution sweeps (non-quick) are not part of `pytest`. Run them explicitly and archive `results/*.json` with the figures.
- Feed-forward before the last measurement is not supported by the engine. It raises `NotImplementedError`, which is not yet covered by a test.

## Added 2026-11-14 (evening): run-policy tests (tests/test_experiments.py)
- `test_budget_adds_imperfections_cumulatively`: updated, 5 stages; the last swaps the ideal decoder for the Steane round.
- `test_reported_runs_use_joint_steane_with_ideal_syndrome_detectors`: every Steane run in every preset is the joint
  round (the separate one only where labelled "separate") and never sets eta_syndrome (the report lists it as not
  modelled). Catches a preset that silently reintroduces the old decoder or eta_s.
- `test_default_operating_point`: DELTA = 0.3 and TELE_DB = 12; every fixed-point sweep (loss, detector, dephasing,
  steane_ancilla, decoder_split) uses exactly 12 dB teleportation and the code Delta = 0.3.

- `tests/test_convergence.py::test_steane_grid_converged` (2026-11-15, E2): the Steane round's grid M = 16 gives the same
  1-F as M = 24 within 3 % (square, 12 dB; measured 1.4 %, hexagonal 0.4 %).
