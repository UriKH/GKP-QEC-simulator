# Extension experiments (results_ext/)

**Status (2026-11-17): not used in the submitted report.** These are side studies from 2026-11-10 that go beyond
the assignment: GKP vs the 4-leg cat under rotations and under shifts (E1, E2), a wrong feed-forward gain (E3), the
exact optimal decoder as a benchmark (SDP), and a convergence check. They were planned for an optional "beyond this
project" section that was never included. The report's single cat number (Section 2.2) comes from
`results/cat_reference.csv`, not from here. Short summary: `docs/EXTENSIONS_SUMMARY.md`.

Scripts are in `scripts_ext/`. Nothing in `results/`, `report/` or `cvqec/` depends on this folder. Every row in the
CSV/JSON files stores its full configuration (code, code parameters, decoder, channel parameters, Fock cutoff N).
The figures (`figures/complementarity.pdf`, `figures/gain_mismatch.pdf`) are not stored; `make_ext_figures.py`
recreates them.

| file | experiment | made by |
|---|---|---|
| `e1.csv/.json` | E1: pure dephasing only, sigma_phi = 0.02-0.40 rad; square/hex GKP, cat, bare qubit; standard + Petz decoders | `run_ext.py e1` |
| `e2.csv/.json` | E2: teleportation channel 2-16 dB (no other imperfection); same codes and decoders | `run_ext.py e2` |
| `e3.csv/.json` | E3: feed-forward gain g = 0.85-1.15 at 10 dB; square GKP, bare qubit | `run_ext.py e3` |
| `optimal.csv/.json` | the optimal (SDP) decoder for the cat and the bare qubit on the E1/E2 grids, with certified bounds | `run_optimal.py` (needs cvxpy) |
| `conv.csv/.json`, `convergence_summary.json` | the hardest E1-E3 points recomputed at N = 180 vs 140 | `run_ext.py conv`, `check_convergence.py` |
| `verify_gain_channel.json`, `verify_gain_convergence.json` | the general-gain channel checked against the full three-mode simulation | `verify_gain_channel.py`, `verify_gain_convergence.py` |
| `key_numbers.json/.txt` | every number quoted in this README | `make_ext_figures.py` |

The rest of this file is the original write-up (2026-11-10). The audit section describes a check against an early
report version; its script (`audit.py`) was removed with the old report folders on 2026-11-17.

## Why

GKP codes (translation symmetry) correct small shifts but are weak against random rotations e^{-i phi n}.
Rotation-symmetric codes such as the 4-leg cat are the opposite (Grimsmo, Combes & Baragiola, PRX 10, 011058, 2020).
The project tested only one side of this: the cat under teleportation (shift) noise. Two problems with that test:
- the cat's parity decoder only undoes loss, not rotations;
- Petz is not a guaranteed "best case".

The extension runs both sides of the comparison, with decoders that give each code its best chance:

| | channel | sweep | codes | decoders |
|---|---|---|---|---|
| **E1** | pure dephasing only (no shifts) | sigma_phi = 0.02 ... 0.40 rad, step 0.02 | square GKP, hex GKP, 4-leg cat, Fock qubit | standard + Petz; optimal (SDP) for cat and Fock |
| **E2** | effective teleportation channel (eta = eta_d = 1, no dephasing) | 2 ... 16 dB, step 1 | same | same; SDP on 2-16 dB (Fock) and 8-16 dB (cat), even dB |
| **E3** | teleportation with feed-forward gain g != sqrt(eta) | g = 0.85 ... 1.15, step 0.01, at 10 dB, eta = eta_d = 1 | square GKP, Fock qubit | standard + Petz |

The codes are compared at equal energy:
- GKP: Delta = 0.3, so n-bar = 5.071 for both lattices.
- Cat: alpha = 2.2521, bisected to the same n-bar with the same method as `presets.code_lineup(include_cat=True)`; n-bar = 5.071.
- Fock qubit {|0>, |1>}: n-bar = 0.5. It is the break-even reference.

Standard decoders: `gkp_ideal` for GKP, `cat_parity` for the cat, and `projection` for the Fock qubit (no correction).

## Changes to the requested design (and why)

1. **An optimal-recovery benchmark was added: an SDP.** In the audit, Petz did worse than plain projection
   for the Fock qubit.
   - At 10 dB teleportation: F_avg 0.875 with Petz vs 0.901 with projection.
   - Under pure dephasing, Petz is exactly 2x worse than doing nothing, because it applies E^dag, i.e. dephases a second time.

   So "Petz for every code" is not a fair best case. `scripts_ext/optimal_recovery.py` solves the exact
   optimal recovery (maximize F_e over all recovery channels) as a semidefinite program (cvxpy + Clarabel). It
   restricts the input space to the support of the channel output and certifies each result with an
   independent dual-feasible point.
   - All 53 SDP points are certified to F_e gap < 8.1e-7, i.e. under 5.4e-7 in F_avg (`optimal.json`, column `F_e_opt_upper`).
   - Neglected support weight: < 3.3e-12.
2. **GKP has no SDP.** Its output support (68-140 Fock directions) makes the SDP too slow here. For GKP we use
   the rigorous Barnum-Knill bracket instead: F_e(Petz) >= F_e(opt)^2. So the optimal GKP infidelity lies between
   (2/3)(1 - sqrt(F_e,Petz)) and the Petz value. That band is shaded in the figure.
3. **The SDP grid is coarser for E2.** The SDP cost grows quickly at low squeezing (cat: d = 59 at 2 dB). So it
   runs on even dB only, from 8 dB for the cat. In E2 the cutoff for keeping a support direction was 1e-11
   (E1: 1e-13), and the neglected weight is added to the certified upper bound.
4. **The Reimpell-Werner iteration was tried and rejected.** It is kept in `optimal_recovery.py` for reference
   only. It stalled at values that depended on the starting point, e.g. 0.753 instead of the true 0.7755 for
   the cat at 10 dB.
5. **E3 needed a new channel.** `cvqec.circuits.teleportation.effective_channel` refuses g != sqrt(eta). The
   general case is derived and implemented in `scripts_ext/gain_channel.py` without changing the library. It
   is verified against the full three-mode Fock engine (see below).

## Audit of the existing results (2026-11-10; historical, the script and its output were removed)

- **Report numbers match the data.** All 53 macros in `report_v2/generated/numbers.tex` and the Pauli table
  were regenerated from `results/*.json` (into `results_ext/audit/`), and they are **identical**.
- **Saved rows match the current code.** Three saved rows were rerun with the current code, with identical
  results (difference 0.0):
  - hexagonal GKP at 12 dB;
  - square GKP in the dephasing sweep;
  - square GKP with Petz at 16 dB.
  E2 also reproduces every shared row of `results/squeezing.json` exactly.
- **HANDOFF section 6 claims, all supported:**
  - Square and hexagonal GKP beat the bare qubit at every dB from 2 to 16: by ×7.1 / ×9.5 at 12 dB, ×25.8 / ×44.9 at 16 dB.
    At 2 dB the margin is thin, a ratio of 1.05 / 1.06. Rectangular GKP also beats it, by at least ×1.01.
  - The Gaussian model is within 0.37% (max, at 2 dB).
  - Hexagonal GKP has p_x = p_y = p_z.
  - Loss 0.038 vs detector 0.076.
  - Asymmetric squeezing: 0.018 / 0.049 / 0.038.
  - Steane: 0.026 at 16 dB, and 0.049 -> 0.017 with better ancillas.
  - Petz 2.5e-4 vs ideal 1.1e-3 at 16 dB.
  - Dephasing crossovers 0.137 / 0.142 / 0.063.
  - Cat vs bare qubit at 10 dB: 0.772 vs 0.901, recomputed.
- **Minor inconsistency, docs only:** `docs/TEST_PLAN.md` says the `truncation_tail` column is "below 1e-8 in every saved row".
  - Only `code_quality.json` has that column, and its maximum is **7.5e-7** at N = 200, the `\maxTail` macro.
  - The other sweeps have no tail column; their convergence was checked by N = 140 -> 170 tests instead.
- **No data problem.** Nothing in results/ needed rerunning.
- **New context for an existing claim.** The "cat 0.77 vs bare 0.90" claim was measured with the parity decoder and projection. With
  optimal recovery for both, at 10 dB, it is F_avg 0.850 (cat) vs 0.915 (bare). So the conclusion that the cat loses to the bare qubit under teleportation noise does not
  depend on the decoder.

## Truncation convergence (conv.json, convergence_summary.json)

The hardest points of each experiment were recomputed at N = 180 (vs 140), for all codes and both decoders:
- E1: sigma = 0.2 and 0.4.
- E2: 2, 10 and 16 dB.
- E3: g = 0.85, 1.0 and 1.15.

Result:
- Largest change in 1 - F_avg: 3.3e-8 absolute, 2.1e-6 relative, at E3 g = 1.15 with Petz, where the amplifier raises photon numbers.
- Everything else changed by < 2e-6 relative.

The SDP was run at N = 140 only. Its input is the converged channel output, and the neglected weight is < 3.3e-12.

## Verification of the general-gain channel (verify_gain_channel.py, verify_gain_convergence.py)

Heisenberg derivation, in course units (EPR from a q-squeezed mode A, r1, and a p-squeezed mode B, r2).
The output is q_out = g q_in + nu_q, p_out = g p_in + nu_p, with

    Var nu_q = [(sqrt(eta) - g)^2 e^{2 r2} + (sqrt(eta) + g)^2 e^{-2 r1}]/4 + (1 - eta)/2 + g^2 (1 - eta_d)/eta_d
    Var nu_p = the same with r1 <-> r2.

As a channel:
- g < 1: loss(g^2) followed by classical shifts of variance Var nu - (1 - g^2)/2.
- g > 1: a quantum-limited amplifier(g^2) followed by shifts of variance Var nu - (g^2 - 1)/2.

The added classical noise is never negative. For eta = eta_d = 1 it equals (g cosh r - sinh r)^2 for g < 1, and
(g sinh r - cosh r)^2 for g > 1. At g = sqrt(eta) it reduces to the library's channel; the numerical difference is 1e-11.

Checks against the full three-mode Fock simulation of the circuit (N = 25 per mode):
- Cases: g = 0.8, 0.9, 1.1 and 1.15, with and without loss, detector inefficiency, asymmetric squeezing and
  dephasing; two input states each.
- Trace distance: <= 7.7e-6 for every case.
- The residual is engine truncation. It falls 6e-5 -> 8e-6 -> 1e-6 for N = 20 -> 25 -> 30.
- The check is sensitive: using the unit-gain channel for the g = 0.9 circuit gives a trace distance of 0.086.

## Key numbers (key_numbers.txt / .json)

1 - F_avg (average infidelity of the logical qubit); "opt" = optimal recovery (SDP).

**E1, pure dephasing**

| sigma_phi | sq. GKP ideal | sq. GKP Petz [opt >=] | hex GKP Petz [opt >=] | cat parity | cat Petz | cat opt | Fock (projection = opt) |
|---|---|---|---|---|---|---|---|
| 0.1 | 0.018 | 0.0039 [0.0020] | 0.0038 [0.0019] | 0.024 | 2.6e-4 | 2.1e-4 | 0.0017 |
| 0.2 | 0.095 | 0.025 [0.013] | 0.031 [0.015] | 0.078 | 0.0041 | 0.0032 | 0.0066 |
| 0.3 | 0.171 | 0.060 [0.031] | 0.089 [0.046] | 0.138 | 0.025 | 0.018 | 0.015 |
| 0.4 | 0.229 | 0.108 [0.056] | 0.159 [0.085] | 0.192 | 0.071 | 0.051 | 0.026 |

- The cat with optimal recovery beats the bare qubit for sigma_phi < 0.27 rad, by up to ×14 at small sigma.
  With Petz the crossover is 0.24 rad.
- GKP beats the bare qubit only for very weak dephasing:
  - with Petz: below 0.060 rad (square) and 0.067 rad (hex);
  - with the ideal decoder: below 0.036 rad (hex). Square GKP with the ideal decoder never beats it on this grid.
- From 0.04 rad up to 0.40 rad, the cat with optimal recovery beats GKP even against GKP's best possible decoder (the lower bound in brackets).
  At 0.02 rad this is not certain: the GKP bound there is 3.2e-6 / 1.2e-6 (square / hex), against 4.7e-6 for the cat.
- The Fock qubit with projection is exactly optimal under pure dephasing: the certified SDP gap is 1e-10.

**E2, teleportation channel**

| dB | sq. GKP ideal | sq. GKP Petz | hex GKP ideal | hex GKP Petz | cat parity | cat Petz | cat opt | Fock proj. | Fock opt |
|---|---|---|---|---|---|---|---|---|---|
| 6 | 0.131 | 0.157 | 0.125 | 0.149 | 0.332 | 0.298 | (not run) | 0.202 | 0.173 |
| 10 | 0.026 | 0.023 | 0.022 | 0.019 | 0.228 | 0.192 | 0.150 | 0.099 | 0.085 |
| 16 | 0.0011 | 2.5e-4 | 6.3e-4 | 1.1e-4 | 0.084 | 0.057 | 0.049 | 0.028 | 0.024 |

- The cat with optimal recovery is worse than the bare qubit with optimal recovery at every dB tested (8-16), by at least ×1.6.
- GKP with the ideal decoder beats the bare qubit at every dB.
- Petz is worse than the ideal GKP decoder below 9 dB (both lattices), and even worse than the bare qubit below about 4 dB.
  Petz is not optimal.

**E3, gain mismatch, square GKP, 10 dB, eta = 1**

- Ideal decoder: 1 - F = 0.0263 at g = 1. The best is g = 0.99, with 0.0260.
  - g = 0.95 / 1.05: 0.034 / 0.043.
  - g = 0.90 / 1.10: 0.070 / 0.091.
  - The error grows faster for g > 1.
- Petz improves monotonically with g over the whole range: 0.049 at g = 0.85, 0.023 at g = 1, 0.016 at g = 1.15 (converged in N).
- Fock qubit with projection: the best is g = 0.92, with 0.083 (0.099 at g = 1).

## Interpretation (plain language)

- **Complementarity confirmed, at equal energy.**
  - Under shifts (E2), GKP is far better than the cat, and the cat is worse than doing nothing, even with the best possible decoder.
  - Under rotations (E1), the order of the two codes flips: at 0.2 rad, the cat with optimal recovery has 3.2e-3, while
    GKP cannot go below 1.3e-2 with any decoder.
- **The decoder matters as much as the code.**
  - The cat with its parity decoder is 25× worse than with optimal recovery at 0.2 rad: the parity decoder only undoes loss.
  - GKP's ideal decoder is about 4× worse than Petz under dephasing, because it assumes shifts, not rotations.
- **The bare qubit is hard to beat under pure dephasing.** It has n-bar = 0.5, and rotation damage grows with
  photon number. So at n-bar about 5, even the rotation-protected cat wins only below about 0.27 rad.
- **E3, why the ideal decoder is best at g slightly below 1.** A gain error scales the whole lattice, so the error
  grows with the distance from the origin, like a rotation. The best gain for the ideal decoder sits slightly
  below 1 because the added noise Var nu is smallest at g = tanh(2r) ≈ 0.98 (10 dB).
- **E3, why Petz improves for g > 1.** The added classical noise (g sinh r - cosh r)^2 vanishes at g = coth r ≈ 1.22.
  There the teleportation is exactly a quantum-limited amplifier. Similarly, at g = tanh r ≈ 0.82 it is a pure loss channel.

  A decoder that knows the gain can exploit this. The ideal GKP decoder, which assumes unit gain, cannot. So
  "gain mismatch is bad" is a statement about decoders that assume g = 1, not about the channel itself.

  This is my reading of the numbers. It is not a claim from a reference, and it should be checked before it
  goes in the report.

## Caveats

- **Petz and the SDP optimum are benchmarks that know the noise exactly; they are not realistic decoders.** A realistic cat
  decoder for rotations would need a phase measurement (Grimsmo et al. 2020). Neither Petz nor the SDP is a
  hardware protocol.
- **GKP's optimum is only bracketed, not computed.** It is somewhere in the shaded band; the band is wide, at most a factor 2 in F_e infidelity.
- **E1 is pure dephasing, with no finite-squeezing noise.** It is a different setting from the report's dephasing sweep (12 dB plus
  dephasing, crossover 0.137 rad). Do not compare their crossovers directly.
- **The cat has one size only.** Only alpha = 2.25 (energy matched) was tested. Larger cats resolve rotations better but suffer
  more from shifts.
- **E3 sweeps the gain only.** It assumes a constant (deterministic) gain error. Its eta = eta_d = 1, and r1 = r2 = 10 dB.
- **The SDP grid is coarser than the other sweeps.** E2: cat 8-16 dB, Fock 2-16 dB, even dB only. E1: full grid.

## Reproduce (~25 min in total on this machine)

```bash
cd scripts_ext
conda run -n qcsim python verify_gain_channel.py                   # ~1 min
conda run -n qcsim python verify_gain_convergence.py               # ~2 min
conda run -n qcsim python run_ext.py e1 e2 e3 conv                 # E1 1 min, E2 5 min, E3 14 min, conv 3 min
# cvxpy is NOT in the qcsim environment; install it outside it (does not modify qcsim):
conda run -n qcsim pip install --target /some/dir cvxpy clarabel   # used: cvxpy 1.9.3, clarabel 0.11.1
conda run -n qcsim env CVXPY_PATH=/some/dir python run_optimal.py  # ~5 min (resumes from optimal.json)
conda run -n qcsim python check_convergence.py
conda run -n qcsim python make_ext_figures.py > ../results_ext/key_numbers.txt
```

## Files

- `e1`, `e2`, `e3`, `conv`, `optimal` (`.csv` + `.json`): data rows with their configuration.
- `convergence_summary.json`, `verify_gain_channel.json`, `verify_gain_convergence.json`: the checks above.
- `key_numbers.json` / `.txt`: every number quoted here.
  - A `nan` crossing means the two curves do not cross on the grid.
  - A "Fock crosses Fock" line is a self-comparison; ignore it.
- `figures/complementarity.pdf` (recreated by `make_ext_figures.py`): (a) E2, shifts; (b) E1, rotations. Same palette and markers as the report.
  - Solid lines: standard decoder.
  - Dotted lines: Petz.
  - Stars: SDP optimum.
  - Shaded: GKP Barnum-Knill band.
- `figures/gain_mismatch.pdf` (recreated by `make_ext_figures.py`): E3.
