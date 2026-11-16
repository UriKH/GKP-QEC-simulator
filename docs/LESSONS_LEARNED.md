# Lessons learned: building a physics simulator and writing a report about it

Project: GKP error correction through CV teleportation (Technion "Modern Quantum Computing", 2026-11-02 to 2026-11-17).
One human (the student, who wrote the report) and one AI agent (Claude Code, which built the simulator, ran the
experiments and reviewed the text). This file collects what worked, what did not, and what to do differently. It is
written as input for a research harness, so each lesson ends with a rule that an agent or a workflow can apply.

Timeline for scale: the simulator (9 milestones, 131 tests) took about 2 days. The report took about 13 days and
15 versioned copies. Most of the effort, and most of the mistakes, were in turning results into a correct, defensible text.

---

## 1. Planning and structure

**1.1 Lock decisions early, in writing, and do not re-open them.**
The plan fixed the circuit, the codes, the error types and the decoders before any code was written ("locked
decisions" in HANDOFF.md). This prevented scope churn across many sessions.
Rule: start every project with a short plan containing locked decisions, milestones and verification checks. Ask the
human to approve it. Re-open a locked decision only when the human asks.

**1.2 Build in milestones that each end with a passing test.**
Conventions, then gates, then codes, then noise, then circuit, then engines, then decoders, then metrics, then sweeps.
Each milestone had physics tests with documented tolerances (squeezed variance, coherent-state teleportation
fidelity, stabilizer expectation e^{-pi Delta^2}). No later bug traced back to an early layer.
Rule: one layer per milestone, each with analytic checks, never "write the whole codebase then test".

**1.3 One source of truth for conventions.**
`cvqec/conventions.py` held hbar, the quadratures and the gate signs, copied from the course book with page numbers.
Sign and factor-of-sqrt2 errors are the classic CV bugs, and this project had none in the code.
Rule: put conventions in one module, cite where they come from, and test them first.

**1.4 Map every requirement of the assignment to a concrete artifact, early.**
Two instructor requirements were only matched at the very end. "Run simulations of the circuit" raised the question
whether the single-mode effective channel counts as the circuit. "A summary of how each code succeeds in correcting
the errors" needed a summary table that did not exist until the last round.
Rule: on day one, write a requirements table (requirement in the original wording, the planned artifact, its
status) and re-read it at every milestone. Quote the original wording and do not paraphrase it.

**1.5 Extensibility through registries paid off.**
New codes, decoders, metrics and experiment presets were each one registry entry. Late additions (joint Steane
decoding, a sensitivity study, a code-width sweep) took minutes, not refactors.

---

## 2. Physics and numerics

**2.1 Derive a cheap model, then verify it against the expensive one.**
The full three-mode Fock simulation is too large at the needed cutoff (N^3 states). We derived that the circuit acts
on Bob's mode as one channel (shift noise, then loss, then dephasing). We verified it against the full simulation at a
small cutoff (trace distance 8.5e-7) and ran all sweeps on the single-mode channel at N = 140.
Rule: when the exact simulation is too costly, derive a reduced model, verify it numerically against the exact one
at a size where both run, and state that check in the report.

**2.2 A small truncation tail does not prove convergence.**
Exponentiating a truncated generator gave wrong gates even when the top Fock level looked empty. The fix was analytic
displacements, plus convergence checks (N = 140 vs 170, and 200 vs 250 for narrow codes).
Rule: test convergence by changing the cutoff and comparing the final metric. Do not rely on tail population alone.

**2.3 Keep an analytic estimate next to every simulated curve.**
A one-line Gaussian model (each quadrature fails when its shift leaves the cell) matched the simulation to three
digits for pure shift errors, and was only approximate for loss and dephasing. It caught a wrong formula in the
report, predicted the best rectangular stretch, and told us where the simulation was really needed.
Rule: plot a closed-form estimate next to the numerics. Agreement validates both; disagreement must be explained.

**2.4 Bracket results with a reference and a benchmark.**
The bare qubit (no encoding) gave the break-even line. The Petz recovery gave a near-optimal upper benchmark, with a
proven bound (Barnum-Knill: error at most twice the optimum).
Rule: every comparison needs a "do nothing" baseline and, if possible, a provable best-case benchmark.

**2.5 Performance traps.** Three-operand `np.einsum` was 20x slower than batched matrix products. Profile before optimising.

---

## 3. Experiment design

**3.1 One imperfection at a time around a fixed operating point, then all at once.**
Each sweep varied one parameter with the rest ideal, so each effect was isolated. Then a combined run, an error
budget (add one imperfection at a time), and a sensitivity table (the combined run at three efficiencies), because
the combined operating point was only illustrative.
Rule: isolate, then combine, then test how much the conclusion depends on the chosen operating point.

**3.2 Choose the operating point for readability and say why.**
12 dB teleportation squeezing was chosen so that its noise is comparable to the code's own width (neither hides the
other) and the error is about 1%, inside what has been demonstrated experimentally. All codes had equal energy.
Rule: document why each default value was chosen. Fair comparisons need an equal resource (here energy).

**3.3 Store every result with its full configuration, and make figures and numbers a separate stage.**
Stage 1 (`run_sweeps.py`, about 12 min) wrote csv and json with all parameters. Stage 2 (`make_figures.py`, about
2 min) made every figure and a file of LaTeX macros. Re-plotting never required re-simulating.
Rule: separate "compute" from "present". Never hand-copy a number from a plot.

**3.4 Every number in the text is a generated macro.**
426 macros (each with a one-line meaning in `docs/results_numbers.md`) meant that after each rerun the report updated
itself, and stale numbers were impossible. It also forced every quoted claim to have a data source.
Cost: the macro list grows large, and a meaning file is needed to navigate it.
Rule: generate all quoted numbers from data. A number with no generator is a red flag.

**3.5 Late physics findings cause reruns. Budget for them.**
The hexagonal code underperformed with a realistic correction round. The reason was that it corrected q and p
separately. Joint decoding fixed it, and that became one of the main results, but it required a full rerun and a
new table. A removed parameter (syndrome-detector efficiency) caused another rerun.
Rule: expect at least one full rerun after the first results. Keep sweeps cheap enough (minutes, not hours) that
this is affordable, and always ask before starting a long run.

---

## 4. Validation (where most real errors were found)

A dedicated validation pass near the end found real problems. Doing it earlier, once per section, would have been better.

**4.1 Re-derive every formula in the text.** Found: a wrong failure-probability formula (wrong quadrature width in the
rectangular-code derivation).

**4.2 Check every qualitative claim against the data, not against memory.** Found:
- "The rectangular code behaves like the bare qubit" when it was slightly better everywhere.
- "Rectangular is best for biased noise" when it is best only for noise along its wide side, and worst for the opposite bias.
- "Exact" for a model that agrees to three digits. "Every fidelity" for a check done on one test input.
- A truncation claim that covered sweeps it was never measured on.

**4.3 Check that the text describes what the code does.** Found three implementation details the report did not state
(perfect readout after the correction round, how the bare qubit is read out, which test state was used).

**4.4 Check citations in the paper itself.** For each citation: download the paper, find the exact sentence or
equation, and record its location (section, equation, page). Cite only when the paper states the exact claim.
WebFetch could not read PDFs. Downloading the PDF and extracting text locally worked.

Rule for all of section 4: run a validation pass per section as soon as it is drafted, with four checks:
formulas, numbers and claims against data, text against code, citations against sources.

---

## 5. Writing with a human author

**5.1 The human writes, the agent proposes, and every proposal is visible and traceable.**
Edits were proposed as a new file version with each change coloured (`\chg{...}`), tagged (`% [C#]`), the original
kept as a comment, and a numbered change list in the reply. The human could accept or reject each one and defend the
text in person. This worked well and is worth keeping.
Costs: leftover colour after accepting; and two recurring LaTeX bugs from the markup itself. A `%` tag glued to a
word ate the following newline ("twothings"). A leading space inside `\chg{ ...}` was swallowed.
Rule: keep traceable proposals, but automate "accept round" (strip markup) and lint for these two patterns, or use
latexdiff instead of hand markup.

**5.2 Calibrate depth to the audience, and ask early.**
The agent repeatedly wrote more derivation than the audience (3rd/4th-year EE students) needed: a thermal-state
appendix, a Bloch-sphere fidelity derivation, a Kraus derivation from the beam splitter, a full Gaussian-integral
derivation. The human cut each one later.
Rule: before writing explanatory text, agree on the reader and the depth ("cite standard results, derive only what
is new"). Prefer a citation that states the exact result over a long derivation.

**5.3 Keep the author's voice.**
The human deliberately kept an informal tone and mentions of which course each concept came from, so the work would
read as their own. The agent should not "polish" this away.

**5.4 Fewer citations, each verified.** The human wanted the shortest reference list that is not plagiarism. Several
citations were dropped (software library, a historical first experiment, a classical bound that the text did not
need). Rule: cite what the text actually uses; drop "nice to mention" references.

---

## 6. Workflow and tooling: what was good and what was bad

**Good**
- A handoff document (HANDOFF.md) plus a per-phase log let work continue across many context resets.
- Persistent memory for user preferences (mentor style, course names, edit workflow).
- Two project skills (rebuild the report; review a section with tracked changes) encoded multi-step procedures.
- Compiling in a scratch copy, then rendering pages to PNG and looking at them. Visual checks caught problems that
  logs did not: swallowed spaces, bold that did not render, overfull lines.
- Archiving the file before each edit round. This made a damaged file easy to restore.

**Bad, and what to do instead**
- **Folder-per-version instead of git.** There were 15 report copies (about 75 MB), and the human also had an Overleaf
  copy that drifted ahead of the local files. It was often unclear which copy was the latest.
  Instead: use git from day one. Make one commit per review round and use branches for proposals. Sync Overleaf through git.
- **Stale context files.** The project instructions kept pointing at an old report folder long after it was replaced.
  They also stated a wrong fact about the user (that they had not learned the Lindblad equation). The skills still
  referenced an old folder. Instead: date-stamp facts, keep one "current state" line that is updated at every
  phase, and verify remembered facts with the human before relying on them.
- **The instruction file grew into a list of gotchas.** It works, but it is long. Instead: keep gotchas in a separate
  file that the instruction file links to, and prune it at each phase.
- **Unverified shell edits.** Two silent failures happened. A heredoc piped into `conda run` did nothing, because
  conda run does not forward stdin. A line-joining script pulled a paragraph into a comment.
  Instead: write edit scripts to a file, assert that each replacement matched exactly once, and diff against the
  archived copy after every edit.
- **Long rounds of small edits.** Many short review rounds each needed an accept, copy, compile cycle.
  Instead: batch related changes, and agree on "this is the last round" criteria.
- **Context loss in long sessions.** Summaries dropped details. The handoff document and the log were what saved
  continuity. Instead: write decisions to files as they are made, not at the end of a phase.

---

## 7. Checklist for the next project

1. Plan with locked decisions, milestones, verification checks, and a requirements table in the original wording.
2. Git from the first file. One commit per milestone or review round.
3. Conventions module, then layer-by-layer milestones, each with analytic tests.
4. A reduced model verified against the exact one. Convergence checked by changing the cutoff.
5. An analytic estimate, a baseline and a benchmark for every result.
6. Compute stage and presentation stage separated. All quoted numbers generated, each with a meaning.
7. Agree on the reader and the depth before writing. Cite rather than derive standard results.
8. A validation pass per section: formulas, data, code, citations.
9. Traceable edit proposals with automated accept and lint.
10. One handoff file and dated facts, updated at every phase. Never trust a stale pointer.
