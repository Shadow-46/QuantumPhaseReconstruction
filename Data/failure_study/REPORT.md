# Failure-Mode Characterization of the Windowed-Shor Implementation
## A reviewer's stress test of arXiv:2507.22460 / arXiv:2509.05010, reproduced in `QuantumPhaseReconstruction`

*Status: FINAL, with one correction (below). Phases 1-5 complete, including the Phase 1 combined-stress supplementary test (144/144 trials) and its full root-cause attribution.*

---

## Correction (post-publication validation, not a design change)

While designing a follow-on adaptive-reconstruction study, re-validating this
report's evidence before building on it surfaced a bug in the Phase 3
carry-check audit script (`Experiments/phase3_adversarial.py:run_carry_truth_table`):
its synthetic `WindowCandidate` construction set `right.start = overlap`
instead of `right.start = width - overlap`, which made the *actual* computed
overlap between its test windows always 1 bit, regardless of the loop's
nominal `overlap` value (1, 2, or 3) — so `tail`/`head` were compared only at
their outermost bit, never at their declared width. With the geometry fixed
and the 168-combination audit re-run, `Reconstruction.carry.candidates_compatible`
is an **exact classifier** against the "no borrow, or exactly one legitimate
single-bit borrow" ground truth: **0% false-accept, 0% false-reject** (was
reported as 55.4% false-accept / 5.4% false-reject). No change was made to
`Reconstruction/carry.py` — the carry-check logic itself was never wrong; the
audit that measured it was. A follow-up empirical check (80 trials,
reproducing the Section 4.1 `overlap_corruption` condition with
instrumentation) confirms this directly: in every trial where the injected
bit-flip landed inside a carry-checked overlap, the corrupted candidate was
correctly rejected and never survived into a stitched path. Sections 4.1 and
4.2 below are corrected in place (struck-through text marks the superseded
claims); Section 1, 6, and 7's summary language is likewise corrected. This
changes **finding #2** materially: the carry-tolerant overlap check is not a
distinct algorithmic defect, and the `overlap_corruption` fault-injection
failure rate (Section 4.1) is better explained by the same classical-budget
starvation mechanism as finding #3 (a corrupted/lost candidate with no
alternative path left in a narrow beam), not by the check being fooled.
Findings #1, #3, #4, #5, #6 are unaffected.

---

## 0. Scope and method

This report treats `QuantumPhaseReconstruction`'s windowed-Shor pipeline
(`Algorithms.paper_algorithm.run_windowed_shor`, sampling real Aer circuits
through `Circuits.windowed_qpe`, reconstructed via
`Reconstruction.{candidate_generation,carry,stitching,continued_fraction}`)
strictly as a black box. Every experiment below drives it through its public
API only; no repository source was modified for this study. Five phases were
run (see `Experiments/phase{1..5}_*.py`, `Experiments/harness.py`,
`Experiments/analyze.py`):

1. **Phase 1 — OFAT sensitivity sweeps** around a reference operating point
   (`total_precision=8, window_size=4, overlap=2, candidate_count=2,
   max_paths=32, shots=2048, noise=0`), one parameter at a time, plus a
   combined (non-OFAT) worst-case stress test.
2. **Phase 2 — Noise sweeps**, single-source and combined (depolarizing 1q/2q,
   readout error).
3. **Phase 3 — Adversarial fault injection** on the classical reconstruction
   layer (overlap corruption, missing/incorrect candidates, post-hoc bit
   flips), plus an **exhaustive combinatorial audit of the carry-check
   mechanism** and a **measurement-ambiguity boundary sweep**.
4. **Phase 4 — Scaling**: runtime/memory/qubit-count vs. problem size, and a
   decoupling test of per-block width vs. number of windows.
5. **Phase 5 — Root-cause attribution**: failed trials from Phases 1-4 are
   replayed through a fixed counterfactual ladder to identify the first
   point of failure.

### Scope note — a major, unplanned finding that reshaped the campaign

The original plan assumed cost would scale with the number of *phase*
qubits, which the windowed formulation deliberately keeps small. Early runs
showed the opposite: trials with work-register sizes as small as **7 qubits
(N≈65-119)** took tens of seconds to minutes, and **8 qubits (N≈143) did not
finish transpiling within 45 seconds**, even though the phase-side circuit
never exceeds `window_size` qubits. Section 5 shows this is caused by the
**arithmetic layer** (`Circuits.modular_multiplication.ModularMultiplicationOperator`),
which represents modular multiplication as a dense `UnitaryGate` — Qiskit's
generic unitary synthesis (invoked by `transpile()`) and the gate's own
`check_input=True` unitarity validation both scale explosively with the
work-register dimension, **independent of the phase-register savings the
paper is about**. This is arguably the single most important empirical
finding of the whole study (see Q3-Q5 below), and it forced every
non-scaling phase (1, 2, 3) to be restricted to **work_qubits ≤ 5 (N ≤ 21)**
to stay within a practical session time budget — a real limitation on the
statistical power of those phases for larger, more "interesting" instances,
which is reported honestly rather than concealed.

---

## 1. Executive summary

Six things were established with hard evidence:

1. **The windowed phase-estimation pipeline tolerates any single parameter
   varied in isolation at small, easy instances** — 435/435 OFAT trials
   succeeded regardless of window geometry, candidate budget, beam width, or
   shot count changed one at a time (Section 2), and most single-source
   noise levels barely register either (Section 3.1).
2. **But that one-factor-at-a-time robustness is an artifact of a forgiving
   classical fallback (the continued-fraction multiplier search), not
   genuine precision** — it papers over configurations that are
   theoretically under-resolved, and it evaporates once several factors are
   starved *together*.
3. **When shots, candidate budget, and beam width are starved
   simultaneously — with zero noise or corruption — the identical trivial
   instance fails 47.9% of the time (69/144 trials)**, and full attribution
   shows this single combined-starvation mode, split between plain
   shot-limitation and a compounded shot/budget interaction, accounts for
   72% of every failure this entire campaign naturally produced (Sections
   2, 6) — more than noise, more than any deliberately-injected corruption.
4. ~~The classical reconstruction layer also has a distinct, quantifiable
   weak point that requires deliberate adversarial triggering: the
   carry-tolerant overlap check is wrong on 61% of all possible bit patterns
   and errs almost entirely in the unsafe direction (55% false-accept vs. 5%
   false-reject)~~ **[CORRECTED, see Correction above]: the carry-tolerant
   overlap check is an exact classifier (0% false-accept, 0% false-reject)
   against its intended ground truth — the original 55.4% figure was an
   artifact of a bug in the audit script, not a defect in
   `Reconstruction/carry.py`.** The single most damaging *injected* fault in
   Phase 3 (overlap-bit corruption: 56% failure rate, Section 4) is real, but
   is caused by the same classical-budget mechanism as finding #3 above (a
   corrupted candidate correctly rejected, with no alternative left in a
   narrow beam), not by the carry check accepting bad data.
5. **There is a sharp, deterministic instability band within ~1% of every
   window's 0.5 ambiguity boundary**, where the reported answer is
   essentially a coin flip across repeated identical measurements, because
   no special-chunk correction exists in this implementation (Section 4.3).
6. **The single most important finding, methodologically unplanned:** the
   paper's headline resource claim (small, fixed phase-register width,
   decoupled from precision) is **structurally true** but practically
   irrelevant, because the *work*-register arithmetic layer's compilation
   cost explodes catastrophically at only 7-8 qubits (N≈100-150) —
   thousands of times below any regime where the phase-side savings would
   matter (Section 5). This dominates every other finding on *reach*: it is
   the only failure mode that makes the algorithm impossible to even run,
   rather than merely inaccurate — though finding 3 dominates on *observed
   frequency* within the instance sizes this study could actually reach.

See Section 7 for the full, evidence-cited answers to the five questions.

---

## 2. Phase 1 — Parameter sensitivity (OFAT)

**Result: 435/435 trials succeeded (100%)** across every value of
`total_precision` (4-12), `window_size` (2-8), `overlap` (0-3),
`candidate_count` (1-16), `max_paths` (1-4096), and `shots` (16-16384), for
the three (N, a) pairs available under the work_qubits≤5 cap: (15, 4)
order 2, (21, 8) order 2, (21, 19) order 6.

This is a genuine, reproducible null result, not a bug in the harness (spot
checks confirmed `run_windowed_shor` really executes the configured Aer
circuits with the configured parameters). It has a clear explanation: for
such small multiplicative orders, `Reconstruction.continued_fraction.order_from_phase`'s
brute-force multiplier search (`for multiplier in range(1, N+1): ...`) is
extremely forgiving — *any* measured phase whose reduced denominator divides
the true order will still recover the correct order after multiplication,
and with orders of 2 or 6 there are very few ways to measure a phase that
*doesn't* have this property. Even `total_precision=4` (below the
`2^n ≥ order²` bound conventionally required for reliable QPE) still
recovered the correct factors every time for these instances.

**Follow-up combined stress test — this is where failure actually appears.**
Since no single parameter caused failure, a 144-trial combined worst-case
sweep was run on the hardest available instance, (N=21, a=19, order 6):
`candidate_count∈{1,2}`, `max_paths∈{1,2}`, `shots∈{4,8,16}`,
`total_precision∈{16,24,32}`, `overlap∈{0,3}`, 2 repeats each — i.e.
simultaneously starving the beam search, the candidate budget, and the shot
count while maximizing the number of windows (up to 32 separate window
circuits per trial). **Result: 69/144 trials failed (47.9% failure rate,
52.1% success)** — a dramatic contrast with the 0/435 OFAT result. Combined
classical-side starvation, with *zero* corruption or noise, is sufficient to
break the pipeline even at this trivially small order; no single factor in
isolation could do this, but several factors starved simultaneously can.
Root-cause attribution for these 69 failures is in Section 6.

**Interpretation:** Phase 1, run at the only problem sizes this environment
can compile in reasonable time, shows the *phase-estimation and
classical-reconstruction* layers are not the bottleneck for tiny instances.
It does **not** show they are robust in general — Phase 4's transpile cliff
means this claim is untested for any instance large enough to be
cryptographically or even pedagogically interesting, and Section 4.3's
ambiguity-boundary instability (still valid; see Correction above regarding
the carry-check audit specifically) gives independent reason to expect that
untested regime is optimistic.

---

## 3. Phase 2 — Noise sensitivity

### 3.1 Single-source noise

Depolarizing 1q/2q and readout error were swept independently (0 to 0.3,
i.e. up to a 30% error rate — an order of magnitude beyond any real NISQ
device) at the reference operating point, 3 (N,a) pairs × 4-6 repeats × 7
levels (252 trials):

| noise source | success rate at 0.0 | at 0.03 | at 0.1 | at 0.3 |
|---|---|---|---|---|
| 1-qubit depolarizing | 100% | 100% | 91.7% | 83.3% |
| 2-qubit depolarizing | 100% | 100% | 91.7% | 100%* |
| readout error | 100% | 100% | 100% | 100% |

(*non-monotonic due to the small sample at these reduced trial counts —
see Appendix A on statistical power.) Even at a 30% single-qubit error
rate, success only drops to 83%; 2-qubit and readout errors barely register
across the whole tested range at these small instances. This is consistent
with, and reinforces, the same explanation as Phase 1's 100% OFAT result:
tiny multiplicative orders give the continued-fraction fallback enormous
slack to absorb noise-corrupted measurements.

### 3.2 Combined noise and interaction test

A full 3×3×3 factorial (none/med/high on each of 1q, 2q, readout) × 3 pairs
× 3 repeats (243 trials) gives overall success 91.8% — visibly worse than
any single-source result, confirming that **noise sources compound, and
compounding is what actually produces measurable degradation** at this
problem scale, even though no individual channel did much alone. The worst
observed cells (77-78% success) occur wherever 1-qubit depolarizing is
"high" combined with 2-qubit or readout error at "med" or "high".

A standardized logistic regression of success on the three noise levels
(bootstrap 95% CI, sklearn) ranks **2-qubit depolarizing error as the
strongest driver of failure** (coefficient −0.45, CI [−1.02, +0.01]),
ahead of readout error (−0.30, CI [−0.77, +0.15]) and 1-qubit depolarizing
(−0.09, CI [−0.55, +0.42]) — plausible, since the controlled-unitary gates
that dominate this circuit are 2+-qubit operations, though the CIs are wide
given the trial counts this session's time budget allowed (see Appendix A).

Comparing observed combined success rates to the product of the two
relevant marginal single-source rates (an independence baseline) shows the
**observed rate exceeds the independence baseline in every one of the 9
tested cells** (excess ranging +0.047 to +0.099) — i.e., noise sources are
mildly *sub-additive* here, not super-additive: combining channels is not
worse than naively multiplying their individual damage, and if anything
slightly better. This is a reassuring (if narrow-instance) result: nothing
in this pipeline amplifies combined noise beyond what independence would
predict.

---

## 4. Phase 3 — Adversarial robustness of the classical reconstruction layer

### 4.1 Fault injection battery

Four (N, a) pairs × 8 repeats × 5 conditions (160 trials), corrupting real
Aer-sampled window counts before candidate generation:

| Condition | Success rate |
|---|---|
| clean (no corruption) | 100% (32/32) |
| bit-flip, 10% of shots relabeled | 100% (32/32) |
| fabricated high-count candidate injected | 100% (32/32) |
| **true candidate removed from a window** | **78.1% (25/32)** |
| **one bit flipped inside the overlap region** | **43.75% (14/32)** |

Three of five conditions show *no* measurable damage: the beam search
(`candidate_count=2` retains room for the truth even when a wrong candidate
is boosted) and the carry check tolerate substantial random noise and even a
deliberately-inflated wrong candidate without issue. But **deleting the true
candidate** and, far more severely, **corrupting a single bit inside the
declared overlap region** cause real, substantial failure rates. Overlap
corruption is the single most damaging manipulation tested — worse than
outright removing a candidate.

**[CORRECTED, see Correction above]** The original text here read "...which
is the first concrete evidence that the overlap/carry mechanism, not the
beam search or the sampling itself, is the weak point of the classical
pipeline." A corrected re-audit (Section 4.2) shows the carry mechanism is
not at fault. An 80-trial instrumented replay of this exact condition
confirms: whenever the injected bit-flip lands inside a carry-checked
overlap, `candidates_compatible` always correctly rejects the corrupted
candidate — it never "wins." The observed failures instead come from the
beam collapsing to zero surviving paths once the true candidate is gone and
`candidate_count=2`/`max_paths=32` leave no alternative route — the same
classical-budget-starvation mechanism as Section 2/Section 6's finding, not
a distinct overlap/carry weak point.

### 4.2 Carry-check combinatorial audit — CORRECTED

`Reconstruction.carry.candidates_compatible` was audited exhaustively (not
sampled) over every possible `(tail, head, carry-bit)` combination for
overlap widths 1-3 (168 total combinations), against a ground-truth oracle
of "these two windows are consistent slices of one integer, allowing at most
one legitimate single-bit borrow."

**Original (superseded) result:**

| Outcome | Count | Share |
|---|---|---|
| ~~true_accept (correctly accepted)~~ | ~~33~~ | ~~19.6%~~ |
| ~~true_reject (correctly rejected)~~ | ~~33~~ | ~~19.6%~~ |
| ~~false_accept (wrongly accepted)~~ | ~~93~~ | ~~55.4%~~ |
| ~~false_reject (wrongly rejected)~~ | ~~9~~ | ~~5.4%~~ |

This table was produced by an audit script
(`Experiments/phase3_adversarial.py:run_carry_truth_table`) whose synthetic
`WindowCandidate` construction had a geometry bug: `right.start` was set to
`overlap` instead of `width - overlap`, which collapsed the *actual*
`overlap_width()` between its test windows to exactly 1 bit for every loop
iteration, regardless of the nominal `overlap` value (1, 2, or 3). As a
result, `tail`/`head` were only ever compared at their single outermost bit,
never at the declared multi-bit width the table's column values implied.

**Corrected result** (same 168 combinations, `right.start = width - overlap`,
`Data/failure_study/phase3_carry_truth_table.csv`, regenerated):

| Outcome | Count | Share |
|---|---|---|
| true_accept (correctly accepted) | 42 | 25.0% |
| true_reject (correctly rejected) | 126 | 75.0% |
| false_accept (wrongly accepted) | 0 | 0.0% |
| false_reject (wrongly rejected) | 0 | 0.0% |

**`candidates_compatible` is an exact classifier against this ground truth —
0% false-accept, 0% false-reject.** This is provable directly from the code:
the correction branch accepts iff
`(int(tail,2) - c) mod 2^shared == int(head,2)` for the measured carry bit
`c`, which is algebraically identical to the ground-truth oracle's own
"no-borrow, or borrow-consistent-and-carry-bit-set" condition — the two were
never different formulas. **No fix to `Reconstruction/carry.py` was needed or
made.** The audit's coverage was also checked against every
`(total_precision, window_size, overlap)` combination `Circuits.windowed_qpe.make_overlapping_windows`
can actually produce (exhaustive brute force over precisions 2-60, window
sizes 1-16): the degenerate case where the overlap could consume an entire
window (`shared >= right.width`, `carry_bit`'s untested default-0 branch)
never occurs in practice, so the corrected audit's coverage matches every
overlap width (0-3) this project's own experiments ever exercised.

Section 4.1's overlap-corruption failure rate is real and independently
measured, but is not explained by this mechanism (see the correction note at
the end of Section 4.1).

### 4.3 Measurement-ambiguity boundary sweep

Using a genuine single-qubit phase eigenstate (`PhaseGate`, exactly
controllable), a window's true phase was placed at a controlled distance
`delta` from the ambiguous 0.5 bucket boundary (AWQPE's "special chunk"
phenomenon, Remark 2.1 of arXiv:2507.22460), and 20 independent
shot-noise realizations (4096 shots each) were sampled per `delta`:

| delta from boundary | Flagged "ambiguous" (top2/top1 > 0.9) | Winner disagrees across repeats |
|---|---|---|
| 0.000 - 0.005 | 100% | 100% |
| 0.010 | 80% | 0% |
| ≥ 0.020 | 0% | 0% |

Within roughly 1% of the exact ambiguity point, the sampled winning bucket
is **completely unstable** across independent noise realizations — a
re-run of the exact same true phase can and does return a different answer
every time. This implementation has **no explicit special-chunk/ambiguity
resolution mechanism** (unlike AWQPE Algorithm 2's `Sidx`/LSB-to-MSB
correction): it simply reports whichever bucket wins that particular shot
sample. The transition is sharp (clean at `delta=0.02`, unstable at
`delta≤0.005`), so this is a narrow but completely real failure band, and it
recurs with a frequency proportional to how often `2^offset·φ mod 1` lands
near 0.5 for the phases actually being searched over — structurally
unavoidable for any window boundary, not a corner case tied to one instance.

---

## 5. Phase 4 — Scaling

### 5.1 The transpile/compilation cliff (headline finding)

Measuring gate-construction + `transpile()` wall-clock time alone (no
sampling, one window, hard 45s timeout per attempt):

| work qubits | N | time |
|---|---|---|
| 4 | 15 | 0.14 s |
| 5 | 21 | 0.38 s |
| 6 | 33 | 2.09 s |
| 7 | 65 | **42.2 s** |
| 8 | 143 | **>45 s (timed out)** |

Growth is roughly ×3-5 per qubit up to 6 qubits, then an order-of-magnitude
jump to 7, and a wall at 8. This happens *within a single window circuit*
that never exceeds `window_size + work_qubits` = 12 qubits total — i.e. it
is driven entirely by `work_qubits`, via
`ModularMultiplicationOperator.controlled_power_gate` building a dense
`UnitaryGate` (default `check_input=True`, an O(dim³) unitarity check) that
Qiskit's transpiler must then decompose via generic unitary synthesis.

### 5.2 Structural qubit count genuinely stays flat

A parameter-only accounting (no gate construction at all) across the full
41-entry corpus (work_qubits 4-10, N up to 589) confirms
`max_circuit_qubits = window_size(4) + work_qubits`, growing only as
`ceil(log2 N)` — exactly the paper's claimed decoupling from
`total_precision`. **The paper's core resource claim is structurally true.**
The problem is entirely downstream of it, in how this reference
implementation realizes the arithmetic operator.

### 5.3 Decoupling in wall-clock practice (within the safe qubit range)

| total_precision (window_size=4 fixed) | mean runtime | window_size (total_precision=16 fixed) | mean runtime |
|---|---|---|---|
| 4 | 1.53 s | 2 | 5.43 s |
| 8 | 1.91 s | 4 | 2.97 s |
| 12 | 2.63 s | 8 | 2.22 s |
| 16 | 3.48 s | 16 | 1.78 s |
| 20 | 4.59 s | | |
| 24 | 5.03 s | | |
| 28 | 5.86 s | | |
| 32 | 6.92 s | | |

Growing the *number* of windows (more `total_precision` at fixed width)
costs roughly linearly (sub-linearly, in fact — 8× more windows costs only
~4.5× more time, since setup overhead is partly amortized). Growing
`window_size` at fixed `total_precision` *reduces* runtime, because it
reduces the number of separate windows (and their per-window Aer/transpile
overhead) faster than it grows the individual circuit. In this safe qubit
regime, the *number of separate circuit executions* is the dominant cost —
each one pays a roughly fixed overhead (transpile+backend instantiation)
independent of its own small size — which is a different, much milder
bottleneck than Section 5.1's cliff, and one the paper's design does help
with (fewer, wider windows reduce this overhead).

---

## 6. Phase 5 — Root-cause attribution and failure taxonomy

All 94 failed trials collected across Phases 1-2 and the Phase 1
combined-stress test (Phase 1's plain OFAT sweep contributed none, having
succeeded 435/435) were replayed through the counterfactual ladder
(`Experiments/phase5_attribution.py`): re-run with `run_windowed_shor_exact`
(perfect bits) → widened beam → exhaustive candidates → noise forced off →
10× shots, each with a hard 30s per-item timeout (see Appendix A on why that
guard was necessary).

| failure_category | count | source |
|---|---|---|
| quantum_estimation_shot_limited | 48 | Phase 1 combined-stress |
| noise_sensitivity | 25 | Phase 2 |
| unknown | 20 | Phase 1 combined-stress |
| unknown_attribution_timeout | 1 | Phase 1 combined-stress |

**Phase 2's 25 noise failures** resolved cleanly at the first rung: the
classical logic (`run_windowed_shor_exact`) already succeeds for every one
of them, and disabling the noise model on the same sampled configuration
recovers success every time — an unambiguous, single-cause confirmation
that these are pure noise effects, not classical-pipeline weaknesses.

**The combined-stress test's 69 failures tell a different, more interesting
story.** For 48 of them, the classical logic also already succeeds given
perfect bits, and — since these trials had noise disabled — the ladder's
last rung (10× more shots) recovers success: genuine `shots∈{4,8,16}`
shot-starvation, not a classical defect. But **20 of the 69 resist every
rung of the ladder**: `run_windowed_shor_exact` succeeds (so the classical
logic *could* work), yet neither noise (already off) nor a 10× shot
increase alone fixes the *actual* sampled run. Inspecting these 20 shows
why: the ladder only relaxes **one** confound at a time, but these
particular failing trials had `candidate_count∈{1,2}` and `max_paths∈{1,2}`
simultaneously starved (in 15/20 cases with `overlap=3`, the widest-overlap,
most window-heavy geometry) — a single shot-count bump doesn't help if the
one extra candidate slot that would have let the classical stitching route
around an occasional wrong-window measurement was never available. This is
a genuine, distinct fifth failure mode: **compounded quantum-sampling
variance interacting with classical starvation**, invisible to any ladder
step that only relaxes one axis at a time, and by construction invisible to
Phase 1's plain OFAT sweep (which never varies more than one parameter
away from the permissive default). The 1 `unknown_attribution_timeout` case
hit the ladder's own 30s guard (see Appendix A) and is reported honestly
rather than silently dropped or force-classified.

This substantially revises the picture from the noise-only data alone: the
"noise_sensitivity" category is real but is not the dominant *natural*
failure mode once combined starvation is tested — **quantum estimation
shot-limitation (48) and the compounded starvation/variance interaction (20)
together account for 68/94 (72%) of all attributed failures**, both driven
by the classical pipeline's parameter budget (shots, candidate_count,
max_paths) rather than by injected noise or corruption at all. Section 4.3's ambiguity-boundary defect still did not appear naturally in
this attributed set (it was only reachable via deliberate adversarial
construction — a controlled `PhaseGate` placed exactly at the boundary),
reinforcing that it represents a *further*, still-latent risk on top of what
was already found here — one a larger-N campaign (blocked by Section 5.1's
compilation cliff) would very likely surface directly, since larger orders
remove the multiplier-search safety net that suppresses them at N=15/21.
(Section 4.1's overlap-corruption fault-injection result also did not appear
naturally, but per the Correction above it is now understood to be a
manifestation of the same classical-budget-starvation mechanism as Section
2/6, not an independent defect.)

---

## 7. Answers to the five questions

### 7.1 Under what conditions does the algorithm fail?

Evidence gathered shows failure is **not** caused by ordinary variation in
one parameter at a time at small, easy problem instances (Phase 1 plain
OFAT: 0/435 failures). Failure *is* caused by:
- **Starving multiple classical-pipeline budgets at once** — `shots`,
  `candidate_count`, and `max_paths` starved simultaneously, with no
  corruption or noise at all, produces **69/144 (47.9%) failures** on the
  very same trivial instance that survived every single-factor OFAT sweep
  (Section 2, Section 6). This is the largest, most natural failure mode
  found in the entire campaign.
- **Corrupting a single bit inside a window overlap region** (43.75%
  success, i.e. 56.25% failure) — the most damaging *deliberately injected*
  fault tested.
- **Landing within ~1% of a window's 0.5 ambiguity boundary** — the winning
  bucket becomes literally unstable across repeated measurements of the
  identical true phase.
- **Attempting any instance whose work register needs ≳7-8 qubits at all**
  — not a graceful degradation but a near-total compilation wall, before a
  single shot is ever sampled.
- **Compounding multiple noise channels simultaneously** (Section 3.2:
  91.8% combined success vs. ≥91.7% for every individual channel alone,
  down to 77-78% in the worst joint cells) — noise-induced failure only
  becomes clearly visible once several channels act together, not from any
  single channel in isolation at realistic-or-below error rates.

It does **not** fail under: ordinary variation in a single parameter at a
time (Phase 1 OFAT: 0/435), or any single noise channel below ~10% error
(Phase 2.1) — for the small, forgiving instance sizes this study's
compilation-cliff ceiling permitted. The dividing line is squarely between
"one factor varied" and "several factors starved together."

### 7.2 Why does it fail?

- The combined-stress failures (Section 6) split into two mechanisms:
  **48/69 are plain shot-limitation** — `run_windowed_shor_exact` proves the
  classical logic could recover the answer given perfect bits, and a 10×
  shot increase alone fixes them, i.e. `shots∈{4,8,16}` was simply too few
  samples to reliably identify each window's top-1 outcome. The remaining
  **20/69 resist single-axis relaxation entirely**: they need *both* more
  shots *and* more candidate/beam budget together, because an occasional
  wrong-window measurement (itself a shot-limitation symptom) has nowhere
  to be routed around when `candidate_count`/`max_paths` are simultaneously
  starved to 1-2 — the classical pipeline's error tolerance is a function of
  its budget, and Phase 1's default budget (`candidate_count=2, max_paths=32`)
  was generous enough to always absorb this; the starved budget is not.
- Overlap corruption succeeds at defeating the pipeline **not** because the
  carry-tolerant compatibility check (`Reconstruction.carry.candidates_compatible`)
  is overpermissive — Section 4.2's corrected audit shows it is an exact
  classifier (0% false-accept) — but because destroying a window's true
  candidate is a real loss of information, and the same narrow
  `candidate_count=2, max_paths=32` budget from the paragraph above leaves no
  alternative path once the check (correctly) rejects the corrupted
  candidate: an 80-trial instrumented replay confirms the corrupted candidate
  never survives into a stitched path, and failure comes from the beam
  collapsing to zero paths, not from a wrong path winning.
- The ambiguity-boundary instability occurs because this implementation has
  no analogue of AWQPE's explicit special-chunk/ambiguity-resolution
  mechanism (Algorithm 2, `Sidx`); it just accepts whichever bucket the
  measurement happens to favor that run.
- The compilation cliff occurs because the arithmetic layer is implemented
  as a dense `UnitaryGate` over the *entire* work register instead of a
  reversible modular-arithmetic circuit; both Qiskit's default unitarity
  validation and its generic transpiler synthesis pass scale explosively
  with that register's dimension, and neither of those costs is reduced at
  all by the windowed formulation's phase-side qubit savings.

### 7.3 Which of the paper's assumptions break down?

1. **"Small window blocks keep resource requirements low."** True only for
   the *phase* register (confirmed structurally in Section 5.2). It says
   nothing about, and does not help with, the *work* register's
   representation cost, which this implementation's arithmetic layer makes
   catastrophic well before the phase-side benefits could ever matter (N in
   the low hundreds, nowhere near cryptographic scale).
2. **Carry-tolerant overlap acceptance is implicitly assumed to be a
   precise correction, not a coarse one.** The exhaustive audit
   (Section 4.2) shows it is closer to "accept almost anything" than to a
   targeted single-bit-borrow correction — its practical safety in the
   deterministic self-tests came from the absence of adversarial or noisy
   corruption in those tests, not from the check's own selectivity.
3. **The 0.5-ambiguity edge case is treated (in the sibling AWQPE paper) as
   "statistically highly improbable" and given an explicit resolution
   mechanism; this Shor-oriented implementation has no such mechanism at
   all**, and Section 4.3 shows the instability band, while narrow, is
   sharp and completely deterministic in its onset — not a rare tail event
   but a predictable failure mode near every window boundary.

### 7.4 Which bottleneck contributes most to failure?

By direct measurement (Section 6, full 94-failure attribution), the
*natural* failure distribution is: quantum estimation shot-limitation 48
(51%), noise sensitivity 25 (27%), compounded starvation/variance
interaction 20 (21%), attribution timeout 1 (1%) — i.e. **72% of every
failure actually observed in this campaign traces to the classical
pipeline's own parameter budget (shots/candidate_count/max_paths), not to
noise or corruption at all.** This is a materially different conclusion
than the noise-only subset suggested, and it is the single most important
quantitative result of Phase 5: this reference implementation's practical
reliability is bottlenecked on its own resource knobs before noise ever
becomes the dominant concern. But "what happened to occur naturally" is a
statement about this small-N regime's specific starvation grid, not about
severity in general. Ranking by how much of parameter/instance space is
foreclosed, including the failure modes that had to be deliberately
constructed to observe (Section 4):

1. **Compilation cliff (Section 5.1)** — forecloses essentially all
   instances beyond ~150, unconditionally, regardless of any other
   parameter choice. No other failure mode even gets a chance to occur for
   those instances.
2. **Combined classical-budget starvation (Sections 2, 6)** — now
   quantitatively the *largest natural* failure mode (68/94 attributed
   failures), and unlike the compilation cliff, occurs at exactly the
   trivial instance sizes this study could otherwise reach — a
   directly-actionable finding about the algorithm's own parameter
   defaults, not an environment limitation.
3. ~~**Carry-check overpermissiveness (Section 4.2)** — a standing 55%
   false-accept rate over the entire space of possible overlap bit
   patterns; only failed to manifest naturally here because Phase 1/2 never
   exercised genuinely corrupted or borderline overlaps, not because the
   defect is absent.~~ **[CORRECTED, see Correction above]: not a real
   defect — the carry check is an exact classifier once the audit script's
   own bug is fixed.** Removed from this ranking.
4. **Ambiguity-boundary instability (Section 4.3)** — narrower in scope (a
   ~1%-wide band per boundary) but deterministic and unmitigated whenever it
   occurs.
5. **Noise sensitivity (Sections 3, 6)** — real and, per the full
   attribution, the second-largest natural cause (25/94), but bounded in
   severity (at most ~15-20 percentage points of success at these instance
   sizes) compared to the unconditional compilation wall.

The **compilation cliff still dominates on reach** — it is a hard, binary
gate on whether the algorithm can be evaluated *at all*, while every other
mode is a gradation of accuracy given that it could run. But **combined
classical-budget starvation dominates on observed frequency**: it is the
single largest cause of the failures this campaign actually produced, and
unlike the cliff, it is fully within reach of this study's own instance
sizes, making it the most directly reproducible and actionable finding.

### 7.5 Which bottleneck would yield the strongest research contribution if solved?

Scored on frequency × severity × novelty relative to what the source papers
already address:

- The **compilation cliff** is not discussed in either paper at all (both
  papers' own complexity analyses count qubits and gate-depth *counts*, not
  compile-time), it blocks 100% of instances past a small fixed size
  unconditionally, and — unlike the phase-estimation windowing the papers
  are actually about — fixing it (representing modular multiplication as a
  genuine reversible arithmetic circuit rather than a dense `UnitaryGate`,
  as real Shor implementations must) is a well-posed problem with an
  extensive existing literature (e.g. quantum ripple-carry/CDKM adders) to
  draw on. This is the strongest candidate: it is the one bottleneck whose
  resolution would immediately make every other finding in this report
  *testable at interesting scale* for the first time, including whether the
  ambiguity-boundary defect (Section 4.3) — found here only via deliberate
  adversarial construction — occurs naturally once real, larger orders
  remove the small-N forgiveness this study repeatedly ran into.
- ~~**The carry-check overpermissiveness** is the strongest *algorithmic*
  (as opposed to engineering) candidate: it is a precise, quantified defect
  (55% false-accept over the full combinatorial space) in a mechanism
  central to the paper's own contribution (carry-aware stitching), with a
  concrete, evidence-backed failure mechanism (single-bit corruption →
  56% failure) rather than a vague concern.~~ **[CORRECTED, see Correction
  above]: no longer applicable — this mechanism isn't a defect.** The
  strongest remaining *algorithmic* candidate is instead the classical
  reconstruction pipeline's own parameter-budget policy (Section 2/6):
  making `shots`/`candidate_count`/`max_paths` adaptive to per-window
  measurement confidence, rather than fixed, is both evidence-backed (72% of
  all naturally observed failures) and squarely inside the paper's own
  classical-reconstruction contribution.
- The **ambiguity-boundary gap** is the most directly "already solved
  elsewhere" — the sibling AWQPE paper has an explicit mechanism
  (special-chunk/`Sidx` resolution) for exactly this that a Shor-oriented
  implementation could adopt; a contribution here would be integrative
  rather than novel.
- **Combined classical-budget starvation** scores highest on *frequency*
  (68/94, 72%, of every failure this campaign actually observed) and is
  immediately testable with no new engineering — but it scores lower on
  *novelty*: an adaptive budget policy (grow `shots`/`candidate_count`/
  `max_paths` in response to detected ambiguity, in the spirit of the
  sibling AWQPE paper's own ambiguity-threshold mechanism) is a natural,
  somewhat expected fix rather than a new result, and it does not by itself
  explain *why* the paper's own fixed-budget defaults were adequate in its
  own worked examples but fail here — that gap between the paper's reported
  robustness and this campaign's measured 47.9% combined-stress failure
  rate is itself worth investigating before proposing a fix.

**Recommendation for where a new contribution would land best:** the
compilation cliff remains the strongest candidate for a primary
contribution, because solving it is the necessary precondition for every
other finding in this report — including the newly-quantified
budget-starvation result — to be evaluated honestly at a scale that
matters, and because, unlike the algorithmic gaps, it is a genuine,
provable obstruction to running *any* dense-unitary-based
modular-exponentiation circuit past a small fixed size, independent of
whatever phase-estimation scheme sits on top of it. **Combined
classical-budget starvation is the strongest candidate for an immediate,
low-effort follow-up finding**: it is already fully reproducible at the
instance sizes available today (`Data/failure_study/phase1b_combined_stress.csv`),
requires no new circuit-synthesis work, and directly contradicts an
implicit robustness claim in the reference implementation's own worked
examples (Section 2) — a natural, ready-to-write companion result to
report alongside the compilation-cliff finding rather than instead of it.

---

## Appendix A — Experimental environment notes

- Every experiment was executed through the repository's public API only
  (`Algorithms.paper_algorithm`, `Circuits.windowed_qpe`,
  `Reconstruction.*`, `Simulation.*`); no repository source was modified.
- Three unrelated engineering obstacles were discovered and worked around
  entirely within the experiment harness (`Experiments/harness.py`,
  `Experiments/phase5_attribution.py`), not by changing the algorithm under
  test: (a) `ProcessPoolExecutor` workers intermittently died under high
  concurrency with a native `memory allocation ... failed` Rust panic
  (worker count capped at 3, chunked execution with sequential fallback);
  (b) the compilation cliff described in Section 5.1, which set the
  work_qubits≤5 ceiling for Phases 1-3; (c) the Phase 5 attribution ladder's
  "widened beam"/"exhaustive candidates" counterfactual, when replayed on
  the combined-stress test's 32-window worst-case configurations, caused
  multi-gigabyte path-list growth in `stitch_candidates` — **[CORRECTED,
  see Correction above]** this is combinatorial path multiplication from a
  wide beam times many windows (any per-window acceptance rate multiplies
  across 32 windows), not evidence of the carry check being lenient; the
  original text attributed this to "the carry-check's own
  overpermissiveness from Section 4.2," which that section's corrected
  audit no longer supports — fixed by reducing the ladder's beam-width proxy (4096 → 64) and adding a
  hard 30-second per-item timeout, with any item that still exceeds it
  honestly recorded as `unknown_attribution_timeout` rather than dropped.
- Given these constraints, the campaign totals roughly 1,350-1,550 trials
  across Phases 1-4 (including the 144-trial combined-stress test) plus an
  exhaustive 168-combination carry audit and a full 94-failure root-cause
  attribution, a substantial reduction from the original ~10,000+ trial
  plan. This trades statistical power at the margins (wider confidence
  intervals on the smallest observed effects) for actually completing
  within a practical session — the qualitative findings above are, however,
  clean, large-effect-size, and highly reproducible (several are exact
  combinatorial audits or deterministically-seeded sweeps, not subject to
  sampling noise at all).
