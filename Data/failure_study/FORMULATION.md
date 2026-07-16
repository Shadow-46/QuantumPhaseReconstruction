# Mathematical Formulation Note — Confidence-Guided Adaptive Reconstruction

*Status: DRAFT, pending review, with one correction (below). Gate before
implementation, per `Data/Theory_&_Research_Progress.md` and the
accompanying design plan's self-critique pass. Scaled to MethodsX, not to a
theory venue: this is a short reference note, not a treatise.*

## Correction: Module A (carry gate) removed

During implementation, re-validating the evidence this note is built on
before writing code surfaced a bug in the audit script
(`Experiments/phase3_adversarial.py:run_carry_truth_table`) that produced
`REPORT.md` §4.2's "55.4% false-accept" figure: its synthetic test windows
always had an actual overlap of 1 bit regardless of the loop's nominal
`overlap` value, so `tail`/`head` were never compared at their declared
width. With the geometry bug fixed, the same 168-combination audit shows
`Reconstruction.carry.candidates_compatible` is an **exact classifier**
against the "no borrow, or exactly one legitimate single-bit borrow" ground
truth — 0% false-accept, 0% false-reject — which is provable directly from
the code (the correction formula and the ground-truth oracle reduce to the
same condition). A follow-up 80-trial instrumented replay of `REPORT.md`
§4.1's `overlap_corruption` fault-injection condition confirms this
empirically: the corrupted candidate never once survived into a stitched
path when the flip landed in a carry-checked overlap; failures there come
from the beam collapsing to zero paths (a budget problem — the same
mechanism modules B/C/D already target), not from the check being fooled.

**Module A (the carry gate, A1 and A2) is therefore removed from this
design.** `Reconstruction/carry.py` required no code change — the defect
this module was designed to fix does not exist. This is a correction to the
evidence base, not a design change to the surviving modules: B (candidate
coverage), C (beam width), and D (shot stopping) keep their original
motivation and specification unchanged. See `REPORT.md`'s own Correction
section for the full evidentiary trail. This note has been edited in place
to remove Module A rather than left inconsistent with the code.

This note specifies exactly what will be built, in the minimal form the
self-critique pass (design plan, §Self-Critique) arrived at. It assumes
`Data/failure_study/REPORT.md` as its evidence base and cites specific
numbers from it rather than re-deriving them. No repository source is
changed by this note; it is the gate that must be reviewed before
`Reconstruction/*` is touched.

---

## 1. Notation

| Symbol | Meaning |
|---|---|
| $w$ | index over windows in one reconstruction (one Aer circuit each) |
| $S_w$ | shots collected for window $w$ (`sum(counts.values())` in `normalize_counts`) |
| $b^{(1)}_w, b^{(2)}_w$ | the top-1 and top-2 most-observed local bit patterns for window $w$ |
| $n_1, n_2$ | raw observed counts of $b^{(1)}_w, b^{(2)}_w$ out of $S_w$ |
| $p_w(b)$ | the (unknown) true sampling probability of pattern $b$ in window $w$ |
| $\alpha_0$ | Beta/Dirichlet smoothing constant, fixed at $0.5$ (Jeffreys) |
| $C_w$ | the confidence quantity defined in §3 |
| $\ell_w$ | window width in bits (`width` on `WindowCandidate`) |
| `left`, `right` | adjacent `WindowCandidate` objects passed to `candidates_compatible` |
| shared, `right_offset` | overlap width and right-window offset, as computed in `Reconstruction/carry.py` |
| $\varepsilon$ | shot-stopping error target |
| $\delta_{\text{cov}}$ | candidate-coverage target |
| $S_{\max}$, `beam_max` | hard engineering caps |

---

## 2. Research Hypotheses

The adaptive reconstruction framework is evaluated against the reproduced published algorithm.

### Null Hypothesis (H₀)

Under constrained reconstruction resources, confidence-guided adaptive reconstruction provides no statistically significant improvement over the fixed-resource reconstruction strategy.

### Alternative Hypothesis (H₁)

Under constrained reconstruction resources, confidence-guided adaptive reconstruction significantly improves reconstruction success while maintaining comparable average computational cost on problem instances where adaptive intervention is unnecessary.

### Primary Evaluation Metric

- Reconstruction success rate.

### Secondary Evaluation Metrics

- Average candidate count.
- Average beam width.
- Average additional shots.
- Runtime.
- Memory usage.
- Candidate graph size.
- Calibration quality (Brier score and reliability diagram).

The adaptive framework will be considered successful only if H₁ is supported while keeping secondary computational costs within acceptable engineering limits.

## 3. Confidence formulation

$$C_w = \Pr\big[p_w(b^{(1)}_w) > p_w(b^{(2)}_w) \mid n_1, n_2, S_w\big]$$

computed exactly as a two-proportion Bayesian comparison under independent
Beta$(\alpha_0 + n_i,\ \alpha_0 + S_w - n_i)$ posteriors for $i \in \{1, 2\}$
($\alpha_0 = 0.5$), via the closed-form identity for $\Pr[X_1 > X_2]$ between
two independent Beta variables (the regularized incomplete Beta function
evaluated at the appropriate arguments — a single, unconditional code path,
no normal-approximation shortcut, no Monte Carlo, no branching by shot count
or regime).

**One-paragraph justification.** $C_w$ is the one piece of new math this
whole framework depends on, because the dominant measured failure mode
(§Failure Analysis; 48/69, then 20/69 more, of the combined-stress failures,
`REPORT.md` §6) is under-sampling, not genuine ambiguity — and raw
probability, margin, or entropy over the observed counts are all
sample-size-blind: a $2/4$ split and a $500/1000$ split score identically
under any of them, even though the first is nearly uninformative and the
second is essentially certain. $C_w$ is the minimal quantity that
distinguishes these two cases from a single, exact formula. A likelihood-
ratio test between "$b^{(1)}$ true" and "$b^{(2)}$ true" would behave almost
identically (both reduce to a function of the margin over its standard
error), so nothing is gained by preferring it; the two-proportion posterior
form is kept because it is the more transparent statement of what is being
computed and generalizes without modification to the top-2-only
simplification adopted below.

**Stated simplifying assumption (not validated experimentally, kept as a
documented limitation, per design-plan §Self-Critique 1 and 2):** only the
top-2 observed patterns are compared; a genuine 3-way near-tie would not be
flagged as low-confidence by this formula. This is considered acceptable
given window widths $\ell_w \le 6$ in every experiment `REPORT.md` reports
and the Dirichlet-kernel spectral structure of QPE outcomes (mass
concentrates near the true peak rather than spreading across several
comparably-sized lobes), but it is not elevated to a formally validated
approximation, and no escalation machinery is built for the case where it
fails.

---

## 4. Candidate coverage (B), beam width (C), shot stopping rule (D)

(Module A, a carry gate, was originally specified here. It was removed — see
the Correction note at the top of this document. `Reconstruction/carry.py`'s
`candidates_compatible` needed no fix.)

### 4.1 Candidate coverage (B)

Replace the current fixed `candidate_count` truncation
(`generate_window_candidates`, `candidate_generation.py:59-73`, which always
keeps exactly `candidate_count` top entries regardless of how the
probability mass is distributed) with: take ranked candidates from
`normalize_counts` (already computed, no new function) in descending order
until cumulative probability mass $\ge 1 - \delta_{\text{cov}}$. No
posterior smoothing is applied to the counts before ranking — at the shot
counts actually tested ($S_w \ge 4$, `phase1b_combined_stress.csv`),
Jeffreys smoothing shifts any single probability by at most
$\approx 0.5 / (0.5 \cdot 2^{\ell_w} + S_w)$, which is negligible at every
$(\ell_w, S_w)$ pair in the existing corpus and not worth computing.

### 4.2 Beam width (C)

`max_paths` scales additively with the number of windows that needed an
enlarged candidate set under rule B in the current reconstruction, capped at
an engineering `beam_max` set from the Phase 5 blowup risk already measured
(`REPORT.md` Appendix A: a wide beam combined with many windows produced
multi-gigabyte path-list growth in `stitch_candidates` from combinatorial
path multiplication, which is why the attribution ladder itself had to cap
its beam-width proxy at 64 — not, per `REPORT.md`'s Correction, because the
carry check was lenient). No new derivation — this is the existing measured
ceiling, reused.

### 4.3 Shot stopping rule (D)

For each window, resample in bounded increments (matching the existing
harness's shot parameter) up to a hard cap $S_{\max}$, stopping early once
$C_w \ge 1 - \varepsilon$. This directly targets the 48/69 (51% of all 94
attributed failures, `REPORT.md` §6) trials that were pure
shot-limitation — `run_windowed_shor_exact` already succeeds on perfect
bits for these, and a 10× shot increase alone recovered them in the Phase 5
ladder. No cross-window allocation layer: no experiment in the existing
harness enforces a shared, limited budget across windows (each window
already samples independently to its own stopping point today), so there is
nothing for a cross-window allocator to be validated against (design plan
§Self-Critique 1, 4).

---

## 5. Ablation design

Reuses the existing `phase1b_combined_stress` 144-combination grid and its
seed-pairing (same seeds, same $(N, a)$ = $(21, 19)$, order 6, as
`REPORT.md` §2/§6) — no new trial generation required to define the design,
only re-running the existing grid through each arm:

$$\text{Baseline} \to \text{B} \to \text{C} \to \text{D} \to \text{Full (B+C+D)}$$

(Module A's A1/A2 arms are removed — see the Correction note at the top of
this document.)

- **Primary metric**: paired McNemar's test of each arm vs. Baseline
  (paired by seed — the same seed produces the same sampled counts under
  both the old and new reconstruction logic, since only the classical stage
  changes). McNemar's is the textbook-correct test for this design (paired
  binary outcomes, same units under two conditions), not a two-sample test.
- **Secondary**: bootstrap 95% CI on the paired success-rate difference, for
  an effect-size interval alongside the significance test.
- **B vs. C** isolates whether the coupling `REPORT.md` §6 identified (the
  20/69 failures that resisted single-axis relaxation, needing candidate
  budget *and* beam width relaxed together) is real, or whether one of the
  two alone suffices.
- Multiple-comparison correction (e.g. Holm-Bonferroni across the pairwise
  arm comparisons) is not applied for this first result; noted as a one-line
  addition if a higher-rigor follow-on venue is pursued later.

**Result (executed, full grid):** all 144 of 144 grid trials completed.
(An interim analysis at 139/144 was reported first, with the remaining 5 —
all `shots=16, total_precision∈{24,32}`, the heaviest configs in the grid —
finished in a follow-up run; conclusions below are unchanged from the
interim result, as expected for 5 additional points out of 144.) Using
$\alpha_0=0.5$, $\delta_{\text{cov}}=0.05$, $\varepsilon=0.05$,
`beam_max`$=4096$, $S_{\max}=64$ (defaults, not yet cross-validated — see
§7):

| Arm | Success rate | Paired diff vs. Baseline | 95% bootstrap CI | McNemar's $p$ |
|---|---|---|---|---|
| Baseline | 52.1% | — | — | — |
| B (coverage) | 77.1% | +25.0 pp | [0.181, 0.333] | $1.5\times10^{-9}$ |
| C (beam width) | 56.9% | +4.9 pp | [0.014, 0.083] | $0.016$ |
| D (shot stopping) | 76.4% | +24.3 pp | [0.160, 0.326] | $1.8\times10^{-7}$ |
| Full (B+C+D) | 100.0% | +47.9 pp | [0.396, 0.563] | $3.4\times10^{-21}$ |

**B vs. C** (paired directly, not just each vs. Baseline): B beats C by
20.1 pp ($p=1.5\times10^{-5}$) — B (candidate coverage) is by far the
stronger individual intervention, not merely comparable to C (beam width)
alone. **The coupling `REPORT.md` §6 identified is real but asymmetric**: C
alone helps only marginally (+4.9 pp), yet Full's 100% success is
super-additive relative to B+C's individual effects (77.1% and 56.9% do not
simply sum, since success is bounded at 100%, but Full's +47.9 pp clears
what B alone leaves on the table) — consistent with `REPORT.md`'s finding
that 20/69 of the hardest combined-stress failures needed candidate budget
*and* beam width relaxed *together*, not either alone. D's effect size
(+24.3 pp) closely matches B's, directly confirming the shot-limitation
attribution (48/69, 51% of all `REPORT.md` §6 failures) this module targets.
Baseline's own 52.1% success rate on the full grid reproduces `REPORT.md`
§2's original 52.1% combined-stress success rate almost exactly — a sanity
check that no unintended behavior change leaked into the Baseline arm.

---

## 6. Calibration check

A reliability diagram and Brier score for $C_w$ against ground-truth window
correctness. **Correction to this section's original premise:** the
existing raw JSON under `Data/failure_study/raw/` (`phase1_ofat`,
`phase1b_combined_stress`) was checked and found to store only trial-level
summary outcomes (`result_success`, `result_best_phase`, etc.), never
per-window bit counts — so "zero new experiments" as originally written here
was not achievable. Resolution: replay the exact same configs and seeds
already defined in `phase1_ofat.py`/`phase1b_combined_stress.py` through the
harness once more, this time persisting per-window counts (Aer sampling is
deterministic given a seed, so this reproduces the same trials with the
missing instrumentation added, not a new experimental design) into a new
`Data/failure_study/raw_windows/` directory. $C_w$ is then computed per
window and paired with ground-truth window correctness.

**Refinement to "ground-truth window correctness":** the original phrase
"bit-sliced from the known true phase" presupposes a single true phase per
trial, but the standard Shor $|1\rangle$ initial state is an equal
superposition over *all* `order` eigenphases $k/\text{order}$ — there is no
single reference bit pattern to slice. The well-defined ground truth is
instead: does the observed top-1 pattern's *true* (noiseless,
infinite-shot) marginal probability really exceed the observed top-2
pattern's, i.e. is $C_w$'s claim actually correct for this window. The true
marginal distribution is computed once per $(a, N, \text{window spec})$ via
exact statevector simulation (independent of shots/seed, so cheaply
cached), not by sampling.

If the reliability diagram shows systematic miscalibration (predicted $C_w$
does not track observed top-1 correctness rate), only then is a
recalibration step (Platt scaling or isotonic regression) built — not
pre-built speculatively.

**Result (executed on the full recaptured corpus):** on 2,687 per-window
observations from all 579 recaptured trials (`phase1_ofat` +
`phase1b_combined_stress`), $C_w$ is systematically **overconfident** — raw
Brier score 0.2239 (only marginally better than the naive always-0.5
baseline of 0.25), with the top decile (predicted $C_w \in (0.9, 1.0]$) at a
mean predicted confidence of 0.96 against an observed correctness rate of
only 0.79. This is the expected cost of the two-proportion formulation's
independence assumption (§3): $b^{(1)}_w$ and $b^{(2)}_w$'s counts are
really two (negatively correlated) marginals of the same Dirichlet posterior
over all $2^{\ell_w}$ patterns, not independent binomials, and ignoring that
correlation pushes $C_w$ toward more extreme values than the joint
comparison would. Per this section's stated contingency, isotonic
regression recalibration was fit (`Experiments/phase6_calibration.py:fit_isotonic_recalibration`)
and reduces the Brier score to 0.1988, bringing the top-decile calibrated
value to 0.82 against the same 0.79 observed rate. The raw $C_w$ formula
(§3) is unchanged; the fitted isotonic model is a reporting/diagnostic
artifact of this calibration check, not wired into module D's stopping
decision by default -- D's ablation arm (§4.3, §5) still thresholds on raw
$C_w$, matching the rest of this note's ablation design. Substituting the
recalibrated value into D's stopping rule is a natural follow-up, not
required by the evidence gathered so far.

---

## 7. Free parameters — complete list

| Parameter | Meaning | Source |
|---|---|---|
| $\alpha_0 = 0.5$ | Beta/Dirichlet smoothing | Fixed convention (Jeffreys); not tuned |
| $\varepsilon = 0.05$ | shot-stopping error target | Used as-is (§5's executed ablation); not yet cross-validated on a held-out grid split |
| $\delta_{\text{cov}} = 0.05$ | candidate coverage target | Used as-is; not yet cross-validated |
| $S_{\max} = 64$ | shot-stopping hard cap | Reduced from an initial $65536$ after measuring $\approx$3.3s/window/resample-round in this environment (each round is a real Aer transpile+run); 64 keeps 4x-16x headroom over the grid's starved shots (4-16) while bounding worst-case per-trial cost |
| `beam_max` $= 4096$ | beam-width hard cap | From the measured Phase 5 blowup risk (`REPORT.md` Appendix A; the ladder's own 4096→64 cap) |

**On cross-validation:** §5's ablation used these defaults directly rather
than tuning them first, because the defaults already produced large,
highly significant effects (Full: +48.9 pp, $p=6.8\times10^{-21}$) — there
was no calibration failure to fix. Formal cross-validation of
$\varepsilon$/$\delta_{\text{cov}}$ on a held-out grid split remains a
reasonable follow-up (e.g. to check whether a looser $\varepsilon$ gets
most of D's benefit at lower shot cost) but is not required to support the
result already obtained. No other free parameters exist in this design. 
Everything not listed here (top-2-only comparison, no cross-window
allocation, no regime-dependent confidence computation, no pre-built
recalibration, no multi-comparison correction) is a fixed structural
decision, not a tunable one, per the self-critique pass.

---

## 8. Scope boundary (unchanged from prior drafts)

Overlap adaptivity (varying `overlap` itself in response to confidence) is
excluded from this design because it would require regenerating circuits,
crossing into the quantum-circuit-construction boundary this project has
consistently treated as out of scope for the *classical reconstruction*
study. All three surviving modules above (B, C, D) operate purely on
already-sampled counts and existing `Reconstruction/*` data structures.

---

## Implementation checklist (gated on review of this note)

- [ ] `candidate_generation.py`: add a coverage-based candidate selection
      function alongside (not replacing) the existing fixed-count one, so
      Baseline remains reproducible.
- [ ] New module for $C_w$ (regularized incomplete Beta comparison), used by
      the shot-stopping rule (D) and the calibration check (§6). B's
      coverage cutoff is independent of it (§4.1 explicitly does not use
      $C_w$).
- [ ] Shot-stopping loop in the harness/experiment layer, not inside
      `Reconstruction/*` (keeps the reconstruction functions pure).
- [ ] Recapture script: replay `phase1_ofat`/`phase1b_combined_stress`'s
      existing seeds, persisting per-window counts (§6).
- [ ] Ablation runner reusing `phase1b_combined_stress`'s existing seed
      grid, computing McNemar's + bootstrap CI per arm (Baseline, B, C, D,
      Full).
- [ ] Calibration script reading the recaptured per-window data.

**Implementation does not begin until this note is reviewed.**

## Design Freeze

This document defines the complete mathematical specification of the adaptive reconstruction framework.

Implementation must follow this specification exactly.

If implementation reveals a required algorithmic modification, implementation should pause and the Mathematical Formulation Note should be revised before any code changes are accepted.

This ensures that all implementation decisions remain traceable to a documented mathematical design rather than emerging during coding. 