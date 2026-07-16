# Correctness Arguments — Confidence-Guided Adaptive Reconstruction

*Written for a MethodsX-style paper. Companion to `THEORETICAL_ANALYSIS.md`
(mathematical foundations) and `COMPLEXITY_ANALYSIS.md` (cost analysis).
Every lemma/proposition below is anchored to a specific function in
`Reconstruction/*` or `Experiments/harness.py`. Claims are labeled:*

- ***(PROVEN)*** — follows deductively from the code as written, stated as
  a theorem/lemma with a complete or code-complete proof.
- ***(HEURISTIC)*** — a design rationale with a plausibility argument, not
  a proof; no counterexample is known, but none is ruled out either.
- ***(EMPIRICAL)*** — established only by measurement (the ablation,
  calibration, or audit studies in `FORMULATION.md`/`REPORT.md`); explicitly
  not derived from the code's structure.

No claim in this document is asserted at a stronger level than the evidence
in `FORMULATION.md`, `REPORT.md`, and `PHASE6_SUMMARY.md` supports.

---

## 1. Preliminaries: what "correctness" means for this pipeline

The reconstruction pipeline is correct on a trial if the phase estimate
$\hat\phi$ it returns satisfies
`recover_order_and_factors(\hat\phi, a, N).factors == ` the true
factorization of $N$. This is an end-to-end, black-box notion of
correctness (`REPORT.md` §0's own scope statement), not a per-component
specification. The lemmas below decompose it into per-stage invariants
sufficient (but, per §1.1, not individually necessary) for end-to-end
correctness.

**Lemma 0 (Sufficient decomposition) *(PROVEN)*.** If (i) every window's
top-ranked retained candidate bit pattern equals the true window bit
pattern, and (ii) `candidates_compatible` accepts exactly the
carry-consistent pairs among those true patterns, and (iii) beam pruning
never discards the resulting fully-true path, and (iv) $2^{\ell}\hat\phi/2^\ell$'s
continued-fraction expansion recovers the true order, then the pipeline
returns the correct factorization.

*Proof.* Under (i), `stitch_candidates`'s per-window candidate groups each
contain the true local bit pattern. Under (ii) (proven exactly in Lemma 1
below), the true patterns are pairwise compatible across every overlap, so
a fully-true path exists among the generated `next_paths` at every level.
Under (iii), this path is not pruned by the `[:max_paths]` truncation, so
it survives to the final `stitched` list; since its score is a product of
observed probabilities and the true pattern is (by hypothesis) top-ranked
in every window, it receives the maximal product score and is therefore
first in `sorted(stitched, ..., reverse=True)`. `recover_order_and_factors`
is then called on it first; under (iv) this succeeds and returns the true
factors, so the outer loop in `reconstruct_adaptive`/`reconstruct_from_counts`
returns on this first iteration. $\blacksquare$

This decomposition is a logical convenience for organizing the lemmas
below, **not** a claim that (i)–(iv) are necessary — a correct answer can
also arise from a non-top-ranked path, a legitimate carry-borrow correction
substituting for exact equality, or a non-first `stitched` candidate
succeeding where an earlier one fails continued-fraction recovery. The
sufficiency framing is retained because it is what the correctness lemmas
below can actually establish component-by-component.

---

## 2. Lemma 1 — the carry-aware compatibility check is an exact classifier

**Statement *(PROVEN)*.** For any two `WindowCandidate` objects `left`,
`right` with shared overlap width $o=$`overlap_width(left, right)` $>0$,
`candidates_compatible(left, right)` returns `True` if and only if the
overlapping bit strings are consistent with "no borrow" or "exactly one
legitimate single-bit borrow with carry indicator $c=$`carry_bit(right, right_offset+shared)`" —
i.e., iff

$$
\mathrm{tail} = \mathrm{head} \quad\text{or}\quad (\mathrm{int}(\mathrm{tail},2) - c) \bmod 2^{o} = \mathrm{int}(\mathrm{head},2).
$$

*Proof.* This is the function's literal implementation
(`candidates_compatible`, `Reconstruction/carry.py:62-76`): it computes
`tail`, `head` from the two candidates' `local_bits` at the overlapping
offsets, returns `True` immediately if `tail == head`, and otherwise
returns the truth value of exactly the stated modular equation with `c =
carry_bit(right, right_offset + shared)`. There is no other branch. The
"ground-truth oracle" used by the corrected audit
(`Experiments/phase3_adversarial.py:run_carry_truth_table`, post-fix) is
defined as "no borrow, or exactly one legitimate single-bit borrow" — the
same predicate, applied to the same tail/head/carry-bit triple — so the
function's return value and the oracle's are the same Boolean function of
`(tail, head, c)` by construction, not merely in agreement on tested
inputs. $\blacksquare$

**Corollary (0% false-accept, 0% false-reject) *(PROVEN, corroborated
EMPIRICALLY)*.** Because the implementation and the ground-truth oracle
are the same formula, `candidates_compatible` cannot disagree with the
oracle on *any* input, tested or not — this is stronger than the 168-combination
exhaustive audit result (`Data/failure_study/phase3_carry_truth_table.csv`:
42 true-accept, 126 true-reject, 0 false-accept, 0 false-reject) by itself
would establish, since a purely empirical exhaustive check over a finite
input space is still a proof only *for that space*. Here, because the two
formulas are algebraically identical (not merely observed to agree), the
result generalizes to every `(tail, head, c)` triple over any overlap
width, not only the audited $o\in\{1,2,3\}$ range. The audit is retained
in the evidence base as independent corroboration (and as the mechanism by
which the *original*, buggy audit's false 55.4%-false-accept claim was
caught and retracted — see `FORMULATION.md`'s Correction note), not as the
sole basis for this lemma.

**What this lemma does *not* say.** It says nothing about whether the
*true* bit patterns ever reach `candidates_compatible` as `left`/`right` in
the first place — that depends on candidate generation (§3) and shot
sufficiency (§4), not on the check itself. `REPORT.md` §4.1's residual
overlap-corruption failure rate (43.75% success under deliberate bit
corruption) is consistent with Lemma 1: the corrupted candidate is
*correctly rejected* by this exact classifier (confirmed by an 80-trial
instrumented replay, `FORMULATION.md`'s Correction note), and the resulting
failures trace to beam collapse (§5 below), not to a compatibility-check
defect.

---

## 3. Invariant — candidate generation always returns a non-empty, correctly-encoded set

**Invariant 1 *(PROVEN)*.** For any non-empty `counts` and valid window
geometry, both `generate_window_candidates` and
`generate_window_candidates_by_coverage` return a list of length $\ge 1$,
each element's `local_bits` a width-$\ell_w$ bit string extracted from a
key of `counts`, and each element's `value` equal to that local pattern
left-shifted into its correct absolute bit position
(`value = local << suffix_width`, `suffix_width = total_precision - start - width`).

*Proof (length $\ge1$).* `generate_window_candidates`: `candidate_count`
is validated `>= 1` (raises `ValueError` otherwise), and `ranked` is a
non-empty list (since `counts` is required non-empty by
`normalize_counts`'s own `ValueError` guard) sliced to
`[:candidate_count]`, which is non-empty whenever `candidate_count>=1` and
`ranked` is non-empty. `generate_window_candidates_by_coverage`: the loop
condition `if cumulative >= target and candidates: break` cannot trigger
on the first iteration (`candidates` is empty), so at least one candidate
is always appended before the loop can exit via the `break`; the loop can
only otherwise exit by exhausting `ranked`, which is non-empty. Either way
$\ge 1$ candidate is returned. $\blacksquare$

*Proof (correct encoding).* Both functions compute `window_bits` as the
last `width` characters of the counts key (or the key itself if already
that length), convert via `int(window_bits, 2)`, and shift by
`total_precision - start - width` — this is a direct transcription of "place
this window's local pattern at absolute bit offset `start`, leaving `start`
low bits (from `suffix_width`) unset, matching the placement `stitch_candidates`
later assumes when merging via `bits[candidate.start + offset] = bit`." No
arithmetic beyond this shift is performed, so the encoding is correct by
construction whenever the input `counts` keys are correctly-widthed
strings — an invariant maintained by `run_windowed_qpe_block`'s own output
format (upstream of these functions, outside their scope of responsibility). $\blacksquare$

This invariant is what makes Lemma 0's step (i) at least *satisfiable* in
principle (the true pattern, if present in `counts`, is correctly encoded
if selected) — it does not itself guarantee the true pattern is selected.

---

## 4. Proposition — Module B's selection rule is minimal-prefix-correct

Restated from `THEORETICAL_ANALYSIS.md` §6.2 for completeness in this
document's proof register.

**Proposition B1 *(PROVEN)*.** `generate_window_candidates_by_coverage`
returns the shortest prefix, under descending-probability ranking, whose
cumulative empirical mass is $\ge 1-\delta_{\mathrm{cov}}$.

*Proof.* See `THEORETICAL_ANALYSIS.md` §6.2, Proposition 5 — restated here:
the loop appends candidates in descending-probability order and checks the
break condition (`cumulative >= target`) only *after* the most recent
append, using the value of `cumulative` from before that append was
included in a prior check — concretely, tracing the loop: at the start of
iteration $i$, `cumulative` holds $\sum_{j<i}\hat p_w(b_{(j)})$; if this
already meets `target` (and at least one candidate exists), the loop exits
*without* including $b_{(i)}$. So the returned set is exactly
$\{b_{(1)},\dots,b_{(m)}\}$ where $m$ is the least index at which the
partial sum up to $m$ (inclusive) first reaches `target` — i.e. the
provably minimal prefix. $\blacksquare$

**What this proposition does *not* say *(explicitly a limitation, not an
oversight)*.** It is silent on whether the true pattern is in this prefix.
A formal statement of that stronger property would require a
finite-sample concentration bound (e.g., relating $\hat p_w$ to $p_w$ via
a Chernoff-type inequality) that is neither present in `FORMULATION.md`
nor implemented in `candidate_generation.py`. Framing $\delta_{\mathrm{cov}}$
as an "error target" in the surrounding documentation is therefore a naming
choice that should not be read as implying a proven error-rate guarantee.

---

## 5. Proposition — beam-search stitching is exact only when unsaturated

**Proposition S1 (Conditional exactness) *(PROVEN)*.** If, at every level
$w=2,\dots,W$ of `stitch_candidates`, the number of *compatible*
next-paths generated (before the `[:max_paths]` truncation) is $\le$
`max_paths`, then `stitch_candidates` returns the true global optimum
(under the product-of-weights score) among all paths constructible from the
given per-window candidate groups — i.e. the beam search is exact
(equivalent to exhaustive search) whenever it is never actually forced to
discard a candidate path by the truncation.

*Proof.* If `len(next_paths) <= max_paths` at every level, the line
`paths = sorted(next_paths, ...)[:max_paths]` is a no-op truncation (the
slice takes everything), so `paths` after each level is exactly the full
set of surviving compatible paths, sorted. By induction over levels, the
final `paths` before the stitching return is exactly all globally
consistent paths through all $W$ window groups, and the returned
`stitched` list — sorted by score — necessarily contains the true maximum
as its first element. $\blacksquare$

**Corollary (Beam search is a lossy heuristic when saturated) *(PROVEN)*.**
Whenever `len(next_paths) > max_paths` at some level, the discarded paths
are permanently unrecoverable (no backtracking; the outer `for group in
windows[1:]` loop only ever extends `paths`, never widens it after
truncation) — if the true path is among the discarded set at any level, it
cannot appear in the final output, regardless of its true-pattern content
at later levels. This is a direct reading of the loop structure (no
exception, no fallback), not an inference.

**Interpretation for modules B and C.** Proposition 6
(`THEORETICAL_ANALYSIS.md` §7.2) shows $P' \ge P$ under module C's cap
hypothesis, which — via Proposition S1's monotonicity in `max_paths` — can
only *reduce or preserve*, never increase, the risk of the discard case in
the Corollary above, relative to baseline, for the *same* candidate
groups. Whether this is *sufficient* to avoid discarding the true path on
any specific trial is an empirical question, answered only by the ablation
(§7 below), not by this proposition alone — Proposition S1 characterizes
*when* the search is exact, it does not prove that condition holds for any
particular configured `max_paths`/`beam_max`.

---

## 6. Lemma — shot-stopping termination and threshold semantics

**Lemma D1 (Termination) *(PROVEN)*.** `_resample_shot_stopping`'s
per-window loop terminates in at most
$\lceil\log_2(S_{\max}/S_w^{(0)})\rceil + O(1)$ rounds. *(Full proof:
`THEORETICAL_ANALYSIS.md` §8.2, Proposition 7 — restated here by
reference, not reproduced, to avoid duplicating the doubling-schedule
argument.)*

**Lemma D2 (Stopping-condition semantics are relative, not absolute)
*(PROVEN, with an EMPIRICAL sharpening)*.** The condition `C_w >= 1 -
epsilon` that halts resampling for a window is evaluated on *raw* $C_w$
(`top1_vs_top2_confidence`), and *(PROVEN)* raw $C_w$ is not recalibrated
anywhere in the call path from `_resample_shot_stopping` (the function
imports and calls only `top1_vs_top2_confidence`; the isotonic
recalibration model lives entirely in
`Experiments/phase6_calibration.py` and is never imported by
`harness.py`). Combined with the *(EMPIRICAL)* calibration-study result
(top predicted-confidence decile: mean $C_w=0.96$, observed correctness
$0.79$; raw Brier score $0.2239$, `FORMULATION.md` §6), it follows that a
window halting resampling at, e.g., $\varepsilon=0.05$ is empirically
correct only $\approx79\%$ of the time within that top-confidence decile,
not $\ge95\%$. **The threshold is a monotone decision rule (Proposition 3,
`THEORETICAL_ANALYSIS.md` §5.2), not a calibrated error-control mechanism**,
and this document deliberately does not restate it as one.

---

## 7. Structural invariant — the framework is a strict generalization of baseline

**Invariant 2 (Reduction to baseline) *(PROVEN)*.** Calling
`reconstruct_adaptive(shor, window, noise_model, use_coverage=False,
use_beam_rule=False, use_shot_stopping=False)` executes a call sequence
identical in structure to `reconstruct_from_counts` composed with
`sample_window_counts`: window counts are sampled once, unmodified
(`use_shot_stopping=False` skips the `_resample_shot_stopping` call
entirely — the `if` guard short-circuits, no partial resampling occurs);
candidate groups are built exclusively via `generate_window_candidates`
with the fixed `window.candidate_count` (the `else` branch of the
`if use_coverage:` conditional); `max_paths` is exactly `window.max_paths`
(the `else` branch of `if use_beam_rule:`). No code path touches B/C/D's
data-dependent logic when all flags are `False`.

*Proof.* Direct inspection of `reconstruct_adaptive`
(`Experiments/harness.py:409-452`): each of the three `if use_*:` blocks
has a complete, mutually exclusive `if/else` with the `else` branch
reproducing the pre-adaptive computation exactly (same function calls, same
arguments, sourced from `window.candidate_count` / `window.max_paths`
rather than any adaptive quantity). $\blacksquare$

**Why this matters for correctness argumentation.** It licenses treating
the *observed* Baseline arm of the ablation (`FORMULATION.md` §5:
$52.1\%$; `PHASE6_SUMMARY.md` §5: $51.1\%$) as a direct, code-level
reproduction of `REPORT.md`'s original $52.1\%$ combined-stress figure —
the near-exact match across three independent measurements (original
study, `FORMULATION.md`'s ablation run, `PHASE6_SUMMARY.md`'s ablation run)
is *consistent with* Invariant 2 (a necessary check, not a proof of the
invariant by itself — an invariant proven from code structure is
corroborated, not established, by matching empirical replication) and
serves as the sanity check `PHASE6_SUMMARY.md` §5 explicitly names it as.

---

## 8. Empirical observations (explicitly not theorems)

These are the results the propositions above provide a partial mechanistic
account for, but do **not** derive. Listed here, separated from the proven
material above, per this document's stated requirement not to blur the two.

1. **Ablation success rates** *(EMPIRICAL, `FORMULATION.md` §5 /
   `PHASE6_SUMMARY.md` §5)*: Baseline $51$–$52\%$, B $76$–$77\%$
   ($+25$ pp, McNemar $p<10^{-8}$), C $56$–$57\%$ ($+5$ pp, $p=0.016$), D
   $76$–$77\%$ ($+24$–$25$ pp, $p<10^{-6}$), Full $100\%$
   ($+48$–$49$ pp, $p<10^{-20}$). No proposition above derives these
   *magnitudes* — Proposition 9 (`THEORETICAL_ANALYSIS.md` §9) proves only
   that each arm's effect is non-negative in expectation from the
   monotone-relaxation structure of the code, not its size.
2. **Super-additivity of the Full arm relative to B and C individually**
   *(EMPIRICAL, `FORMULATION.md` §5)*: consistent with, but not proven by,
   the qualitative coupling hypothesis in `REPORT.md` §6 (20/69 failures
   needing both budget axes relaxed together) — no formal joint-failure
   model is derived in this framework to predict the super-additivity
   quantitatively.
3. **$C_w$ miscalibration** *(EMPIRICAL, `FORMULATION.md` §6)*: raw Brier
   $0.2239$ vs. naive-$0.5$ baseline $0.25$; isotonic recalibration reduces
   Brier to $0.1988$. This *is* explained mechanistically by (A2)'s
   falsity (`THEORETICAL_ANALYSIS.md` §3, §5.5) but the specific numeric
   improvement is a fitted, data-dependent result, not a derived one.
4. **Carry-check corrected audit** *(EMPIRICAL, corroborating Lemma 1)*:
   168/168 combinations classified correctly, 0% false-accept/false-reject
   — Lemma 1 explains *why* this must hold for all inputs, not only the
   168 audited ones; the audit is retained as an independent empirical
   check that the code matches its own specification, and as the record of
   how the original (buggy-audit-derived) 55.4% false-accept claim was
   caught.
5. **Baseline reproduction across three measurements** (`REPORT.md` §2,
   `FORMULATION.md` §5, `PHASE6_SUMMARY.md` §5: $52.1\%$, $52.1\%$,
   $51.1\%$) *(EMPIRICAL)*: corroborates, but does not prove, Invariant 2.

---

## 9. Summary table

| Claim | Status | Location |
|---|---|---|
| `candidates_compatible` is an exact classifier | **PROVEN** (algebraic identity) + corroborated by exhaustive audit | §2, Lemma 1 |
| Candidate generation returns $\ge1$ correctly-encoded candidate | **PROVEN** | §3, Invariant 1 |
| Module B returns the minimal covering prefix | **PROVEN** | §4, Proposition B1 |
| Module B's prefix contains the *true* pattern | **not claimed / no guarantee exists** | §4 |
| Beam stitching is exact when unsaturated | **PROVEN** (conditional) | §5, Proposition S1 |
| Beam stitching can permanently lose the true path when saturated | **PROVEN** | §5, Corollary |
| Shot-stopping loop terminates | **PROVEN**, bound $O(\log(S_{\max}/S_w^{(0)}))$ | §6, Lemma D1 |
| Shot-stopping threshold is a calibrated error rate | **false** — monotone but uncalibrated | §6, Lemma D2 |
| All-flags-off reproduces baseline exactly | **PROVEN** (code structure) | §7, Invariant 2 |
| B/C/D each improve success rate, and by how much | **EMPIRICAL** only | §8.1 |
| Full arm reaches 100% on the evaluated grid | **EMPIRICAL** only | §8.1 |
| $C_w$ is measurably overconfident | **EMPIRICAL** only (mechanistically consistent with (A2)'s falsity) | §8.3 |
