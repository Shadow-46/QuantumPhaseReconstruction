# Confidence-Guided Adaptive Reconstruction for Windowed Quantum Phase Estimation

## A Compiled Research Report: Novelty, Mathematical Models, Proofs, and Experimental Validation

*This document compiles and cross-references the project's existing research
artifacts — `README.md`, `Data/Theory_&_Research_Progress.md`, and the five
companion notes in `Data/failure_study/` (`REPORT.md`, `FORMULATION.md`,
`THEORETICAL_ANALYSIS.md`, `COMPLEXITY_ANALYSIS.md`,
`CORRECTNESS_ARGUMENTS.md`, `PHASE6_SUMMARY.md`) — into a single narrative.
It introduces no new claims, formulas, or numbers; every result below is
sourced from, and cross-referenced back to, those documents and the code
they describe. Evidentiary labels are preserved throughout:
**(PROVEN/T)** — proven from stated assumptions; **(HEURISTIC/H)** —
plausible, not proven; **(EMPIRICAL/E)** — established only by measurement.*

---

## Table of Contents

1. [Overview](#1-overview)
2. [Background & Novelty](#2-background--novelty)
3. [Mathematical Models](#3-mathematical-models)
4. [Proofs & Correctness Arguments](#4-proofs--correctness-arguments)
5. [Complexity Analysis](#5-complexity-analysis)
6. [Experimental Validation](#6-experimental-validation)
7. [Limitations](#7-limitations)
8. [Conclusion & Future Work](#8-conclusion--future-work)

---

## 1. Overview

`QuantumPhaseReconstruction` is a research-grade Python/Qiskit platform that
reproduces and extends **windowed quantum phase estimation (WQPE)** for
Shor's algorithm, following arXiv:2507.22460 and arXiv:2509.05010. The
codebase is organized around a hard boundary between two stages
(`README.md`):

- **Quantum sampling** (`Circuits/`, `Simulation/`) — modular
  multiplication, inverse QFT, standard QPE, and windowed measurement
  circuits, run on Aer with optional noise models.
- **Classical reconstruction** (`Reconstruction/`) — window candidates,
  carry/overlap checks, beam-search stitching, and continued-fraction order
  recovery, operating purely on already-sampled measurement counts.

The project's own stated research objective (`Data/Theory_&_Research_Progress.md`)
is explicit about which side of that boundary it targets:

> "The goal of this research is to improve the classical reconstruction
> stage of Windowed Quantum Phase Estimation (WQPE) used inside the Modular
> Shor's Algorithm. The objective is NOT to modify the quantum circuit
> itself. Instead, the objective is to improve the classical reconstruction
> pipeline after quantum measurements have been obtained."

The project moved through four stages, each documented in a dedicated
source file this report draws from:

1. **Reproduction** — a from-scratch Qiskit implementation of the published
   windowed algorithm (§2).
2. **Failure-mode characterization** — a black-box stress-testing campaign
   (`REPORT.md`, Phases 1–5) that found the published algorithm's apparent
   robustness was misleading, and that a large majority of natural failures
   traced to a fixed, non-adaptive classical resource budget (§6).
3. **Mathematical formulation** of a fix — a confidence-guided adaptive
   reconstruction framework (`FORMULATION.md`, `THEORETICAL_ANALYSIS.md`;
   §3–§4).
4. **Implementation and validation** (Phase 6, `PHASE6_SUMMARY.md`) — the
   framework was built, ablated, and calibrated against ground truth (§6).

A distinctive feature of this project, preserved intact throughout this
report rather than smoothed over, is that **one of its own headline findings
was later retracted**: an initially reported 55.4% false-accept rate in the
carry-check mechanism was traced to a bug in the *audit script* that
measured it, not a defect in the mechanism itself. The correction is
documented at the source (`REPORT.md`'s and `FORMULATION.md`'s own
"Correction" sections) and is carried through this compiled report exactly
as the source material states it — see §6.5.

---

## 2. Background & Novelty

### 2.1 Standard QPE and Shor's reduction

For a unitary $U$ with eigenstate $|u\rangle$, $U|u\rangle = e^{2\pi i
\varphi}|u\rangle$. Quantum phase estimation (QPE) estimates $\varphi$ by
preparing an $n$-qubit phase register in $|+\rangle^{\otimes n}$, applying
$U^{2^j}$ controlled by phase qubit $j$ ($j = 0,\dots,n-1$), and applying
the inverse QFT; measurement yields an $n$-bit integer $y$ with $\varphi
\approx y/2^n$.

Shor's order-finding reduction uses $U_a|y\rangle = |ay \bmod N\rangle$,
whose eigenphases are $\varphi_s = s/r$ for multiplicative order
$r = \mathrm{ord}_N(a)$. Starting QPE from $|1\rangle$ (an equal
superposition of the eigenstates $|u_s\rangle$) and applying continued
fractions to the measured phase recovers $r$; factors then follow from
$\gcd(a^{r/2} \pm 1, N)$.

`Circuits/modular_multiplication.py` builds the $2^{n_w}\times 2^{n_w}$
permutation matrix for $U_a$ ($n_w = \lceil\log_2 N\rceil$ work qubits).
`controlled_power_gate(exponent)` computes $U_a^{2^k}$ not by literal matrix
powering but by exponentiating the *scalar* first — `pow(a, exponent, N)` —
then building the permutation for that single reduced base, which is exact
since $U_a^k = U_{a^k \bmod N}$. `Circuits/iqft.py` implements the inverse
QFT explicitly (controlled phase rotations + Hadamards), verified against
the analytic matrix $(\mathrm{IQFT})_{row,col} = \omega^{row\cdot
col}/\sqrt{2^n}$. `Circuits/qpe.py` assembles the standard, full-register
version of this circuit for reference and testing.

### 2.2 The windowed-circuit novelty (quantum side, reproduced from the source papers)

`Circuits/windowed_qpe.py` implements the paper's core resource-reduction
idea: instead of one $n$-qubit phase register, the phase interval is
partitioned into (possibly overlapping) contiguous **windows**.
`make_overlapping_windows(total_precision n, window_size w, overlap o)`
produces windows advancing by `step = w - o`. `build_windowed_qpe_circuit`
builds a **block-local** circuit that allocates only `spec.width` phase
qubits — never $n$ — with phase qubit $j$ within the block controlled by
$U_a^{2^{\text{start}+j}}$ via the same modular-exponent-reduction trick as
above, followed by that block's own `width`-qubit IQFT. Circuit width and
depth are therefore decoupled from the total requested precision $n$ — the
paper's headline resource claim, which §6.4's structural measurement
confirms is *true* in isolation.

### 2.3 The classical reconstruction novelty (also reproduced from the source papers)

Because each window is measured independently, the classical stage must
reassemble a full-precision phase estimate from $W$ partial, possibly
ambiguous local measurements:

- **Candidate generation** (`Reconstruction/candidate_generation.py`) ranks
  each window's observed bit patterns by empirical probability and keeps
  the top `candidate_count`, encoding each as a `WindowCandidate` at its
  correct absolute bit position.
- **Carry-aware overlap compatibility** (`Reconstruction/carry.py`) checks
  whether two adjacent windows' overlapping bits are consistent, tolerating
  a single legitimate carry/borrow across the block boundary — this is
  "Algorithm 3" of arXiv:2509.05010.
- **Beam-search stitching** (`Reconstruction/stitching.py`) extends
  surviving paths across all $W$ windows, scoring by the product of
  candidate probabilities and keeping only the top `max_paths` at each
  level (a bounded, not exhaustive, search).
- **Continued-fraction order recovery** (`Reconstruction/continued_fraction.py`)
  applies the standard Shor post-processing step to each stitched candidate
  phase in ranked order until one yields a valid, non-trivial factorization.

### 2.4 How the three `Algorithms/` variants relate

- `standard_shor.py` — a classical/textbook reference: exact
  order-finding by brute-force search, bypassing quantum sampling for its
  actual result path (it can also build a single monolithic full-width QPE
  circuit as a reference artifact, but does not use it for factoring).
- `paper_algorithm.py` — the actual windowed-QPE + carry-aware stitching +
  candidate reconstruction pipeline; the genuine reproduction of the source
  papers, sampling real Aer circuits with no shortcuts.
- `window_policy.py` (formerly `adaptive_algorithm.py`) — a thin, deterministic policy layer that derives
  `WindowConfig` (window width, overlap) from requested precision via a
  fixed formula; it does **not** touch the reconstruction pipeline itself
  and is explicitly scoped as a parameter-selection extension point, not
  the confidence-driven adaptivity described next.

### 2.5 The project's own contribution: confidence-guided adaptive reconstruction

A literature-review pass recorded in `Data/Theory_&_Research_Progress.md`
concluded:

> "Multi-candidate reconstruction already exists. Carry-aware stitching
> already exists. Beam-search-like pruning already exists. Static overlap
> already exists. Adaptive overlap has not been proposed. **Continuous
> confidence-guided adaptive reconstruction appears to be novel.**"

This is the project's explicit novelty claim — not against standard Shor's
algorithm or QPE directly, but against the published Modular-Shor/AWQPE
literature and the existing family of Bayesian/ML/DP/beam/graph-search
reconstruction methods it surveyed. The claim is scoped narrowly: what is
new is *making the classical reconstruction pipeline's own resource
allocation (candidate count, beam width, shot count) a function of a
computed per-window confidence statistic*, rather than fixed constants —
while deliberately leaving the quantum circuit and its resource profile
untouched (§2.2, §2.3 are reproduced, not modified).

This confidence-guided framework — three surviving modules, B (candidate
coverage), C (beam width), D (shot stopping) — is specified formally in §3,
proven correct to the extent the source documents establish in §4, analyzed
for cost in §5, and validated experimentally in §6. It was motivated
directly by the failure-mode campaign summarized in §6: the dominant
*natural* failure mode found there was combined classical-resource
starvation (§6.2), not noise or an algorithmic defect, which is exactly
what B/C/D target.

---

## 3. Mathematical Models

*Source: `FORMULATION.md` (specification) and `THEORETICAL_ANALYSIS.md`
(formal treatment). Notation is unified across both.*

### 3.1 Formal problem definition

An instance is a modulus $N$, base $a$ coprime to $N$ with order $r =
\mathrm{ord}_N(a)$, and a windowed configuration $(\ell, s, o)$ — total
precision $\ell$, window width $s$, overlap $o$ — from which
`make_overlapping_windows` deterministically derives $W$ window
specifications, each covering bit positions $[\mathrm{start}_w,
\mathrm{start}_w + \mathrm{width}_w)$ of an $\ell$-bit phase register, with
$\mathrm{start}_{w+1} = \mathrm{start}_w + \mathrm{width}_w - o$.

Each window $w$ is sampled independently, producing a count dictionary
$\mathrm{counts}_w : \{0,1\}^{\mathrm{width}_w} \to \mathbb Z_{\ge0}$ with
$S_w = \sum_b \mathrm{counts}_w(b)$ shots.

**The reconstruction problem**: given $\{\mathrm{counts}_w\}_{w=1}^W$,
produce $\hat\phi \in [0,1)$, $2^\ell\hat\phi \in \mathbb Z$, such that
continued-fraction order recovery applied to $\hat\phi$ yields the true
order $r$ and a non-trivial factorization of $N$.

The baseline (published) pipeline solves this by: (i) truncating each
window's ranked counts to a fixed `candidate_count` $k$; (ii) beam-searching
globally consistent paths under a fixed `max_paths` bound $P$, via the
carry-aware overlap check; (iii) running continued-fraction order search on
each stitched phase, in ranked order, until one succeeds. The three
adaptive modules below replace the fixed constants $k$, $P$, and a fixed
$S_w$ with data-dependent quantities computed from the already-sampled
counts, **without altering the structure of steps (ii)–(iii)** — no circuit
is ever regenerated by any adaptive module, so the framework is a strict
*post-processing* refinement, not a new algorithm.

### 3.2 Notation

| Symbol | Meaning |
|---|---|
| $w$ | window index, $1 \le w \le W$ |
| $S_w$ | shots collected for window $w$ |
| $b^{(1)}_w, b^{(2)}_w$ | top-1, top-2 observed bit patterns for window $w$ |
| $n_1, n_2$ | raw observed counts of $b^{(1)}_w, b^{(2)}_w$ |
| $p_w(b)$ | unknown true sampling probability of pattern $b$ |
| $\alpha_0$ | Beta/Dirichlet smoothing constant, fixed at $0.5$ (Jeffreys) |
| $C_w$ | the confidence statistic (§3.4) |
| $\ell_w$ | window width in bits |
| $K_w$ | number of *distinct* observed patterns in window $w$, $K_w \le \min(2^{\ell_w}, S_w)$ |
| $k$ | fixed candidate count (baseline) |
| $m_w$ | data-dependent candidate count under Module B |
| $P$ / $P'$ | beam width, baseline / adaptive (Module C) |
| $\delta_{\mathrm{cov}}$ | candidate-coverage slack target |
| $\varepsilon$ | shot-stopping error target |
| $S_{\max}$, `beam_max` | hard engineering caps |
| $o$ (`shared`) | overlap width between adjacent windows |
| $R_w$ | number of resampling rounds Module D executes on window $w$ |

### 3.3 Candidate probability model

The empirical distribution $\hat p_w(b) = \mathrm{counts}_w(b)/S_w$ is the
MLE of $p_w(b)$; no smoothing is applied to the *ranking* used by candidate
generation. Smoothing enters only inside $C_w$'s posterior model, and its
effect on ranking is bounded:

$$
|\tilde p_w(b) - \hat p_w(b)| \;\le\; \frac{\alpha_0}{\alpha_0\cdot 2^{\ell_w} + S_w} = O\!\left(\frac{1}{2^{\ell_w}}\right) \text{ for } S_w \gg 2^{\ell_w}
$$

At the smallest tested shot count ($S_w=4$) and $\ell_w \le 6$, this bound
is $\approx 0.015$ — judged negligible, which is why the implementation
ranks by raw $\hat p_w$ rather than the smoothed $\tilde p_w$. (Proof in §4.3.)

### 3.4 The confidence statistic $C_w$ — the framework's one new piece of math

$$
C_w \;=\; \Pr\big[p_w(b^{(1)}_w) > p_w(b^{(2)}_w) \mid n_1, n_2, S_w\big]
$$

computed as an exact two-proportion Bayesian comparison under independent
posteriors $X_i \sim \mathrm{Beta}(\alpha_0+n_i,\ \alpha_0+S_w-n_i)$,
$i\in\{1,2\}$, $\alpha_0=0.5$:

$$
C_w = \Pr[X_1 > X_2] = \int_0^1 f_{X_1}(x)\, F_{X_2}(x)\, dx
$$

evaluated by adaptive quadrature (`scipy.integrate.quad`) — "a single,
unconditional code path, no normal-approximation shortcut, no Monte Carlo,
no branching by shot count or regime" (`FORMULATION.md` §3).

**Why this specific statistic.** Raw probability, margin, or entropy over
observed counts are sample-size-blind: a $2/4$ split and a $500/1000$ split
score identically under any of them, even though the first is nearly
uninformative and the second nearly certain. $C_w$ is presented as the
minimal quantity that distinguishes these two cases from a single exact
formula.

**Why no closed form is used.** For general Beta parameters, $\Pr[X_1>X_2]$
has no elementary closed form (it reduces to a Gauss hypergeometric
${}_3F_2$ at 1). A known finite-sum identity exists when one shape
parameter is a positive *integer*, but here $b_i = 0.5 + (S_w - n_i)$ is
always a half-integer, so that identity does not apply without a separate
derivation. The implementation instead evaluates the integral directly by
quadrature — exact up to floating-point/solver tolerance, no case-specific
formula required.

**Stated simplifying assumption (not experimentally validated, a
documented limitation):** only the top-2 observed patterns are compared; a
genuine 3-way near-tie is not flagged as low-confidence. Accepted given
window widths $\ell_w \le 6$ in every experiment reported and the
Dirichlet-kernel spectral concentration of QPE outcomes.

### 3.5 Module B — candidate coverage

Replaces the fixed `candidate_count` truncation with a coverage rule: take
ranked candidates in descending order until cumulative probability mass
reaches a target:

$$
\sum_{i=1}^{m_w} \hat p_w(b_{(i)}) \;\ge\; 1-\delta_{\mathrm{cov}}, \qquad m_w \ge 1
$$

No posterior smoothing is applied before ranking (§3.3's bound shows this
is negligible at the tested shot counts).

### 3.6 Module C — beam width

$$
\mathrm{expanded} = \#\Big\{ w : \sum_{i=1}^{k} \hat p_w(b_{(i)}) < 1-\delta_{\mathrm{cov}} \Big\}
\qquad
P' = \min\big(\texttt{beam\_max},\ P\cdot(1+\mathrm{expanded})\big)
$$

`expanded` is computed against the fixed baseline $k$, independent of
whether Module B is enabled — this is what keeps B and C separately
ablatable.

### 3.7 Module D — shot stopping

For each window, resample in bounded increments up to a hard cap
$S_{\max}$, stopping early once $C_w \ge 1-\varepsilon$; each round adds
$\min(S_w, S_{\max}-S_w)$ shots from the same operator/eigenstate and
merges counts. This directly targets shot-limited failures (§6.2, §6.5).

### 3.8 A retracted fourth module

The original design also specified "Module A," a proposed fix to the carry
check, motivated by an initially reported 55.4% false-accept rate. During
implementation, re-validating that evidence before writing code surfaced a
bug in the *audit script* that produced the 55.4% figure (§6.5 gives the
full story). With the bug fixed, the carry check was found to already be an
exact classifier — **Module A was therefore removed from the design before
any code was written**, and `Reconstruction/carry.py` required no change.

### 3.9 Assumptions

Stated explicitly because every proposition in §4 is conditional on them,
and one is known to be measurably false:

- **(A1) Per-window multinomial sampling.** $\mathrm{counts}_w \sim
  \mathrm{Multinomial}(S_w, p_w(\cdot))$, independent across windows.
  Undisputed by any experiment.
- **(A2) Independence of the top-2 marginals** (used only by $C_w$). Treats
  $b^{(1)}_w, b^{(2)}_w$'s counts as independent Binomials rather than
  correlated coordinates of one Dirichlet posterior. **Proven false** by
  the calibration study (§6.6): raw Brier score $0.2239$ vs. a naive
  always-$0.5$ baseline of $0.25$. Retained anyway as a deliberate,
  documented scope decision.
- **(A3) Top-2 sufficiency.** Not experimentally validated; accepted given
  $\ell_w \le 6$ in the evaluated corpus.
- **(A4) No cross-window coupling.** Each module allocates budget per
  window independently (Module C aggregates a per-window *indicator*
  across windows into one scalar, but does not otherwise model inter-window
  dependence).
- **(A5) Static circuit geometry.** Window start/width/overlap are fixed at
  configuration time and never revised by confidence — all adaptivity
  operates strictly on already-sampled classical data; overlap adaptivity
  is explicitly out of scope (it would require regenerating circuits,
  crossing the classical/quantum boundary the project maintains).

### 3.10 Free parameters (complete list)

| Parameter | Value | Source |
|---|---|---|
| $\alpha_0$ | $0.5$ | Fixed convention (Jeffreys); not tuned |
| $\varepsilon$ | $0.05$ | Used as-is in the executed ablation; not cross-validated |
| $\delta_{\mathrm{cov}}$ | $0.05$ | Used as-is; not cross-validated |
| $S_{\max}$ | $64$ | Reduced from an initial $65536$ after measuring $\approx3.3\,$s/round; keeps 4×–16× headroom over the grid's starved shots (4–16) |
| `beam_max` | $4096$ | From the measured Phase 5 blowup risk (§6.5) |

No other free parameters exist in the design; everything else (top-2-only
comparison, no cross-window allocation, no regime-dependent confidence
computation, no pre-built recalibration) is a fixed structural decision.

---

## 4. Proofs & Correctness Arguments

*Source: `CORRECTNESS_ARGUMENTS.md` and `THEORETICAL_ANALYSIS.md`. Claims
are labeled exactly as the source documents label them — upgrading a
heuristic or empirical claim to a proof would misrepresent the evidence.*

### 4.1 What "correctness" means for this pipeline

The pipeline is correct on a trial if its returned phase estimate
$\hat\phi$ satisfies `recover_order_and_factors(\hat\phi, a, N).factors ==`
the true factorization of $N$ — an end-to-end, black-box notion.

**Lemma 0 (Sufficient decomposition) *(PROVEN)*.** If (i) every window's
top-ranked retained candidate equals the true window bit pattern, (ii) the
carry check accepts exactly the carry-consistent pairs among those true
patterns, (iii) beam pruning never discards the resulting fully-true path,
and (iv) continued-fraction expansion recovers the true order, then the
pipeline returns the correct factorization.

*Proof.* Under (i)–(ii), the true patterns are pairwise compatible across
every overlap, so a fully-true path exists at every stitching level. Under
(iii) it survives pruning; since the true pattern is top-ranked in every
window, it receives the maximal product score and is stitched first. Under
(iv), continued-fraction recovery on it succeeds first, so the outer loop
returns on this first iteration. $\blacksquare$

This decomposition is a logical convenience, **not** a claim that (i)–(iv)
are necessary — a correct answer can also arise from a non-top-ranked path,
a legitimate carry-borrow substituting for exact equality, or a later
`stitched` candidate succeeding where an earlier one fails.

### 4.2 The carry-check is an exact classifier

**Lemma 1 *(PROVEN)*.** For overlap width $o = \mathrm{overlap\_width}(\mathrm{left},\mathrm{right}) > 0$,
`candidates_compatible(left, right)` returns `True` iff

$$
\mathrm{tail} = \mathrm{head} \quad\text{or}\quad (\mathrm{int}(\mathrm{tail},2) - c) \bmod 2^{o} = \mathrm{int}(\mathrm{head},2),
$$

where $c$ is the carry bit just past the shared overlap.

*Proof.* This is the function's literal implementation
(`Reconstruction/carry.py:62-76`): compute `tail`, `head`; return `True`
immediately on exact equality, else return the stated modular equation.
There is no other branch. The corrected audit's "ground-truth oracle" is
defined by the *same* predicate applied to the same `(tail, head, c)`
triple — so the implementation and the oracle are the same Boolean function
by construction, not merely observed to agree. $\blacksquare$

**Corollary (0% false-accept, 0% false-reject) *(PROVEN, corroborated
EMPIRICALLY)*.** Because the two formulas are algebraically identical, this
generalizes to *every* `(tail, head, c)` triple over any overlap width, not
only the audited range — stronger than the 168-combination exhaustive audit
alone would establish (42 true-accept, 126 true-reject, 0 false-accept, 0
false-reject; see §6.5 for how this superseded the original 55.4% figure).

**What this lemma does not say.** It says nothing about whether the *true*
bit patterns ever reach the check in the first place — that depends on
candidate generation and shot sufficiency, not the check itself. The
overlap-corruption experiment's 43.75% success rate (§6.3) is *consistent*
with this lemma: the corrupted candidate is correctly rejected, and the
resulting failures trace to beam collapse (budget starvation), not a
compatibility-check defect.

### 4.3 Candidate generation is always non-empty and correctly encoded

**Invariant 1 *(PROVEN)*.** For any non-empty `counts`, both candidate
generation functions return $\ge 1$ candidate, each correctly encoded at
its absolute bit position (`value = local << suffix_width`).

*Proof (length ≥ 1).* Fixed-count generation validates `candidate_count >=
1` and slices a non-empty ranked list. Coverage-based generation's break
condition (`cumulative >= target and candidates`) cannot trigger on the
first iteration (`candidates` is empty), so at least one candidate is
always appended first. $\blacksquare$

*Proof (encoding).* Both functions shift the local pattern by
`total_precision - start - width`, matching the placement `stitch_candidates`
assumes when merging — direct transcription, no further arithmetic.
$\blacksquare$

### 4.4 Module B is minimal-prefix-correct

**Proposition B1 / Proposition 5 *(PROVEN)*.** `generate_window_candidates_by_coverage`
returns the shortest descending-probability prefix whose cumulative
empirical mass is $\ge 1-\delta_{\mathrm{cov}}$.

*Proof.* The loop appends candidates in descending order and checks the
break condition only after the append; since $\hat p_w \ge 0$ and the list
is sorted descending, partial sums are non-decreasing in $m$, so the first
$m$ meeting the target is the unique minimal such $m$. $\blacksquare$

**What this does *not* prove *(explicit limitation)*.** Minimality of the
*empirical* coverage set says nothing about whether the *true* top-1
pattern is among the selected candidates — that would require a
finite-sample concentration bound (e.g. Chernoff or DKW-type), which is
neither derived nor implemented. $\delta_{\mathrm{cov}}$ controls sample
coverage, not a calibrated population error rate.

### 4.5 Beam stitching is exact only when unsaturated

**Proposition S1 (Conditional exactness) *(PROVEN)*.** If at every level of
`stitch_candidates` the number of compatible next-paths generated (before
truncation) is $\le$ `max_paths`, the search returns the true global
optimum (equivalent to exhaustive search).

*Proof.* If `len(next_paths) <= max_paths` at every level, the truncation
`[:max_paths]` is a no-op, so `paths` remains exactly the full surviving
set at each level; by induction, the final list contains all globally
consistent paths, sorted, with the true maximum first. $\blacksquare$

**Corollary (lossy when saturated) *(PROVEN)*.** Whenever the count exceeds
`max_paths` at some level, discarded paths are permanently unrecoverable —
no backtracking exists in the loop structure.

**Interpretation for B/C.** Module C's Proposition 6 ($P' \ge P$ whenever
`beam_max` $\ge P$) can only reduce, never increase, the risk of this
discard case relative to baseline — but *whether* this suffices on any
given trial is an empirical question (§6.5), not something Proposition S1
alone answers.

### 4.6 Properties of $C_w$

**Proposition 1 (Symmetry) *(PROVEN)*.** $C_w(n_1,n_2,S_w) + C_w(n_2,n_1,S_w) = 1$,
since $X_1, X_2$ independent and continuous $\Rightarrow \Pr[X_1{=}X_2]=0$.

**Proposition 2 (Boundary value) *(PROVEN)*.** $C_w(n,n,S_w) = 1/2$,
immediate from Proposition 1.

**Proposition 3 (Monotonicity) *(PROVEN)*.** $C_w$ is non-decreasing in
$n_1$ for fixed $n_2, S_w$ (and correspondingly non-increasing in $n_2$),
by the Beta family's monotone-likelihood-ratio property and preservation of
stochastic dominance under $\Pr[\cdot > X_2]$.

**Proposition 4 (Asymptotic consistency) *(HEURISTIC, proof sketch)*.** As
$n_1,n_2,S_w \to \infty$ with $n_1/S_w \to p_1 > p_2 \leftarrow n_2/S_w$,
$C_w \to 1$ in probability — via the Bernstein–von Mises normal
approximation of the Beta posterior and Slutsky's theorem. Labeled (H), not
(T), because the normal approximation used is itself asymptotic.

### 4.7 Shot-stopping: proven termination, explicitly unproven calibration

**Lemma D1 (Termination) *(PROVEN)*.** The resampling loop for window $w$
terminates in at most $R_w \le \lceil \log_2(S_{\max}/S_w^{(0)}) \rceil$
rounds.

*Proof.* While $S_w \le S_{\max}/2$, the added amount equals $S_w$ itself,
so $S_w$ at least doubles each round; once $S_w > S_{\max}/2$, it is
strictly increasing and bounded above by $S_{\max}$, forcing termination
within one further round beyond the doubling phase. $\blacksquare$ At the
configured defaults ($S_w^{(0)}=4$, $S_{\max}=64$), this is exactly 4
rounds — matching `config.py`'s own comment.

**Lemma D2 (Stopping condition is relative, not absolute) *(PROVEN, with an
EMPIRICAL sharpening)*.** The stopping check uses *raw*, uncalibrated
$C_w$ — the isotonic recalibration model is never imported by the
harness. Combined with the calibration result (top decile: mean predicted
$C_w=0.96$, observed correctness $0.79$; §6.6), a window halting at
$\varepsilon=0.05$ is empirically correct only $\approx79\%$ of the time,
not $\ge95\%$. **$\varepsilon$ is a tuning parameter for a monotone
heuristic, not a statistically guaranteed error tolerance**, despite its
name — this does not contradict Module D's large measured effect (§6.5): a
miscalibrated but monotonically informative statistic can still be a highly
effective *relative* ranking signal.

### 4.8 The framework is a strict generalization of baseline

**Invariant 2 / Proposition 8 (Reduction to baseline) *(PROVEN)*.** Calling
the adaptive reconstruction entry point with all three flags `False`
executes a call sequence structurally identical to the pre-adaptive
baseline pipeline.

*Proof.* Direct inspection: each of the three `if use_*:` blocks has a
complete `if/else`, and every `else` branch reproduces the pre-adaptive
computation exactly (same calls, same arguments, sourced from the fixed
`WindowConfig` fields). $\blacksquare$

This licenses treating the observed Baseline arm of the ablation (§6.5:
52.1% / 51.1% across independent runs) as a direct, code-level reproduction
of the original failure study's 52.1% figure — corroborating, not proving,
the invariant, but a meaningful sanity check.

**Proposition 9 (Monotone non-restrictiveness) *(PROVEN, conditional)*.**
Under each module's stated hypotheses, B, C, and D can each only *relax* a
resource constraint relative to baseline, never tighten one — explaining
why no ablation arm underperforms baseline, though effect *sizes* remain
empirical.

### 4.9 Summary table

| Claim | Status | Location |
|---|---|---|
| Carry check is an exact classifier | **PROVEN** (algebraic identity) + corroborated by exhaustive audit | §4.2 |
| Candidate generation returns $\ge1$ correctly-encoded candidate | **PROVEN** | §4.3 |
| Module B returns the minimal covering prefix | **PROVEN** | §4.4 |
| Module B's prefix contains the *true* pattern | **not claimed / no guarantee exists** | §4.4 |
| Beam stitching is exact when unsaturated | **PROVEN** (conditional) | §4.5 |
| Beam stitching can permanently lose the true path when saturated | **PROVEN** | §4.5 |
| $C_w$ symmetry, boundary value, monotonicity | **PROVEN** | §4.6 |
| $C_w$ asymptotic consistency | **HEURISTIC** (proof sketch) | §4.6 |
| Shot-stopping loop terminates | **PROVEN**, bound $O(\log(S_{\max}/S_w^{(0)}))$ | §4.7 |
| Shot-stopping threshold is a calibrated error rate | **false** — monotone but uncalibrated | §4.7 |
| All-flags-off reproduces baseline exactly | **PROVEN** (code structure) | §4.8 |
| B/C/D each improve success rate, and by how much | **EMPIRICAL** only | §6.5 |
| $C_w$ is measurably overconfident | **EMPIRICAL** only (mechanistically consistent with (A2)'s falsity) | §6.6 |

---

## 5. Complexity Analysis

*Source: `COMPLEXITY_ANALYSIS.md`. All derivations are against the actual
code, not a hypothetical alternative implementation.*

### 5.1 Baseline complexity

$$
T_{\text{baseline}} = O\Big(\underbrace{\textstyle\sum_w K_w\log K_w}_{\text{candidate gen}} \;+\; \underbrace{W\cdot P\cdot k\cdot(o+\log(Pk))}_{\text{stitching}} \;+\; \underbrace{P\cdot N\log N\cdot M(\log N)}_{\text{order recovery}}\Big)
$$

$$
S_{\text{baseline}} = O\big(\textstyle\sum_w K_w + P\cdot k\cdot W\big)
$$

The order-recovery term is the one place baseline complexity is not
polynomial in circuit parameters ($W,k,P$) but linear in the arithmetic
parameter $N$ directly — the search is a brute-force multiplier scan, not
an order-finding shortcut. Empirically (`REPORT.md` §2), for the small
orders tested this search is extremely forgiving and exits early — not a
formally derived bound, but a reported **(E)** observation.

### 5.2 Per-module cost

- **Module B** — same asymptotic class as baseline candidate generation
  ($O(K_w \log K_w)$ time); the truncation constant becomes data-dependent
  ($m_w$) instead of fixed ($k$). Worst case equals baseline's worst case
  (a flat distribution forces $m_w \to K_w$); best case ($m_w=1$ on a
  near-deterministic window) is strictly cheaper. B's measured benefit is a
  *statistical*, not computational-complexity, improvement — its cost is
  paid downstream in stitching.
- **Module C** — beam width becomes a bounded random variable $P' \in [P,
  \texttt{beam\_max}]$; stitching cost can grow by up to a factor of
  `beam_max`$/P$ (up to 128× at evaluated defaults) relative to baseline,
  paid only on trials where `expanded > 0`.
- **Module D** — adds $O(\log_2(S_{\max}/S_w^{(0)}))$ additional circuit
  executions per window in the worst case (Lemma D1, §4.7), each costing
  one real Aer transpile+run ($\approx3.3\,$s measured) — **the single most
  expensive operation in the entire adaptive framework in absolute
  wall-clock terms**, more expensive per unit than any purely classical
  step in B or C.

### 5.3 Aggregate complexity

$$
T_{\text{full}} = O\Big(\textstyle\sum_w K_w\log K_w \;+\; W\cdot P'\cdot m_w\cdot(o+\log(P'm_w)) \;+\; P'\cdot N\log N\cdot M(\log N) \;+\; W\log_2(S_{\max}/S_{\min}^{(0)})\cdot T_{\mathrm{sample}}\Big)
$$

$$
S_{\text{full}} = O\Big(\textstyle\sum_w m_w \;+\; P'\cdot m_w\cdot W \;+\; W\cdot S_{\max}\Big)
$$

Every term is bounded by an explicit engineering cap carried over from the
evidence base ($m_w \le K_w$, $P' \le$`beam_max`$=4096$, $R_w \le
\lceil\log_2(S_{\max}/S_w^{(0)})\rceil$, $S_{\max}=64$). The framework is
**worst-case bounded but not worst-case cheap**: no module introduces a new
complexity class relative to baseline; the cost is a bounded multiplicative
inflation, with the bound set by the three free parameters in §3.10.

### 5.4 What actually dominates in practice (E)

By direct wall-clock reasoning: circuit execution ($\approx3.3\,$s/round)
is orders of magnitude more expensive per operation than any purely
classical step (sorting a $\le16$-entry dict, a handful of $O(o{=}3)$
compatibility checks, a bounded-iteration quadrature call). **Module D's
resampling is therefore the dominant wall-clock cost whenever it
triggers** — even though B and C carry the larger asymptotic-complexity
footnotes on paper (exponential-in-$\ell_w$ worst case, combinatorial beam
blowup risk).

---

## 6. Experimental Validation

*Source: `REPORT.md` (Phases 1–5) and `PHASE6_SUMMARY.md` (Phase 6). All
experiments drive the pipeline through its public API only ("black box");
no repository source was modified for Phases 1–5.*

### 6.0 Method and scope

Six phases were run:

1. **OFAT sensitivity sweeps** around a reference point (`total_precision=8,
   window_size=4, overlap=2, candidate_count=2, max_paths=32, shots=2048,
   noise=0`), plus a combined worst-case stress test.
2. **Noise sweeps** — single-source and combined (depolarizing 1q/2q,
   readout error).
3. **Adversarial fault injection** on the classical layer, plus an
   exhaustive carry-check audit and a measurement-ambiguity boundary sweep.
4. **Scaling** — runtime/memory/qubit-count vs. problem size.
5. **Root-cause attribution** — failed trials replayed through a fixed
   counterfactual ladder.
6. **Confidence-guided adaptive reconstruction** — implementation, ablation,
   and calibration of Modules B/C/D.

An unplanned scope-shaping finding governs everything else: work-register
sizes as small as **7 qubits (N≈65–119)** took tens of seconds to minutes
to *compile*, and **8 qubits (N≈143) did not finish transpiling within
45s** — even though the phase-side circuit never exceeds `window_size`
qubits. This forced Phases 1–3 to stay at work_qubits ≤ 5 (N ≤ 21), a real
limitation on statistical power for larger instances, reported honestly
rather than concealed.

### 6.1 Phase 1 — parameter sensitivity

**435/435 OFAT trials succeeded (100%)** across every value of
`total_precision` (4–12), `window_size` (2–8), `overlap` (0–3),
`candidate_count` (1–16), `max_paths` (1–4096), and `shots` (16–16384), for
three small (N,a) pairs. This is a genuine null result, explained by the
continued-fraction multiplier search being extremely forgiving at small
orders (2 or 6) — it is not evidence of general robustness.

**Combined-stress follow-up (the pivotal experiment):** simultaneously
starving `candidate_count∈{1,2}`, `max_paths∈{1,2}`, `shots∈{4,8,16}` while
maximizing `total_precision∈{16,24,32}` and `overlap∈{0,3}` on the same
trivial instance (N=21, a=19, order 6), 144 trials: **69/144 failed (47.9%
failure, 52.1% success)** — with zero noise or corruption. No single factor
in isolation could cause this; several starved together can.

### 6.2 Phase 2 — noise sensitivity

Single-source sweep (252 trials, 0–30% error): even 30% single-qubit
depolarizing error only drops success to 83.3%; 2-qubit depolarizing and
readout error barely register up to 30% at this instance size. Combined
3×3×3 factorial (243 trials): **91.8% overall**, worst cells 77–78%.
Logistic regression ranks 2-qubit depolarizing as the strongest driver
(coefficient −0.45). Observed combined rates *exceed* the independence
baseline in every tested cell — noise sources are mildly **sub-additive**
here, not super-additive.

### 6.3 Phase 3 — adversarial robustness

| Condition | Success rate |
|---|---|
| clean (no corruption) | 100% (32/32) |
| bit-flip, 10% of shots relabeled | 100% (32/32) |
| fabricated high-count candidate injected | 100% (32/32) |
| true candidate removed from a window | 78.1% (25/32) |
| **bit flipped inside the overlap region** | **43.75% (14/32)** |

Overlap corruption is the single most damaging manipulation tested. An
80-trial instrumented replay confirms the corrupted candidate is *always*
correctly rejected by the carry check (§4.2, §6.5) — the resulting failures
come from beam collapse to zero surviving paths, not the check being
fooled.

Ambiguity-boundary sweep (single-qubit `PhaseGate`, controlled distance
from the 0.5 bucket boundary): within delta ≤0.005 of the boundary, the
winning bucket is **completely unstable** across independent repeats
(100% winner-disagreement); by delta≥0.02, fully stable. This implementation
has no special-chunk/ambiguity-resolution mechanism (unlike AWQPE Algorithm
2's `Sidx` correction).

### 6.4 Phase 4 — scaling and the compilation cliff

| work qubits | N | transpile time |
|---|---|---|
| 4 | 15 | 0.14 s |
| 5 | 21 | 0.38 s |
| 6 | 33 | 2.09 s |
| 7 | 65 | **42.2 s** |
| 8 | 143 | **>45 s (timeout)** |

This happens within a single window circuit that never exceeds `window_size
+ work_qubits` = 12 qubits total — the compilation cliff is driven entirely
by `work_qubits`, via `ModularMultiplicationOperator` building a dense
`UnitaryGate` whose default unitarity check and Qiskit's generic unitary
synthesis both scale explosively with work-register dimension.

A structural, gate-construction-free accounting confirms `max_circuit_qubits
= window_size + work_qubits`, growing only as $\lceil\log_2 N\rceil$ —
**the paper's phase-register decoupling claim is structurally true.** The
problem is entirely downstream, in how this reference implementation
realizes the arithmetic operator. Runtime vs. `total_precision` grows
sub-linearly (1.53s → 6.92s, precision 4→32); runtime vs. `window_size`
actually *decreases* (5.43s → 1.78s, size 2→16), since fewer, wider windows
amortize per-circuit overhead.

### 6.5 Phase 5 — root-cause attribution, and the self-correction

All 94 failed trials from Phases 1–2 were replayed through a counterfactual
ladder (perfect bits → widened beam → exhaustive candidates → noise off →
10× shots):

| failure category | count | share |
|---|---|---|
| quantum_estimation_shot_limited | 48 | 51% |
| noise_sensitivity | 25 | 27% |
| unknown (resists every ladder rung) | 20 | 21% |
| unknown_attribution_timeout | 1 | 1% |

**72% of all naturally observed failures trace to classical-pipeline budget
starvation** (shots + shot/budget interaction), not noise. The 20 "unknown"
cases needed *both* more shots *and* more candidate/beam budget
simultaneously (15/20 at `overlap=3`) — a distinct fifth failure mode
invisible to any single-axis relaxation, and invisible to Phase 1's plain
OFAT sweep by construction.

**The self-correction, in full.** The original Phase 3 carry-check audit
reported: 55.4% false-accept, 5.4% false-reject — a headline finding that
the carry-tolerant overlap check was badly overpermissive. While designing
the follow-on adaptive-reconstruction study, re-validating this evidence
before building on it surfaced a bug in the audit script
(`Experiments/phase3_adversarial.py:run_carry_truth_table`): its synthetic
`WindowCandidate` construction set `right.start = overlap` instead of
`right.start = width - overlap`, which pinned the *actual* computed overlap
to exactly 1 bit regardless of the loop's nominal `overlap` value — `tail`
and `head` were never compared at their declared width. With the geometry
bug fixed and the 168-combination audit re-run:

| Outcome | Original (superseded) | Corrected |
|---|---|---|
| true_accept | 19.6% | 25.0% |
| true_reject | 19.6% | 75.0% |
| false_accept | **55.4%** | **0.0%** |
| false_reject | 5.4% | 0.0% |

`Reconstruction/carry.py` required **no code change** — the audit that
measured it was wrong, not the mechanism (see §4.2's proof that this result
generalizes beyond the audited cases). A follow-up 80-trial instrumented
replay of the overlap-corruption fault-injection condition confirmed this
directly: the corrupted candidate never once survived into a stitched path
when the flip landed in a carry-checked overlap.

Three independent checks establish this is a genuine correction, not a
judgment call: (1) **mathematical** — the carry-correction formula and the
audit's own ground-truth oracle are algebraically the same formula, so the
original test could only ever have detected its own bug, never a defect in
the formula; (2) **empirical re-run** — 0%/0% after fixing the geometry;
(3) **instrumented replay** — 80/80 trials show correct rejection. The
correction was made *before* Module A (a proposed carry-check fix) was
implemented, retiring that module from the design (§3.8) — the original
`REPORT.md` and `FORMULATION.md` were edited in place with visible
strikethrough/correction markup rather than silently rewritten, preserving
an auditable trail of what was wrong and why.

**Ranking of bottlenecks (`REPORT.md` §7.4):** the compilation cliff
dominates on *reach* (a hard, unconditional gate on whether the algorithm
can run at all, past ~150); combined classical-budget starvation dominates
on *observed frequency* (68/94 = 72% of every failure the campaign actually
produced) and, unlike the cliff, is fully within reach of this study's own
instance sizes — making it the most directly actionable finding, and the
one the Phase 6 adaptive framework (B/C/D) was built to address.

### 6.6 Phase 6 — implementing and validating the fix

**Ablation** (`FORMULATION.md` §5 / `PHASE6_SUMMARY.md` §5), reusing the
144-trial combined-stress grid, paired by seed, McNemar's test + bootstrap
95% CI:

| Arm | Success rate | Paired diff vs. Baseline | McNemar's $p$ |
|---|---|---|---|
| Baseline | 52.1% (51.1% in the independent Phase 6 run) | — | — |
| B (coverage) | 77.1% (76.3%) | +25.0 pp | $1.5\times10^{-9}$ |
| C (beam width) | 56.9% (56.1%) | +4.9 pp | $0.016$ |
| D (shot stopping) | 76.4% (76.3%) | +24.3 pp | $1.8\times10^{-7}$ |
| **Full (B+C+D)** | **100.0%** | **+47.9–48.9 pp** | $3.4\times10^{-21}$ (or $6.8\times10^{-21}$) |

B beats C by 20.1 pp ($p=1.5\times10^{-5}$) — candidate coverage is by far
the stronger individual intervention. The coupling identified in §6.5's
attribution (20/69 failures needing both axes relaxed together) is real but
asymmetric: C alone helps only marginally, yet Full's 100% is
super-additive relative to B and C's individual effects. D's effect size
closely matches B's, directly confirming the shot-limitation attribution
(48/69, 51% of §6.5's failures). Baseline's own reproduction of the
original 52.1% figure across three independent measurements is a sanity
check that no unintended behavior change leaked into the baseline arm.

**Calibration check.** On 2,687 per-window observations from 579
recaptured trials, $C_w$ is systematically **overconfident**: raw Brier
score **0.2239** (barely better than the naive always-0.5 baseline of
0.25); the top decile (predicted $C_w \in (0.9,1.0]$, mean 0.96) had
observed correctness of only **0.79**. This is the expected cost of the
independence assumption (A2, §3.9) being false: $b^{(1)}_w$ and
$b^{(2)}_w$'s counts are really correlated Dirichlet marginals, not
independent binomials. Isotonic regression recalibration reduces Brier
score to **0.1988** (top-decile calibrated value 0.82 vs. 0.79 observed).
The raw $C_w$ formula is unchanged in the shipped pipeline; the
recalibration model is a diagnostic artifact, not wired into Module D's
stopping decision by default.

### 6.7 Answers to the five research questions (`REPORT.md` §7, condensed)

1. **Under what conditions does it fail?** Not from ordinary single-factor
   variation at small instances (0/435). It fails from: starving multiple
   classical budgets at once (47.9%); corrupting a bit inside an overlap
   region (56.25%); landing within ~1% of a 0.5 ambiguity boundary; any
   instance whose work register needs ≳7–8 qubits at all; compounding
   multiple noise channels together.
2. **Why does it fail?** Combined-stress failures split into plain
   shot-limitation (48/69) and a compounded shot/budget interaction (20/69)
   with nowhere to route around an occasional wrong-window measurement.
   Overlap corruption fails not because the carry check is overpermissive
   (it is not — §4.2) but because losing a window's true candidate leaves
   no alternative path in a narrow beam. The ambiguity instability occurs
   because no special-chunk resolution mechanism exists. The compilation
   cliff occurs because the arithmetic layer is a dense `UnitaryGate`
   rather than a reversible modular-arithmetic circuit.
3. **Which of the paper's assumptions break down?** "Small window blocks
   keep resource requirements low" is true only for the phase register, not
   the work register. The 0.5-ambiguity edge case, treated as improbable in
   the sibling AWQPE paper and given an explicit resolution mechanism
   there, has no analogue in this implementation.
4. **Which bottleneck contributes most?** By natural-failure frequency:
   classical-budget starvation (72%). By reach: the compilation cliff,
   unconditionally blocking essentially all instances beyond ~150.
5. **Which would yield the strongest contribution if solved?** The
   compilation cliff is the strongest *engineering* candidate (undiscussed
   in either source paper, blocks 100% of larger instances, and has an
   existing literature — reversible ripple-carry/CDKM adders — to draw on).
   Combined classical-budget starvation is the strongest candidate for an
   immediate, low-effort *algorithmic* contribution — and is exactly what
   the confidence-guided adaptive framework (§3, §4, this section) was
   built to address, with results now in hand.

---

## 7. Limitations

*Consolidated from `THEORETICAL_ANALYSIS.md` §10.*

1. **(A2)'s independence assumption is measurably false** (§3.9, §6.6) —
   raw Brier $0.2239$, "the single most consequential gap between the model
   and reality in the whole framework."
2. **No PAC-style coverage guarantee for Module B** (§4.4) — $\delta_{\mathrm{cov}}$
   bounds sample coverage, not population/true-pattern coverage.
3. **$\varepsilon$ is not a calibrated error rate** (§4.7) — Module D's
   stopping rule is a useful, monotone, but uncalibrated ranking signal.
4. **Module C's beam-scaling rule is heuristic** (§3.6) — `expanded` is a
   simple per-window ambiguity count, not a derived sufficient beam width.
5. **Top-2-only simplification (A3)** is unvalidated and structurally blind
   to genuine 3-way ties.
6. **No cross-window resource coupling (A4)** — no formal analysis of
   joint, budget-constrained allocation across windows was attempted; no
   experiment in the harness exercises a shared cross-window budget.
7. **All experimental validation of the adaptive framework was conducted at
   small multiplicative order** ($r=6$, $N=21$), within the
   compilation-cliff-imposed work-register ceiling — the continued-fraction
   stage's own forgiveness at small order (§6.1) is a known confound for
   *any* claim about classical-stage robustness, adaptive or not; measured
   gains are not yet validated at an order large enough to remove that
   confound.
8. **Overlap geometry is not adaptive (A5)** — none of the propositions in
   §4 say anything about the 0.5 ambiguity boundary instability (§6.3), a
   distinct, unaddressed failure mode this framework does not target.

---

## 8. Conclusion & Future Work

**Headline conclusions.**

1. The windowed formulation's phase-register resource claim — small,
   fixed-width circuits decoupled from total precision — is **structurally
   true** (§5.1, §6.4) but practically overshadowed by a work-register
   compilation cliff this reference implementation's arithmetic layer
   introduces, which is not addressed by anything in this project's scope.
2. **Combined classical-resource starvation, not noise or an algorithmic
   defect, is the dominant natural failure mode** at the instance sizes
   this study could reach (72% of all observed failures, §6.5) — this is
   the finding that motivated the entire confidence-guided adaptive
   reconstruction framework.
3. That framework — a per-window Bayesian confidence statistic $C_w$
   driving adaptive candidate coverage (B), beam width (C), and shot
   stopping (D) — is proven to be a strict generalization of the baseline
   pipeline (§4.8) with proven termination and correctness properties where
   they exist, and clearly labeled heuristic/empirical claims where they do
   not (§4.9's summary table).
4. On the evaluated grid, **the combination of all three modules eliminates
   every measured combined-stress failure** (52.1%→100%, $p\approx6.8\times10^{-21}$),
   with B (candidate coverage) and D (shot stopping) each independently
   responsible for roughly half the baseline failure rate.
5. A significant self-correction — the carry-check false-accept claim,
   originally 55.4%, traced to an audit-script bug and corrected to
   0%/0% — is preserved and documented rather than smoothed over,
   including its downstream consequence: Module A was retracted from the
   design before implementation began.

**Open items** (from `PHASE6_SUMMARY.md` §7):

1. Finish or formally close out the small number of remaining ablation
   trials cut short for time (results were already stable at the point they
   were stopped).
2. Decide whether the isotonic recalibration layer ships in the live
   pipeline or remains a diagnostic-only artifact — currently, shot-stopping
   decisions in Module D are made on the known-overconfident raw $C_w$, not
   the recalibrated value.
3. Re-check the Phase 1–2 sweep scripts for the same "parameter doesn't
   actually vary the thing being tested" bug pattern that caused the
   carry-check retraction, since only the Phase 3 script was confirmed
   fixed.
4. Push the corrected reports and new Phase 6 code to review — the
   carry-check retraction in particular should be visible to anyone who
   previously relied on the original 55.4% figure.
5. Scale the ablation beyond its current sample size if targeting
   publication, particularly to tighten confidence intervals on Module C's
   comparatively small effect size.
6. Most substantively: validate the adaptive framework's measured gains at
   a multiplicative order large enough to remove the small-order
   continued-fraction forgiveness confound (Limitation 7, §7) — which is
   blocked on solving the compilation cliff (Limitation/Finding 1 above)
   first.
