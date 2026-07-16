# Theoretical Analysis — Confidence-Guided Adaptive Reconstruction

*Companion note to `FORMULATION.md` (design specification), `REPORT.md`
(evidence base), and `PHASE6_SUMMARY.md` (execution record). This document
adds no new algorithmic content and specifies no behavior not already frozen
in those three files; it develops the mathematics underlying the modules
they describe, as implemented in `Reconstruction/confidence.py`,
`Reconstruction/candidate_generation.py`, `Reconstruction/carry.py`,
`Reconstruction/stitching.py`, and wired in `Experiments/harness.py`. Every
symbol below resolves to a specific function or line range; nothing here is
aspirational.*

*Status labels used throughout: **(T)** theorem/proposition proven from the
stated assumptions; **(H)** heuristic argument, plausible but not proven;
**(E)** empirical observation, established by measurement, not derivation.*

---

## 1. Formal problem definition

An instance is a modulus $N$, base $a$ coprime to $N$ with multiplicative
order $r = \mathrm{ord}_N(a)$, and a windowed phase-estimation configuration
$(\ell, s, o)$ — total precision $\ell$ (`total_precision`), window width $s$
(`window_size`), overlap $o$ (`overlap`) — from which
`Circuits.windowed_qpe.make_overlapping_windows` deterministically derives a
sequence of $W$ window specifications $\mathrm{spec}_1,\dots,\mathrm{spec}_W$,
each covering bit positions $[\mathrm{start}_w, \mathrm{start}_w+\mathrm{width}_w)$
of an $\ell$-bit phase register, with $\mathrm{start}_{w+1} = \mathrm{start}_w + \mathrm{width}_w - o$.

For eigenstate $|1\rangle$ of $U_a: |x\rangle \mapsto |ax \bmod N\rangle$, each
window $w$ is sampled independently by a real Aer circuit
(`run_windowed_qpe_block`), producing a count dictionary
$\mathrm{counts}_w : \{0,1\}^{\mathrm{width}_w} \to \mathbb{Z}_{\ge 0}$ with
$S_w = \sum_b \mathrm{counts}_w(b)$ shots.

**The reconstruction problem** is: given $\{\mathrm{counts}_w\}_{w=1}^W$,
produce a full-precision phase estimate $\hat\phi \in [0,1)$, $2^\ell\hat\phi \in \mathbb{Z}$,
such that `Reconstruction.continued_fraction.recover_order_and_factors`
applied to $\hat\phi$ yields the true order $r$ and a non-trivial
factorization of $N$.

The published (baseline) pipeline solves this by: (i) truncating each
window's ranked counts to a fixed `candidate_count` $k$
(`generate_window_candidates`); (ii) beam-searching globally consistent
paths through the window candidates under a fixed `max_paths` bound $P$,
using the carry-aware overlap check (`stitch_candidates`,
`candidates_compatible`); (iii) running the continued-fraction order search
on each stitched phase, in ranked order, until one yields valid factors
(`recover_order_and_factors`).

The three surviving adaptive modules (**B**, **C**, **D**; module A was
retracted — see `FORMULATION.md`'s Correction) replace fixed constants $k$
and $P$, and a fixed $S_w$, with data-dependent quantities computed from the
*already-sampled* counts, without altering steps (ii)–(iii) structurally.
No circuit is ever regenerated or re-specified by any of the three modules
(`FORMULATION.md` §8, scope boundary); this fact is used below to establish
that the adaptive framework is a strict *post-processing* refinement of the
classical stage, not a new algorithm.

---

## 2. Notation

Reused from `FORMULATION.md` §1, restated for self-containedness:

| Symbol | Meaning | Code anchor |
|---|---|---|
| $w$ | window index, $1 \le w \le W$ | `make_overlapping_windows` |
| $S_w$ | shots collected for window $w$ | `sum(counts_w.values())` |
| $b^{(1)}_w, b^{(2)}_w$ | top-1, top-2 observed bit patterns in window $w$ | `normalize_counts` ranking |
| $n_1, n_2$ | raw counts of $b^{(1)}_w, b^{(2)}_w$ | — |
| $p_w(b)$ | unknown true sampling probability of pattern $b$ in window $w$ | — |
| $\alpha_0$ | Beta/Dirichlet smoothing constant, fixed at $0.5$ (Jeffreys) | `top1_vs_top2_confidence` |
| $C_w$ | confidence statistic, §5 below | `top1_vs_top2_confidence` |
| $\ell_w$ | window width in bits | `WindowCandidate.width` |
| $K_w$ | number of *distinct* observed bit patterns in window $w$, $K_w \le \min(2^{\ell_w}, S_w)$ | `len(counts_w)` |
| $k$ | fixed candidate count (baseline) | `WindowConfig.candidate_count` |
| $m_w$ | data-dependent candidate count under module B | `generate_window_candidates_by_coverage` |
| $P$ | beam width / `max_paths` (baseline, fixed) | `WindowConfig.max_paths` |
| $P'$ | adaptive beam width under module C | `reconstruct_adaptive` |
| $\delta_{\mathrm{cov}}$ | candidate-coverage slack target | `AdaptiveReconstructionConfig.delta_cov` |
| $\varepsilon$ | shot-stopping error target | `AdaptiveReconstructionConfig.epsilon` |
| $S_{\max}$ | shot-stopping hard cap | `AdaptiveReconstructionConfig.s_max` |
| `beam_max` | beam-width hard cap | `AdaptiveReconstructionConfig.beam_max` |
| $o$ (`shared`) | overlap width between adjacent windows | `overlap_width` |
| $R_w$ | number of resampling rounds module D executes on window $w$ | `_resample_shot_stopping` |

---

## 3. Assumptions

Made explicit because every proposition below is conditional on them, and
several are known — from the calibration study (`FORMULATION.md` §6) — to be
literally false at the measured level of precision. They are retained
because they are the assumptions the *implementation* encodes, not because
they are asserted to hold exactly.

**(A1) Per-window multinomial sampling.** For fixed true window
probabilities $p_w(\cdot)$, $\mathrm{counts}_w \sim \mathrm{Multinomial}(S_w, p_w(\cdot))$,
independently across windows. This is the standard sampling model for
projective measurement counts from a fixed circuit and is not itself
disputed by any experiment in `REPORT.md`; it is the substrate every
downstream statistic is built on.

**(A2) Independence of the top-2 marginals (used only by $C_w$).** $C_w$
treats $b^{(1)}_w$'s and $b^{(2)}_w$'s counts as two independent Binomial
observations rather than two (negatively correlated) coordinates of one
$K_w$-way multinomial / Dirichlet posterior. **This assumption is proven
false in exactly the sense the calibration study measures** (raw Brier
score $0.2239$ vs. a naive always-$0.5$ baseline of $0.25$; top-decile
predicted confidence $0.96$ against observed correctness $0.79$,
`FORMULATION.md` §6). It is retained in the shipped $C_w$ (module D
thresholds on raw $C_w$, not the isotonic-recalibrated value) as a
deliberate, documented scope decision, not an oversight. All propositions
in §7 that depend on (A2) are marked accordingly, and §9 restates this as a
first-class limitation rather than a footnote.

**(A3) Top-2 sufficiency.** Reconstruction ambiguity for a window is fully
captured by comparing its two highest-count patterns; a genuine three-way
near-tie is not distinguished from a clean two-way split by $C_w$. Not
validated experimentally (`FORMULATION.md` §3); accepted as reasonable
given $\ell_w \le 6$ in the evaluated corpus and the Dirichlet-kernel
spectral concentration of QPE outcomes, but not proven.

**(A4) No cross-window coupling in the adaptive rules.** $C_w$, module B's
coverage cutoff, and module D's stopping decision are all computed
per-window from that window's own counts alone; no joint distribution
across windows is modeled or exploited, and no shared resource budget is
enforced across windows (`FORMULATION.md` §4.3). Module C is the one
exception — it aggregates a per-window *indicator* (needed-expansion) across
windows into a single scalar, but does not otherwise model inter-window
statistical dependence.

**(A5) Static circuit geometry.** $\mathrm{start}_w, \mathrm{width}_w, o$ are
fixed by the initial `WindowConfig` and never revised in response to
confidence (`FORMULATION.md` §8, explicit scope exclusion). All adaptivity
operates strictly on already-sampled classical data.

---

## 4. Candidate probability model

For window $w$, the empirical distribution
$\hat p_w(b) = \mathrm{counts}_w(b)/S_w$ (`normalize_counts`) is the
maximum-likelihood estimate of $p_w(b)$ under (A1), with no smoothing
applied to the *ranking* used by candidate generation (`FORMULATION.md`
§4.1: "No posterior smoothing is applied to the counts before ranking").
Smoothing enters only inside $C_w$'s posterior model (§5), and its effect on
the *ranking* itself is bounded and shown negligible at the tested shot
counts:

**(T, restated from `FORMULATION.md` §4.1).** For $\alpha_0 = 0.5$ and a
window with $2^{\ell_w}$ possible patterns, the Jeffreys-smoothed estimate
$\tilde p_w(b) = (\alpha_0+\mathrm{counts}_w(b))/(2^{\ell_w}\alpha_0+S_w)$
differs from $\hat p_w(b)$ by at most

$$
|\tilde p_w(b) - \hat p_w(b)| \;\le\; \frac{\alpha_0}{\alpha_0\cdot 2^{\ell_w} + S_w} = O\!\left(\frac{1}{2^{\ell_w}}\right) \text{ for } S_w \gg 2^{\ell_w},
$$

*Proof.* $\tilde p_w(b) - \hat p_w(b) = \dfrac{\alpha_0 S_w - \alpha_0 \cdot 2^{\ell_w}\,\mathrm{counts}_w(b)}{S_w(\alpha_0 2^{\ell_w}+S_w)}$;
bounding $\mathrm{counts}_w(b) \ge 0$ and $\mathrm{counts}_w(b) \le S_w$ gives
the numerator magnitude at most $\alpha_0 S_w$, and dividing by the
denominator yields the stated bound. $\blacksquare$

At the smallest tested shot count ($S_w = 4$, `phase1b_combined_stress.csv`)
and $\ell_w \le 6$, this bound is $\le 0.5/(0.5\cdot 64+4) \approx 0.015$,
which `FORMULATION.md` judges "negligible" — this is a quantitative
restatement of that judgment, not a new claim, and it is why the
implementation ranks by raw $\hat p_w$ rather than $\tilde p_w$.

---

## 5. Confidence model

### 5.1 Definition

$$
C_w \;=\; \Pr\big[X_1 > X_2\big], \qquad
X_1 \sim \mathrm{Beta}(\alpha_0+n_1,\ \alpha_0+S_w-n_1), \quad
X_2 \sim \mathrm{Beta}(\alpha_0+n_2,\ \alpha_0+S_w-n_2), \quad X_1 \perp X_2,
$$

exactly `top1_vs_top2_confidence(n1, n2, s_w, alpha0=0.5)`, computed as

$$
C_w = \int_0^1 f_{X_1}(x)\, F_{X_2}(x)\, dx
$$

(`integrand = beta.pdf(x, a1, b1) * beta.cdf(x, a2, b2)`, evaluated by
adaptive quadrature — no closed form is used in the implementation; see §5.4
on why none is available in general for these parameters).

### 5.2 Basic properties (T)

**Proposition 1 (Symmetry).** $C_w(n_1,n_2,S_w) + C_w(n_2,n_1,S_w) = 1$.

*Proof.* $C_w(n_2,n_1,S_w) = \Pr[X_2 > X_1]$ under the same joint law (only
the labels swap). Since $X_1, X_2$ are independent continuous random
variables, $\Pr[X_1 = X_2] = 0$, so $\Pr[X_1>X_2]+\Pr[X_2>X_1] = 1$.
$\blacksquare$

(Verified computationally by `confidence.py:self_test`'s
`c_fwd + c_rev == 1` check; the proof above is why that test is expected to
hold exactly, not merely approximately, up to quadrature error.)

**Proposition 2 (Boundary value).** $C_w(n,n,S_w) = 1/2$ for any $n, S_w$.

*Proof.* Immediate from Proposition 1 with $n_1=n_2=n$: $2C_w(n,n,S_w)=1$. $\blacksquare$

**Proposition 3 (Monotonicity in $n_1$).** For fixed $n_2, S_w$, $C_w$ is
non-decreasing in $n_1$ on $\{0,\dots,S_w\}$.

*Proof sketch.* $\mathrm{Beta}(a,b)$ is stochastically non-decreasing in $a$
for fixed $b$ and stochastically non-increasing in $b$ for fixed $a$ (a
standard total-positivity property of the Beta family — the family has
monotone likelihood ratio in $a$). Increasing $n_1$ by 1 increases
$a_1=\alpha_0+n_1$ by 1 *and* decreases $b_1=\alpha_0+S_w-n_1$ by 1
simultaneously, so both effects push $X_1$ upward in the stochastic-order
sense: $X_1^{(n_1+1)} \succeq_{\mathrm{st}} X_1^{(n_1)}$. Stochastic
dominance is preserved under $\Pr[\cdot > X_2]$ for $X_2$ held fixed, so
$C_w(n_1+1,n_2,S_w) \ge C_w(n_1,n_2,S_w)$. $\blacksquare$

By Proposition 1, $C_w$ is correspondingly non-increasing in $n_2$ for fixed
$n_1, S_w$.

### 5.3 Asymptotic consistency (H, with a proof sketch under stated regularity)

**Proposition 4 (Consistency).** Fix $p_1 > p_2 > 0$ and let
$n_1, n_2, S_w \to \infty$ with $n_1/S_w \to p_1$, $n_2/S_w \to p_2$ (e.g.
under (A1) with true window probabilities $p_1, p_2$ for the top two
patterns and $S_w \to \infty$). Then $C_w \to 1$ in probability.

*Proof sketch.* By the Bernstein–von Mises property of the Beta posterior,
$X_i \big| n_i, S_w \;\approx\; \mathcal{N}\!\big(n_i/S_w,\ \hat p_i(1-\hat p_i)/S_w\big)$
for large $S_w$, uniformly on compacta away from $\{0,1\}$ (standard
Beta-normal approximation; the $\alpha_0$ smoothing term is $O(1/S_w)$ and
vanishes in the limit). Then $X_1 - X_2$ is asymptotically
$\mathcal N(p_1-p_2,\ O(1/S_w))$, whose variance $\to 0$ while the mean is
fixed at $p_1-p_2>0$; by Slutsky's theorem,
$\Pr[X_1-X_2>0] \to \mathbf 1\{p_1-p_2>0\} = 1$. $\blacksquare$

This is the formal counterpart of `FORMULATION.md` §3's stated motivation
("raw probability... [is] sample-size-blind... $C_w$ is the minimal
quantity that distinguishes these two cases") and of the self-test's
empirical check (`c_large > 0.999` at $n_1{:}n_2{:}S_w = 600{:}400{:}1000$).
It is labeled (H) rather than (T) because the Beta-normal approximation
used is itself asymptotic, not exact for finite $S_w$; no finite-sample
rate (e.g. a Berry–Esseen-type bound) is derived or needed by anything
downstream.

### 5.4 On the absence of a closed form (H — a documented implementation choice, not a defect)

For general real $a_1,b_1,a_2,b_2>0$, $\Pr[X_1>X_2]$ for independent Beta
variables has no elementary closed form; it is expressible via a Gauss
hypergeometric function ${}_3F_2$ evaluated at $1$ (a standard but
non-elementary special-function identity). A finite closed-form *sum*
exists in the special case where $b_1$ (or $b_2$) is a positive **integer**
— the identity used by, e.g., Evan Miller's Beta-comparison formula for
two-proportion A/B tests — but here $b_i = \alpha_0 + S_w - n_i =
0.5 + (S_w-n_i)$ is always a half-integer (since $\alpha_0=0.5$ and
$S_w, n_i \in \mathbb Z$), so that finite-sum identity does not apply
without a separate derivation for half-integer parameters. The
implementation instead evaluates the exact integral definition directly by
adaptive quadrature (`scipy.integrate.quad`), which is exact up to floating
point and solver tolerance and requires no case-specific formula — consistent
with `FORMULATION.md` §3's description ("a single, unconditional code path,
no normal-approximation shortcut, no Monte Carlo, no branching by shot
count or regime"). No claim is made here that a closed form is
unobtainable in principle for half-integer parameters (a hypergeometric
reduction likely exists); only that the implementation does not use one,
and correctness does not depend on one existing.

### 5.5 What $C_w$ is *not* (limitation, stated early because it recurs throughout §7–9)

$C_w$ is a posterior probability under (A1)+(A2) that pattern 1's *true*
marginal exceeds pattern 2's, **conditional on (A2) holding**. It is not:
(i) a frequentist confidence level with a guaranteed miscoverage rate — no
frequentist coverage guarantee is proven or claimed anywhere in
`FORMULATION.md`; (ii) calibrated in the sense of matching observed
correctness rates — it is measurably *not*, per §6's Brier-score result,
which is the direct empirical consequence of (A2) being false. Module D's
threshold $1-\varepsilon$ is therefore a *decision threshold on an
uncalibrated statistic*, not a statement that windows crossing it are
correct with probability $\ge 1-\varepsilon$. This distinction is developed
formally in §8.3.

---

## 6. Adaptive candidate coverage (Module B)

### 6.1 Definition

`generate_window_candidates_by_coverage` ranks patterns by $\hat p_w$
descending and returns the shortest prefix $b_{(1)},\dots,b_{(m_w)}$ such
that

$$
\sum_{i=1}^{m_w} \hat p_w(b_{(i)}) \;\ge\; 1-\delta_{\mathrm{cov}}, \qquad m_w \ge 1.
$$

### 6.2 Correctness of the selection rule (T)

**Proposition 5 (Minimality).** $m_w$ as computed by
`generate_window_candidates_by_coverage` equals
$\min\{ m \ge 1 : \sum_{i=1}^m \hat p_w(b_{(i)}) \ge 1-\delta_{\mathrm{cov}}\}$,
i.e. the coverage rule returns the *shortest* prefix meeting the target.

*Proof.* The implementation iterates the descending-ranked list, appending
each candidate and checking `cumulative >= target and candidates` as the
break condition *before* appending the next one — i.e., it stops as soon as
the cumulative sum after the previous append meets the target, never
appending an unneeded element. Since $\hat p_w \ge 0$ and the list is sorted
descending, the partial sums $\sum_{i=1}^m \hat p_w(b_{(i)})$ are
non-decreasing in $m$, so the first $m$ at which the target is met is the
unique minimal such $m$ (existence is guaranteed since
$\sum_i \hat p_w(b_{(i)}) = 1 \ge 1-\delta_{\mathrm{cov}}$ for any
$\delta_{\mathrm{cov}}>0$). $\blacksquare$

This is a statement about the *observed* distribution only — a purely
combinatorial fact about prefix sums of a sorted list, true regardless of
how well $\hat p_w$ approximates $p_w$.

### 6.3 What Module B does *not* prove (limitation)

Minimality of the *empirical* coverage set says nothing about whether the
*true* top-1 pattern is among the $m_w$ selected candidates. A formal
version of that stronger claim would require a concentration inequality
(e.g., a multiplicative Chernoff or Dvoretzky–Kiefer–Wolfowitz-type bound)
relating $\hat p_w$ to $p_w$ as a function of $S_w$, translated into a
statement like "$\Pr[\text{true top-1} \notin \text{selected set}] \le
g(\delta_{\mathrm{cov}}, S_w)$ for some $g$." **No such bound is derived or
implemented.** $\delta_{\mathrm{cov}}$ controls coverage of the *sample*,
not a calibrated error rate on the *population*; the two coincide only in
the $S_w \to \infty$ limit (informally, by the same consistency argument as
Proposition 4). The observed effect of module B (+25.0 pp success,
`FORMULATION.md` §5) is accordingly reported in §9 and
`CORRECTNESS_ARGUMENTS.md` as an **empirical (E)** result, not a
consequence of Proposition 5.

### 6.4 Interaction with stitching

Module B changes the *number* of candidates per window from a constant $k$
to a data-dependent $m_w \le K_w \le \min(2^{\ell_w}, S_w)$, propagated
unchanged into `stitch_candidates`. Since `candidates_compatible` (the
carry-aware overlap check) is applied identically regardless of which
candidate-generation function produced its inputs, **Module B cannot change
the correctness of the overlap check** — only the population of candidates
it is applied to. This separation (candidate population vs. compatibility
test) is what makes B and C independently ablatable per `FORMULATION.md`
§5's design.

---

## 7. Adaptive beam width (Module C)

### 7.1 Definition

$$
\mathrm{expanded} = \#\Big\{ w : \sum_{i=1}^{k} \hat p_w(b_{(i)}) < 1-\delta_{\mathrm{cov}} \Big\}
\qquad\text{(`_windows_needing_expansion`, using the fixed baseline }k\text{, independent of whether B is enabled)}
$$

$$
P' = \min\big(\texttt{beam\_max},\ P\cdot(1+\mathrm{expanded})\big).
$$

Note $\mathrm{expanded}$ is computed by re-ranking each window's counts
against the *fixed* `candidate_count` $k$, not against module B's own
data-dependent $m_w$ — this is what keeps B and C measurable as separately
switchable arms in the same ablation (`FORMULATION.md` §5: "computed
independently of whether module B's coverage-based candidate generation is
actually enabled").

### 7.2 Conservativity (T)

**Proposition 6 (Beam width is never more restrictive than baseline).**
Whenever $\texttt{beam\_max} \ge P$ (true for every configuration in
`FORMULATION.md`'s ablation: $\texttt{beam\_max}=4096 \gg P \le 32$),
$P' \ge P$.

*Proof.* $\mathrm{expanded} \ge 0$ so $P\cdot(1+\mathrm{expanded}) \ge P$.
$\min(\texttt{beam\_max}, P\cdot(1+\mathrm{expanded})) \ge \min(\texttt{beam\_max}, P) = P$
under the stated hypothesis $\texttt{beam\_max}\ge P$. $\blacksquare$

**Corollary.** Under (T)-Proposition 6's hypothesis, module C's beam-search
truncation can only *retain* more candidate paths than baseline at each
level of `stitch_candidates`, never fewer — the true globally-consistent
path is never pruned by C where it would have survived under baseline (see
§8.2 for the general pruning-optimality statement this depends on).

### 7.3 What Module C does *not* prove (limitation)

$\mathrm{expanded}$ is a simple count of ambiguous windows, not a bound on
the number of globally-consistent paths that actually exist — that quantity
depends on the joint structure of all $W$ windows' candidate sets and the
carry check, which is combinatorial and not captured by a per-window scalar
(see §8.2 for the growth bound, and `COMPLEXITY_ANALYSIS.md` §4 for its
consequences). $P'$ is therefore a **heuristic (H)** budget-scaling rule
justified by the qualitative intuition in `PHASE6_SUMMARY.md` ("if many
windows needed extra candidates... the search space is more ambiguous than
usual"), not a derived sufficient beam width for exact recovery. No claim
is made, nor was one specified in `FORMULATION.md`, that $P'$ guarantees the
true path survives beam pruning.

---

## 8. Adaptive shot stopping (Module D)

### 8.1 Definition

For each window $w$, `_resample_shot_stopping` repeats: compute
$C_w$ from current counts; if $C_w \ge 1-\varepsilon$, stop; else draw
$\min(S_w, S_{\max}-S_w)$ additional shots from the *same* operator/eigenstate
(`ModularMultiplicationOperator`, `run_windowed_qpe_block`) and merge counts;
repeat until either the threshold is met or no further shots can be added
($S_w \ge S_{\max}$).

### 8.2 Termination (T)

**Proposition 7 (Bounded rounds).** The resampling loop for window $w$
terminates after at most
$R_w \le \lceil \log_2(S_{\max}/S_w^{(0)}) \rceil$ rounds, where
$S_w^{(0)}$ is the initial shot count.

*Proof.* Each round adds $\min(S_w, S_{\max}-S_w)$ shots, where $S_w$ is the
*current* total. While $S_w \le S_{\max}/2$, the added amount is $S_w$
itself (since $S_{\max}-S_w \ge S_w$), so $S_w$ at least doubles each such
round; hence there are at most $\lceil\log_2(S_{\max}/S_w^{(0)})\rceil$
doubling rounds before $S_w > S_{\max}/2$. Once $S_w > S_{\max}/2$, the
added amount is $S_{\max}-S_w < S_w$, i.e. strictly less than a doubling,
and $S_w$ is strictly increasing and bounded above by $S_{\max}$, so the
loop must reach $S_w = S_{\max}$ (at which point `extra_shots <= 0` forces
exit) in at most one further round beyond the doubling phase. The total
round count is therefore $O(\log_2(S_{\max}/S_w^{(0)}))$. $\blacksquare$

This proves the loop always halts — it is a genuine termination guarantee,
not merely an engineering cap — and matches the configured defaults exactly:
`config.py`'s comment ("4 doublings from shots=4") is the $S_w^{(0)}=4$,
$S_{\max}=64$ instance of this proposition, $\lceil\log_2(64/4)\rceil = 4$.

### 8.3 What Module D does *not* prove (limitation, sharpened from §5.5)

Because $C_w$ is not calibrated (§5.5, §9), the stopping condition
$C_w \ge 1-\varepsilon$ does **not** imply
$\Pr[\text{top-1 pattern is truly correct}] \ge 1-\varepsilon$. The
calibration study's own numbers make this precise: in the top predicted-confidence
decile $C_w \in (0.9,1.0]$ (mean predicted $0.96$), the *observed*
correctness rate was $0.79$ (`FORMULATION.md` §6) — i.e. windows that
satisfy D's stopping rule at $\varepsilon=0.05$ (threshold $0.95$, inside
this decile) are, empirically, wrong roughly $21\%$ of the time, not $5\%$.
**$\varepsilon$ is therefore a tuning parameter for a heuristic stopping
rule, not a statistically guaranteed error tolerance**, despite its name.
This does not contradict module D's large measured effect (+24.3–25.2 pp
success, `FORMULATION.md` §5 / `PHASE6_SUMMARY.md` §5): a miscalibrated but
monotonically-informative statistic (Proposition 3, 4) can still be a highly
effective *relative* signal for "resample this window" even though its
absolute value should not be read as a probability. The distinction between
"useful ranking signal" and "calibrated probability" is exactly the
distinction the isotonic-recalibration fit (`Experiments/phase6_calibration.py`)
was built to close, and which module D deliberately does not use
(`FORMULATION.md` §6, `_resample_shot_stopping`'s own docstring).

---

## 9. Theoretical properties of the full framework

**Proposition 8 (Conservativity / baseline reproduction) (T).** With
`use_coverage=use_beam_rule=use_shot_stopping=False`,
`reconstruct_adaptive` executes the identical sequence of operations as the
baseline path (`generate_window_candidates` with the fixed `candidate_count`,
`stitch_candidates` with the fixed `window.max_paths`, no resampling) —
verifiable directly from `reconstruct_adaptive`'s branching, all three
`if use_*:` blocks are no-ops when their flags are `False`, leaving exactly
baseline's calls. This is why `FORMULATION.md` §5 can validly compare
Baseline against B/C/D/Full as *arms of one function*, not as separately
implemented pipelines, and why the reproduced $52.1\%$ (design-note ablation)
/ $51.1\%$ (`PHASE6_SUMMARY.md` execution) baseline success rates against
`REPORT.md`'s original $52.1\%$ serve as a valid sanity check on this
proposition (an empirical corroboration of a structural code fact, not a
statistical result in its own right).

**Proposition 9 (Monotone non-restrictiveness of B, C, D individually) (T,
conditional).** Under the stated per-module hypotheses (§6.2 for B — trivially,
since $m_w \ge 1$ can only add candidates relative to no floor being enforced
by B beyond what coverage requires; §7.2 Proposition 6 for C; §8.2
Proposition 7 for D's bounded termination), each module can only *relax* a
resource constraint relative to baseline (more candidates, wider beam, more
shots), never tighten one. No module removes information or candidates that
baseline would have kept. This is a structural property of the code, and it
is the reason the observed effect sizes in `FORMULATION.md` §5 are all
positive (no arm underperforms baseline) — though the *sizes* of those
effects (§9's remaining content) are empirical, not derivable from
monotonicity alone, since a monotone relaxation could in principle have zero
or negligible effect if the relaxed resource were never the binding
constraint.

**Convergence discussion (H).** As $S_{\max}\to\infty$, $\delta_{\mathrm{cov}}\to 0$,
and `beam_max` $\to\infty$ jointly (i.e., all three adaptive budgets
unbounded), the framework's classical reconstruction stage converges to an
*exhaustive* search: module B's coverage cutoff admits every observed
pattern, module D drives every window's $C_w \to 1$ by Proposition 4
(assuming a true probability gap $p_1>p_2$ exists, e.g. from decisive quantum
amplitude structure — see (A1)), and module C's beam width no longer prunes
any consistent path. In this joint limit, and only in this limit, the
classical stage becomes complete relative to the sampled/samplable
information: recovery fails only if no consistent stitched path exists at
all (a genuine information-theoretic deficiency, e.g., all windows agreeing
on a wrong bit pattern by chance) rather than from budget starvation. This
is presented as a *qualitative limit statement*, not a rate-of-convergence
result — no finite-$S_{\max}$, finite-`beam_max` guarantee is derived, and
the ablation's own finite-cap results ($S_{\max}=64$, `beam_max`$=4096$)
are reported as the empirical operating point actually validated, not as an
approximation to this limit with a known error bound.

---

## 10. Limitations (consolidated)

1. **(A2) independence assumption is measurably false** (§3, §5.5) —
   Brier $0.2239$ raw, the direct empirical cost of treating two correlated
   Dirichlet marginals as independent binomials. This is the single most
   consequential gap between the model and reality in the whole framework.
2. **No PAC-style coverage guarantee for module B** (§6.3): $\delta_{\mathrm{cov}}$
   bounds sample coverage, not population/true-pattern coverage.
3. **$\varepsilon$ is not a calibrated error rate** (§8.3): module D's
   stopping rule is a useful, monotone, but uncalibrated ranking signal, not
   a statistical guarantee, and the shipped configuration deliberately does
   not apply the recalibration that would partially close this gap.
4. **Module C's beam-scaling rule is heuristic** (§7.3): $\mathrm{expanded}$
   is a simple per-window ambiguity count, not a derived sufficient beam
   width; no exactness guarantee is claimed or provable from the rule alone.
5. **Top-2-only simplification (A3)** is unvalidated for the general case
   and known to be structurally blind to genuine 3-way ties.
6. **No cross-window resource coupling (A4)**: each module allocates budget
   per-window independently; no formal analysis of joint, budget-constrained
   allocation across windows is attempted (`FORMULATION.md` §4.3 notes this
   was out of scope because no experiment in the harness exercises a shared
   cross-window budget).
7. **All experimental validation of this framework (`FORMULATION.md` §5, §6)
   was conducted at small multiplicative order** ($r=6$, $N=21$), within the
   compilation-cliff-imposed work-register ceiling documented in `REPORT.md`
   §0/§5. The continued-fraction stage's own forgiveness at small order
   (`REPORT.md` §2) is a known confound for *any* claim about classical-stage
   robustness, adaptive or not — this framework's measured gains are not yet
   validated at an order large enough to remove that confound.
8. **Overlap geometry is not adaptive** (A5, `FORMULATION.md` §8): none of
   the propositions above say anything about instability near a window's
   0.5 ambiguity boundary (`REPORT.md` §4.3), a distinct, unaddressed
   failure mode this framework does not target.
