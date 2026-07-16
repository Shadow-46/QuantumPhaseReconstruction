# Complexity Analysis — Confidence-Guided Adaptive Reconstruction

*Companion to `THEORETICAL_ANALYSIS.md` and `CORRECTNESS_ARGUMENTS.md`.
Analyzes the computational complexity of the published (baseline)
reconstruction pipeline and of modules B (candidate coverage), C (beam
width), and D (shot stopping), as implemented in
`Reconstruction/candidate_generation.py`, `Reconstruction/carry.py`,
`Reconstruction/stitching.py`, `Reconstruction/continued_fraction.py`, and
wired in `Experiments/harness.py::reconstruct_adaptive`. All derivations are
against the actual code; no complexity claim below describes a hypothetical
alternative implementation.*

---

## 0. Common quantities and cost primitives

| Symbol | Meaning |
|---|---|
| $W$ | number of windows (`len(make_overlapping_windows(...))`) |
| $\ell_w$ | width of window $w$ in bits |
| $K_w$ | number of *distinct* observed bit patterns in window $w$'s counts, $K_w \le \min(2^{\ell_w}, S_w)$ |
| $k$ | fixed `candidate_count` (baseline) |
| $m_w$ | data-dependent candidate count kept for window $w$ under module B, $1 \le m_w \le K_w$ |
| $P$ | fixed `max_paths` (baseline beam width) |
| $P'$ | adaptive beam width under module C |
| $o$ (`shared`) | overlap width in bits between adjacent windows, $o \le \min(\ell_w,\ell_{w+1})$; $o \le 3$ in every experiment `REPORT.md`/`FORMULATION.md` reports |
| $N$ | RSA-style modulus |
| $R_w$ | number of resampling rounds module D performs on window $w$ |
| $T_{\mathrm{sample}}$ | wall-clock cost of one Aer transpile+run for one window's resample round (external, non-classical cost) |

Two cost primitives recur and are stated once:

- **Sorting a window's counts.** `normalize_counts` builds a `dict` of size
  $K_w$ (one pass, $O(K_w)$); ranking it (`sorted(..., key=..., reverse=True)`)
  costs $O(K_w \log K_w)$ comparisons. This same sort is re-executed,
  independently, by `generate_window_candidates`, `generate_window_candidates_by_coverage`,
  and `_windows_needing_expansion` — each call site re-sorts from scratch;
  none reuses another's sorted result (an identifiable, though minor,
  redundancy — see §4.4).
- **A single carry-aware compatibility check.** `candidates_compatible`
  computes `overlap_width` ($O(1)$), slices two bit strings of length $o$
  and converts them to integers ($O(o)$ under a unit-cost integer model, or
  $O(o)$ word operations), and performs one modular subtraction — total
  $O(o)$, and since $o \le 3$ throughout the evaluated corpus, effectively
  $O(1)$ in practice but stated as $O(o)$ for generality.

---

## 1. The published (baseline) algorithm

Baseline = `reconstruct_adaptive` with all three flags `False`, structurally
identical to `reconstruct_from_counts` (Proposition 8,
`THEORETICAL_ANALYSIS.md` §9).

### 1.1 Candidate generation — `generate_window_candidates`

- **Time:** $O(K_w \log K_w)$ per window (sort), plus $O(k)$ to build the
  returned list $\Rightarrow O(K_w\log K_w + k)$ per window,
  $O\!\big(\sum_w K_w \log K_w\big)$ total across $W$ windows.
- **Space:** $O(K_w)$ transient (the normalized dict) $+ O(k)$ retained
  (the returned candidate list) per window.
- **Worst case:** $K_w = 2^{\ell_w}$ (every possible pattern observed at
  least once) $\Rightarrow O(2^{\ell_w}\ell_w)$ per window — exponential in
  window width, but window widths in the evaluated corpus are $\le 8$
  (`REPORT.md` §2, `window_size` swept $2$–$8$), so $2^{\ell_w}$ stays
  small ($\le 256$) in practice; this is a property of $\ell_w$ being kept
  deliberately small by the windowed formulation itself, not of the
  reconstruction algorithm's efficiency.
- **Expected case:** since $K_w \le S_w$ always, and the evaluated
  combined-stress grid uses `shots`$\in\{4,8,16\}$ (`FORMULATION.md` §5),
  $K_w \le 16$ in the hardest tested configurations, making this stage cheap
  in absolute terms regardless of $\ell_w$.

### 1.2 Carry-aware compatibility — `candidates_compatible`

$O(o)$ per call (§0), invoked once per (path, candidate) pair inside
stitching — folded into §1.3's cost, not separately bottleneck-relevant.

### 1.3 Beam stitching — `stitch_candidates`

For each of the $W-1$ transitions between adjacent window-candidate groups,
with $|paths|$ the current surviving path count (initially $|windows_1| = k$):

- **Per-level expansion:** for each of $\le P$ surviving paths and each of
  $k$ candidates in the next group, one $O(o)$ compatibility check $\Rightarrow$
  $O(P\cdot k\cdot o)$ candidate pairs generated (before pruning).
- **Per-level sort + truncate:** sorting up to $P\cdot k$ generated
  next-paths costs $O(Pk\log(Pk))$, then truncating to $P$.
- **Total across $W-1$ levels:** $O\!\big(W\cdot P\cdot k\cdot(o + \log(Pk))\big)$ time.
- **Space:** each surviving path stores references to up to $W$ candidates
  $\Rightarrow O(P\cdot W)$ for the path list at any level, plus $O(P)$ for
  the final `StitchedPhase` list. The *transient* per-level candidate list
  before pruning is $O(P\cdot k)$ path objects, each $O(W)$ in size in the
  worst case (a path near the end of the window sequence) — worst-case
  transient space $O(P\cdot k\cdot W)$.
- **Worst case (beam never binds, i.e. every pair compatible and $P\ge k^{W-1}$):**
  path count grows to $k^{W}$ before any pruning — the beam width $P$ exists
  specifically to cap this combinatorial blowup; with $P$ fixed, per-level
  cost is bounded as above regardless of $W$'s growth, at the cost of
  potentially pruning the correct path (a completeness/optimality trade,
  see `CORRECTNESS_ARGUMENTS.md` §3).

### 1.4 Order/factor recovery — `recover_order_and_factors` / `order_from_phase`

- **Time:** `order_from_phase` loops `multiplier in range(1, N+1)`, each
  iteration performing one `pow(a, candidate, N)` modular exponentiation —
  $O(\log N)$ modular multiplications via Python's built-in square-and-multiply,
  each multiplication $O(M(\log N))$ under a chosen multiplication cost
  model ($M(b) = O(b^2)$ schoolbook or $O(b\log b)$ with fast multiplication;
  Python's CPython uses schoolbook-ish multiplication for the moduli sizes
  here). Total per call: $O(N\log N \cdot M(\log N))$ worst case (loop never
  exits early), i.e. **linear in $N$ times poly-log factors** — this is the
  one place baseline complexity is not polynomial in the *circuit*
  parameters ($W,k,P$) but in the *arithmetic* parameter $N$ directly,
  because the search is a brute-force multiplier scan, not an
  order-finding shortcut.
- **Expected case:** `REPORT.md` §2 documents that for the small orders
  tested ($r\in\{2,6\}$), the search "is extremely forgiving" and exits on
  the first or an early multiplier for the great majority of measured
  phases — the *effective* number of iterations before success is small and
  order-dependent, not $\Theta(N)$, though no formal expected-iteration-count
  bound is derived anywhere in the evidence base; this is reported as an
  **(E)** observation, not a theorem.
- **Outer loop cost:** `reconstruct_adaptive`/`reconstruct_from_counts`
  calls `recover_order_and_factors` once per stitched candidate, in ranked
  order, stopping at the first success — up to $P$ calls in the worst case
  (no candidate recovers valid factors) $\Rightarrow O(P\cdot N\log N\cdot M(\log N))$
  worst case for this stage.

### 1.5 Baseline total

$$
T_{\text{baseline}} = O\Big(\underbrace{\textstyle\sum_w K_w\log K_w}_{\text{candidate gen}} \;+\; \underbrace{W\cdot P\cdot k\cdot(o+\log(Pk))}_{\text{stitching}} \;+\; \underbrace{P\cdot N\log N\cdot M(\log N)}_{\text{order recovery}}\Big)
$$

$$
S_{\text{baseline}} = O\big(\textstyle\sum_w K_w + P\cdot k\cdot W\big)
$$

For the evaluated regime ($W\le 32$, $k\le 2$, $P\le 32$, $o\le3$,
$N\le 21$ in the combined-stress grid; up to $N=589$ elsewhere in the
corpus per `REPORT.md` §0), the stitching and order-recovery terms dominate
in wall-clock terms not because of their asymptotic class but because of
the real quantum-circuit sampling cost that produces $\mathrm{counts}_w$ in
the first place (external to this classical-complexity analysis; see
`REPORT.md` §5.1's compilation-cliff finding, which is a property of
circuit construction, not of the classical reconstruction code analyzed
here).

---

## 2. Module B — candidate coverage

`generate_window_candidates_by_coverage` replaces the fixed truncation at
$k$ with a data-dependent truncation at $m_w$ (Proposition 5,
`THEORETICAL_ANALYSIS.md` §6.2).

- **Time:** identical sort $O(K_w\log K_w)$, plus a single linear scan with
  early exit $O(m_w) \le O(K_w)$ — **same asymptotic class as baseline's
  candidate generation**, with the truncation constant now $m_w$ instead of
  $k$.
- **Space:** $O(m_w)$ retained per window instead of $O(k)$.
- **Worst case:** a perfectly flat empirical distribution forces
  $m_w \to K_w$ (every pattern needed to reach $1-\delta_{\mathrm{cov}}$
  coverage) — identical to baseline's worst case with `candidate_count` set
  to $2^{\ell_w}$, i.e. module B's *worst case* candidate-generation cost
  equals baseline's worst case, not better.
- **Best case:** a near-deterministic window collapses to $m_w=1$
  (`candidate_generation.py:self_test`'s `confident_counts` example),
  strictly *cheaper* per-window than baseline's fixed $k\ge2$ in the
  evaluated defaults.
- **Expected case (E):** ambiguous windows in the evaluated corpus needed
  "multiple candidates" (self-test, `ambiguous_counts`); no aggregate
  distribution of $m_w$ across the 144-trial ablation grid is published in
  `FORMULATION.md`/`PHASE6_SUMMARY.md`, so no numeric expected-case figure
  is asserted here beyond the qualitative "adapts up or down" behavior the
  code guarantees by construction.
- **Comparison with published implementation:** same time/space complexity
  *class* per window ($O(K_w\log K_w)$ time, $O(m_w)$ vs. $O(k)$ space); the
  actual number kept becomes data-dependent rather than constant, which is
  the entire mechanism of B's effect and is **not** a complexity
  improvement or regression in the classical Big-O sense — B's measured
  benefit (+25.0 pp success, `FORMULATION.md` §5) is a *statistical*, not
  *computational-complexity*, improvement, and its downstream cost is paid
  in §3 (larger $m_w$ propagating into stitching).

---

## 3. Module C — beam width

### 3.1 `_windows_needing_expansion`

- **Time:** for each window, sorts `normalize_counts(counts).values()`
  ($O(K_w\log K_w)$, re-sorting counts already sorted elsewhere in the same
  call — see §4.4) and sums the top $k$ ($O(k)$) $\Rightarrow O(K_w\log K_w)$
  per window, $O(\sum_w K_w\log K_w)$ total — **the same order as one full
  pass of baseline's own candidate generation**, effectively doubling the
  counts-sorting work when module C is enabled alongside B, since B, C, and
  (if D is also enabled) D's resampling loop's own $C_w$ computation are all
  independent re-derivations from the same underlying counts.
- **Space:** $O(1)$ beyond the transient sorted list.

### 3.2 Beam-width formula and its downstream cost

$P' = \min(\texttt{beam\_max},\ P\cdot(1+\mathrm{expanded}))$ is $O(1)$ to
compute (§3.1 dominates). Its effect is entirely on §1.3's stitching cost,
substituting $P'$ for $P$:

$$
T_{\text{stitch, C}} = O\big(W\cdot P'\cdot k'\cdot(o+\log(P'k'))\big), \qquad k' = k \text{ or } m_w \text{ (depending on whether B is also enabled)}.
$$

- **Worst case:** $\mathrm{expanded}=W$ (every window ambiguous) gives
  $P' = \min(\texttt{beam\_max}, P\cdot(1+W))$; at the evaluated defaults
  ($P=32$, `beam_max`$=4096$, $W$ up to $32$ per `FORMULATION.md` §5's
  grid), $P\cdot(1+W) = 32\cdot33 = 1056 < 4096$ — the cap does not bind at
  this specific grid's worst case, but *would* bind for larger $W$
  (`beam_max` becomes the effective ceiling once $W > \texttt{beam\_max}/P - 1 \approx 127$
  at these defaults). This is exactly the mechanism `REPORT.md` Appendix A
  documents as having produced "multi-gigabyte path-list growth" before the
  ladder's own beam-width proxy was capped at 64 — an empirical (E)
  instance of this worst case being reached in practice (in the Phase 5
  attribution ladder, not in module C itself, but via the identical
  mechanism: beam width $\times$ candidates-per-window $\times$ windows,
  uncapped, multiplies combinatorially).
- **Multiplicative overhead vs. baseline:** since $P' \le \texttt{beam\_max}/P$
  times larger than $P$ in the worst case, stitching cost can grow by up to
  a factor of $\texttt{beam\_max}/P$ (here, up to $128\times$ at the
  evaluated defaults) relative to baseline — this is the *quantified price*
  of module C's relaxation, paid only on windows/trials where
  $\mathrm{expanded}>0$; trials with no ambiguous windows incur zero
  overhead ($P'=P$).
- **Space:** scales identically to time's multiplicative factor:
  $O(P'\cdot k'\cdot W)$ transient, $O(P')$ retained.

### 3.3 Comparison with published implementation

Baseline's $P$ is a compile-time constant with $O(1)$ computation and a
fixed worst-case stitching cost. Module C keeps the *same* stitching
algorithm and asymptotic *shape* (linear in $W$, linear in beam width,
$\log$ in beam-width$\times$candidates for the sort) but makes the beam
width itself a bounded random variable $P' \in [P, \texttt{beam\_max}]$
depending on the data — the complexity class is unchanged, but the
realized constant is no longer fixed at design time, and its ceiling
(`beam_max`) is an engineering cap carried over unchanged from the Phase 5
blowup measurement (`FORMULATION.md` §4.2), not re-derived here.

---

## 4. Module D — shot stopping

### 4.1 Per-round cost

Each round of `_resample_shot_stopping` for window $w$: one
$O(1)$-ish confidence evaluation `top1_vs_top2_confidence` (adaptive
quadrature with a bounded iteration count independent of $n_1,n_2,S_w$
magnitude — `scipy.integrate.quad`'s adaptive subdivision is governed by
the integrand's smoothness and the requested tolerance, not by the Beta
parameters' size, so this is treated as $O(1)$ per call for complexity
purposes, though not literally constant-time at the level of floating-point
operation counts), plus one real circuit resample
(`run_windowed_qpe_block`, cost $T_{\mathrm{sample}}$ — external to this
classical analysis, dominated by Aer transpile+execution, measured at
$\approx 3.3\,\mathrm{s}$/round in this environment per `FORMULATION.md`
§7), plus $O(\min(K_w, S_{\max}))$ to merge the new counts dict.

### 4.2 Total rounds and cost

By Proposition 7 (`THEORETICAL_ANALYSIS.md` §8.2),
$R_w \le \lceil\log_2(S_{\max}/S_w^{(0)})\rceil$. Total added classical
work across all windows:

$$
T_{\text{D, classical}} = O\Big(\sum_w R_w\Big) = O\big(W\log_2(S_{\max}/S_{\min}^{(0)})\big),
$$

and total added *wall-clock* (dominated by external sampling, not classical
computation):

$$
T_{\text{D, sampled}} = O\Big(\sum_w R_w\Big)\cdot T_{\mathrm{sample}} = O\big(W\log_2(S_{\max}/S_{\min}^{(0)})\big)\cdot T_{\mathrm{sample}}.
$$

- **Space:** $O(S_{\max})$ worst case per window for the merged counts dict
  (bounded, since counts can have at most $\min(2^{\ell_w}, S_{\max})$
  distinct keys) — no growth beyond what a single window's final,
  fully-resampled counts dict would need.
- **Worst case:** every window starts at the smallest tested shot count
  ($S_w^{(0)}=4$) and never crosses the confidence threshold before
  $S_{\max}$ — $R_w = \lceil\log_2(64/4)\rceil = 4$ rounds per window, the
  bound realized exactly at the configured defaults
  (`config.py`'s own comment: "4 doublings from shots=4"), giving
  $T_{\text{D,sampled}} = O(4W)\cdot T_{\mathrm{sample}}$, i.e. up to
  $5\times$ the circuit-sampling cost of a single baseline pass per window
  (original sample + 4 resample rounds) in the worst case.
- **Why $S_{\max}=64$ specifically (E, from `config.py`):** chosen after
  measuring $\approx 3.3\,\mathrm s$/round in this environment; a single
  D-arm trial at the original, unbounded ($S_{\max}=65536$, $14$ doubling
  rounds) setting took $\approx120\,\mathrm s$ — an empirically-driven
  engineering cap on $R_w$'s worst case, not a value derived from the
  complexity bound itself (the bound in §4.2 explains *why* capping
  $S_{\max}$ linearly caps $R_w$ logarithmically, but the specific value
  64 was chosen by measurement, not derivation).

### 4.3 Comparison with published implementation

Baseline samples each window exactly once, at a fixed $S_w$, with zero
classical resampling logic — $O(1)$ circuit executions per window. Module D
adds $O(\log_2(S_{\max}/S_w^{(0)}))$ *additional* circuit executions per
window in the worst case, each costing one real transpile+run
($T_{\mathrm{sample}}$), which is the single most expensive operation in
the entire adaptive framework in absolute wall-clock terms (§4.1) — more
expensive, per unit, than any of the purely classical operations in §1–3,
none of which touch the quantum simulator at all. This asymmetry is why
`FORMULATION.md` §7 caps $S_{\max}$ far below its original value: the
classical stages (B, C) are cheap by comparison and were not similarly
capped for performance reasons (their caps — $\delta_{\mathrm{cov}}$,
`beam_max` — are set for statistical/blowup reasons, not wall-clock ones).

---

## 5. The complete adaptive framework

### 5.1 Aggregate complexity

$$
T_{\text{full}} = O\Big(\underbrace{\textstyle\sum_w K_w\log K_w}_{\text{B, and C's redundant re-sort, §2, §3.1}} \;+\; \underbrace{W\cdot P'\cdot m_w\cdot(o+\log(P'm_w))}_{\text{stitching, now with adaptive }P', m_w\text{ , §3.2}} \;+\; \underbrace{P'\cdot N\log N\cdot M(\log N)}_{\text{order recovery, now over up to }P'\text{ candidates}} \;+\; \underbrace{W\log_2(S_{\max}/S_{\min}^{(0)})\cdot T_{\mathrm{sample}}}_{\text{D's resampling, §4.2}}\Big)
$$

$$
S_{\text{full}} = O\Big(\textstyle\sum_w m_w \;+\; P'\cdot m_w\cdot W \;+\; W\cdot S_{\max}\Big)
$$

**Every term is bounded by an explicit engineering cap** carried over
unchanged from the evidence base, not newly derived here: $m_w \le K_w \le
\min(2^{\ell_w}, S_w)$, $P' \le \texttt{beam\_max}=4096$, $R_w \le
\lceil\log_2(S_{\max}/S_w^{(0)})\rceil$ with $S_{\max}=64$. The framework is
therefore **worst-case bounded but not worst-case cheap**: relative to
baseline, the realized cost on any given trial is a random variable whose
*distribution* depends on how ambiguous the sampled counts turn out to be,
with a hard ceiling set by the caps rather than by any property of the
input instance.

### 5.2 Which term actually dominates, in practice (E)

By direct wall-clock reasoning from the measurements cited in
`FORMULATION.md` §7 and `REPORT.md` §5: circuit execution
($T_{\mathrm{sample}} \approx 3.3\,\mathrm s$/round, and the underlying
per-window base sampling cost, itself dominated by transpile time per
`REPORT.md` §5.1) is orders of magnitude more expensive per operation than
any purely classical step in §1–3 (sorting a $\le 16$-entry dict, a handful
of $O(o{=}3)$ compatibility checks, a bounded-iteration quadrature call).
**Module D's resampling is therefore the dominant wall-clock cost of the
adaptive framework whenever it triggers**, not modules B or C, even though
B and C are the ones with the larger asymptotic-complexity footnotes
(exponential-in-$\ell_w$ worst case for candidate generation, combinatorial
beam blowup risk). This is an empirical ranking of practical cost, not a
restatement of the asymptotic bounds in §5.1 (which do not, by themselves,
rank the terms without plugging in the measured constants).

### 5.3 Comparison summary table

| Stage | Baseline time | Baseline space | Adaptive time | Adaptive space | Adaptive worst-case multiplier vs. baseline |
|---|---|---|---|---|---|
| Candidate generation | $O(K_w\log K_w)$/window | $O(k)$/window | $O(K_w\log K_w)$/window (B) | $O(m_w)$/window | up to $K_w/k$ (B alone can only *grow* the kept set) |
| Beam stitching | $O(W P k(o+\log(Pk)))$ | $O(PkW)$ | $O(W P' m_w(o+\log(P'm_w)))$ (C, and B via $m_w$) | $O(P'm_wW)$ | up to $(\texttt{beam\_max}/P)\times(K_w/k)$ jointly |
| Order recovery | $O(P N\log N\, M(\log N))$ | $O(1)$ extra | $O(P' N\log N\, M(\log N))$ | $O(1)$ extra | up to $\texttt{beam\_max}/P$ (more stitched candidates to try) |
| Sampling | $O(1)$ circuit run/window | — | $O(1+R_w)$ circuit runs/window (D) | $O(S_{\max})$/window counts | up to $1+\lceil\log_2(S_{\max}/S_w^{(0)})\rceil$ ($5\times$ at evaluated defaults) |

All four rows keep the **same asymptotic functional form** as baseline
(polynomial in $W, P, k, N$, with the same $\log N$ order-recovery term
untouched by any module); the adaptive framework's cost is a **bounded
multiplicative inflation** of baseline's cost, with the bound set by
$\delta_{\mathrm{cov}}$ (indirectly, via $m_w \le K_w$), `beam_max`, and
$S_{\max}$ — exactly the three free parameters `FORMULATION.md` §7 lists as
the framework's complete tunable surface. No module introduces a new
complexity class (e.g. nothing becomes exponential in $W$ that was not
already potentially exponential in $W$ at baseline, absent the $P$/$P'$ cap).
