# Boundary-Local Adaptive Overlap (BLAO): design proposal, revision 2 (FOR REVIEW; nothing implemented)

- **Status.** DESIGN ONLY. No BLAO code, no pilot, no dev or test runs. P5-P10b are frozen and untouched.
- **Revision history.**
  - Revision 1 (2026-09-28).
  - Revision 2 (2026-09-29): the conditional-rescue claim is qualified; BLAO is split into two separate research questions; exact signal definitions; a local-decoder analysis; an n = 16 worked example; falsification criteria; the minimal pilot.

---

## 0. What the P6 follow-up does and does NOT establish

**Withdrawn as evidence for BLAO.** Revision 1 cited the P6 follow-up: "overlap placed at the responsible boundary j* rescues 91-100% of failures". That figure is **not evidence that boundary-local placement works.**
- In the analysed variant (ext_v1_gated_A2 with awqpe_ext), the trigger `lowerhalf_mask` and the decoder rule in `awqpe_ext_decode` require the **same** predicate: the corrected lower chunks j+1..B equal exactly 10...0.
- Overlap therefore lands at j* only in trials where the decoder is also able to use it.
- The 91-100% is a rescue rate **conditional on that shared predicate**, a selected subset. It says nothing about overlap placed at j* for any other reason, under either decoder.
- RESULTS_LOG carries the same qualification.

**What the P6 follow-up does establish.**
- (i) At n = 16, the frozen trigger essentially never places overlap at upper boundaries.
- (ii) The shot-fixable failures there are mostly higher-order errors.

Neither says a boundary-local method would rescue them. That is what BLAO would test.

## 1. Notation

**Chunks, boundaries and the boundary variable.**
- Chunks i = 0..B−1 are MSB-first, with widths m_i, offsets k_i = Σ_{l<i} m_l, and M_i = 2^(m_i).
- Boundary j (0 ≤ j ≤ B−2) separates chunk j from chunk j+1. It has K_j = k_{j+1} bits above it.
- Let x_j = 2^(K_j) φ, p_j = ⌊x_j⌋ mod M_j, and r_j = frac(x_j) ∈ [0,1), the remainder below boundary j.
- Borrow bit: β_j = 1[r_j ≥ 1/2].

**Block law** (paper, Lemma 3.3).
- Chunk i's block sees δ_i = frac(2^(k_i) φ) and returns y with P(y|φ) = K_(M_i)(δ_i − y/M_i).
- The chunk-(j+1) block sees δ_(j+1) = r_j **exactly**.
- The chunk-j block sees M_j δ_j = (p_j mod M_j) + r_j, so it straddles p_j and p_j + 1 when r_j ≈ 1/2.

**Target.** The best n-bit numerator N* = ⌊2^n φ + 1/2⌋. Its chunk j equals p_j unless the lower part rounds up with a carry.

**Overlap actions** (unchanged from P6).
- O-ext_j: chunk j re-measured with width m_j + v, v = 1. It reads round(2^v x_j) mod 2^(m_j+v) in the limit.
- O-bridge_j: a width-m_j block at offset K_j − s, s = 1.
- A batch is S0 shots, with at most A overlap batches per trial.
- The two mechanisms are never mixed within an arm.

**Observables at decision time t.**
- Counts N_i(y) of every chunk block, and N_o(y) of every overlap block already measured.
- The shot ledger.
- The D2 grid log-likelihood ℓ_t(g) on φ_g = g/G, which is a deterministic function of the counts.
- Nothing else. In particular, not φ, not j*, and not the stratum label.

---

## 2. Two separate research questions: do not conflate

| | **Track A: BLAO-D2** | **Track B: boundary-local extended decoder ("awqpe_bl")** |
|---|---|---|
| Question | Does choosing overlap *geometry* by a boundary-local criterion improve likelihood (D2) reconstruction at equal U-query cost? | Can a D1-type decoder consume an O-ext block at boundary j using only local information, without the whole-lower-part 10...0 predicate? |
| What changes | Only the trigger (where and when overlap is placed). The decoder (D2) is unchanged and uses every block. | The decoder rule. The trigger is secondary. |
| Nature | The cleanest test of the core BLAO hypothesis (measurement geometry). | A new decoder research problem, labelled non-paper; never "AWQPE". |
| Comparators | Frozen P6 D2 variants (bridge_s1_A2, ext_v1_eig_A2), which already place overlap by *global* information gain at any boundary, plus the P5 eig-shots baseline. | awqpe_ext + lowerhalf_gated (frozen), and faithful AWQPE eps_safe with P5 eig shots. |
| Mechanisms | O-ext and O-bridge. | O-ext only (no D1 stitching rule for bridge blocks, as in P6). |
| Claims | Its own primary family. | Its own primary family, and only if its deterministic validation (§4.4) passes. |

**Rules against conflation.**
- A Track-A result says nothing about Track B, and vice versa.
- A Track-B decoder is never evaluated with a Track-A trigger's success used as evidence, or the reverse.
- No arm mixes D2 decoding with awqpe_bl.

---

## 3. Track A: BLAO-D2, exact signal definitions

**Common protocol.**
- The decision process is P6's B4-style arm: information-gain shots on chunks, plus at most A overlap batches.
- At each step where an overlap batch is budget-eligible, compute a score s_j for each internal boundary j = 0..B−2 from the observables at time t.
- Candidate set: C_t = { j : s_j > τ }. If C_t is empty, the step is a chunk-shot step chosen by the frozen P5 eig_cell signal. Otherwise overlap is placed at j_t = argmax_(j∈C_t) s_j, with ties broken by seeded jitter.
- Every signal has exactly one scalar constant, chosen on dev only.

### BL-1: Local ambiguity (counts only)

For block i, let t1_i and t2_i be its two most frequent outcomes (seeded jitter for ties), and ρ_i = N_i(t2_i) / N_i(t1_i) ∈ [0,1].

  U_j = ρ_j · 1[ t2_j ≡ t1_j ± 1 (mod M_j) ]
    (the chunk-j block split between neighbouring values, i.e. r_j near 1/2)

  L_j = ρ_(j+1) · 1[ {t1_(j+1), t2_(j+1)} = {M_(j+1)/2 − 1, M_(j+1)/2} ]
    (the chunk-(j+1) block split across its MSB flip, i.e. β_j itself undecided)

  s_j^(BL1) = max(U_j, L_j),  τ = τ1.

- **Observable.** The counts of blocks j and j+1 only.
- **Boundary selection.** The largest local split ratio.
- **Why it is not a hard-coded predicate.** It is a continuous statistic of two blocks' counts with one threshold. It does not require any exact decoded value of any chunk.

### BL-2: Adjacent-window consistency (counts plus the known kernel)

Two independent local posteriors on the same bit β_j, each with a flat prior and a 64-point grid per unit:

- **Upper window (chunk j).** Write M_j δ_j = t1_j + u, with u ∈ [−1/2, 1/2).
  - L_up(u) = Π_y K_(M_j)((t1_j + u − y)/M_j)^(N_j(y)).
  - Since r_j = u mod 1, β_j = 1 ⇔ u < 0.
  - P_up = P(u < 0 | counts of block j).
- **Lower window (chunk j+1).** On r ∈ [0,1): L_low(r) = Π_y K_(M_(j+1))(r − y/M_(j+1))^(N_(j+1)(y)).
  - P_low = P(r ≥ 1/2 | counts of block j+1).

  s_j^(BL2) = P_up (1 − P_low) + (1 − P_up) P_low,  τ = τ2.

This is the probability that the two windows' independent readings of β_j disagree. It is high when either is uncertain, and highest when both are confident but contradict each other.

- **Observable.** The counts of blocks j and j+1, and the kernel.
- **Boundary selection.** The largest disagreement probability.
- **Known weakness, to be measured in the pilot, not assumed.** The upper window's information about r_j comes only from sub-bin peak shape. It may be weak at small m_j or S0, pushing P_up toward 1/2 and making s^(BL2) ≈ 1/2 almost regardless of P_low.

### BL-3: Boundary-local expected information gain (the D2 posterior)

- Let π_t(g) ∝ exp ℓ_t(g). Label each grid point b_j(g) = 1[frac(2^(K_j) φ_g) ≥ 1/2]; these labels are precomputed once per (G, K_j).
- For the overlap block O_j at boundary j, with outcome law p_(O_j)(y|g), define the one-shot mutual information with the boundary bit:

  I_j = H( Σ_g π_t(g) p_(O_j)(·|g) ) − Σ_(b=0,1) π_t(b_j = b) · H( Σ_(g: b_j(g)=b) π_t(g | b_j = b) p_(O_j)(·|g) ).

- This is computed on the posterior support (tail mass ≤ 10^−6), as in eig_cell_support.

**Two-stage selection.** Let E_t* = max_i EIG_cell(chunk i), the frozen P5 signal.
- The score is s_j^(BL3) = I_j / E_t*. Overlap is placed at argmax_j s_j if s_j > κ; otherwise the step is a chunk-shot step.
- The ratio is compared with a threshold, not merged into one argmax. The two informations have different targets, which makes that comparison explicit and tunable (one constant, κ).

**Relation to the frozen comparator.** ext_v1_eig_A2 and bridge_s1_A2 score overlap blocks by information about the whole n-bit cell. BL-3 scores each overlap block only for the bit it exists to resolve.
- **Track A therefore tests a specific hypothesis:** a boundary-local target beats a global-cell target. It is not just "local beats lowerhalf_gated".

**D2 decoding.** Unchanged. The MAP of the grid posterior over all blocks, with each shot entering the likelihood once.

**Firewall (all three).**
- The signals are functions of counts and the likelihood, implemented under `research/awqpe/overlap/`. The AST leakage test covers that package.
- A truth-raises oracle test will be added.
- j* and the strata are analysis-only.

---

## 4. Track B: can a decoder use O-ext information locally?

### 4.1 Where the current rule needs global information

- Algorithm 2 decides chunk j's borrow from the MSB of the corrected chunk j+1. It fails when the chunk-j block rounds the "wrong" way while unflagged, or when the lower part is exactly 10...0 (the V1 floor).
- `awqpe_ext` consults the widened block only when chunks j+1..B are exactly 10...0 (2^−(n−K_j) of phases in the limit).
- The reason recorded in its docstring: substituting the widened reading unconditionally raised failures to 16-54%, because the widened reading's own rounding can carry.

### 4.2 Two lemmas (infinite-shot limit)

**L1 (widened floor).**
- With x_j = q + r_j, where q = ⌊x_j⌋, the limit reading of O-ext_j is w = round(2^v x_j) = 2^v q + round(2^v r_j).
- Hence ⌊w / 2^v⌋ = q exactly when round(2^v r_j) < 2^v, i.e. when

  r_j < 1 − 2^−(v+1)  (v = 1: r_j < 3/4),

  and q + 1 otherwise.
- No condition on chunks j+2..B enters.

**L2 (target).**
- The best n-bit chunk j equals q mod M_j unless the lower n − K_j bits round up with a carry, which requires r_j ≥ 1 − 2^−(n−K_j+1).
- Since n − K_j ≥ 2 for every internal boundary (all widths are > 1), that threshold is ≥ 7/8 > 3/4.
- Hence **for r_j < 3/4 the substitution chunk j := ⌊w/2^v⌋ mod M_j is exactly correct**, whatever the lower chunks are.

**Local certificate for r_j < 3/4.**
- The chunk-(j+1) block measures r_j directly: in the limit it reads y = round(M_(j+1) r_j) mod M_(j+1).
- Reading y ∈ Y_safe = { ⌈M/4⌉, …, ⌈3M/4⌉ − 1 − μ } (with M = M_(j+1)) implies r_j < (y + 1/2)/M < 3/4.
- It also excludes the wrap reading y = 0, which could mean r_j ≈ 1.
- For M = 16 and μ = 1: Y_safe = {4, …, 10}, so r_j < 10.5/16 ≈ 0.656.
- The lower edge ⌈M/4⌉ restricts substitution to where it matters: r_j away from 0, where the chunk-j block can straddle two values.

### 4.3 The candidate rule (awqpe_bl), stated for review, not implemented

Algorithm 2 runs LSB to MSB, unchanged, except at boundary j (upper chunk j). **If** an O-ext_j block has shots **and** the raw top outcome of the chunk-(j+1) block lies in Y_safe, **then** set chunk j := ⌊t1(O-ext_j) / 2^v⌋ mod M_j with no further borrow. Otherwise apply Algorithm 2's rule.
- **Inputs used at boundary j:** the O-ext_j block, the chunk-(j+1) block's raw top outcome, and the already-decoded state that Algorithm 2 carries. Chunks j+2..B are never read.
- **In the limit it is exact.** Proven for r_j < 3/4 by L1 + L2. Otherwise it reduces to Algorithm 2.
- **With finite shots it is NOT known to help.** It errs if (a) the O-ext top outcome is off by one bin, or (b) the chunk-(j+1) reading falls inside Y_safe while r_j ≥ 3/4. Event (b) needs an error of more than μ + 1/2 bins.
- **Both error probabilities are local and computable from the kernel.** The question is whether, on realistic phases, the rescues outnumber the new errors.
- **Open issue: μ.** μ ∈ {0, 1} is to be fixed by the deterministic validation before any stochastic data. At M_(j+1) = 4 (m = 2), Y_safe = {1, 2} with μ = 0 leaves only a half-bin margin, so the rule may be unsafe for 2-bit chunks. That is a reason to restrict Track B to m ≥ 3, to be decided at validation.
- **Honest status.** Local decoding is *possible in the infinite-shot limit* (L1 + L2 + certificate). Whether it is *beneficial at finite shots* is exactly what Track B would test.

### 4.4 Deterministic validation (must pass before any Track-B pilot)

1. With no O-ext blocks: equal to faithful AWQPE, bit for bit.
2. On every phase where chunks j+1..B are exactly 10...0: equal to awqpe_ext. This case lies inside Y_safe, so awqpe_bl strictly generalises the frozen extension.
3. Infinite-shot limit on the V1 failure set: never worse than awqpe_ext.
4. Infinite-shot limit on a dense generic grid: zero new failures (the carry check of L1 and L2).
5. Locality: permuting or resampling the counts of chunks j+2..B does not change the decision at boundary j (unit test).
6. Exact computation of the finite-shot error terms (a) and (b) from the kernel for m ∈ {3,4} and S ∈ {4,8,16}, reported before the pilot.

---

## 5. Worked example: n = 16, [4,4,4,4]

**The phase.** Constructed for illustration; it is not from any P6 pool. The numbers come from the existing decoders (awqpe_vectorised, awqpe_ext_decode, lowerhalf_mask, likelihood_decode): 200,000 trials at S = 16 shots per block, eps = 0.9, and 2,000 trials for D2. This is not a pilot.

  2^16 φ = 5·4096 + 2051.3,  φ = 0.343800354…
  N* = 22531 = 0101 1000 0000 0011  → chunks (5, 8, 0, 3).

**Boundary 0 (upper).** x_0 = 16φ = 5.50081, so q = 5 and r_0 = 0.50081 (just above 1/2).

| Block | Reads | Most likely outcomes (probability) |
|---|---|---|
| chunk 0 (k = 0) | 16δ = 5.5008 | 6 (0.408), **5 (0.405)**, 7 (0.046). Split almost exactly between 5 and 6. |
| chunk 1 (k = 4) | 16δ = 8.0129 | 8 (0.9995). Sharp: β_0 = 1 is locally certain. |
| chunk 2 (k = 8) | 16δ = 0.206 | 0 (0.868) |
| chunk 3 (k = 12) | 16δ = 3.300 | 3 (0.738) |
| O-ext_0 (k = 0, m = 5) | 32φ = 11.0016 | 11 (0.99999). Sharp. |

**Why the frozen rule fails.**
- The lower part below boundary 0 is 1000 0000 0011, not 1000 0000 0000. Below boundary 1 it is 0000 0011, not 1000 0000. Below boundary 2 it is 0011, not 1000.
- So lowerhalf_gated fires at **0.0%** of trials at every boundary.
- Even if an O-ext_0 block were measured, awqpe_ext would ignore it. Its output was identical to faithful AWQPE in 200,000/200,000 trials.

**Faithful AWQPE (eps 0.9) success is 55.8%.** The failure is entirely one case:

| Case | Probability | Algorithm 2 result | Success given the case |
|---|---|---|---|
| A: t1 = 6, unflagged | 0.450 | 6 − 1 (borrow, since chunk-1 MSB = 1) = 5 | 0.999 |
| **B: t1 = 5, unflagged** | **0.440** | **5 − 1 = 4** | **0.000** (error 4096 LSB, a higher-order failure at j* = 0) |
| C: flagged | 0.109 | min(5,6) = 5, no borrow | 0.994 |

**What a genuinely local decoder needs, and has, here.**
- (i) β_0 from the adjacent chunk-1 block, which reads 8 with probability 0.9995.
- (ii) An unambiguous reading of chunk 0's integer part, which O-ext_0 gives: ⌊11/2⌋ = 5 in 100% of simulated trials.
- The certificate holds: the chunk-1 reading 8 lies in Y_safe = {4..10}, so r_0 < 3/4. By L1 + L2, chunk 0 = 5 is then exact.
- Chunks 2 and 3 are **not needed**. The frozen rule's requirement that they be 0000 and (at the n-bit level) 0000 is exactly what excludes this phase.

**D2 contrast (Track A context).** D2 decodes this phase correctly in 99.4% of trials with no overlap at all. Its likelihood uses chunk 0's counts at 5 *and* 6 jointly with chunk 1's sharp reading.
- So this failure is **specific to the faithful decoder** and is a Track-B example.
- Track A's targets are D2's residual failures: posterior mass split between hypotheses that differ at an upper boundary, for example φ vs φ ± 2^−K_j under low shot counts.
- There an overlap block measuring the bits around that boundary can discriminate the modes: O-bridge_0, at offset 3, reads frac(8φ), which depends on q's LSB.
- Whether local targeting does this better than the frozen global-cell eig is the Track-A question.

---

## 6. Falsification criteria (declared before any data)

BLAO, or the affected track, is **rejected** if any of the following holds.

**F1: Pool dependence.** Any result that holds on the already-inspected P6 phases (dev or test) but not on the fresh pools of §7 is void. The P6 phases are never used for BLAO selection, tuning or claims.

**F2: No equal-cost gain (Track A).** On the fresh dev pool, no BL signal beats both frozen P6 D2 variants on success at equal U-queries (upper bracket, D-022 matching, phase-cluster CI including 0).
- The equal-shot result is reported but cannot rescue the claim.

**F3: No equal-cost gain (Track B).** awqpe_bl (with any trigger) does not beat awqpe_ext + lowerhalf_gated at equal U-queries (upper bracket) on the fresh dev pool.
- Or it fails any item of §4.4 at validation, which rejects it before any stochastic run.

**F4: Hidden hard-coded predicate.** A trigger's firing decisions are essentially an exact decoded-value predicate. Operationally, on the pilot or dev pool, Cohen's κ ≥ 0.9 between its fire/no-fire decisions (per boundary and step) and either:
- (a) `lowerhalf_mask`, or
- (b) "the decoded chunk j+1 equals M/2", or
- (c) "raw chunk j+1 ∈ Y_safe".

Such a trigger is treated as a decoder-specific predicate, not a boundary-local statistic, and is dropped. (awqpe_bl's own Y_safe certificate is a *decoder* rule derived from L1 + L2, fixed before data. It is not a trigger and is not subject to F4.)

**F5: Fresh-pool failure.** The Track-A dev winner shows no gain on the held-out BLAO test pool, using the same primary comparison: phase-cluster CI including 0, or Holm p ≥ 0.05. The same applies to Track B.

**F6: Placement without payoff.** A trigger raises placement at j* (analysis-only label) without raising success at equal U. This is reported as a negative result. It shows the P6 conditional rescue rate does not transfer.

**F7: Budget dilution.** The overlap budget is exhausted (≥ 95% of trials use all A batches) while success at equal U does not improve. The signal fires indiscriminately.

---

## 7. Minimal pilot (DESIGN ONLY; not implemented, not run)

**Purpose.** Feasibility and instrumentation only: firing rates, placement, runtime, memory, matching validity. **No selection, tuning or claim** is made from pilot data.

**Fresh phase pool.**
- `make_phase_table` with a new master seed, 20261010, reserved only for BLAO; the split is "dev". This gives a seed domain disjoint from every earlier phase.
- Before running, an assertion checks that no pilot φ equals any φ in the P6 dev or test phase tables (or the P7-P9 tables) to within 2^−40. The phases are logged.
- The later BLAO dev and test pools use further new master seeds (20261011 dev, 20261012 test). The P6 test selections and tables are never reused.

**Common settings.** These are the P6 cell settings, not re-tuned:
- partitions [4,4,4] (n = 12) and [4,4,4,4] (n = 16);
- strata S1_final_half, S2_boundary, S4_uniform;
- 4 phases per stratum and 3 replicates, giving 36 trials per partition per arm;
- S0 = 4, dS = S0, budget multiplier r = 2, A = 2.

**Provisional constants.** τ1 = 0.5, τ2 = 0.25, κ = 1, μ = 1. They are fixed a priori and not tuned on the pilot; dev will select them later on its own pool.

**Arms.**

| Track | Arm | Decoder | Overlap trigger | Mechanism |
|---|---|---|---|---|
| (ref) | B2: P5 eig shots, no overlap | D2 / faithful eps_safe | none | none |
| A | frozen comparator | D2 | ext_v1_eig_A2 (global eig) | O-ext |
| A | frozen comparator | D2 | bridge_s1_A2 (global eig) | O-bridge |
| A | BL-1 | D2 | s^(BL1) > τ1 | O-ext and O-bridge (separate arms) |
| A | BL-2 | D2 | s^(BL2) > τ2 | O-ext and O-bridge (separate arms) |
| A | BL-3 | D2 | s^(BL3) > κ | O-ext and O-bridge (separate arms) |
| B | frozen comparator | awqpe_ext | lowerhalf_gated | O-ext |
| B | awqpe_bl | awqpe_bl | lowerhalf_gated | O-ext (isolates the decoder change) |
| B | awqpe_bl | awqpe_bl | BL-1 | O-ext |
| B | awqpe_bl | awqpe_bl | BL-2 | O-ext |

- That is 11 overlap arms plus the reference.
- Each overlap arm is paired with its equal-U lower and upper brackets (D-022) against B2.
- Track B arms run only after §4.4 passes.

**Pilot outputs (instrumentation, no inference).**
- Firing rate per boundary j and step.
- The distribution of the placement boundary.
- The placement rate at j* (analysis-only label).
- Overlap batches used, and the fraction of trials exhausting A (for F7).
- κ-agreement with the predicates of F4.
- The upper-window P_up distribution (the BL-2 weakness).
- Equal-U bracket validity rates.
- Runtime and peak memory per shard.
- Invariant checks: the leakage test, bit-identity of the frozen comparators with the P6 code paths, and awqpe_bl = awqpe_ext on the 10...0 subset.

**Gate from pilot to dev.** It passes if all of the following hold:
- no crashes;
- the invariants pass;
- no trigger has a firing rate of 0% or 100% at every boundary;
- runtime is within budget.

Performance numbers from the pilot are not used for any decision.

**Estimated cost.** 2 partitions × 36 trials × 12 arms, plus brackets, at the P6 per-trial cost: under 10 minutes on the existing CPU runner.

---

## 8. Explicitly not done

- No BLAO or awqpe_bl code.
- No pilot, dev or test runs.
- No change to P5-P10b artefacts or selections.
- No CUDA or system-level installs.
- No claim that BLAO works.
