# Boundary-Local Adaptive Overlap (BLAO): design proposal (FOR REVIEW; nothing implemented)

- **Status.** DESIGN ONLY. No code, no pilot, no dev or test runs.
- **Frozen phases.** P5-P10b are frozen and must not be modified. BLAO, if approved, gets its own cycle: deterministic validation, pilot, dev, freeze, then a held-out test on a NEW phase pool (see section 6).

## 0. Motivation, and a correction to how the P6 follow-up should be read

**What the P6 follow-up found** (RESULTS_LOG, "why extended-decoder overlap rescue declines at n = 16"):
- In the frozen variant `ext_v1_gated_A2` with `awqpe_ext`, overlap is almost never placed at the boundary responsible for the error at upper boundaries: 0% at boundaries 0 and 1 when n = 16.
- When it is placed there, 91-100% of those failures are rescued.

**Correction: that rescue rate is conditional and must not be assumed to transfer.** Both the trigger (`overlap/triggers.py::lowerhalf_mask`) and the decoder's substitution rule (`decode/awqpe_overlap.py::awqpe_ext_decode`) require the same predicate: the already-corrected lower part (chunks j+1..B) equals exactly 10...0.
- So "placed at j*" in that analysis implies "the decoder's rule could act".
- The 91-100% is a rescue rate *given that predicate*. It does not show that overlap placed at j* for another reason would be used, or would help, under the current decoder.
- With `awqpe_ext` unchanged, a boundary-local trigger that fires where the predicate fails would place overlap blocks the decoder ignores.
- **BLAO therefore has two separable parts:** (i) where to place overlap (a trigger), and (ii) for D1-type decoding, a boundary-local rule for using it (a new decoder).
- For D2 only (i) is needed, because the likelihood uses every block. D2 overlap rescue did not decline with n in P6.

## 1. Notation (as in the code)

**Chunks and boundaries.**
- Chunks i = 0..B−1 are MSB-first, with widths m_i, offsets k_i = Σ_{l<i} m_l, and M_i = 2^(m_i).
- Boundary j (j = 0..B−2) lies between chunk j and chunk j+1, with K_j = k_{j+1} bits above it.

**Block law** (Lemma 3.3 of the paper; `model/kernel.py`).
- Chunk i's block sees δ_i = frac(2^(k_i) φ) and returns y with probability K_(M_i)(δ_i − y/M_i).

**The boundary variable.**
- r_j = frac(2^(K_j) φ) is the remainder below boundary j, and note r_j = δ_(j+1).
- The borrow bit is β_j = 1[r_j ≥ 1/2]. Algorithm 2 decides it from the MSB of the corrected chunk j+1. Chunk j's raw estimate rounds 2^(m_j) δ_j, so it is off by one exactly when r_j ≥ 1/2 and the borrow is not applied.

**Observables at any time.**
- Counts N_i(y) per chunk block, and N_ov(y) for placed overlap blocks.
- Top-two outcomes t1_i, t2_i and ratio ρ_i = N_i(t2_i)/N_i(t1_i).
- D1 raw and corrected chunks and flags (`awqpe_vectorised`, eps 0.9 or eps_safe).
- The D2 grid log-likelihood ℓ(φ_g), g = 0..G−1.

**Overlap actions** (unchanged from P6, `overlap/candidates.py`).
- O-ext at j: chunk j re-measured with width m_j + v, where v = 1.
- O-bridge at j: a width-m_j block at offset K_j − s, where s = 1.
- A batch is S0 shots. At most A overlap batches per trial.

## 2. Candidate signals

### BL-1: Local ambiguity

**Definition.** Using only chunks j and j+1:
- Upper term: U_j = ρ_j · 1[t1_j, t2_j adjacent mod M_j]. This is the chunk-j block itself being split between two neighbouring values, the signature of r_j near a rounding point of chunk j.
- Lower term: L_j = ρ_(j+1) · 1[{t1_(j+1), t2_(j+1)} = {M_(j+1)/2 − 1, M_(j+1)/2}]. This is the chunk-(j+1) block split exactly across its MSB flip, i.e. the borrow bit β_j itself undecided.
- The score is S1_j = max(U_j, L_j).
- **Bayesian variant, to be decided on dev.** Replace each ratio with the Beta-posterior probability that the second outcome's true probability exceeds 0.5 × the first's. This is less noisy at small counts.

**Observables.** The counts of blocks j and j+1 only. No other chunk, no posterior.

**Decision rule.** Among boundaries with S1_j > τ1, place overlap at argmax_j S1_j, subject to the budget A. τ1 is a single threshold chosen on dev only.

**Overlap action.** O-ext (v = 1) or O-bridge (s = 1) at j. The mechanism is fixed per arm and never mixed.

**Decoder requirements.**
- D2: none.
- D1-type: the new boundary-local rule of section 3 is required; `awqpe_ext` would ignore these blocks.

**Cost.** O(Σ_i M_i) per trial per step (top-two by partial sort). Negligible.

**Firewall.** It is a function of counts only. It lives in `overlap/`, so it falls under the AST leakage test.

**Difference from lowerhalf_gated.**
- lowerhalf_gated needs every lower chunk j+1..B to decode to exactly 10...0: a global and very restrictive condition, which at boundary 0 of [4,4,4,4] means 12 exact bits.
- BL-1 looks only at the two chunks adjacent to the boundary.
- It also reacts to statistical ambiguity (split counts) rather than to one exact decoded value.

### BL-2: Adjacent-window consistency

**Definition.** Both neighbouring blocks carry information about the same bit β_j:
- The upper block, through the sub-bin shape of its peak: the relative counts at t1_j and its neighbours locate 2^(m_j) δ_j within the bin, and hence r_j at coarse resolution.
- The lower block, directly, through the MSB of 2^(m_(j+1)) δ_(j+1).

Each window gets its own local 1-D posterior with a flat prior:
- **Upper.** On u ∈ [−1/2, 1/2), where 2^(m_j) δ_j = t1_j + u, with likelihood Π_y K_(M_j)((t1_j + u − y)/M_j)^(N_j(y)). Here P_up = P(β_j = 1 | upper) = P(u < 0 | upper). Rounding up from below means the remainder is above 1/2.
- **Lower.** On δ_(j+1) ∈ [0, 1) with its own kernel likelihood, giving P_low = P(δ_(j+1) ≥ 1/2 | lower).

The score is the probability that the two windows disagree on β_j under independence:

  S2_j = P_up (1 − P_low) + (1 − P_up) P_low.

It is high when either window is uncertain, and highest when both are confident but disagree, which signals an error in one of them.

**Observables.** The counts of blocks j and j+1, and the kernel (known model). No global posterior.

**Decision rule.** Among boundaries with S2_j > τ2, place overlap at argmax_j S2_j, within A. τ2 is chosen on dev.

**Overlap action.** O-ext or O-bridge at j. O-bridge is the natural action for BL-2, because a bridge block straddles the boundary and measures β_j together with the adjacent bits of both chunks.

**Decoder requirements.** As BL-1: none for D2; D1-type needs section 3.

**Cost.**
- Two 1-D posteriors on sub-grids of about 64 points per boundary: O(B · (M_j + M_(j+1)) · 64) per trial per step. Negligible.
- The upper-window posterior uses only the neighbourhood of t1_j (for example ±2 outcomes), to stay cheap and local.

**Firewall.** Counts plus the known kernel only; in `overlap/`.

**Caveat, to verify in the deterministic validation.** The upper block's sub-bin information about r_j is weak for small m_j and few shots, so P_up may sit near 1/2 almost always. BL-2 could then collapse to a function of P_low alone. The validation must report the distribution of P_up and the fraction of decisions it changes relative to P_low alone.

**Difference from lowerhalf_gated.** It uses statistical evidence from both sides of the boundary and never requires exact decoded values of chunks below j+1.

### BL-3: Boundary-local expected information gain

**Definition.** Under the current D2 grid posterior π(g) ∝ exp ℓ(g), the target is the local variable β_j(g) = 1[frac(2^(K_j) φ_g) ≥ 1/2]. For the overlap block O_j at boundary j, with outcome Y ~ p(y | φ_g):

  BL3_j = I(β_j ; Y) = H(Σ_g π(g) p(·|g)) − Σ_(b∈{0,1}) π(β_j = b) H(Σ_(g: β_j(g)=b) π(g|b) p(·|g)).

Optionally report it per U-query, BL3_j / q(O_j).
- **Closest existing signal.** The P5/P6/P7 `eig_cell` policy scores the mutual information of a block with the global n-bit cell (or the tau-cell in P7). That pools evidence about all boundaries into one target.
- **What BL-3 does instead.** It scores each overlap candidate only for the bit it exists to resolve. The chunk-shot actions keep the global signal, so shots and overlap compete on commensurate terms only through the decision rule below.

**Observables.** The D2 grid log-likelihood and the known kernel.

**Decision rule** (a two-stage choice that keeps allocation unchanged):
- (a) Choose the best chunk action by the frozen `eig_cell` signal.
- (b) Take overlap at j* = argmax_j BL3_j instead, if BL3_(j*) > κ · EIG_cell(best chunk) and the overlap budget allows. κ is one scalar chosen on dev.
- **Rejected alternative.** Folding BL3 into a single argmax with eig_cell was rejected, because the two informations have different targets and are not directly comparable.

**Overlap action.** O-ext or O-bridge at j*.

**Decoder requirements.** For D2, none; BL-3 is the natural D2 trigger. For D1-type, section 3.

**Cost.**
- For each candidate: a posterior-predictive over M_ov outcomes on the posterior support (K_s points), O(T · K_s · M_ov). This is the same order as eig_cell_support, which was feasible at n = 16 in P6/P7.
- The β_j labels per grid point are precomputed once per (G, K_j).

**Firewall.** It uses the posterior only, not truth. It lives in `overlap/`, next to or reusing `allocation/eig_cached.py` machinery.

**Difference from lowerhalf_gated.** It has no structural precondition at all. It places overlap wherever the posterior says a boundary bit is uncertain *and* the overlap block would resolve it.

## 3. Decoder requirement for D1-type decoding: "awqpe_bl" (a new, labelled non-paper decoder)

`awqpe_ext` substitutes the widened reading at boundary j only when the lower part is exactly 10...0. Its docstring records why unconditional substitution is harmful: the widened reading's own rounding can carry into chunk j's bits with probability about 2^−(v+1), which gave 16-54% failures in P6 development. A boundary-local rule must therefore say *when* the widened reading is more trustworthy than Algorithm 2's borrow.

**Candidate rule, to be specified exactly and validated before any pilot.** At boundary j, with a widened block for chunk j present:
- Let w = t1_ext. Its extra bit e = w mod 2^v (v = 1) is an independent reading of β_j, with its own carry risk when r_j is near 1 − 2^−(v+1).
- Let b = the MSB of the corrected chunk j+1, which is Algorithm 2's reading of β_j.
- If e and b agree: keep Algorithm 2's result.
- If they disagree AND the local evidence favours the widened block, substitute chunk j := floor(w / 2^v) mod 2^(m_j). The favouring condition is BL-2's P_low near 1/2, or chunk j+1's split across its MSB (BL-1's L_j above threshold).
- Otherwise keep Algorithm 2's result.

**Required properties** (a deterministic validation suite, like the A-E suite for awqpe_ext):
1. With no overlap blocks, it equals faithful AWQPE exactly.
2. When the lower part is 10...0, it equals awqpe_ext exactly, so it strictly generalises the frozen extension.
3. In the infinite-shot limit, on the V1 failure set, it is never worse than awqpe_ext.
4. On a generic dense grid in the infinite-shot limit, it introduces no new failures. This is the carry-risk check.

**Labelling.** It must always be labelled "adaptive overlap + boundary-local extended decoder", never AWQPE.

**Scope.** O-bridge stays D2-only, as in P6. No stitching rule is invented for bridge blocks under D1.

## 4. Summary comparison

| | lowerhalf_gated (frozen P6) | BL-1 | BL-2 | BL-3 |
|---|---|---|---|---|
| Information used | exact corrected values of chunks j+1..B | counts of chunks j, j+1 | counts of chunks j, j+1 plus the kernel | the global D2 posterior |
| Scope | global (whole lower part) | local | local | local target, global evidence |
| Fires at upper boundaries for n = 16 | almost never (0% at j* in P6 test) | yes, if split counts | yes, if the windows disagree or are uncertain | yes, if β_j is uncertain and resolvable |
| Tuned constants | none | τ1 | τ2 | κ |
| Cost per step | one D1 decode | O(ΣM) | O(B·M·64) | O(T·K_s·M_ov) per candidate |
| D2 needs a new decoder | no | no | no | no |
| D1 needs a new decoder | no (awqpe_ext) | yes (awqpe_bl) | yes | yes |

## 5. Hypotheses and how each could fail (not assumed)

- **H1 (placement).** BL triggers place overlap at the responsible boundary more often than lowerhalf_gated at n = 12 and 16.
  - Falsified if the dev placement rate at j* (analysis-only truth label) is not higher.
- **H2 (D2 value).** With D2, a BL trigger beats the frozen P6 D2 variants at equal U (upper bracket, D-022 matching).
  - It may well fail. P6 D2 overlap gains were small (+0.56 to +0.60 pts), and BL-3 may be nearly equivalent to the existing `eig` candidate policy.
- **H3 (D1 value).** awqpe_bl with a BL trigger beats awqpe_ext + lowerhalf_gated at equal U, especially at n = 16.
  - It fails if the carry risk of section 3 outweighs the gain, or if local signals fire too often and waste the A budget.
- **Risk: budget dilution.** Local triggers fire more often. Under a fixed A and equal-U accounting they can lose, even with better placement. This is why the equal-U upper bracket stays primary.

## 6. Evaluation protocol sketch (for approval; not started)

1. **Deterministic validation.** For awqpe_bl: properties 1-4 of section 3. For BL-1/2/3: unit tests on constructed count vectors; the leakage AST test; and a truth-raises oracle test.
2. **Pilot.** At about 1% scale, report the firing rate per boundary, the placement rate at j*, runtime and memory.
3. **Dev.** Select one signal and its constant per (mechanism, decoder) by a rule predeclared before dev data.
   - Primary: B4-style BLAO arm vs B2m_upper at equal U.
   - Comparators: the frozen P6 variants run on the same dev phases.
4. **Freeze.** A new D-record with the config hash.
5. **Held-out test on a NEW phase pool.** The P6 test phases were inspected post hoc (the n = 16 analysis), so they are no longer clean.
   - BLAO's test must use a fresh seed domain: a new master seed and split tag, e.g. `test_blao`.
   - It must not reuse or modify the P6 held-out selections or tables.
   - The frozen P6 variants are re-run on the new pool only as comparators.

## 7. What is explicitly NOT being done now

- No BLAO or awqpe_bl code.
- No experiments.
- No change to P5-P10b artefacts.
- No claim that BLAO works.
