# Decisions log

Entries are append-only. A decision that changes later gets a new entry that references the old one; the old entry is not edited.

---

### D-001 (2026-09-27): Separate branch, worktree and package
- **Decision.**
  - New branch `research/awqpe-unified-adaptive`, created from tag `paper-a-v1.0-frozen` (90bc134).
  - It is checked out as a git worktree at `../QPR-awqpe`, because the main checkout holds uncommitted Paper A/B edits.
  - All new code lives in `research/awqpe/`.
  - Old modules are never edited.
- **Why.** This preserves the frozen Paper A and keeps the user's in-progress manuscript edits untouched.
- **Rejected.** Putting new folders at the repo root: `experiments/` would collide with `Experiments/` on case-insensitive Windows.

### D-002 (2026-09-27): Separate environment `.venv-awqpe`
- **Decision.**
  - Create a Windows venv from Python 3.14 with the Paper A pins (qiskit 2.4.1, qiskit-aer 0.17.2, numpy 2.4.6, scipy 1.17.1, pandas 3.0.3, matplotlib 3.11.0), plus pyarrow 25.0.1, PyYAML 6.0.3 and pytest 9.1.1.
  - Pin these in `research/requirements-awqpe.txt`.
- **Why.** The repository `.venv` is a WSL venv without Aer. `C:\Python314` is Paper A's reproducibility environment and must not change.

### D-003 (2026-09-27): Kernel implemented in Fejér form
- **Decision.** Evaluate K_M(θ) = M⁻² Σ_{|d|<M} (M−|d|) cos(2πdθ).
- **Why.**
  - It is algebraically identical to the paper's sin² ratio: the tests find agreement to 1e-13, and agreement with the Qiskit statevector to 7e-15.
  - It has no removable singularity.
  - It gives exact K′ and K″ for Fisher information.
  - Gaussian phase jitter becomes an exact per-frequency damping.
- **Not an approximation.** No approximation of the paper's model is introduced.

### D-004 (2026-09-27): Interpretations where arXiv v3 is ambiguous
- **(a) Ties** are broken uniformly at random, using U(0, 0.5) jitter drawn from a seeded generator.
- **(b) Eq. 2.1** is silent on non-wrap pairs that contain 0 or n−1. These use the ordinary minimum.
- **(c) Flag indexing.** Algorithm 2 line 17 reads chunk j's flag as `A[j−1]` in 0-indexed Python.
- **(d) Borrow source.** The borrow reads chunk j+1's most-significant bit after chunk j+1's own correction (a verbatim in-place loop).
- **(e) Final block.** If the final block is ambiguous, its flag is set and b_ml stays t1.

### D-005 (2026-09-27): Two decoders as a crossed factor
- **Decision.**
  - D1 = faithful AWQPE (Algorithms 1+2).
  - D2 = grid-likelihood MAP over φ (G = 2^(n+4)) using all counts of all blocks.
  - Every allocation or overlap comparison is reported within each decoder.
- **Why.** Otherwise a gain from better decoding could be attributed to adaptive allocation. D2 is also the principled way to fuse overlapping blocks.
- **Status.** D2 is not AWQPE and is never labelled as such.

### D-006 (2026-09-27): Primary budget = total shots (user decision)
- **Decision.** Total shots is the primary budget. Controlled-U queries, Σ shots·2^k(2^m−1), are always co-reported. Cost-aware policy variants are included.

### D-007 (2026-09-27): Special-chunk ablation decoder
- **Decision.** `awqpe_ablate_special` disables only the S_idx borrow suppression.
- **Why.** It was added to test causally whether the special-chunk rule explains the infinite-shot failures (it does not; see RESULTS_LOG P2a).
- **Status.** Diagnostic only; never called AWQPE.

### D-008 (2026-09-27): Tolerance metric for baseline summaries
- **Decision.** Baseline summaries report P(|φ̂−φ|_circ ≤ 2⁻ⁿ) and exclude stratum S1 from exact-best-n-bit summaries.
- **Why.** In S1 the best n-bit value itself is within 5% of a tie, so "exact" is ill-defined there. Tolerance 2⁻ⁿ accepts either neighbouring grid point.
- **Scope.** The full tolerance study (τ grid) is Phase 8.

### D-009 (2026-09-27): Dev/test split by seed domain
- **Decision.**
  - Phases for tuning come from `split=dev`.
  - Held-out claims come from `split=test`, which is generated from a disjoint seed domain.
  - The paper's example phases (S5) appear only in dev.
- **Rule.** Phases 2–4 are characterisation runs and use dev only. Test is first touched after policy constants are frozen (Phase 5+).

### D-010 (2026-09-27): Independent verification implementation
- **Decision.** `research/awqpe/verification/awqpe_strings.py` is a second, string-based implementation of Algorithms 1–2 with switchable readings of every ambiguity.
- **Status.** Verification only. The faithful baseline (`baseline/awqpe.py`) is unchanged, and no reading is substituted into it.

### D-011 (2026-09-27): Safe-epsilon baseline, a separate labelled condition
- **Rule, fixed a priori with no tuning on any split:**

  ε_safe(partition) = floor_{0.01}(min_j ε*(k_j, m_j)) − 0.02

  where ε*(k, m) = K_M((½ + 2^−(k+1))/M) / K_M((½ − 2^−(k+1))/M), k = bits below boundary j, and m = width of the block above it.
- **Examples.** [4,4] → 0.75; [3,2,3] → 0.60; [2,2,2,2] → 0.35; [6,6] → 0.91.
- **Verification.** Zero infinite-shot failures at ε_safe for all 13 swept partitions (`verify_safe_epsilon.csv`).
- **Usage.** Always reported *alongside* the faithful ε=0.9 condition, never replacing it. Its finite-shot behaviour (more flags from sampling noise) is measured, not assumed.

### D-012 (2026-09-27): Terminology for the ε floor
- **Decision.** Call it a "potential decoder limitation under the arXiv-v3 implementation" until the published version is checked. Do not call it an "AWQPE flaw".

### D-013 (2026-09-27): Circuit-level verification uses explicit controlled matrices
- **Decision.** `general_unitary_block` applies controlled-U^(2^(k+p)) as an explicit (2·dim)×(2·dim) `UnitaryGate` on [control, *targets] rather than `UnitaryGate(U^p).control(1)`.
- **Why.** Identical mathematics, without Qiskit's synthesis of large controlled unitaries: about 70× faster and more accurate (1e-15 against 1e-13).
- **Status.** This is a verification construction. Depth-realistic circuits for noise are a Phase 10 concern.

### D-014 (2026-09-27): Mixture success metrics
- **Decision.** For non-eigenstate inputs, report P(hit a nonzero eigenphase s/r) and P(output = 0) separately.
- **Why.** "Hit any eigenphase" is inflated by the trivial s = 0 component, which the min(t1, t2) rule selects under ties.

### D-015 (2026-09-27): Signal definitions fixed a priori
- **Boundary risk.** Uses a half-width of 0.15 outcome units.
- **Local posterior grid.** 2^(m+6) points, uniform prior on δ_b.
- **Global information gain.** eig_phi and eig_cell are the exact mutual information (bits) of one extra shot, taken from the D2 grid posterior (G = 2^(n+4)).
- **Legacy C_w.** Imported unchanged from the frozen Paper A code.
- **Status.** None of these were tuned.

### D-016 (2026-09-27): P4 analysis plan, declared before the full P4 results were inspected
- **Signal orientation ("higher = more need for shots").** This is fixed in advance:
  - Oriented as-is: ratio, entropy, local_post_sd, chunk_risk, boundary_risk, eig_phi, eig_cell and fisher_phi, i.e. a Fisher-driven allocator gives shots where Fisher information is largest.
  - Negated: c1, margin, p_chunk_correct and legacy_cw.
  - global_risk is trial-level; it is used for stopping, not for choosing a block.
- **Units.** Unit = (trial, block) pair within a decoder and cell (partition, S0, ΔS).
- **Metrics:**
  - Spearman ρ between the oriented signal and Δerr = err0 − err1.
  - AUROC for "fix" (tol0 = False and tol1 = True, τ = 2⁻ⁿ) over all pairs, and within initially failing trials only (the "which block to fix" question).
  - Greedy one-step policy value: P(tol1) when the batch goes to argmax(signal) in each trial, against a random block (mean over blocks) and the oracle-best block.
  - Calibration of p_chunk_correct against t1_correct, in 10 bins.
- **Breakdowns.** By stratum, block position, n and m.
- **Inference.** Phase-cluster bootstrap CIs (1000 resamples) for the policy-value differences.
- **Status.** P4 is dev-split characterisation. Its conclusions choose which signals go into P5; they are not final claims.

### D-017 (2026-09-27): Signals carried into P5 (chosen from P4 evidence)
- **Primary adaptive signal:** eig_cell.
- **Simple comparator:** entropy, representing the interchangeable count/top-two family.
- **Legacy comparator:** C_w.
- **Negative control:** Fisher information. Keep it in P5 to test whether one-step harm persists under multi-step allocation, and because P1-F2 predicts it becomes phase-dependent under noise (P9).
- **Reference arms:** uniform (budget-matched, tuned on dev) and oracle (reference only, `is_oracle=True`).
- **Decoders stay separate arms:** uniform+D1 vs adaptive+D1, and uniform+D2 vs adaptive+D2.
- **Deferred.** The unified {shots, overlap, stop} controller is not built until P5 and P6 show which actions pay.

### D-018 (2026-09-27): P5 design and protocol, declared BEFORE the dev run
- **Budget.** Every trial gets S0 shots per block up front. A further (r−1)·S0·B shots follow in batches of ΔS, one block per decision. The total r·S0·B is identical for every policy.
  - Uniform is round-robin and ends with exactly r·S0 shots per block.
  - Every block owns a pre-generated outcome stream (common random numbers), so policies are paired shot for shot.
- **Factors.**
  - S0 ∈ {2, 4, 8}, r ∈ {2, 4}, ΔS ∈ {1, S0}.
  - Partitions [4,4], [3,2,3], [2,2,2,2], [4,4,4], [3,3,3,3] and [4,4,4,4] (n = 8, 12, 16).
  - Decoder is a crossed factor: every trajectory is scored under AWQPE ε=0.9, AWQPE ε_safe and D2.
- **ΔS = 4·S0 is dropped.** With r ∈ {2, 4} the extra budget per block is only S0 or 3·S0, so ΔS = 4·S0 allows at most (r−1)·B/4 decisions, often zero or one. That is not a sequential policy.
- **Policies.** uniform, eig_cell (primary), entropy, ratio, legacy_cw, fisher (negative control) and random.
  - greedy_oracle[decoder] is a myopic one-step oracle that uses realized counterfactuals. It is an analysis reference, tagged is_oracle.
- **Fixed constants, not tuned.**
  - EIG support: posterior tail 1e-6, at most 8192 grid points. The dropped mass is checked on dev.
  - Likelihood grid refine: 3 for n = 8 and 12, 2 for n = 16.
  - Tie-breaking uses a seeded random choice.
- **What dev may decide.** Only:
  1. feasibility (cells whose runtime is prohibitive may be removed from the test grid, with the reason logged);
  2. test sample size (phases per stratum and replicates), from dev variance;
  3. confirmation that EIG support truncation is negligible.

  No policy, signal, threshold or ΔS is chosen from dev accuracy. Both ΔS values stay as factors in test.
- **Test (held-out split). The primary family is 3 comparisons: eig_cell vs uniform within each decoder.**
  - Endpoint: final P(|φ̂−φ| ≤ 2⁻ⁿ), with equal weight on every test cell.
  - Paired difference per phase, averaged over cells and replicates. The CI is a phase-cluster bootstrap (2000 resamples).
  - The p-value is a two-sided sign-flip permutation test on per-phase mean differences, with Holm correction across the 3 decoders.
  - TOST margin ±1 percentage point for "no difference".
- **Secondary** (reported, no multiplicity claims):
  - every other policy vs uniform;
  - regimes: stratum, n, partition, S0, r, ΔS and decoder-limited status;
  - efficiency (adaptive − uniform)/(oracle − uniform);
  - RMSE, median error and exact-n-bit rate;
  - U-queries and the shot distribution;
  - decision analysis: chose-best rate, regret, and the Spearman correlation of signal vs realized gain.
- **Decoder-limited trials** (the decoder fails even in the infinite-shot limit) are labelled `limit_correct = False` and analysed separately from shot-fixable trials.

### D-019 (2026-09-27): P5 configuration FROZEN for the held-out test
- **Test config.** `research/configs/p5_adaptive_test.yaml`
  - sha256 b4494f0b…331f7d
  - resolved config_hash 55b4a056e37c1ac1
- **Code at freeze.**
  - `allocation/policies.py` sha256 19de6cff…c04b
  - `run_adaptive_shots.py` sha256 2063f512…2a59f
  - `summarize_p5.py` sha256 7f5cfbb8…01cc
  - The commit that contains this entry is the freeze commit.
- **Identical to dev** (D-018 grid) except `split: test`, the experiment id, the status text, and removal of the pilot overrides. The test split draws phases from a disjoint seed domain and excludes the paper's example phases.
- **What dev decided (only what D-018 allows):**
  1. Feasibility: the whole grid ran in about 2 h on dev, so no cell is dropped.
  2. Sample size: kept at 8 phases per stratum × 5 replicates. Dev phase-cluster 95% CIs on the primary differences were about ±0.5 points, adequate to detect effects of ~1 point.
  3. EIG support truncation: fixed at tail 1e-6 with at most 8192 points, unchanged.
- **Not chosen from dev accuracy:** no policy, signal, ΔS or constant. Both ΔS values remain factors.
- **Analysis addition before test** (descriptive; not part of the primary family): the decision table gains `p_chose_best_when_informative`, the chose-best rate restricted to decisions where some block's next batch changes tolerance success or error. Without it the unconditional rate is dominated by ties.
- **Primary family unchanged from D-018:** eig_cell vs uniform within each decoder, with Holm across the 3 decoders.

### D-020 (2026-09-27): Phase 6 design; P5 artefacts frozen
- **P5 frozen.** P5 is frozen at commit 6fd3a1d with config hash 55b4a056e37c1ac1. P6 code, configs and results are new files: `run_adaptive_overlap.py`, `overlap/`, `decode/awqpe_overlap.py`, `configs/p6_*`, `results/p6_*`. P5 files are only imported, never modified.
- **Unchanged from P5:**
  - the same phase generator and dev/test split;
  - the same CRN outcome-stream strategy. Chunk streams are identical across every arm and variant of a shard; overlap blocks get their own streams.
  - the same three decoders, the same τ = 2⁻ⁿ primary metric, and the same phase-cluster bootstrap with Holm correction;
  - the same leakage firewall. The overlap triggers live in `overlap/`, which the AST test covers.
- **Methodological constraint: faithful AWQPE cannot use overlap blocks.** Algorithms 1–2 have no input for extra blocks.
  - **O-ext under D1:** handled by a separately labelled NON-PAPER decoder, `awqpe_ext`. When Algorithm 2 reaches boundary j and the already-corrected lower part is exactly 10…0 (the observable ε-floor situation, V1), chunk j is set to floor(t1_ext / 2^v); otherwise Algorithm 2 is unchanged.
    - Unconditional substitution was tried first and **rejected**: the widened reading's own rounding carries into chunk j about 2^−(v+1) of the time, which raised infinite-shot failures to 16–54%.
    - The conditional rule gives 0% infinite-shot failures on all partitions tested, and is bit-identical to faithful AWQPE when no widened block exists (tests).
  - **O-bridge:** has no D1 interpretation without new stitching logic, so it is evaluated with D2 only. Faithful D1 ignores the block; those rows are reported as "not usable by D1", not as bridge results.
- **Arms.**
  - B1 uniform; B2 P5 information gain (eig).
  - B3 overlap alone: A overlap batches chosen by the trigger, then uniform chunk shots.
  - B4 eig + overlap: the trigger chooses among chunks and at most A overlap batches.
  - B1u and B2u: uniform or eig on chunks with a per-trial U-query budget equal to what B3 or B4 consumed (equal quantum work; any overshoot favours the baseline).
  - Rescue continuations of B2 by A batches: greedy-oracle shots, greedy-oracle overlap, practical eig shots, practical trigger overlap.
- **Triggers.**
  - `eig`: the P5 one-shot cell information.
  - `lowerhalf`: overlap first at boundaries whose decoded lower part is 10…0, with eig as the fallback.
- **Cost accounting.** Every row carries total shots, chunk shots, overlap shots and U-queries.
- **Dev decides** (after the pilot, before the full dev and test runs): trigger, v, bridge shift, A_max and budgets. Nothing is tuned on test.

### D-021 (2026-09-27): P6 pilot configuration
- **Config.** `research/configs/p6_overlap_pilot.yaml` (dev split only).
- **Grid.**
  - Partitions: n=8 ([4,4], [3,2,3], [2,2,2,2]) and n=12 ([4,4,4], [3,3,3,3]).
  - S0 ∈ {2, 4}, r = 2, ΔS = S0.
  - Strata S2_boundary, S1_final_half and S4_uniform, 8 phases each × 4 replicates.
- **Variants (10).**
  - ext: v ∈ {1, 2}, trigger ∈ {eig, lowerhalf}, A ∈ {1, 2}.
  - bridge: shift ∈ {half, 1}, A ∈ {1, 2}.
- **Status.** The pilot is exploratory. Its results choose the variants for the full dev run; they are not claims.

### D-022 (2026-09-27): Corrected equal-U-query matching (supersedes the `uq` baselines of the P6 pilot)
- **Problem.** The pilot's B1u/B2u stopped at the first batch reaching the overlap arm's U cost. Because least-significant-block batches are expensive, they overshot by about 17% on average, which favoured the baselines. The pilot results stay as recorded and are labelled preliminary.
- **Procedure** (`baseline_path` / `match_to_target` in `run_adaptive_overlap.py`):
  1. For each shard, run the baseline chunk-only policy (uniform for B1, P5 eig for B2) on the same CRN streams. Use the same decision rule and tie-breaking seed as the nominal B1/B2 arms. Decode after every batch and record cumulative U-queries, until every trial's U reaches the largest target.
     - Decisions never depend on the target, so one path serves every comparison.
     - At the nominal budget the path reproduces B1/B2 exactly (tested).
  2. For each trial and overlap arm (target = that arm's realised U-queries), take the baseline state at the **last path point with U ≤ target (lower)** and at the **first point with U ≥ target (upper)**.
     - Both are recorded as rows, with `u_queries`, `u_target` and `u_diff_vs_target`.
     - The upper bracket is invalid (NaN, `match_valid = False`) if the path cap is reached; this is reported.
  3. **Primary equal-U comparison, conservative for the claim being made:**
     - a superiority claim for overlap is tested against the **upper** bracket (baseline with ≥ the same quantum work);
     - an inferiority claim is tested against the **lower** bracket;
     - both brackets and the U differences are always reported.
- **No truth** enters the path, the matching or the decisions.

### D-023 (2026-09-27): The `awqpe_ext` decoder extension, validated deterministically
- **Validation suite.** `research/tests/test_awqpe_ext_validation.py`, 22 cases:
  - (A) no overlap block, or blocks with no shots → bit-identical to faithful AWQPE;
  - (B) overlap present but trigger false → identical;
  - (C) trigger true → only bits at or above the least-significant triggered boundary can change;
  - (D) infinite shots → every recorded ε=0.9 failure phase (V1 failure records) is repaired, and every already-correct phase is untouched;
  - (E) no evaluator, simulator-truth or oracle import in the decoder or trigger code, and no phase argument in the decoder signature.
- **Status.** This validates the decoder extension, not P6 performance.
- **Terminology (binding).** `awqpe_ext` is a NEW EXPERIMENTAL DECODER EXTENSION.
  - Results using it are "adaptive overlap + extended decoder", never "overlap improves AWQPE".
  - O-bridge stays D2-only; no D1 stitching rule is invented.

### D-024 (2026-09-27): P6 pilot relabelled as a stress test
- **Relabel.** The P6 pilot (`p6_adaptive_overlap_pilot/20260927T152618Z`) used ambiguity-enriched strata only. It is a stress test and pilot; it is not used to select P6 variants.
- **Evidence.** The safe-ε diagnostic on the predeclared strata (P6D-1) shows that its positive `awqpe_ext` B4 − B2 does not carry over.

### D-025 (2026-09-27): Dropping O-ext v = 2 (mechanism and cost, not pilot accuracy); proposed, pending review
- **Theory.** v = 1 already removes the infinite-shot floor, with 0 failures on every partition (D-023 D). In the danger band, 2^v·δ rounds without carry for any v ≥ 1, so v > 1 adds no resolving power for the targeted situation.
- **Cost.** An O-ext block's U per shot is 2^k(2^(m+v) − 1), about 2× for v = 2. On dev (P6D-1), v = 2 added about 28% more U-queries per trial than v = 1 (2129 against 1665) with no larger immediate gain.
- **Approved by the user (2026-09-27).**

### D-026 (2026-09-27): P6 full DEV protocol and selection rule, declared BEFORE the dev run
- **Grid** (`research/configs/p6_overlap_dev.yaml`):
  - The P5 partitions [4,4], [3,2,3], [2,2,2,2], [4,4,4], [3,3,3,3] and [4,4,4,4] (n = 8, 12, 16).
  - All predeclared strata S0–S5, 8 phases per stratum (S5 = the paper's six phases), 5 replicates.
  - S0 ∈ {2, 4, 8}, r ∈ {2, 4}. The batch size is fixed at ΔS = S0.
- **Excluded: 1-shot batches (ΔS = 1) and S0 = 1.**
  - Each overlap action is one batch. With ΔS = 1 an overlap action would be a single shot of the overlap block, and the number of decisions per trial would grow by ×S0, multiplying compute.
  - S0 = 1 gives zero-information initial windows (P2c: all decoders equal at 1 shot per block).
  - This is a compute and design decision, not an accuracy-based one.
- **Arms.**
  - Uniform (B1), P5 information gain (B2), overlap-only (B3), information gain + overlap (B4).
  - Lower and upper U-matched B1/B2 baselines (D-022).
  - Rescue diagnostics: practical extra eig batch and practical trigger-overlap batch, plus greedy-oracle shot and overlap batches for awqpe, awqpe_eps_safe, awqpe_ext, awqpe_eps_safe_ext and likelihood.
- **Variants.**
  - O-ext v = 1 with triggers {eig, lowerhalf, lowerhalf_gated} × A ∈ {1, 2}.
  - O-bridge (D2 only) with shift ∈ {1, half} × A ∈ {1, 2}, eig trigger.
- **`lowerhalf_gated` (dev candidate, approved).** Overlap candidates are eligible only at boundaries whose current decoded lower part (faithful ε=0.9 chunk decode) is exactly 10…0.
  - In the overlap-only arm, a batch with no eligible boundary goes to the next uniform chunk, keeping shots equal.
  - In the rescue diagnostic it takes no action (no fallback), so "overlap rescue" is never contaminated by shots.
- **Performance-only changes (results bit-identical, tested):**
  - Chunk-only arms (B1, B2, extra-shot rescues, oracle-shot rescues, U-matching paths) are computed once per shard rather than once per variant. They are identical across variants: same chunk streams, seeds and arm keys. Verified against the frozen 9af6716 runner: 2,628/2,628 final rows and all action rows identical.
  - The P6 runner uses `allocation/eig_cached.py`, a memoised line-for-line copy of the frozen P5 signal. It is asserted bit-identical to `policies.eig_cell_support`, and the P5 file is untouched.
- **Analysis** (`summarize_p6.py`, tested):
  - Equal-U comparisons use only the D-022 brackets.
  - The **upper bracket is PRIMARY** for positive overlap claims; the lower bracket is a sensitivity analysis.
  - Invalid upper brackets are excluded and counted.
  - Legacy `*Umatched*` arms are refused.
- **Selection rule (fixed now, before any dev result):**
  - For each (mechanism, decoder) where the decoder can use the mechanism (ext: awqpe_ext, awqpe_eps_safe_ext, likelihood; bridge: likelihood), select the variant with the largest dev difference **B4 − B2m_upper** (equal U, conservative), pooled over all dev cells with equal cell weight via phase clusters.
  - Report its phase-cluster 95% CI and the Holm-adjusted p (across variants within that mechanism and decoder).
  - The lower bracket is reported as a sensitivity analysis. Lower mean U breaks ties only.
  - Faithful decoders are never selected (they cannot use overlap).
  - A selected variant whose CI includes 0 is still recorded, with the CI. Whether to take it to test is decided at user review.
- **Claims kept separate:**
  - (A) D2: overlap geometry adds information.
  - (B) O-ext + `awqpe_ext`: a new extended decoder can exploit O-ext.
  - (C) Resource allocation: overlap versus extra shots at equal U-queries.


### D-027 (2026-09-28): P6 variants FROZEN for the held-out test (selection by the D-026 rule on dev)
- **Dev run.** `p6_adaptive_overlap_dev/20260927T174305Z`: 216/216 shards, 0 failed.
- **Selected** (primary = B4 - B2m_upper at equal U, phase-cluster CI, Holm across variants within mechanism x decoder):
  - O-bridge, D2 likelihood: **bridge_s1_A2**, +0.75 pts [0.50, 1.00], Holm p < 1e-4 (runner-up bridge_half_A2, +0.59).
  - O-ext, D2 likelihood: **ext_v1_eig_A2**, +0.51 [0.31, 0.70], Holm p < 1e-4 (runner-up ext_v1_lowerhalf_A1, +0.36).
  - O-ext + awqpe_ext (eps 0.9, extended decoder): **ext_v1_gated_A2**, +6.86 [5.12, 8.74], Holm p < 1e-4 (runner-up ext_v1_gated_A1, +6.85; A2 also has the lower U).
  - O-ext + awqpe_eps_safe_ext (extended decoder): **ext_v1_gated_A1**, +2.73 [1.81, 3.72], Holm p < 1e-4 (runner-up ext_v1_gated_A2, +2.72).
- **Test config.** `research/configs/p6_overlap_test.yaml`: the D-026 grid on the test split, restricted to these four variants. S5 (paper phases) never appears in test. config_hash a5ffa3c9e917e95f, sha256 7f7ed4d6b1f2a292.
- **Frozen selection file.** `research/configs/p6_frozen_selection.csv`, sha256 7a5e30b2e1b57b20.
- **Test primary family.** The 4 (mechanism, decoder, variant) triples above, comparison B4 vs B2m_upper (equal U, conservative), Holm across all 4. `summarize_p6.py --tag test --frozen-selection` refuses to re-select.
- **Secondary.** Equal-shot and lower-bracket comparisons, overlap-only vs uniform, the rescue table and strata. No multiplicity claims.
- **Test is run exactly once.** No parameter may change after test data exist.


### D-028 (2026-09-28): Phase 8 phase-tolerance design, declared before any P8 data
- **Question (Mode B).** For a declared tolerance tau and alpha = 0.05, what resources give realised coverage P(|phi_hat - phi|_circle <= tau) >= 1 - alpha? Do posterior-credible stopping and information-guided allocation reduce them relative to a fair fixed-budget baseline?
- **Tolerance grid.** tau = 2^-j, j in {3, 4, 6, 8, 10, 12, 14, 16}, with j <= n.
- **tau-matched truncation.** tau_bits(tau) = ceil(log2 1/tau) - 1, because rounding to K bits has error <= 2^-(K+1). The "_trunc" arms use the shortest MSB-first block prefix with at least tau_bits(tau) bits. Every arm is run on both the full partition and the prefix, so the saving from stopping is kept separate from the trivial saving of measuring fewer (and U-expensive) low-order bits.
- **Arms.** All arms share per-block common-random-number streams.
  - fixed_S (S in {4, 8, 16, 32, 64, 128} per block), decoded by D2 (grid MAP) and by D1 = faithful AWQPE at eps_safe. A one-block prefix uses eps 0.9, where eps is inert.
  - stop_uniform: S0 = 4 per block, then round-robin batches of 4.
  - stop_eig: S0 = 4, then batches of 4 to the block maximising the frozen P5 information-gain signal. The cell resolution is set to tau_bits(tau) rather than n. This is the only change, and it is declared here.
  - The cap is 128 shots per block. Trials that reach the cap without stopping are kept (stopped = False) and scored.
- **Stopping rule** (`stopping/rules.py`, inside the leakage firewall). Stop when the uniform-prior grid posterior mass within tau of the MAP is >= 1 - alpha. alpha is fixed at 0.05 and is not tuned.
- **Dev calibration (the only tuned quantity).** For each (partition, tau, full|trunc, decoder), S* = the smallest fixed S whose dev coverage (pooled over strata S0-S4, equal weight) is >= 0.95. This fixed-S* arm is the "tau-matched fixed" baseline. If no S in the grid reaches coverage, the cell has no baseline and is reported as such.
- **Test primary family** (Holm, predeclared). Partitions [4,4], [4,4,4], [4,4,4,4]; tau in {2^-3, 2^-(n/2), 2^-n}. That gives 9 cells, each comparing stop_eig_trunc (or stop_eig where the prefix is the full partition) against fixed_S*_trunc with D2.
  - Endpoint: the paired difference in mean total shots per trial, with a phase-cluster bootstrap CI and a sign-flip permutation p.
  - A cell counts as a positive result only if (i) the difference favours stopping after Holm and (ii) stop_eig realised test coverage is >= 0.95 at the point estimate. Coverage is also reported with its cluster CI.
  - U-queries are co-reported for every comparison.
- **Secondary (no multiplicity claims):**
  - stop_eig vs stop_uniform (allocation effect);
  - full vs trunc (truncation effect);
  - D1 fixed S* vs D2 fixed S*;
  - coverage per stratum, including under-coverage in S1/S2;
  - calibration of credible mass vs realised coverage.
- **Mode A** (tolerance success vs budget) is a re-analysis of the frozen P5 test run at several tau. No new data; it is labelled as a re-analysis.
- **Pipeline.** Pilot, then the dev grid (`research/configs/p8_tolerance_dev.yaml`), then a freeze record with the S* table and test config hash, then one test run.


### D-029 (2026-09-28): P8 FROZEN for the held-out test
- **Dev run.** `p8_tolerance_dev/20260927T220813Z`: 25/25 shards, 0 failed, 40 phase clusters x 4 replicates per partition.
- **Frozen S* table.** `research/configs/p8_frozen_sstar.csv` (sha256 79838c7bc0465356): the smallest fixed S per block with pooled dev coverage >= 0.95, per (partition, tau, full|trunc, decoder). Nothing else was tuned. alpha, S0, cap, the grids and the stopping rule are unchanged from D-028.
- **Test config.** `research/configs/p8_tolerance_test.yaml`: the D-028 design on the test split, 12 phases per stratum (S0-S4), 4 replicates. config_hash 248abb85674d7dd9.
- **Primary family.** 9 cells, Holm, as declared in D-028.
  - Where dev already shows the stopping arm spending more shots than the fixed S* baseline (coarse tau, where S* = S0 = 4), the cell stays in the family. It was declared before the data; it is not dropped.
  - `summarize_p8 --tag test` refuses to run without the frozen S* table.
- **Test is run exactly once.**
