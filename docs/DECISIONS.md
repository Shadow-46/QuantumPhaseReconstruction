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
