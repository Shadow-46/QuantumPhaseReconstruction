# Results log

Append-only. Each entry names its run directory under `research/results/`. Large shards there are git-ignored; summaries are committed. "Preliminary" means dev-split and/or pilot scale, not a final claim.

---

## P1: Dirichlet-kernel and Fisher validation (2026-09-27)
- **Objective.** Validate the measurement model before any experiment.
- **Configuration.** Unit tests `research/tests/test_kernel.py`: m = 1–8, offsets 0–7; statevector cross-check for m ∈ {2,3,4} and k ∈ {0,2,5}.
- **Results.**
  - **Model accuracy.**
    - Normalisation, symmetry, grid-exactness and Lemma 3.3 bounds all hold.
    - The halfway top-two probabilities are equal and sum to ≥ 8/π².
    - Kernel versus Qiskit statevector: max |Δ| = 7e-15.
    - Fejér versus sin² ratio: 7e-16.
  - **Fisher information agrees.** Analytic Fisher matches finite differences to 3e-9 relative.
  - **Finding F-1: ideal block Fisher information is phase-independent.** For the ideal kernel the per-shot classical Fisher information of a block equals the QFI of its control register for every φ, including grid points via the double-zero limit:

    I_b(φ) = 4^k · (4π²/3)(M²−1)

    Numerically CFI/QFI = 1 ± 5e-10 over 20,001 phases for m = 2, 3, 4, 6.
  - **Finding F-2: under noise, Fisher information depends strongly on φ.**
    - It is lowest at grid-aligned phases. With depolarisation or readout error it is exactly 0 there, because all derivatives vanish.
    - It is highest near half-way points.
    - With jitter σ = 0.02 and m = 4, CFI/QFI ranges 0.20–0.33.
- **Interpretation.**
  - A "Fisher-driven" allocation, evaluated at the true or estimated phase, cannot distinguish ambiguous from unambiguous blocks in the noiseless case. It reduces to a fixed schedule set by (k, m).
  - Under noise, Fisher information is largest exactly where AWQPE's discrete decision is hardest, and zero where it is easiest. Fisher information measures local sensitivity, not bit-decision reliability.
  - These are hypotheses for P4 to test empirically, not conclusions about allocation value.
- **Supports or refutes.** Refutes the prior "Fisher information will flag ambiguous windows" in the ideal model.
- **Next.** P4 compares Fisher information with ambiguity and posterior signals as predictors of the marginal value of shots.

## P2a: Faithful AWQPE in the infinite-shot limit (2026-09-27)
- **Run.** `p2a_infinite_shot_limit/20260926T204034Z` (config `research/configs/p2_infinite_shot_limit.yaml`).
  - The earlier run `20260926T203851Z`, which has no ablation column, is superseded but kept.
- **Objective.** Test Theorem 3.7 with no sampling noise: counts equal the exact kernel probabilities.
- **Configuration.**
  - Eight partitions: [4,4], [2,2,2,2], [3,2,3], [3,3,2], [2,3,3], [4,4,4], [3,3,3,3], [2×6].
  - Dense grid 2^(n+6) at offset 0 (contains dyadic and exact half-way points) and at offset 0.37.
  - ε ∈ {0.5, 0.7, 0.8, 0.9, 0.95}.
- **Results.** P(estimate ≠ best n-bit), where the final-block width is the width of the last partition entry:

  | final block width | ε = 0.5 | ε = 0.9 | ε = 0.95 |
  |---|---|---|---|
  | 2 | 3.4% | 12.8% | 14.4% |
  | 3 | 0% | 4.9–5.1% | 5.7–6.4% |
  | 4 | 0% | 1.8% | 2.5% |

  - **Mechanism, verified.** Every failure has some boundary whose lower bits round, or tie, to exactly 10…0 while the upper block is unflagged, i.e. ε ≥ ε*(k, m).
  - **Thresholds.** ε*(k) ≈ 0.36, 0.61, 0.78, 0.88 and 0.94 for k = 2…6, nearly independent of m. They predict exactly where failures start. For example, [4,4] shows 0% at ε = 0.7 and failures at ε = 0.8, against ε* = 0.779.
  - **The only other failures** are measure-zero exact (n+1)-bit ties on the offset-0 grid.
  - **Ablation.** Disabling the special-chunk rule leaves the failure rates essentially unchanged; it only moves which side of the tie fails.
  - **Failure size.** Median failure error ≈ 7.7 × 2⁻ⁿ. These are chunk-level (catastrophic) errors, not last-bit errors.
- **Interpretation.** Under our reading of v3, AWQPE with the suggested ε ≈ 0.9 has an irreducible failure floor unless each boundary has at least 6 lower bits. ε is a correctness-critical parameter, not only a false-positive trade-off.
- **Supports or refutes.** Refutes "failures are exceptionally rare" for practical small blocks at ε = 0.9. It supports Remark 3.4's claim that low ε is safe, and shows low ε is *required*.
- **Caveat.** Based on arXiv v3; confirm against the published version.
- **Next.** Include ε as a factor everywhere. Consider a "geometry-safe ε" baseline, ε < min_boundaries ε*(k), as a stronger AWQPE-faithful reference.

## P2b: Finite-shot AWQPE baseline Monte Carlo, preliminary (2026-09-27)
- **Run.** `p2b_awqpe_baseline_mc/20260926T204053Z` (config `p2_baseline_mc.yaml`, dev split).
- **Configuration.**
  - Partitions [4,4], [2,2,2,2], [3,2,3], [3,3,2].
  - Shots per block ∈ {16, 32, 64, 128, 256, 1024, 10240}.
  - 206 dev phases (40 per stratum S0–S4, plus 6 paper phases), with 100 replicates each: 20,600 trials per cell.
  - Counts are shared across ε, decoders and the ablation.
- **Results.** P(|err| ≤ 2⁻ⁿ), excluding S1:

  | Configuration | AWQPE, ε=0.9 | AWQPE, ε=0.5 | D2 |
  |---|---|---|---|
  | [4,4], 16 shots/block | 86.6% | 96.2% | 99.6% |
  | [4,4], 256 shots/block | 93.6% | 100.0% | 100.0% |
  | [4,4], 10240 shots/block | 94.5% | 100.0% | 100.0% |
  | [2,2,2,2], 16 shots/block | 71.9% | 90.0% | 98.8% |
  | [2,2,2,2], 10240 shots/block | 88.8% | 99.1% | 100.0% |

  - Most of the ε=0.9 loss is in S2 (boundary stratum): 73% at [3,2,3] with 256 shots, against 99.9% at ε=0.5.
  - **Leakage checks on D2.**
    - Shuffling counts across trials drops D2 to chance (0.88% against 0.81% expected).
    - Shuffling one block drops it to 29%.
    - An independent brute-force MLE agrees on 100% of trials.
- **Interpretation.**
  - At these scales most of AWQPE's error is **decoder-induced**, not shot-limited: more shots cannot remove the ε-floor.
  - A joint-likelihood decoder on the same counts nearly eliminates it with 16 shots per block.
  - For the adaptive programme this means:
    1. Adaptive-allocation claims must be made within a decoder (D-005).
    2. Against faithful AWQPE at ε = 0.9, any method that incidentally changes decoding would look spuriously good.
    3. The headroom left for adaptive *shot* allocation under D2 is small in the noiseless eigenstate setting at n = 8. Adaptive allocation may matter mainly at lower budgets, larger n, or under noise and mixtures. That is to be measured, not assumed.
- **Preliminary.** Dev split, n = 8 only, ideal model. Confidence intervals are pending the analysis script.
- **Next.**
  - P3: initial-state and mixture sensitivity.
  - P4: information signals and marginal-value-of-shots study.
  - Lower-budget regime (1–16 shots per block).
  - n = 12 and 16.
