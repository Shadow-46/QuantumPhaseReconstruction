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

## V1: Independent verification of the ε=0.9 infinite-shot floor (2026-09-27), PRELIMINARY (arXiv-v3 reading)
- **Status.** This is a **potential decoder limitation under our arXiv-v3 implementation**, not a claim about the published algorithm. Earlier entries (P2a) used stronger wording ("refutes"). This entry supersedes that wording; the earlier entries are left unedited per the append-only rule.
- **Run.** `verify_epsilon_floor/20260926T210310Z` (config `research/configs/verify_epsilon.yaml`). 11 shards failed with out-of-memory errors at m=6. They were re-run via `--resume --rerun-failed` after making the kernel memory-linear (identical values; tests unchanged).
- **Tables.** `research/analysis/tables/verify_*.csv`, including every failure record in `verify_failure_records.csv.gz`. Each record carries: phi, partition, ε, affected block, the block's full probability vector, t1/t2 and their probabilities, C2/C1, flag, raw/decoded/true bits, true nearest and reconstructed phase, and circular/absolute error.
- **Results.**
  1. **Independent implementation.** A string-based Algorithm 1–2 decoder written from the pseudocode (`verification/awqpe_strings.py`) agrees with the primary decoder on **262,144/262,144** infinite-shot decodes: 4 partitions × 2 grids × 8 ε.
  2. **Dense ε sweep.** ε from 0.30 to 0.99 in steps of 0.01, over 13 partitions with n = 8–12 and block widths 2–6. On the generic grid (offset 0.37), the exact per-phase theoretical predicate (`verification/theory.py`) gives **0 false positives and 0 false negatives** over 330,904 observed failures (faithful decoder) and 340,264 (special-chunk ablation). On the offset-0 grid, the only mispredictions are exact (n+1)-bit ties, which are measure-zero and excluded from the predicate by construction.
  3. **Threshold.** The empirical failure onset lies in (ε*, ε*+0.01] for every partition, where ε*_min = min over boundaries of ε*(k,m). Checked against the ε* computed from the Dirichlet-kernel probabilities, e.g. [4,4]: ε* = 0.7792, onset 0.79; [6,6]: ε* = 0.9394, onset 0.95. [6,6] and [2,4,6] show 0% failures at ε=0.9, as predicted.
  4. **All plausible readings.** All 216 combinations of the ambiguity readings were tested: flag index, min-mod (3 forms), borrow source, final block, special chunk (3 forms) and ties (3 forms).
     - At ε=0.9 **every reading keeps a failure floor**: per partition, the best reading for [4,4] still fails 1.56% and for [3,2,3] 4.9%.
     - **108/216 readings reproduce all six §6 examples: exactly those that read the borrow MSB from the already-corrected lower chunk (our default).** Every "raw" reading breaks at least one example.
  5. **Failure anatomy at ε=0.9.** In every failure the responsible upper block is unflagged. Errors span 4–455 LSB (median 4–32 LSB by partition).
  6. The paper's worked examples still reproduce (tests `test_verification.py`, `test_awqpe_baseline.py`).
- **Interpretation.** The failure floor is a property of Algorithms 1–2 as written in v3 under every reading we could construct. It is fully explained by one mechanism: the excluded case δ̂_k = 0.5 of Lemma 3.2, left unflagged when ε ≥ ε*(k, m).
- **Next.** Confirm against the published version. The faithful ε=0.9 baseline is kept unchanged. The ε_safe baseline (D-011) is added as a separately labelled condition.

## P3: Initial-state / eigenphase sensitivity (2026-09-27), PRELIMINARY (dev phases, n=8)
- **Run.** `p3_state_sensitivity/20260926T211750Z` (config `p3_state_sensitivity.yaml`). Tables: `research/analysis/tables/p3_*.csv`.
- **Objective.** Establish why the Paper A |1⟩/U_a experiments were not a clean single-eigenphase AWQPE experiment. Separate state preparation, eigenphase, grid location and resulting ambiguity.
- **Configuration.**
  - 142 cases in six families:
    - A: phase gate, 53 stratified phases plus 1/4, 1/6, 1/10.
    - B: 3-qubit diagonal U, 8 basis eigenstates.
    - C: V·D·V† with a Haar-random V, 8 eigenvectors.
    - D: U_a eigenstates |u_s⟩ for (15,2), (21,2), (33,2), (15,4), (21,19), (15,13) and (21,8).
    - E: |1⟩ on the same U_a instances (the Paper A input).
    - F: controlled mixtures, w0 ∈ {1, .9, .75, .5}.
  - Partitions [4,4] and [3,2,3]; shots per block {16, 256, 2048}; 200 replicates.
  - Decoders: AWQPE at ε=0.9, AWQPE at ε_safe (D-011), and D2 (which assumes a single eigenphase, so it is misspecified for mixtures).
- **Verification.**
  - All 111 claimed eigenstates have ‖Uψ − e^{2πiφ}ψ‖ ≤ 2.7e-15.
  - Spectral decomposition of U_a|1⟩ gives weights exactly 1/r on s/r.
  - For every case, real Qiskit circuits (StatePreparation of ψ, controlled U^(2^(k+p)), IQFT) match the mixture-of-kernels prediction to ≤ 2.1e-8.
- **Results.**
  1. **State-preparation route is irrelevant for exact eigenstates.** B (diagonal) and C (V·D·V†) with the same eigenphases differ by ≤ 0.02 in success, which is sampling noise.
  2. **For eigenstates, grid location drives ambiguity and error.**

     | Grid location | mean true top-two ratio | AWQPE ε=0.9 @2048 shots | AWQPE ε_safe @2048 | D2 @2048 |
     |---|---|---|---|---|
     | generic | 0.22 | 1.00 | 1.00 | 1.00 |
     | dyadic_n | 0.28 | 1.00 | 1.00 | 1.00 |
     | final_half | 0.81 | 0.89 | 0.98 | 1.00 |
     | boundary_hard | 0.87 | 0.87 | 0.99 | 1.00 |

     At 16 shots per block, boundary_hard gives 0.61 / 0.75 / 0.99 for the three decoders.
  3. **U_a eigenstates |u_s⟩ work for every instance, dyadic or not.** With 2048 shots, P(|err| ≤ 2⁻ⁿ) = 1.00 with zero ambiguity flags for r = 2, 4, 6 and 10.
  4. **The Paper A input |1⟩ (uniform mixture over s/r) behaves very differently.**
     - Every block's true top-two ratio is exactly 1.0, so AWQPE flags 84–100% of trials as ambiguous, regardless of the phase grid.
     - AWQPE's min(t1, t2) rule then collapses most outputs to 0, the trivial s = 0 component: P(output = 0) = 0.44–1.00.
     - The probability of returning a *useful* nonzero s/r:

       | r | phases | AWQPE | D2 |
       |---|---|---|---|
       | 2 | dyadic | 0.00–0.02 | 0.44–0.54 |
       | 4 | dyadic | 0.47–0.56 | 0.73–0.77 |
       | 6 | non-dyadic | 0.00–0.10 | ≈ 0 |
       | 10 | non-dyadic | 0.04–0.18 | 0–0.25 |

     - For non-dyadic r, independent blocks collapse onto different eigencomponents, so the decoded chunks do not belong to any single phase (cross-block incoherence).
  5. **Mixtures with a dominant component (w0 ≥ 0.75) are nearly harmless.** Success is 1.00 at 2048 shots; degradation begins at w0 = 0.5.
- **Interpretation (answers the P3 question).**
  - The old setup's behaviour was caused **primarily by the non-eigenstate input**: a uniform mixture over s/r creates exact top-two ties in every block, irrespective of phase-grid location.
  - **Phase-grid alignment modulates the damage**: for dyadic s/r (r = 2, 4) later blocks still agree across components; for non-dyadic r they do not.
  - So the answer is **both, interacting**, with the mixture as the dominant cause.
  - Consequence: Paper A's top-two / C_w evidence mostly measured mixture ties, not phase ambiguity. This supports Alok Shukla's point: with a proper eigenstate, the top-two structure reflects grid position and becomes informative.
- **Supports or refutes.** Supports the hypothesis that the initial state, not only the decoder, invalidated the old experiments as AWQPE evidence.
- **Next.**
  - P4 information study on eigenstate inputs, with both D1 conditions and D2.
  - Mixtures return in P9 as a robustness factor.
