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

## P4: Do window signals predict the marginal value of extra shots? (2026-09-27), PRELIMINARY (dev phases, ideal model)
- **Run.** `p4_information_study/20260926T212553Z` (config `p4_information.yaml`; analysis plan D-016 declared before inspection). Tables: `research/analysis/tables/p4_*.csv`.
- **Run incidents, no data lost.**
  - The first attempt lost all 180 shards when an out-of-memory worker death broke the process pool. It was resumed with `--rerun-failed` after two fixes: smaller posterior chunks, and the runner now rebuilds broken pools.
  - A final Windows file-lock on `progress.json` was repaired by a no-op resume. All 180 shards are complete.
- **Design.**
  - Partitions [4,4], [3,2,3], [2,2,2,2] (n=8) and [4,4,4], [3,3,3,3] (n=12).
  - S0 ∈ {4, 16, 64} shots per block; ΔS ∈ {S0, 4·S0}.
  - 206 dev phases × 30 replicates. Each block receives its own pre-drawn extra batch (CRN).
  - Improvement is measured within each decoder: AWQPE ε=0.9, AWQPE ε_safe, and D2.
- **Headroom.** One batch on the best block (oracle) versus a random block.
  - At S0=4, the oracle minus random is 5–14 points in every decoder, largest in the S2 boundary stratum (13.8–16.8 points).
  - At S0 ≥ 16, D2 is already ≥ 99% and has no headroom.
  - **Only 29–49% of AWQPE@0.9 failures are fixable by one extra batch on any block.** The rest is the ε floor (V1). The figures are 47–79% for ε_safe and 72–100% for D2.
- **Which signal predicts benefit?** Greedy one-step policy: put ΔS on the block that maximises the signal.

  | Decoder | Signal | Gain vs random block at S0=4, ΔS=4 (pts) |
  |---|---|---|
  | D2 | eig_cell | +3.1 to +7.2 (95% phase-cluster CIs exclude 0) |
  | D2 | count signals (entropy, c1, margin, ratio, chunk/boundary risk, legacy C_w) | +1.9 to +5.0 |
  | D1 ε_safe | eig_cell | +3.0 to +9.1 |
  | D1 ε_safe | entropy | +2.1 to +7.5 |
  | D1 ε=0.9 | eig_cell | +3.0 to +6.6 |
  | any | Fisher information | −0.4 to −3.2 (CIs exclude 0 for D2) |

  - With ΔS = 16 the gains grow, e.g. ε_safe on [2,2,2,2]: eig_cell +15.4, entropy +12.9.
  - Mean share of oracle headroom captured across all dS=S0 cells:
    - D2 + eig_cell: 78%. D2 + eig_phi: 63%.
    - Best count signal: 49–55% under D2, 26–35% under D1.
    - Fisher: −18% to −46%, i.e. worse than random.
  - Within-failing-trial AUROC for "this block's batch fixes the trial":
    - D2: eig_cell 0.90, eig_phi 0.83.
    - D1: count signals 0.79–0.86.
    - Fisher: 0.25–0.38 (anti-predictive).
- **Mechanism for Fisher.** In the ideal model Fisher is 4^k·const (F-1), so it always chooses the least-significant block. That block has the lowest fix rate (1.1–1.9% at S0=4) against 4–9% for upper blocks.
- **Legacy C_w.** It performs like the plain top-two ratio or margin; it adds nothing beyond simple count statistics.
- **Calibration.** Block-local P(chunk correct) is roughly calibrated. It is overconfident by 5–6 points in the 0.65–0.87 range and underconfident below 0.3.
- **Interpretation (answers G.1 and G.2 for this regime).**
  1. Window evidence *does* predict where an extra batch helps, well above chance.
  2. The most useful signal is a decision-focused expected information gain: the mutual information of an extra shot with the n-bit cell, under the global posterior. It is the best signal for every decoder.
  3. Top-two and count statistics carry about 60% of that value and are interchangeable.
  4. **Fisher information is not a useful allocation signal here; it is harmful.**
  5. Benefits live in the low-shot regime (S0 ≈ 4) and at boundary-hard phases.
  6. For faithful AWQPE@0.9, most failures cannot be fixed by shots at all, so there overlap or a decoder change is required (candidate outcome 4).
- **Caveats.** One-step, single-batch counterfactuals only; ideal model; dev split. Multi-step allocation under a fixed budget (P5) must show whether these one-step gains accumulate.

## P2c: Uniform-shot baseline in the low-shot regime, n = 8/12/16 (2026-09-27), PRELIMINARY (dev split)
- **Run.** `p2c_lowshot_scaling/20260926T220629Z` (config `p2c_lowshot_scaling.yaml`). Table: `research/analysis/tables/p2c_lowshot_scaling.csv` (phase-cluster bootstrap CIs).
- **Run incidents.** The [8,8] partition (m=8, 256 outcomes) exceeded memory in the dense likelihood table. A sparse observed-outcome likelihood path was added, equal to the dense one to 2e-11 on finite entries with 100% MAP agreement.
  - Orphaned workers from the first attempt kept running and starved the re-run. They were stopped, and the runner now takes a per-run lock and clears stale failure files.
  - The [8,8] rows at 8 and 16 shots per block come from the re-run.
- **Design.**
  - Partitions [4,4], [3,2,3], [2,2,2,2], [4,4,4], [3,3,3,3], [6,6], [4,4,4,4] and [8,8].
  - Shots per block ∈ {1, 2, 4, 8, 16}; 206 dev phases × 50 replicates.
  - Decoders on identical counts: AWQPE ε=0.9, AWQPE ε_safe, and D2.
- **Results.** P(|err| ≤ 2⁻ⁿ), all strata:
  - **At 1 shot per block** every decoder is equal (39–66%). There is no information to exploit, and success is set by the number of blocks.
  - **With [4,4,4,4] at n=16:**

    | shots/block | AWQPE ε=0.9 | AWQPE ε_safe | D2 |
    |---|---|---|---|
    | 4 | 68.4% | 68.4% | 79.7% |
    | 16 | 81.3% | 86.9% | 99.4% |

  - **AWQPE@0.9 saturates** at 72–90% by 16 shots per block, well below D2.
  - **ε_safe is not uniformly better.** It helps at ≥ 8 shots but can be *worse* at 4 shots ([2,2,2,2]: 60.1% against 63.6%): with few shots, low ε raises more spurious flags. So the finite-shot ε trade-off is real.
  - **Success depends mainly on the block structure, not on n.** [2,2,2,2] (n=8), [3,3,3,3] (n=12) and [4,4,4,4] (n=16) behave alike at equal shots per block. Fewer, wider blocks do better at equal shots per block, though they cost more U-queries and deeper circuits.
  - **Boundary-hard phases (S2) at 16 shots:** AWQPE@0.9 reaches 59–66%, ε_safe 66–92%, D2 98–99%.
- **Interpretation.** The low-shot regime (2–8 shots per block) is where allocation can matter for every decoder: D2 still has 5–20 points of headroom there, and D1 has more. This regime, at n = 8–16, is the target for P5.

## P5: Multi-step adaptive shot allocation, DEV results (2026-09-27), PRELIMINARY / DEV
- **Hypothesis.** Under a fixed total shot budget, an information-guided sequential allocation raises P(|φ̂−φ| ≤ 2⁻ⁿ) relative to uniform, within each decoder. The effect should be concentrated in shot-fixable ambiguity regimes. Fisher information (negative control) should not help.
- **Configuration.**
  - D-018 grid on the dev split: 6 partitions (n = 8, 12, 16); S0 ∈ {2, 4, 8}; r ∈ {2, 4}; ΔS ∈ {1, S0}.
  - 46 phases × 5 replicates per partition-cell.
  - 7 policies plus 3 decoder-specific greedy oracles.
  - Common-random-number streams. Every trajectory is scored under AWQPE ε=0.9, AWQPE ε_safe and D2.
- **Run.** `p5_adaptive_shots_dev/20260927T0702*`. The analysis used 416/432 shards. The last 16 (n=16, S0=8) finished afterwards; the test analysis uses the full grid.
- **Results.** Paired difference vs uniform in P(tol), phase-cluster bootstrap 95% CI (276 clusters), in percentage points:

  | policy | AWQPE ε=0.9 | AWQPE ε_safe | D2 |
  |---|---|---|---|
  | eig_cell (primary) | **+2.75** [2.26, 3.25] | **+4.51** [3.96, 5.13] | **+2.85** [2.48, 3.24] |
  | legacy_cw | +2.14 [1.55, 2.67] | +2.46 [1.84, 3.06] | +0.45 [0.06, 0.94] |
  | ratio | +1.44 | +1.17 | −0.87 |
  | entropy | +0.60 | +0.21 (ns) | −1.68 |
  | random | −1.41 | −1.16 | −1.33 |
  | fisher (negative control) | **−8.26** | **−12.65** | **−13.86** |
  | greedy oracle (reference) | +9.49 | +7.07 | +2.62 |

- **Efficiency vs the greedy oracle** (eig_cell): 0.29 for ε=0.9, 0.64 for ε_safe, and 1.09 for D2. Above 1 is possible because the oracle is myopic.
- **Regimes** (eig_cell − uniform, points; order: ε=0.9 / ε_safe / D2):
  - **Decoder-limited trials** (AWQPE@0.9 fails even with infinite shots): +0.15. Shot allocation cannot fix the decoder floor.
  - **By S0:** S0=2 gives +4.4 / +7.8 / +5.9; S0=8 gives +1.2 / +1.2 / +0.3.
  - **By stratum:** S2_boundary gives +4.0 / +7.0 / +4.3, the largest of all strata.
  - **By n:** n=16 gives +4.4 / +6.2 / +4.2, versus n=8 +2.4 / +4.5 / +2.5.
  - **By budget:** r=2 gives +3.7 / +6.0 / +4.0; r=4 gives +1.7 / +2.8 / +1.5.
- **Decision analysis.**
  - eig_cell has the lowest error regret of the practical policies: 3.2 LSB for D2, against 9.7 for uniform and 58.9 for Fisher.
  - It has the largest error reduction per extra shot: about 9.0–9.4 LSB per shot, against 7.5–8.2 for uniform and ≈0 for Fisher.
  - Raw signal-vs-realized Spearman is small (0.04–0.07) because most (decision, block) pairs have zero realized gain. The conditional chose-best metric was added for test (D-019).
- **Interpretation (dev only).** The information-gain policy beats uniform for every decoder. The advantage concentrates where theory and P4 predicted: low per-block shots, boundary-hard phases, larger n, and tight budgets. It vanishes for decoder-limited failures. Fisher-driven allocation is strongly harmful. Count-based signals are weak and can hurt under D2.
- **Next.** Held-out test with the frozen configuration (D-019).

## P5: Multi-step adaptive shot allocation, HELD-OUT TEST (2026-09-27)
- **Configuration.** Frozen under D-019 (`research/configs/p5_adaptive_test.yaml`, config_hash 55b4a056e37c1ac1, freeze commit 1c29d1d). The grid is the same as dev: 6 partitions (n = 8, 12, 16), S0 ∈ {2, 4, 8}, r ∈ {2, 4}, ΔS ∈ {1, S0}.
  - 40 held-out phases per partition (disjoint seed domain, no paper phases) × 5 replicates: **240 phase clusters**.
  - 7 policies plus 3 greedy oracles, with common-random-number streams. Each trajectory is scored under all 3 decoders.
- **Run.** `p5_adaptive_shots_test/20260927T084715Z`: 360/360 shards, 0 failed, about 2.5 h. No code or parameter changed after the freeze.
- **Tables.** `research/analysis/tables/p5_test_*.csv`.
- **Primary family** (eig_cell vs uniform, P(|φ̂−φ| ≤ 2⁻ⁿ), paired per phase, 2000-resample phase-cluster bootstrap, 20,000-draw sign-flip test, Holm across 3 decoders):

  | decoder | diff (pts) | 95% CI | p (Holm) | greedy oracle − uniform | efficiency vs oracle |
  |---|---|---|---|---|---|
  | AWQPE ε=0.9 | **+3.06** | [2.47, 3.67] | < 1e-4 | +9.98 | 0.31 |
  | AWQPE ε_safe | **+4.91** | [4.22, 5.60] | < 1e-4 | +7.42 | 0.66 |
  | D2 likelihood | **+2.85** | [2.46, 3.26] | < 1e-4 | +2.67 | 1.06 |

  All three primary hypotheses are supported. No p-value reached the resolution of the permutation test (0/20,000). The TOST ±1 pt equivalence does not hold for any primary comparison.
- **Secondary** (vs uniform, points; order: ε=0.9 / ε_safe / D2):
  - legacy_cw: +3.06 / +3.37 / +0.81.
  - ratio: +2.14 / +1.99 / −0.44 (ns under D2).
  - entropy: +1.23 / +1.10 / −1.22.
  - random: −0.82 / −0.56 / −1.03.
  - **Fisher (negative control): −8.78 / −13.67 / −13.78.**
- **Regimes** (eig_cell − uniform, points; order: ε=0.9 / ε_safe / D2):

  | Regime | Levels |
  |---|---|
  | S0 (shots per block up front) | 2: **+5.0 / +8.6 / +6.1**; 4: +3.1 / +4.8 / +2.3; 8: +1.0 / +1.4 / +0.1 |
  | n | 8: +2.4 / +4.2 / +2.6; 12: +2.9 / +4.9 / +2.8; 16: **+5.2 / +7.3 / +3.7** |
  | r (budget multiplier) | 2: +4.0 / +6.4 / +4.3; 4: +2.1 / +3.5 / +1.4 |
  | stratum | S2_boundary: +3.8 / **+8.0** / +4.2; S1_final_half is the smallest: +1.9 / +2.2 / +1.8 |

  - **Decoder-limited AWQPE@0.9 trials** (fail even with infinite shots): only 600 of 14,400 trial pairs from a few phases. Estimates there are unstable (dev +0.2, test +2.3) and are not interpreted beyond "no larger than in shot-fixable trials".
- **Decision analysis** (does P4's one-step signal survive sequentially?):
  - On informative decisions (some block's next batch changes the outcome), eig_cell picks the realized-best block 56–61% of the time, against 33–36% for uniform and 29–33% for random. Fisher picks it 10–31% of the time.
  - Mean error regret per decision (LSB): eig_cell 3.7 (D2) to 14.3 (ε=0.9); uniform 9.2–22.3; Fisher 57–70.
  - Error reduction per extra shot: eig_cell 8.5–9.6 LSB, the highest of the practical policies; uniform 7.0–8.3; Fisher ≈ 0.
  - Pooled within-decision rank correlation between signal and realized gain is small everywhere (≤ 0.06), because most (decision, block) pairs have zero realized gain.
- **Dev versus test.** Test reproduces dev: eig_cell +2.75 / +4.51 / +2.85 on dev against +3.06 / +4.91 / +2.85 on test.
- **Interpretation.**
  1. Observable window evidence, used as the expected information gain about the final n-bit answer, allocates a fixed global shot budget measurably better than uniform, for every decoder.
  2. The effect is concentrated in identifiable regimes: very few shots per block, tight budgets, larger n, and boundary-hard phases. It is near zero once blocks already have about 8 shots each.
  3. Against the faithful ε=0.9 decoder, adaptive allocation recovers only about 31% of the greedy-oracle headroom, because much of that headroom concerns which failures shots can or cannot fix. With ε_safe it recovers 66%. With D2 it matches or beats the myopic oracle.
  4. Simple count signals are decoder-dependent: they help under AWQPE and are neutral or harmful under D2. Legacy C_w is as good as eig_cell only under ε=0.9.
  5. **Fisher-information allocation is consistently and strongly harmful (−9 to −14 points).** It always targets the least-significant block.
- **Limitations.**
  - Ideal Dirichlet model, eigenstate input, no noise (P9 pending), no Qiskit circuit-level validation of the allocation (P10).
  - The primary endpoint is tolerance 2⁻ⁿ only. Other tolerances are for P8.
  - The oracle is myopic, so "efficiency" above 1 is possible and does not mean optimality.
  - The EIG signal uses the ideal model's posterior; its robustness to model mismatch under noise is untested.
  - Absolute gains are modest (about 3–5 points pooled; up to about 8.6 in the best regime).
- **Next decision (for the user).** Shots help only where uncertainty is shot-fixable. For AWQPE@0.9, most of the remaining oracle headroom lies in outcomes that shots cannot fix. That is the evidence base for P6 (adaptive overlap). P6 is not started, pending review.

### P5 dev addendum (2026-09-27)
- The 6 dev shards that failed with out-of-memory errors were re-run and completed (432/432). The failures were caused by running the dev analysis concurrently with the workers.
- On the full dev grid, eig_cell − uniform is **+2.68 / +4.35 / +2.71** points (ε=0.9 / ε_safe / D2), against +2.75 / +4.51 / +2.85 on 416/432 shards. Fisher is −8.10 / −12.63 / −13.42.
- `p5_dev_*.csv` now reflects the full grid. The test was frozen before these shards finished, and nothing about the test depended on them (D-019).

## P6: Adaptive overlap, PILOT (2026-09-27), PILOT / DEV (not a claim)
- **Run.** `p6_adaptive_overlap_pilot/20260927T152618Z` (config `p6_overlap_pilot.yaml`, design D-020/D-021, commit 327148b). Tables: `research/analysis/tables/p6_pilot_*.csv`.
- **Grid (dev split).**
  - Partitions [4,4], [3,2,3], [2,2,2,2] (n=8) and [4,4,4], [3,3,3,3] (n=12).
  - S0 ∈ {2, 4}, r = 2, ΔS = S0.
  - Strata S2_boundary, S1_final_half and S4_uniform: 24 phases × 4 replicates per partition-cell, 120 phase clusters.
  - 10 overlap variants.
- **Sanity checks.**
  - Overlap rescues 0 failures under faithful AWQPE, which cannot use overlap blocks (as required).
  - Every shot-budget arm uses identical total shots.
  - The no-overlap baselines are identical across variants (CRN).
- **Rescue classification** (P5-policy B2 failures, one extra batch; greedy-oracle capability; counts):

  | decoder | subset | P5 failures | shots only | overlap only | both | neither |
  |---|---|---|---|---|---|---|
  | awqpe_ext (ε=0.9, ext v1) | decoder-limited | 36 | 3 | 21–22 | 5–6 | 6 |
  | awqpe_ext (ε=0.9, ext v1) | shot-fixable class | 208 | 33–45 | 68–71 | 14–17 | 78–90 |
  | awqpe_eps_safe_ext | all (none decoder-limited) | 191 | 48–53 | 54–55 | 27–28 | 56–60 |
  | likelihood D2 (ext or bridge) | all | 77 | 1–8 | 15–26 | 45–59 | 1–9 |

  Practical policies (trigger-chosen overlap vs eig-chosen shots) show the same pattern. For decoder-limited ε=0.9 failures, the lowerhalf-triggered O-ext rescues 24–25 + 3–4 of 36, against 1–2 + 3–4 for extra shots.
- **Practical arms** (P(tol), ext_v1_lowerhalf_A1 / bridge_s1_A2):
  - D2: B2 (P5 eig) 92.0%; B4 (eig + overlap, equal shots) 98.0% / 98.5%; B2u (eig, U-matched) 96.5% / 94.9%.
  - awqpe_ext: B2 74.6%; B4 79.8%; B2u 75.9%.
- **Paired differences** (points; phase-cluster 95% CI; Holm across the 10 variants within a decoder and comparison):
  - **D2, B4 − B2u (equal U):** +1.2 to +3.7 for every variant (CIs exclude 0 except ext_v2 at the edge). Bridge s1 A2: +3.65 [2.2, 5.3], Holm p < 1e-4.
  - **D2, B3 − B1 (overlap alone vs uniform):** +7.5 to +9.7, while B3 uses 23–40% *fewer* U-queries than B1.
  - **awqpe_ext, B4 − B2:** +0.3 to +5.2 at equal shots, +3.9 at equal U (lowerhalf A1). The CI is [−0.6, 8.3]; the pilot is underpowered here.
  - **awqpe_eps_safe_ext, B4 − B2:** −4.5 to +1.0 at equal shots and **−4 to −10 at equal U**, i.e. harmful. B3 − B1 is positive (+2.7 to +4.8).
  - Faithful awqpe / awqpe_eps_safe with overlap arms: lower, because overlap shots are wasted on a decoder that ignores them (expected; not an overlap result).
- **Overlap actions (B4).**
  - Most go to the most-significant boundary: roughly 50% at boundary 1, 30% at boundary 2 and 15% at boundary 3.
  - Immediate effect on the estimate:

    | decoder | improved | worsened |
    |---|---|---|
    | D2 | 24–37% | 5–11% |
    | awqpe_ext | 9–14% | 1–2% |

- **Methodological concerns found.**
  1. **Equal-U matching overshoots.** B2u stops at the first batch reaching the target, and least-significant-block batches are expensive, so B2u received on average about 17% *more* U-queries than B4. The pilot's equal-U results are therefore conservative for overlap. For dev/test: match to the closest cumulative U and report both the under- and over-matched baselines as a bracket.
  2. **B4 uses about 8–12% more U than B2 at equal shots**, and B3 uses 25–40% less than B1. Equal-shot comparisons therefore do not imply equal quantum work; both accountings must be reported (they are).
  3. **The O-ext result under safe-ε is negative** despite positive capability. Hypothesis to check on dev before any test: at safe ε the ambiguity flag already resolves most boundary cases, so overlap batches mostly displace useful chunk shots. This needs a diagnostic, not a tuning fix.
  4. **Pilot phases are ambiguity-enriched** (no S0/S3 strata). Effect sizes will shrink on the predeclared full stratum set, which the full dev/test must use.
  5. **The D1 overlap results depend on the non-paper `awqpe_ext` decoder** (D-020). They are overlap-plus-decoder results and must be labelled as such. Bridge has no D1 counterpart.
  6. **v=2 and A=2 add cost without clear benefit** over v=1 and A=1 in this pilot (dev decision).
- **Assessment.** Worth taking to the full dev grid for D2 (both mechanisms) and for awqpe_ext with the lowerhalf trigger. Pending user review; no full grid started.
