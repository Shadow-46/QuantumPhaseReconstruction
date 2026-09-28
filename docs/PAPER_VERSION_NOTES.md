# AWQPE reference: version used and discrepancies

## Reference
A. Shukla and P. Vedula, "Toward Practical Quantum Phase Estimation: A Modular, Scalable, and Adaptive Approach", *Advanced Quantum Technologies* 9(3), e00683 (2026), first published 20 March 2026, doi:10.1002/qute.202500683.

## Which version the code implements
- **The published (Wiley) version: not accessed.** On 2026-09-27 the article page and PDF (`advanced.onlinelibrary.wiley.com/doi/{full,pdf}/10.1002/qute.202500683`) returned HTTP 403 behind a bot-verification challenge. The challenge was not bypassed.
- **arXiv:2507.22460v3 (15 Nov 2025): read in full, all 21 pages.** Its journal-ref is the AQT article above. All of the algorithm (Algorithms 1–2, Eq. 2.1, Lemmas 3.1–3.3, Remarks 2.1–3.6) is implemented from v3.
- arXiv v1 (30 Jul 2025) and v2 (7 Aug 2025) were not used.

## Known or suspected differences between v3 and the published version
1. **Noise.**
   - A search-engine snippet of the published abstract says the simulations show "robustness relative to standard QPE under representative noise models".
   - The v3 abstract says only "accuracy and robustness", and v3 has **no noise section**. Its sections are: Introduction, AWQPE, Proof of Correctness, Complexity, Advantages, Numerical Simulation, Conclusion.
   - So the published version very likely contains a noise study that v3 lacks.
   - **Action:** when the user supplies the published PDF (target `docs/refs/`), diff it against v3 before designing Phase 9 (noise). Then update `research/awqpe/baseline/awqpe.py` if Algorithms 1–2 changed.
2. The title differs by one word: "Toward" (published) versus "Towards" (arXiv).

## Ambiguities in v3 and how the code resolves them
See `DECISIONS.md` D-004. In brief:
- **Tie-breaking.** The paper picks randomly among ties; we use a seeded random choice.
- **Eq. 2.1 gap.** The equation is silent on non-wrap pairs that contain 0 or n−1. We use the ordinary minimum.
- **Flag indexing.** Algorithm 2 indexes `A[j]` from 1 while the text treats `A` as 0-indexed. We read chunk j's flag as `A[j-1]`.
- **Borrow source.** The borrow reads the most-significant bit of chunk j+1 after that chunk's own correction, as written.
- **Final block.** When the final block is ambiguous, its flag is set but t1 is kept.
- **Unspecified shots.** N_shots is not given. The §6.1 counts sum to 10,240 per block.
- **Unspecified simulation.** The "million-case Dirichlet-kernel simulation" gives no phase distribution, shot count or block sizes, so it could not be reproduced as specified.

## Errata found while reproducing v3
- **Table 2, case 4 (φ = 1/√2).** The table lists m = [3,3,3,3,3,3,3,3,3], which is 27 bits. The raw and final binary strings it prints have 30 bits. With ten 3-bit blocks our implementation reproduces both strings and the printed decimal 0.7071067811921239 exactly. With nine blocks the first 27 raw bits match. The partition in the table is a typo: it should list 10 blocks.
- The other five §6 examples reproduce bit for bit in the infinite-shot limit: 0.8203125 with [3,2,3], 0.3 with [2,2], π/6 with [3,2,2,3], 0.671875 with [4,4] and sin(π/12) with [5,6,7,4]. The printed §6.1 counts also decode to the stated result. See `research/tests/test_awqpe_baseline.py`.

## Potential decoder limitation under the arXiv-v3 implementation (PRELIMINARY; see RESULTS_LOG P2a and V1)
These are findings from this project (see `RESULTS_LOG.md`, P2a/P2b). They should be confirmed against the published version and ideally with the authors before anything is claimed publicly.
- **The deterministic failure mode.** Theorem 3.7 excludes the case δ̂_k = 0.5, where the remaining k bits round to exactly 10…0. Remark 2.1 describes the related case as "exceptionally rare" when the *true* fractional part is exactly 0.5.
- **The excluded case is not rare.** The *rounded* case occurs with probability ≈ 2^−k at a boundary with k lower bits, e.g. 12.5% for a final block of width 3. In that case the lower chunks carry no information about which way the upper block rounded.
- **When it resolves correctly.** Algorithm 1 resolves it correctly only if the upper block is flagged ambiguous. Flagging requires ε < ε*(k, m) = K_M((½+e)/M) / K_M((½−e)/M), with e = 2^−(k+1). This is approximately 0.36, 0.61, 0.78, 0.88 and 0.94 for k = 2…6.
- **Consequence at the suggested ε.** At the paper's suggested ε ≈ 0.9, AWQPE is therefore not exact even with infinitely many shots unless every boundary has at least 6 lower bits.
- **Measured failure rates at ε = 0.9 (dense grid, infinite shots):**

  | final block width | failure rate |
  |---|---|
  | 2 | 12.8% |
  | 3 | 4.9% |
  | 4 | 1.8% |

- **The special-chunk rule does not fix it.** Disabling the rule (ablation) moves failures between the two sides of the tie without reducing them.


Independent verification (RESULTS_LOG V1) found the following:
- A second implementation agrees on 262,144/262,144 infinite-shot decodes.
- The theoretical predicate is exact on the generic grid (0 FP/FN).
- All 216 readings of the ambiguities keep the floor at ε=0.9.
- The worked examples select the in-place borrow reading.


## Update 2026-09-28: free arXiv PDF obtained and checked for noise models
- **File.** `docs/refs/Shukla_Vedula_2026_AWQPE_arXiv.pdf`, from https://arxiv.org/pdf/2507.22460. sha256 4bb53982c02fc3104c94e38a54d51bfe55771218b32ab5b5f4f5679c37186064.
  - It is git-ignored, as a third-party PDF under `docs/refs/*.pdf`. The hash identifies it.
- **Version.** The stamp says v3 (15 Nov 2025): the same version as before and the latest on arXiv. There are no algorithmic changes to check.
- **Noise.** v3 defines **no** noise model; the only stochastic model is shot sampling of the ideal kernel. The full passage-by-passage comparison with `model/noise.py` and P9 is in `docs/PAPER_NOISE_COMPARISON.md`.
- **Still open.** The Wiley published version remains unaccessed. Whether it adds noise models is unknown.
