# Unified Information-Guided Adaptive AWQPE: results summary (as of 2026-09-28)

This is a one-page index of conclusions. The authoritative record is `docs/RESULTS_LOG.md` (append-only, chronological) and `docs/DECISIONS.md` (D-001 to D-035). Every claim below is labelled with its evidence status.

## Scope and conventions

- **Model.** Ideal Dirichlet-kernel model unless marked. Qiskit Aer reproduces it (P10a/P10b).
- **Algorithm reference.** AWQPE is taken from arXiv v3. The published noise models have not been checked because the published PDF has not been supplied (docs/refs/ is empty).
- **Decoders.**
  - **D1 is faithful AWQPE** (Algorithms 1+2) at eps 0.9 and at eps_safe.
  - **D2 is a grid-likelihood MAP decoder, not AWQPE.**
  - **awqpe_ext is a non-paper extension** ("adaptive overlap + extended decoder"). It is never described as faithful AWQPE.
- **Claims kept separate.** (A) D2 information; (B) extended decoder; (C) resource allocation.
- **Statistics.** Phase-cluster bootstrap CIs, sign-flip tests, and Holm correction within each predeclared primary family. Each design was frozen before its single held-out test run.

## Held-out conclusions

| Phase | Question | Held-out result | Status |
|---|---|---|---|
| V1 | Is there an eps = 0.9 failure floor? | Yes, exactly predicted: failure when the lower bits round to 10...0 and the upper chunk is unflagged. It is independent of 216 readings of the paper's ambiguities. | verified (exhaustive) |
| P3 | Initial state | The old repository's \|1> input on U_a is an eigenphase mixture. Faithful AWQPE collapses on exact ties; eigenstates give 100%. | DEV (deterministic) |
| P5 | Global shot allocation | Information-gain (eig_cell) allocation beats uniform: +2.9 to +4.9 pts success, Holm p < 1e-4. Fisher information is harmful (-8.8 to -13.8). | **HELD-OUT TEST** |
| P6 | Adaptive overlap vs extra shots at equal U (upper bracket) | D2: bridge +0.60, ext +0.56 pts. awqpe_ext: +6.60 (eps 0.9), +2.47 (eps_safe). All Holm-significant. Overlap rescues 122-131 of 183 decoder-limited failures; shots rescue 14-22. | **HELD-OUT TEST** |
| P6 follow-up | Why extended-decoder rescue falls at n = 16 | Placement, not efficacy. Overlap at the responsible boundary rescues 91-100% at every n. The whole-lower-part 10...0 gate rarely fires at upper boundaries (0% at n = 16). | post hoc, descriptive |
| P8 | Tolerance stopping (Mode B) | Posterior-credible stopping keeps coverage >= 0.95 in all 96 cells. It saves 33-67% of shots vs a dev-tuned fixed budget at fine tau (4/9 primary cells). At coarse tau it costs 0.02-0.25 extra shots. 2/28 dev-tuned fixed baselines under-cover on test. | **HELD-OUT TEST** |
| P8 Mode A | Tolerance success vs budget | eig_cell beats uniform at every tau in {2^-3, 2^-(n/2), 2^-n}, all 27 cells positive; the gain grows as tau tightens. | re-analysis of P5 test |
| P7 | Unified shots + overlap + stop controller | Beats the same controller without overlap in 11/12 cells. The saving is small: 0.1-6.5% of shots, largest at the finest tau. Most of the gain over uniform + stop (up to 5.6 shots) comes from allocation; overlap adds up to 1.35. ext is U-neutral; bridge costs 0.7-12% more U. | **HELD-OUT TEST** |
| P9 | Noise (analytic channels) | Noise-aware D2 beats faithful AWQPE eps_safe in 12/12 cells (+6.7 to +20.4 pts). Modelling the noise matters most for stopping: an ideal-model posterior under-covers (down to 0.53 at strong jitter); the noise-aware one holds >= 0.95 (primary cells). | **HELD-OUT TEST** (stylised channels) |
| P10a | Aer vs kernel, ideal | Distributions match (480 blocks, KS p 0.72; 200k-shot recheck, TV 0.001). End-to-end is equivalent. A post hoc sign pattern did not replicate. | validation |
| P10b | Aer noise vs analytic laws | Readout and jitter match exactly. Gate-level depolarisation is only approximately global: 43-46% of the deviation is outside the global form, and lambda tracks the block's largest power better than total q. | validation |

## What the evidence supports (short form)

1. **Allocation.** Allocation guided by expected information gain is the largest and most robust adaptive gain. Fisher information is the wrong signal for this task.
2. **Overlap.** Adaptive overlap adds real but modest information on top of allocation and stopping. Its large effects appear only in combination with an extended decoder, which is a new decoder and not AWQPE.
3. **Stopping.** A posterior-credible stopping rule is calibrated in the ideal model. It is also more robust out of sample than a dev-tuned fixed budget. Under noise, it stays calibrated only if the likelihood models the noise.
4. **Faithful AWQPE.** Faithful AWQPE is dominated by likelihood decoding at every budget, precision and noise setting tested.

## Least robust claims (read with care)

- **Borderline P7 cells** (Holm p 0.025-0.046): ext at [4,4], tau 2^-4; bridge at [4,4,4], tau 2^-6; bridge at [4,4,4,4], tau 2^-8. The P10a replication showed that single-run phase-cluster CIs for small effects can be optimistic.
- **P9 depolarisation numbers.** These come from a stylised global channel (see P10b). Only the qualitative conclusions carry over to circuit-level noise.

## Not done / open

- **Published noise models.** Needs the published PDF placed in `docs/refs/`. The plan is to diff it against v3 and run a P9b if its noise models differ.
- **GPU Aer benchmark.** Needs installing qiskit-aer-gpu in WSL, which was deferred rather than done unattended. The circuits used are 3-7 qubits, where CPU is adequate.
- **A boundary-local extended decoder.** This could remove the n = 16 placement limitation. It would be a new decoder needing its own dev/test cycle.
