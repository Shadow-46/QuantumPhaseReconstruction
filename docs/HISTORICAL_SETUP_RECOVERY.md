# Historical setup recovery: Paper A experiments

Recovered on 2026-09-27 from the repository at the frozen commit
`90bc134811b09e12c123a5599128fb02372371cb` (tag `paper-a-v1.0-frozen`).

**Status: recovered from code; nothing below is from memory.** Every claim cites the file that establishes it.

## Unitary
Modular multiplication U_a|y⟩ = |a·y mod N⟩ on ⌈log₂N⌉ work qubits. States y ≥ N are fixed points so the map is a permutation. It is implemented as a dense `UnitaryGate` (`Circuits/modular_multiplication.py:73-115`). The controlled power U_a^(2^e) is built as the permutation for a^(2^e) mod N.

## Initial state: the computational basis state |1⟩, not an eigenstate
- Every Aer-sampled window circuit uses `eigenstate=1`, which applies X to work qubit 0 and so prepares |1⟩. See `Algorithms/paper_algorithm.py:83`, `Algorithms/adaptive_reconstruction.py:110`, `Experiments/phase6_calibration.py:53`, `Experiments/phase4_scaling.py:46,182` and `Circuits/windowed_qpe.py:98-100`.
- |1⟩ is the uniform superposition (1/√r)Σ_s |u_s⟩ of the r eigenstates of U_a, with eigenphases s/r, where r is the multiplicative order of a mod N. The repo says this itself in `Data/failure_study/FORMULATION.md:263-265` and `Algorithms/paper_algorithm.py:136-139`. `Data/failure_study/THEORETICAL_ANALYSIS.md:30` calls |1⟩ an "eigenstate", which is incorrect.
- **Consequence for windowed QPE.** Each window is a separate circuit execution, so each shot of each window collapses independently onto a random s. Different windows therefore do not observe the same phase. Every window's outcome distribution is the mixture Σ_s (1/r)·p(y | s/r). That suits Shor order-finding, where any s/r is useful, but it is not the single-eigenphase setting AWQPE is designed and proved for.

## Instances (N, a)
- The corpus has 41 pairs (`Data/failure_study/corpus.json`), built by `Experiments/harness.py:99-122`.
- The experiments that were actually run were limited to at most 5 work qubits: (15,4) r=2, (15,13) r=4, (21,8) r=2 and (21,19) r=6 (`harness.representative_subset`).
- The Phase 7 held-out instances were (15,2) r=4, (21,8) and (21,2) r=6 (`Experiments/phase7_generalization.py:70-74`).
- For r ∈ {2,4} every eigenphase s/r is dyadic, so windows never show genuine top-two ambiguity; outcomes are deterministic up to the mixture. Only r=6 (N=21) produced non-dyadic phases. This is consistent with `paperA_synthetic_grid.py:171-172`, where "resolvable" marginals came only from (21,19).

## Default configuration (`Experiments/harness.py:148-160`, `config.py`)
Precision n=8, window width 4, overlap 2, candidate_count 2, max_paths 32, 2048 shots, no noise, seed 314159. **All windows of one trial shared the same Aer `seed_simulator`** (`Algorithms/paper_algorithm.py:87`).

## Reconstruction (not AWQPE)
- The method uses overlapping windows with top-k candidates and carry-aware beam stitching, followed by continued fractions for the order. See `Reconstruction/candidate_generation.py`, `carry.py:62-76`, `stitching.py:46-69` and `continued_fraction.py`.
- There is no top-two ratio test, no LSB→MSB borrow correction and no special chunk. The repo acknowledges this at `Data/failure_study/REPORT.md:349-351`.
- The 0.9 ratio appears only as a diagnostic in `Experiments/phase3_adversarial.py`, attributed there to the paper.

## Other findings
- **Terminology.** "Dirichlet" in Paper A refers to a Dirichlet posterior over multinomial counts (`Reconstruction/confidence.py:69-102`, a C_w diagnostic). The repo contains no Dirichlet-kernel measurement model; every "true" distribution came from an exact statevector.
- **The only true-eigenstate experiment** was the Phase 3 ambiguity sweep: a PhaseGate with target |1⟩, width 4, 4096 shots × 20 repeats (`Experiments/phase3_adversarial.py:209-263`).
- **Seed reproducibility.** `phase1_ofat.py:78` and the Phase 3 battery derive seeds with `hash()` of tuples that contain strings, which depends on PYTHONHASHSEED. The recorded seeds still allow replay, but re-deriving them does not reproduce them.
- **Environment.** The pins in `requirements.txt` (qiskit 2.4.1, qiskit-aer 0.17.2, numpy 2.4.6, scipy 1.17.1) match `C:\Python314`, not the repository `.venv`. The `.venv` is a WSL venv with qiskit 2.5.0 and no Aer.

## Implications for the new work
1. Paper A's "top-two / C_w" evidence was measured on mixture distributions from a non-eigenstate input, mostly with dyadic phases. It is not evidence about AWQPE windows on an eigenstate.
2. The new track (`research/awqpe`) uses exact eigenstates by default. It studies non-eigenstate inputs deliberately, as a controlled factor (Phase 3), including this historical |1⟩/U_a case.
