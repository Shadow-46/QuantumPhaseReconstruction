# AWQPE implementation notes

How `research/awqpe` maps onto Shukla & Vedula (arXiv:2507.22460v3; see `PAPER_VERSION_NOTES.md`).

| Paper element | Code | Fidelity |
|---|---|---|
| Block i: m_i controls, powers U^(2^(k+p)), p = 0..m_i−1, with k = Σ_{j<i} m_j (Alg. 1 l.4–12) | `blocks/geometry.py::partition_blocks`; `circuits/qiskit_blocks.py::phase_gate_block` | exact |
| P(j\|δ) = sin²(2^m πθ)/(2^{2m} sin²πθ), θ = δ − j/2^m (Sec. 3.2) | `model/kernel.py` (Fejér form, D-003) | exact (tested vs statevector to 1e-14) |
| Eigenstate \|u⟩ on the target (Alg. 1 l.6); the figures use X on one qubit | `phase_gate_block`: U = P(2πφ), \|1⟩ | exact eigenstate |
| Top-two t1*, t2*, random tie-break (Alg. 1 l.14–15) | `baseline/awqpe.py` | exact (D-004a) |
| Ambiguity C(t2*)/C(t1*) > ε, ε ≈ 0.9 (Alg. 1 l.18) | same | exact; ε is swept |
| min(t1*, t2*) mod 2^m if not the final block (l.20–21, Eq. 2.1) | `modular_min` | exact, plus D-004b for the gap |
| Special chunk: rightmost non-zero chunk = 10…0 (Alg. 2 l.3–14) | `ambiguity_resolution` | exact |
| LSB→MSB borrow using the next chunk's MSB, suppressed if ambiguous or special (Alg. 2 l.15–23) | same | exact (D-004c,d) |
| Target: best n-bit ⌊2ⁿφ + ½⌋ (Thm 3.7) | `evaluation/metrics.py::best_nbit_integer` | exact |

**Not from the paper (extensions, always labelled):**
- the likelihood decoder (D2);
- the analytic noise channels;
- overlap blocks (`role="ext"`/`"bridge"`);
- adaptive allocation, overlap and stopping;
- the special-chunk ablation.

**Qiskit conventions.**
- Control p is measured into classical bit p. The counts key's integer is the kernel outcome y.
- Controlled powers are synthesised as a single CPhase(2π·frac(2^(k+p)φ)), which is exact.
- Circuit depth for noise modelling is handled separately; the repeated-U construction will be used for Phase 10 noise cross-checks.
