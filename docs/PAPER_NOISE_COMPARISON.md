# AWQPE paper vs this project's noise models

## Source examined

- **File.** `docs/refs/Shukla_Vedula_2026_AWQPE_arXiv.pdf`, downloaded 2026-09-28 from https://arxiv.org/pdf/2507.22460. sha256 4bb53982c02fc3104c94e38a54d51bfe55771218b32ab5b5f4f5679c37186064.
- **Version.** The page stamp reads "arXiv:2507.22460v3 [quant-ph] 15 Nov 2025", so this is the latest arXiv version and the same v3 the implementation was built from (PAPER_VERSION_NOTES.md). All 21 pages were read.
- **Published version.** The Wiley paywalled version (Adv. Quantum Technol. 9(3) e00683) was **not** accessed and is not assumed.
- **Scope of every statement below.** Everything here is about arXiv v3 only. The published abstract reportedly mentions "representative noise models". That content, if it exists, is not in v3 and remains unverified.

## A. Noise models actually present in the paper (arXiv v3)

**No noise channel is defined anywhere in v3.** v3 has no noise parameter, no noisy simulation and no noise result. The only stochastic model is **finite-shot sampling of the ideal outcome distribution**. Every passage that touches noise or errors:

| Where (page) | What it says | Is it a noise model? |
|---|---|---|
| Abstract, p1 | Standard QPE's deep circuits challenge NISQ devices; AWQPE is "well-suited for near-term quantum platforms". | No: motivation only. |
| §1, p1-2 | Iterative PE is "susceptible to error accumulation"; RPE (refs 12-15) uses statistics "to mitigate noise". | No: literature context. |
| §3.2, Lemma 3.3, p8-12 | Outcome law of one block: P(j\|δ) = sin²(2^m πθ) / (2^{2m} sin²(πθ)), θ = δ − j/2^m (the squared Dirichlet kernel). Counts over N_shots are i.i.d. draws, with Hoeffding bounds, e.g. N_shots ≥ (2/ΔP²_min) ln((2^m − 2)/ε₁) and N_shots ≥ (1+ε)² ln(1/ε₂) / (2 p²_{T1} Δ²_R). | **Only shot (multinomial sampling) noise of the ideal kernel.** |
| Remarks 3.5-3.6, p12-13 | Bernstein and Sanov bounds as tighter alternatives for the same sampling problem. | Shot noise only. |
| §4.1-4.2 and Table 1, p14-15 | Controlled-U^(2^j) "is typically synthesized by applying the controlled-U operation 2^j times". Block depth is O(2^(p_start+m_i−1) · C_d(U)), and total applications of U are 2^(Σ_{j<i} m_j)(2^(m_i) − 1). Deep circuits are "highly susceptible to decoherence". | No: a resource count with a qualitative decoherence remark. No error rate is given. |
| Remark 2.2, p6 | Re-run with U′ = e^(i2πΔφ) U to cross-check special chunks. | No: a deliberate deterministic phase shift, not noise. |
| Remark 2.3, p6 | Entropy or KL divergence as alternative ambiguity criteria. | No. |
| §5, p16-17 | "Resilience to local gate errors": faulty block outputs "can be selectively discarded or remeasured"; "statistical filtering or targeted repetition". | No: qualitative claims with no model or experiment. |
| §6, p17-20 | Qiskit simulations (ideal; §6.1 counts sum to 10,240 per block) and a "direct Dirichlet kernel based simulation" over a million cases. | Ideal plus shot noise only. |
| §7 Conclusion, p20 | The post-processing "not only resolves ambiguities ... but also corrects for measurement errors ... in the presence of noise". | No: an unsupported claim, since no noise model or experiment in v3 backs it. |

**Conclusion for A.**
- The paper's only error model is multinomial shot noise on the ideal squared-Dirichlet kernel. This project implements it exactly: `model/kernel.py` in Fejér form, identical to the paper's sin² ratio, with multinomial or inverse-CDF sampling.
- It is validated against Qiskit Aer (P10a).
- It is used in every phase from P1 to P8.

## B. This project's own analytic extensions (`research/awqpe/model/noise.py`)

These are phenomenological channels **defined by this project**, acting on a block's outcome distribution p(y|φ), y ∈ {0..M−1}, M = 2^m, block offset k. The code applies them in this order:

1. **Control-phase jitter.** Each shot sees θ + ξ with ξ ~ N(0, s²). In Fejér form, K_M(θ) = M⁻² Σ_{|d|<M} (M − |d|) cos 2πdθ. Jitter multiplies frequency d by g_d = exp(−2π²d²s²), exactly.
   - Parameters: `jitter_sigma` σ and `jitter_mode`.
   - "fixed" means s = σ for every block. "scaled" means s = σ·2^k, i.e. jitter on φ itself amplified by the block's power.
   - Verified against Monte Carlo averaging (test_kernel) and against Aer circuit averaging (P10b, max |z| 2.9).
2. **Global depolarisation.** p → (1 − λ)p + λ/M.
   - Parameters: `depol_gamma` γ and `depol_mode`.
   - "per_query": λ = 1 − exp(−γ q), with q = 2^k(2^m − 1) controlled-U applications per shot. "fixed": λ = γ.
3. **Readout confusion.** An independent per-bit channel with P(1|0) = `readout_p01` and P(0|1) = `readout_p10`, applied as a tensor product over the m measured bits. It is exact vs Aer ReadoutError (P10b).

**The single point of contact with the paper.** The per-query cost q = 2^k(2^m − 1) is exactly the paper's Table 1 "Total applications of U" for a block (with k = Σ_{j<i} m_j). It also matches the paper's own statement that U^(2^j) is synthesised by repetition.
- The paper supplies the *cost count*.
- The *error law* built on it (exponential white-noise mixing per application) is this project's assumption.
- P10b shows gate-level depolarisation is only approximately of that form. About 43-46% of the deviation is outside the global form, and λ tracks the block's largest power better than total q.

## C. What P9 used that the paper does not support

| P9 element | Supported by v3? |
|---|---|
| Readout p ∈ {0.01, 0.03, 0.05} | No. v3 has no readout-error model; its Conclusion only asserts robustness to "measurement errors". |
| Per-query depolarisation γ ∈ {1e-6, 1e-5, 1e-4} | No. The cost basis q is from the paper (Table 1); the error law and all γ values are ours. |
| Fixed control-phase jitter σ ∈ {0.01, 0.02, 0.04} | No. |
| Combined setting (readout 0.02, γ 1e-5, σ 0.02) | No. |
| Noise-aware likelihood decoding (D2 with the channel in the likelihood) | No. D2 is not AWQPE. |
| Posterior-credible stopping under noise | No. v3 §5 mentions "early stopping conditions" only in passing. |
| Faithful AWQPE (D1) evaluated under the channels above | The *algorithm* is the paper's; the *noise conditions* are ours. |

## Consequences

- **P9 stands as recorded.** It is a study of this project's stylised channels (already labelled so in D-030 and in RESULTS_LOG). No P9 statement is a claim about the paper's noise results, and none should be read as one. Adding the PDF changes nothing in P9. **P9 is not re-run.**
- **What P9 does say about v3's claims.** The Conclusion's claim that the post-processing "corrects for measurement errors" is not evaluated in v3. P9 shows only that, under this project's readout channel, faithful AWQPE at eps_safe loses accuracy at fixed shots and is outperformed by noise-aware likelihood decoding. That is a finding about our channel, not a refutation of a model the paper never specifies.
- **Open.** If the published version contains noise models, they can only be compared once that text is legitimately available, via the user's institutional access. The plan (a P9b with the paper's channels, its own dev/test cycle) is unchanged.
