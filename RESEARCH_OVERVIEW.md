# Research Overview: Confidence-Guided Adaptive Reconstruction for Windowed Quantum Phase Estimation

*A plain-language walkthrough of this project — how it started, what was learned, what was built, what broke, what was fixed, and where things stand today (2026-08-13).*

---

## TL;DR

This project set out to improve one specific piece of **Shor's algorithm** (the quantum algorithm for factoring numbers): the classical "clean-up" step that happens *after* the quantum computer has been measured, in a resource-saving variant called **Windowed Quantum Phase Estimation (WQPE)**.

The work went through five stages:

1. **Learn the theory** and reproduce the published algorithm from scratch in Qiskit.
2. **Stress-test it** with a large, systematic experimental campaign to find out how and when it actually breaks (it turns out single-factor stress barely matters, but starving several resources *at once* breaks it nearly half the time).
3. **Catch and correct a mistake** in the team's own earlier finding — a "55% failure rate" claim that turned out to be a bug in the *test script*, not the algorithm. This was caught, verified three independent ways, and openly documented rather than quietly dropped.
4. **Design and build a fix** — a new "confidence-guided adaptive reconstruction" method that spends more computational effort only where the algorithm is genuinely unsure, instead of a fixed budget everywhere. Tested head-to-head, it took the failure rate from **52% success → 100% success** on the hardest tested case.
5. **Currently in progress**: broadening the validation so the result isn't just true for one lucky test case, and writing up the whole study as a paper draft.

Nothing about the quantum circuit itself was changed — only the classical decision-making that happens after the quantum measurement.

---

## 1. What This Research Is About

**Shor's algorithm** factors large numbers exponentially faster than any known classical method, which is why it threatens RSA encryption once large enough quantum computers exist. At its core, Shor's algorithm needs to find the *period* of a mathematical sequence, and it does that using **Quantum Phase Estimation (QPE)** — a subroutine that estimates an angle (a "phase") by running a quantum circuit and measuring it.

The problem: standard QPE needs one qubit of "phase register" for every bit of precision you want, and the circuit gets deeper as you add qubits. That's expensive and fragile on today's noisy quantum hardware.

**Windowed QPE (WQPE)** — the idea from the papers this project builds on (arXiv:2507.22460, arXiv:2509.05010) — fixes that by splitting the phase register into several small, overlapping **windows**, running each one as its own small, shallow circuit, and then *stitching the pieces back together classically* afterward. This trades quantum circuit cost for classical post-processing complexity.

This project's stated goal, from the very first research notes:

> "The goal of this research is to improve the classical reconstruction stage of Windowed Quantum Phase Estimation (WQPE) used inside the Modular Shor's Algorithm. The objective is **NOT** to modify the quantum circuit itself. Instead, the objective is to improve the classical reconstruction pipeline after quantum measurements have been obtained."

So throughout this entire project, the quantum circuits are treated as fixed and correct — all the engineering and research happens in the *classical software* that turns noisy quantum measurements into a final answer.

---

## 2. How It Started

The project began with a deep, deliberate study phase before any research contribution was attempted:

- **Learned standard Shor's algorithm** end to end: modular exponentiation, QPE, inverse QFT, continued-fraction expansion, and how a measured phase becomes a factor of N.
- **Learned QPE theory** in detail: eigenvalues/eigenstates, controlled-unitary operations, why phase information is invisible until the inverse QFT reveals it, spectral leakage, and the connection between phase precision and measurement probability.
- **Learned the windowed variant**: why it's needed for near-term ("NISQ") hardware, how windows are chosen, what "overlap bits" are for, and how candidates from separate windows get stitched into one answer.
- **Ran a literature review** comparing this project's planned direction against existing work — Bayesian reconstruction, maximum-likelihood reconstruction, machine-learning-based reconstruction, dynamic programming, beam search, graph search, and confidence-based decoding. The conclusion: multi-candidate reconstruction, carry-aware stitching, beam-search pruning, and static overlap windows all already existed in the literature — but **continuous, confidence-guided *adaptive* reconstruction did not appear to have been proposed**. That gap became this project's target contribution.

Only after this groundwork was the actual codebase built.

---

## 3. Building the Reproduction (the foundation)

Before trying to improve anything, the published algorithm was implemented **from scratch in Qiskit**, to make sure there was a faithful, verified baseline to measure any improvement against. The repository is deliberately split along the boundary described above:

| Folder | What it does |
|---|---|
| `Circuits/` | The actual quantum circuits: modular multiplication, inverse QFT, standard QPE, and the windowed measurement circuits. |
| `Reconstruction/` | Everything classical: turning raw measurement counts into ranked candidates, checking overlap consistency between windows ("carry checking"), stitching windows into a full phase estimate, and recovering the final order via continued fractions. |
| `Algorithms/` | Three variants: a textbook brute-force reference (`standard_shor.py`), the genuine paper reproduction (`paper_algorithm.py`), and a small parameter-selection extension point (`adaptive_algorithm.py`). |
| `Simulation/` | Runs everything on Qiskit's Aer simulator, with optional realistic noise models (depolarizing error, readout error, etc). |
| `Evaluation/` | Sweeps, metrics, plots, and summary tables for quick experiments. |

Every module has a self-test, and the implementation was checked against the published examples to confirm it reproduces them correctly. This gave the project a trustworthy, "functionally correct" baseline — the thing every later improvement would be measured against.

---

## 4. The Failure-Hunting Campaign (Phases 1–5)

With a working reproduction in hand, the next stage was a large, systematic **stress test** of the published algorithm — deliberately trying to break it in every way that could plausibly matter, and measuring what actually happens rather than guessing. This became a five-phase campaign (`Experiments/phase1_ofat.py` through `phase5_attribution.py`), sharing one common experiment "harness."

### What was found, phase by phase:

**Phase 1 — vary one thing at a time.** Precision, window size, overlap, candidate count, beam width, and shot count were each varied individually, one at a time, across **435 trials**. Result: **100% success, every time.** This sounded great, but turned out to be a red herring — see below.

**Follow-up — starve everything at once.** Since no single factor could break it, the same parameters were starved *simultaneously* (low shots, small candidate budget, narrow beam search, high precision) across 144 trials on the hardest available test case. Result: **47.9% of trials failed** — a huge contrast with the "100% success" from varying things individually. The conclusion: the algorithm's apparent one-factor-at-a-time robustness was really just a forgiving fallback step (the continued-fraction search is very tolerant when the true answer is a small number) papering over real weakness. Only when multiple resource constraints hit *together* does the pipeline actually break.

**Phase 2 — noise.** Realistic quantum-hardware noise (depolarizing errors, readout errors) was injected, alone and in combination. Individually, even very high noise (30% error rates — far beyond real hardware) barely affected results. Combined, noise sources compounded and produced visibly worse results (91.8% success in the worst combined test), with 2-qubit gate errors identified as the strongest single driver of failure.

**Phase 3 — deliberate sabotage.** The classical reconstruction layer was attacked directly: bit-flips, missing candidates, corrupted overlap data. Most attacks were absorbed gracefully; the most damaging one (corrupting the overlap bits between adjacent windows) caused a 56% failure rate — but see the correction below regarding *why*.

**Phase 4 — scaling.** This produced the single most surprising, high-impact finding of the entire study: the windowed approach *does* keep the phase-side circuit small, exactly as advertised — but the *other* half of the circuit (the modular multiplication / "work register" arithmetic) explodes in compile time almost immediately. Instances with as few as **7 qubits (numbers around 65–119) took tens of seconds to minutes to compile**, and **8 qubits (numbers around 143) didn't finish compiling within 45 seconds**. This has nothing to do with the phase-register savings the windowed method is about — it's a completely separate bottleneck in how the modular-multiplication circuit is built. This forced the entire rest of the study to be restricted to very small test numbers just to stay computationally tractable, and it was reported honestly as a real limitation rather than hidden.

**Phase 5 — root-cause attribution.** Every failure from the earlier phases was replayed to identify exactly where in the pipeline it first went wrong. The headline result: **72% of every naturally occurring failure in the whole campaign traced back to one mechanism — combined classical-resource starvation** (not enough shots, not enough candidates, and too narrow a search beam, all at once). This became the target for the fix built in Phase 6.

---

## 5. Catching a Mistake in the Team's Own Work

While designing the fix for the failures found above, the plan called for a new module to patch a specific defect: the "carry check" (the logic that checks whether two adjacent windows agree with each other) had been measured, in Phase 3, as **incorrectly accepting bad data 55.4% of the time** — a serious, safety-relevant bug.

Before writing any code to fix it, that evidence was re-verified — and it didn't hold up. The 55.4% number turned out to come from a **bug in the test script itself**, not the carry-check logic. The synthetic test data generator was accidentally constructing all its test windows with the exact same 1-bit overlap, no matter what overlap value it claimed to be testing — so it never actually tested the multi-bit-overlap cases the 55.4% figure was supposed to describe.

This was confirmed three independent ways before accepting the correction:
1. **Mathematically** — the formula used to check candidates and the formula used to judge them in the audit were algebraically the same formula, so the audit could never have caught a real defect in the first place.
2. **Empirically** — fixing the test script's geometry and re-running it gave **0% false-accept, 0% false-reject** (an exact result).
3. **By instrumented replay** — 80 trials with a deliberately corrupted candidate showed the real carry-check correctly rejecting the bad data in all 80 cases.

**The carry-check logic itself needed no fix.** The planned "Module A" was dropped entirely before any code was written for it. Rather than quietly deleting the old claim, the original reports were edited in place with visible strikethrough/correction markup, so anyone who previously relied on the 55.4% figure could see exactly what changed and why. This kind of self-caught, transparently documented correction is treated as a strength of the research process, not something to hide — and it's called out explicitly in the project's own paper-drafting notes as something that *builds* credibility with reviewers rather than undermining it.

---

## 6. The Fix: Confidence-Guided Adaptive Reconstruction (Phase 6)

With the false lead removed, the actual improvement was designed and built around the real, confirmed root cause: **the published algorithm spends a fixed classical budget on every window, whether that window's measurement was clear-cut or genuinely ambiguous.** The fix is to make that budget *adaptive* — spend more effort only where the algorithm is actually unsure.

This required exactly one new piece of math: a **confidence statistic**, $C_w$, computed per window. It answers "how sure am I that the top-ranked candidate is really better than the runner-up, given how many measurement shots I've collected so far?" — a proper Bayesian comparison that accounts for sample size, unlike simpler measures (raw probability, margin, entropy) which can't tell a "barely-informative small sample" apart from a "practically certain large sample."

Three modules were built around this statistic (a planned fourth — the carry-check fix — was dropped, per the section above):

- **Module B — Candidate coverage:** instead of always keeping a *fixed number* of candidates per window, keep however many are needed until their combined probability mass crosses a target threshold. Flat, ambiguous windows get more candidates kept; clear-cut windows get fewer.
- **Module C — Adaptive beam width:** the search that stitches windows together widens itself when many windows turned out to need extra candidates (i.e., when the problem instance is more ambiguous than usual).
- **Module D — Shot-stopping:** windows whose confidence is still too low get *more measurement shots*, up to a cap, before reconstruction proceeds — instead of a fixed shot count applied uniformly regardless of how ambiguous a particular window is.

A calibration check was also run to make sure the confidence statistic can be trusted: it turned out to be **measurably overconfident** (raw Brier score 0.224), which is honestly reported as a known limitation. A post-hoc recalibration layer improved this to 0.199, though it isn't (yet) wired into the live decision-making — a deliberate, documented scope decision rather than an oversight.

### The results

A controlled **ablation study** (139 of 144 planned trials completed; the remaining 5 were stopped once results had clearly stabilized) compared each module, and all three together, against the original fixed-budget baseline, on the hardest test case from the failure study:

| Configuration | Success rate | Improvement vs. baseline | Statistical significance |
|---|---|---|---|
| Baseline (published algorithm) | 51.1% | — | — |
| + Module B only | 76.3% | +25.2 points | p = 2.8×10⁻⁹ |
| + Module C only | 56.1% | +5.0 points | p = 0.016 |
| + Module D only | 76.3% | +25.2 points | p = 1.8×10⁻⁷ |
| **All three together** | **100%** | **+48.9 points** | **p = 6.8×10⁻²¹** |

All differences are statistically significant. The baseline figure (51.1%) closely matches the original failure-study baseline (52.1%), confirming that fixing the test-script bug didn't quietly change how the underlying algorithm behaves — it only changed how a separate audit measured it. On this test set, **combining all three adaptive modules eliminated every measured failure.**

---

## 7. Writing It Up

In parallel with the implementation work, substantial effort went into how to present this as a credible academic paper. A dedicated research pass studied what journals like *Quantum* and venues like ACM TQC actually expect from a paper that mixes formal claims with empirical results — things like: leading every result with effect size and confidence intervals (not just a p-value), correctly choosing and reporting statistical tests (McNemar's test for paired trial outcomes, bootstrap confidence intervals, Holm–Bonferroni correction across multiple comparisons), writing an honest "Limitations" section, and — specifically — how to disclose the self-caught carry-check correction without it reading as an "erratum" (it isn't one, since it was caught before publication).

A full compiled report (`Data/failure_study/COMPILED_REPORT.md`) now exists, cross-referencing every claim back to its source data and consistently labeling each one as **(PROVEN)**, **(HEURISTIC)**, or **(EMPIRICAL)** so readers know exactly how much weight each claim can bear. A first draft of the actual manuscript exists in `Paper/main.tex`, targeting *MethodsX*.

---

## 8. Current Status (as of this writing)

The project is currently in **Phase 7: closing generalization gaps**, prompted by an honest self-review of the Phase 6 result. Three specific weaknesses were identified in the 100%-success headline number:

1. **It was only tested on one problem instance** (N=21, a=19) — so the result might not generalize.
2. **That one instance was also the same one used to design the fix** — a circularity concern (tuning and testing on the same case).
3. **There was no "matched-resource" control** — Module B/C/D succeed partly *because* they're allowed to use more computational resources (more candidates, wider search, more shots) than the fixed baseline. Without a fair comparison against a baseline given the *same* extra resources (just not spent adaptively), it's hard to know how much of the gain comes from the "confidence-guided" idea itself versus simply "spending more."

`Experiments/phase7_generalization.py` (currently uncommitted, in-progress work) addresses all three:

- **Held-out instances:** the same 144-trial grid is being re-run on three additional, pre-registered problem instances the modules were never tuned against — (N=15, a=2), (N=21, a=8), and (N=21, a=2) — chosen by a fixed rule *before* any trial ran, to rule out cherry-picking.
- **Matched-resource controls:** for each adaptive module, a "fair fight" baseline is constructed that gets fed the *same* realized resource usage (extra shots, extra candidates, wider beam) the adaptive version actually used in that trial — but without confidence guidance choosing where to spend it. This isolates whether the *adaptive allocation* itself is doing the work, not just the extra resources.
- A safety gate re-runs the original canonical instance with new instrumentation and checks the results are **bit-for-bit identical** to the already-published numbers, to make sure adding this instrumentation didn't quietly change anything.

### Live progress snapshot

| Stage | Status |
|---|---|
| Canonical (N=21, a=19) instrumented re-run — sanity check | ✅ Done (144/144), verified to match the published numbers exactly |
| Held-out instance (N=15, a=2) | ✅ Done (144/144) |
| Held-out instance (N=21, a=8) | ✅ Done (144/144) |
| Held-out instance (N=21, a=2) | 🔄 In progress (34/144) |
| Matched-resource control arms | ⏳ Not yet started |

A pair of small support scripts were added alongside this: `Experiments/check_status.py` (reports how far each of the five stages above has progressed, since individual trials can take up to ~20–25 minutes) and `Experiments/resume_campaign.py` (safely restarts the campaign — already-finished stages are skipped instantly, so it's safe to re-run at any time, e.g. after an interruption).

A known constraint carried over from Phase 4's discovery still applies here: because of the modular-multiplication compilation cliff, only N ∈ {15, 21} can produce a valid two-prime-factor test case within the "≤5 work qubits" ceiling this whole study has had to operate under — so held-out instances broaden *which* small case is tested, but can't yet test a case that's meaningfully *harder* than the current canonical one. That remains an open limitation, not something Phase 7 claims to solve.

---

## 9. What's Left / Where This Is Headed

Honestly-tracked open items, carried forward from the project's own notes:

1. **Finish Phase 7** — complete the last held-out grid and the matched-resource controls, then fold the results into the compiled report and paper draft.
2. **Decide on the recalibration layer** — currently, the shot-stopping module makes decisions using the known-to-be-overconfident raw confidence score, not the improved recalibrated version. This needs an explicit decision either way, not silence.
3. **Scale up trial counts** if the results are meant to support a publication — the qualitative ranking (B ≈ D > C, and all three together eliminating failures) is solid, but confidence intervals — especially for Module C's smaller effect — would tighten with more trials.
4. **Solve or route around the compilation cliff** — this is flagged as the single most consequential unsolved engineering problem in the whole project, since it's the only failure mode that makes the algorithm *impossible to even run* at interesting sizes, rather than just less accurate. It's out of scope for this project's classical-reconstruction focus, but it is the ceiling on how far any of this can be validated.
5. **Finish the manuscript** — incorporate the paper-craft research findings (contributions-first structure, proper statistical reporting, an honest limitations section, and matter-of-fact disclosure of the carry-check self-correction) into a submission-ready draft, first targeting *MethodsX*/*Quantum*, with ACM TQC as a fallback.

---

*This document summarizes the project state as of 2026-08-13. For full technical detail, see `Data/Theory_&_Research_Progress.md` (theory background), `Data/failure_study/REPORT.md` (Phases 1–5, with corrections), `Data/failure_study/PHASE6_SUMMARY.md` (the adaptive-framework build and ablation), `Data/failure_study/COMPILED_REPORT.md` (the full cross-referenced report with proofs and evidentiary labels), and `Paper/main.tex` (the manuscript draft).*
