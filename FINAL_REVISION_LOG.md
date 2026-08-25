# Final Revision Log

**Manuscript:** *Confidence-Guided Adaptive Reconstruction for Windowed Quantum Phase Estimation*
**Target venue:** MethodsX (Elsevier)
**Scope of this revision:** manuscript only. No algorithm, module, implementation, experiment, or
numerical result was changed. The four data figures were re-rendered from the frozen CSVs; no
plotted value differs from the value it had before.

---

## 1. The claim audit that drove the revision

Every major claim was classified before editing. Categories: **A** mathematically derived,
**B** formally proven under stated assumptions, **C** implementation invariant, **D**
experimentally validated, **E** engineering observation, **F** hypothesis, **G** unsupported or
overstated. Only the entries that changed the manuscript are listed.

| Claim (as previously worded) | Was | Should be | Action taken |
|---|---|---|---|
| "The framework **provably** reduces to the published pipeline" | stated as B | **C** — a property of one program, established by branch inspection | Converted `proposition` → `invariant`, new `\INVAR` label, retitled "Equivalence to the published baseline", proof environment relabelled "Justification (code inspection, not a mathematical proof)" |
| "The carry check is an exact classifier" | B | **C** | Relabelled `\INVAR`; commentary now says the object constrained is the shipped predicate |
| "$C_w$ … Bayesian **confidence statistic**" | reads as A/B | **A conditional on a false assumption** | Renamed throughout to **posterior comparison statistic**; added an explicit "What $C_w$ is, stated precisely" paragraph naming (A2) as false and disclaiming calibrated-probability / coverage / error-rate readings |
| "no module introduces a new complexity class" | stated unconditionally | **A, but only for classical work** | Qualified: true for classical asymptotics with $m_w,P'$ bounded; explicitly wrong unrestricted, because module D changes the number of circuit executions |
| Module C presented alongside B and D | implied D | **F/E at best** | Repositioned everywhere (see §4) |
| "Success … 100%" as headline | D, but regime-bound | **D, narrower** | Headline is now the four-instance range 97.9–100% |
| Adaptive allocation is what produces the gain | implied D | **not supported** | Reframed to automatic budget *selection* (see §3) |
| Isotonic recalibration "improves Brier to 0.1988" | read as D | **D, in-sample only** | Now stated as an in-sample fit bounding what recalibration could achieve |
| Validation covers the failure landscape | implied | **G — validation is noiseless** | New Limitation 3 and a new design paragraph (see §5) |
| "Dynamic-programming and graph-search reconstruction give exactness…" | read as a literature claim, uncited | **G — uncited** | Reworded to an algorithmic-limit statement; explicitly says no named implementation was benchmarked |
| "To the best of our knowledge … has not been proposed" | F | **F, but weakly evidenced** | Replaced with a survey-outcome statement plus a positioning table |

---

## 2. Sections substantially rewritten

| Section | What changed |
|---|---|
| **Abstract** | Rewritten to the four-part structure: problem → method → results with the matched-control qualification → scope and limitations. Now states the held-out range, all three headroom gains, the three null controls, the noiseless scope, and the Brier score. |
| **Highlights** | Rewritten; a fourth bullet added carrying the matched-control finding. |
| **Introduction** | New third paragraph stating up front what the validation does and does not establish. Contributions 2 and 3 reworded to match the downgraded claims. |
| **Related work / novelty** (§2.3) | New comparison table (Table 2) positioning the work against standard QPE, windowed/modular WQPE, and Bayesian amplitude/phase estimation across six dimensions. Closing positioning paragraph rewritten; the DP/graph-search paragraph de-claimed. |
| **§2.4 posterior comparison statistic** | Retitled. Model construction spelled out (Binomial marginals, Jeffreys priors, resulting Beta posteriors). Two equations merged into one. New precision paragraph on what the statistic is and is not. |
| **§2.6 Module C** | New "Evidential status" paragraph at the point of definition. |
| **§2.7 retracted module** | Shortened by roughly a third; quantitative detail deferred to Appendix C. |
| **§2.10 Properties** | Label legend expanded to four categories with definitions. Two results moved from `\PROVEN` to `\INVAR`. |
| **§2.11 Computational cost** | Fully restructured into (a) classical asymptotics, (b) bounded overhead from the caps, (c) practical wall-clock. New explicit statement that module D's extra sampling is not free merely because $S_{\max}$ bounds it, with the realised 9.3 → 54.9 shots/window inflation quoted. |
| **§3.1 Design** | New paragraph: all three validation stages are noiseless, why, and what that excludes. New paragraph qualifying the word "pre-registered" (repository artifact, not third-party timestamp). |
| **§3.3 Ablation results** | Canonical table deleted and merged into the per-instance table. Module C now described as a null result rather than a small positive one. |
| **§3.5 Matched-resource controls** | New "How exact is the matching?" paragraph quantifying that B's control was slightly *over*-matched (k = 4.29 vs realised $m_w$ = 4.25), D's slightly *under*-matched (43.70 vs 43.82 extra shots/window, a 0.26% shortfall in the direction favouring D), and C's exact on every trial. |
| **§3.6 Calibration** | In-sample caveat added; closing paragraph rewritten to state the two consequences for an adopter. |
| **§3.7 Baseline reproduction** | Now states that the original campaign drove a *different* entry point (`Algorithms/paper_algorithm.py`), making the agreement cross-implementation evidence rather than a self-comparison. |
| **Limitations** | Reordered and expanded to 14 items; new noiseless-validation item placed third; two modelling simplifications merged; in-sample recalibration gap added; simulation-only item rewritten to name hardware-specific exposures. |
| **Reproducibility** | Trimmed by roughly half; artifact map moved to a new Appendix B. |
| **Figures 2–5** | All four re-rendered (see §6). |

---

## 3. The central claim, before and after

**Before.** The abstract attributed the improvement to confidence guidance and led with
52.1% → 100%.

**After.** The manuscript states, in the abstract, the introduction, the related-work
positioning, the results, the discussion and the limitations, that:

> the evidence supports **automatic, instance-adaptive selection of the reconstruction budget** —
> removing the need to choose that budget correctly in advance — **rather than a demonstrated
> superiority of confidence-guided placement over an equivalent uniform allocation.**

This is the reading the matched-resource controls actually support: three of four
adaptive-vs-control comparisons are null (p = 1.00, 1.00, 0.67), and the controls are oracles that
could only be constructed after observing the adaptive rule's own choices.

---

## 4. Modules repositioned according to evidence

| Module | Evidence | How it is now presented |
|---|---|---|
| **B** candidate coverage | Largest, most stable effect; improves every instance with headroom; ≤2 regressions anywhere | Strong empirical support; recommended first adoption |
| **D** shot stopping | Large where shots bind, nil where they do not; most two-sided discordant counts (6,41), (8,50), (8,8) | Strong support *conditional on headroom*; flagged as the least stable per-trial and the dominant wall-clock cost |
| **C** beam width | +4.9 pp on the design grid; **fails Holm correction**; null on two of three held-out instances; identical to its control on every trial | Explicitly *not* independently validated. Framed as a combinatorial safeguard retained on the strength of Prop. 6 (cannot narrow the beam) and the 20 compounded failures. Recommended against enabling alone. |

Module C was **not removed** — it remains part of the completed framework and the Full arm. Only
the interpretation changed.

---

## 5. Claims weakened, corrected, or newly disclosed

1. **"Provably reduces to baseline"** → implementation invariant, verified trial-by-trial.
2. **"100% success"** → 97.9–100% across four instances; the canonical 100% is explicitly labelled
   the more fragile of the two numbers.
3. **"No new complexity class"** → true of classical work only; module D's additional circuit
   executions are called out as a resource the baseline does not spend.
4. **$C_w$ as confidence** → posterior comparison statistic; calibrated-probability reading
   explicitly disclaimed.
5. **Recalibration Brier 0.1988** → disclosed as an in-sample fit with no held-out split.
6. **Noiseless validation** — newly disclosed. Every trial in the validation section was sampled
   without a noise model, while 27% of the motivating campaign's failures were noise-attributed.
   The manuscript now makes no claim in either direction about behaviour under noise, and
   separates B/C (limited exposure) from D (plausibly adverse: resampling can sharpen confidence
   in a wrong pattern under a noise floor).
7. **"Pre-registered"** → qualified as a repository artifact preceding the data, not a
   third-party-timestamped registration.
8. **Module C significance** → the enlarged Holm family (9 comparisons rather than 5) means C no
   longer survives correction; this is stated in the statistics section, the results, the
   discussion and the limitations.
9. **Matched-control exactness** → the two approximate controls' residual errors are quantified
   and their direction stated.

---

## 6. Figures and tables

**New figure pipeline.** `Experiments/paper_figures.py` renders Figures 2–5 from frozen CSVs. It
runs no trial and asserts the calibration Brier scores still equal 0.2239 / 0.1988 before writing.
The original campaign plots remain in `plots/` untouched.

| Figure | Problem in the previous version | Fix |
|---|---|---|
| 2 failure attribution | Raw identifiers as rotated tick labels; no grouping; the 72% message invisible | Horizontal bars, human-readable labels, colour-coded by mechanism, explicit bracket annotating "classical budget starvation: 72%" |
| 3 held-out | Single crowded 20-bar panel; rotated value labels; ceiling effect unexplained | Four small multiples, one per instance; per-panel baseline and ceiling reference lines; ceiling instance annotated in place; collapsed to a single label where all arms tie |
| 4 matched-resource | Central conclusion not stated visually; baseline label collided with bars | "not significant"/"significant" annotations with p-values over each pair; baseline moved into the legend |
| 5 calibration | Two separate diagrams; overconfidence not emphasised; stopping threshold unmarked | Single two-panel figure; overconfidence gap shaded; marker area ∝ bin population; module D's $1-\varepsilon = 0.95$ threshold marked; the 0.96-predicted / 0.79-observed point called out; in-sample nature of the right panel stated in the caption |

**Tables.** The canonical ablation table was deleted and merged into the per-instance table
(removing a genuine duplication). The novelty table was added. The Holm table was already updated
to the family of nine. The artifact map moved to Appendix B. The cost table gained a rule
separating the classical rows from the sampling row.

**Float placement was broken and is fixed.** Before this revision, Tables 3–5 and the appendix
tables were deferred by LaTeX to pages 45–49, up to 20 pages after the text discussing them. Added
`placeins` with section barriers, relaxed float fractions, and changed 12 floats from `[t]` to
`[htbp]`. The held-out table now lands 2 pages after its text; the matched-resource table 4.

---

## 7. Material moved to appendix or supplementary

| Material | From | To |
|---|---|---|
| Carry-audit before/after counts and the geometry-bug detail | §2.7 main text | Appendix C (already there; main text shortened to point at it) |
| Script-to-artifact map (Table) | Reproducibility section | New Appendix B |
| Canonical-only ablation table | §3.3 | Merged into Table 3 (§3.4) |
| Full proofs | — | Appendix A (unchanged; statements remain inline) |
| Realised resource distributions, calibration bins, per-cell grids | — | Appendix D (unchanged) |

**Length.** Redundancy was cut (duplicate table, halved reproducibility section, shortened
retraction narrative, merged limitations) but required content was added in the same pass
(novelty table, restructured complexity section, noise limitation, matching-exactness paragraph,
label taxonomy). Net source length is approximately unchanged at ~2,300 lines; the honest summary
is that **redundancy fell and evidential content rose**, not that the manuscript got shorter. No
evidence needed to defend a central claim was removed.

---

## 8. References

- `kitaev1995`, `modularshor2025`, `cuccaro2004` converted from mislabelled `@article`-with-arXiv
  entries to proper `@misc` preprint entries with `eprint`/`archivePrefix`.
- `wqpe2025` now cites the journal version of record (*Advanced Quantum Technologies*, 2026) with
  the preprint identifier in a note, per "prefer the published version".
- `qiskit2024` converted to `@misc` (it was triggering an empty-journal BibTeX warning).
- **Added** `vandervaart1998` (*Asymptotic Statistics*) and `lehmannromano2005` (*Testing
  Statistical Hypotheses*), filling two genuine citation gaps: the Bernstein–von Mises step in the
  consistency proof sketch, and the monotone-likelihood-ratio property of the Beta family used in
  the monotonicity proof. Both are cited at the point of use, in the main text and the appendix.
- BibTeX now runs with zero warnings.

---

## 9. Verification performed

- Every numeric claim recomputed from the frozen CSVs and matched verbatim against `main.tex`:
  per-instance success rates (4 instances × 5 arms), discordant counts, exact-binomial p-values,
  bootstrap intervals, Brier scores, realised resource table, matched-control exactness figures,
  attribution counts. **No discrepancies.**
- The instrumented canonical grid confirmed bit-identical to the published ablation on all 144
  trials in all five arms.
- Compile: 0 errors, 0 undefined references, 0 undefined citations, 0 BibTeX warnings; worst
  overfull box 15.8 pt.
- Terminology sweep: no occurrence of "calibrated probability" except in explicit negation; no
  occurrence of "provably reduces", "confidence-calibrated", or hype vocabulary.

---

## 10. Limitations that remain unresolved because they require new experiments

These are stated in the manuscript and are **not** fixable by writing:

1. **Instance ceiling.** N ∈ {15, 21}, r ∈ {2, 4, 6}, forced by the measured compilation cliff.
   Removing it requires replacing the dense-unitary arithmetic layer with a reversible
   modular-arithmetic circuit.
2. **No noise validation.** Requires re-running the five-arm grid with the noise models already in
   the codebase — compute only, but a new experiment.
3. **No fully matched control for the Full arm.** The existing control matches shots only.
   Constructing a three-way matched control is a direct extension of the released instrumentation.
4. **No instance with concentrated ambiguity.** The prediction that confidence-guided placement
   would separate from uniform placement where ambiguity is concentrated in a minority of windows
   is stated as a prediction, not a result, and no reachable instance exhibits it.
5. **Hyperparameters not cross-validated.** ε and δ_cov were used at defaults.
6. **Recalibration not validated out of sample.** Requires a held-out calibration split.
7. **Simulation only.** No hardware access.
8. **Run-to-run variability uncharacterised.** Two repeats per grid cell; McNemar speaks to one
   evaluation set.
