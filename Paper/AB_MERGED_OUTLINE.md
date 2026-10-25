# Merged Paper A + B: section-by-section outline

**Status:** planning draft, 25 Oct 2026. Nothing has been rewritten yet.
**Target:** one full research article for a Q1 journal. First choice: IEEE Transactions on Quantum Engineering. Alternatives: Physical Review A, then Quantum Information Processing. Check each against the institution's Q1 list before choosing.
**Rules:**
- Do not edit the frozen Paper A (tag `paper-a-v1.0-frozen`). Build the merged paper as a new manuscript on its own branch.
- No new experiments. Every number comes from A or B as they stand.
- No text or figures shared with Paper C. Cite C as an arXiv preprint only.

**Length budget:** about 10,000 words of main text, plus one combined supplement. The sources hold about 22,000 words of main text (A about 12,600, including MethodsX boilerplate; B about 9,400). Roughly half must move to the supplement or be cut.

---

## Working title

*Budget starvation, not logic defects: failure attribution and adaptive classical reconstruction in windowed quantum phase estimation*

Shorter alternative: *Why windowed phase estimation fails, and adaptive reconstruction-budget selection as a remedy*

## Storyline in one paragraph

Windowed QPE saves qubits in Shor's algorithm by shifting work to a classical reconstruction stage. That stage uses fixed budgets.
1. **Diagnose (B Part I).** One-factor stress never breaks the pipeline: 0 failures in 435 trials. Starving the budgets jointly fails 47.9% of trials with no noise. Of the failures, 51% come from shot-budget starvation, 27% from noise and 0% from logic defects.
2. **Treat (A method).** A per-window posterior statistic C_w drives coverage-based retention and confidence-gated shot stopping. With every module switched off, the method reduces exactly to the published pipeline.
3. **Test (B Part II and A validation).** Success rises from 52.1% to 100% on the 144-cell grid, and to 97.9–100% on three held-out instances. At window level, retention keeps 0.841 of the true mass against 0.571.
4. **Be honest about why.** Matched-resource controls and the window-level replay show the gain comes from *choosing* the budget automatically, not from *placing* it better.

Two evidence levels must be kept separate and named explicitly: **trial level** (B's end-to-end grid) and **window level** (A's 3,354 windows checked against exact noiseless distributions).

---

## Front matter (about 400 words)

| Item | Source | Action |
|---|---|---|
| Abstract (≤250 words) | A and B abstracts | New text. Order: failures → attribution → method → end-to-end gain → held-out results → matched controls → window-level mechanism → scope. Keep "noiseless evaluation, N ∈ {15, 21}" in the last line. |
| Highlights / graphical abstract | A front matter | Only if the journal asks. Reuse `graphical_abstract.tex` after redrawing it for the merged story. |
| Specifications table | A (MethodsX only) | **Drop.** Move code/data availability to the back matter. |
| Keywords | both | windowed QPE; Shor's algorithm; classical post-processing; adaptive resource allocation; failure attribution; reproducibility |

## 1. Introduction (about 900 words)

**Sources:** B §1 Introduction, A Introduction and *Contributions*.

- Problem: windowed QPE trades qubits for classical work, and the classical stage has fixed constants nobody stress-tests.
- Gap: there are no controlled failure studies of this stage, and no adaptive alternative tested against matched-resource controls.
- **Contributions** (merge A's and B's lists into one list of four or five items):
  1. A controlled failure characterisation with counterfactual attribution.
  2. C_w and a modular adaptive reconstruction that reduces exactly to the baseline.
  3. A held-out, pre-specified end-to-end evaluation.
  4. Matched-resource controls and a window-level replay that isolate *selection* from *placement*.
  5. A full public artifact.
- Close with a roadmap of the paper.

## 2. Background and related work (about 900 words)

**Sources:** B §2, A "The original method and the three constants it fixes", A "Relation to existing adaptive and Bayesian approaches".

- 2.1 Windowed QPE in modular Shor. Name the three fixed constants (retained candidates k, beam width, shots).
- 2.2 Related work: adaptive and Bayesian phase estimation, beam search and decoding. State precisely what is new. Take this text from A, which is the more careful version.
- **Merge note:** this is where A and B overlap most. Keep a single notation table and define the symbols once.

## 3. Failure characterisation (about 2,000 words) — *from B Part I*

**Sources:** B §§4–10. Condense each subsection to its finding plus one figure or table.

- 3.1 Experimental design (B §4). Instances, grids, seeds, simulator.
- 3.2 One factor at a time never fails: 0 / 435 (B §5). One short paragraph.
- 3.3 Joint budget starvation: 47.9% failure without noise (B §6). **This is the key motivating figure.**
- 3.4 Noise: depolarising and readout effects combine approximately independently (B §7). Put the heatmap here and the full factorial in the supplement.
- 3.5 Direct fault injection into the classical stage (B §8). Condense it to a validity check that the failure detectors work.
- 3.6 Compilation cliff (B §9). Explain why the study is limited to N ∈ {15, 21}. Give this one paragraph and reuse it in Limitations.
- 3.7 **Failure attribution** (B §10): 94 failures, 51% shot starvation, 27% noise, 21% unattributed, 0% logic. This is the bridge into §4: starvation, not logic, causes the failures, so budgets are what should adapt.
- **Move to the supplement:** B §11 "A Self-Caught Audit Error". Keep a two-sentence note in §5 or in Limitations. Its transparency is a strength, but it takes space.

## 4. Adaptive reconstruction method (about 2,000 words) — *from A*

**Sources:** A "Method details".

- 4.1 Why fixed budgets fail. Now one short paragraph that points back to §3.7 (A "Why fixed budgets can cause failures").
- 4.2 The posterior comparison statistic C_w: definition, what it is, why it was chosen over simpler statistics. Numerical evaluation goes to the supplement.
- 4.3 Module B: coverage-based candidate retention.
- 4.4 Module C: ambiguity-scaled beam width. **Documented, not recommended alone.** Keep it brief.
- 4.5 Module D: confidence-gated shot stopping.
- 4.6 The complete procedure. Algorithm box. Invariant: with all modules off, it reduces to the published pipeline.
- 4.7 Theoretical properties. State Propositions 1–7 and Lemmas 1–2 in one compact table. Proofs go to the supplement (already in A's supplement §S1).
- 4.8 Computational cost: asymptotic cost, bounded overhead from the caps, wall-clock time. Shorten to one paragraph and one table.
- **Move to the supplement:** "A retracted module and an unchanged carry check", "Assumptions and parameters" (keep a one-line summary in the main text), "Protocol for reuse" (or turn it into a short "Practical recommendations" box in §7).

## 5. Experimental setup and statistics (about 800 words)

**Sources:** B §§13–15, A "Statistical methods", B §19 (baseline reproduction).

- 5.1 Evaluation protocol: the 144-cell grid, the arms, held-out instances specified in advance.
- 5.2 Baseline reproduction as a validity check (B §19). Place it **before** the results, to show the baseline is faithful.
- 5.3 Statistical methods. **Unify A's and B's methods into one set:** paired tests, effect sizes, confidence intervals, multiple-comparison family, and the p-value convention.
- 5.4 Window-level ground truth (A "Ground truth at the window level"): 3,354 windows against exact noiseless distributions, and the shared corpus.

## 6. Results (about 2,500 words)

Order: trial level first (the headline), then window level (the mechanism).

- 6.1 **Ablation, trial level** (B §16): published pipeline 52.1% → full adaptive 100%. One main figure.
- 6.2 **Held-out generalisation** (B §17): 97.9–100% on three held-out instances, pre-specified.
- 6.3 **Matched-resource controls** (B §18): no module shows a statistically distinguishable advantage once resources are matched. This is the central honest result and needs one clear table.
- 6.4 **Window-level retention** (A "Retention results"): 0.841 against 0.571. Insensitivity to budget: 0.737–0.986 across k, against 0.29–1.00. Benefit = robustness to the setting, not superiority over a well-chosen fixed k.
- 6.5 **Window-level shot stopping** (A "Module D at the window level"): not distinguishable from uniform spending at the shipped cap, and the cap does not explain the null result. The constructed "when would placement matter" probe goes to the supplement.
- 6.6 **Calibration** (A "Calibration"): AUC 0.640 means C_w ranks windows usefully. Brier 0.2239 means it is not a calibrated probability. Bin-level data goes to the supplement.

## 7. Discussion (about 600 words)

**Sources:** B §20, A discussion items 1–3.

- Reconcile the two levels. The end-to-end gain (§6.1–6.2) and the null placement results (§6.3, §6.5) are consistent: the gain comes from **automatic budget selection**.
- What this means for practitioners: a short "Practical recommendations" box condensed from A's "Protocol for reuse". Enable B and D, leave C off, and use the stated caps.
- Relation to adaptive and Bayesian QPE. Mention in one neutral sentence that later work (Paper C, arXiv) studies decoding for AWQPE. No results from C.

## 8. Limitations (about 600 words)

**Sources:** B §21 (1,109 words, so cut it about in half), A discussion items 4–5.

- Simulation only. The evaluation is noiseless. N ∈ {15, 21} because of the compilation cliff.
- Module C is not validated.
- C_w is not calibrated. Assumption (A2), independence, has an exact cost (A "The exact cost of the independence assumption").
- 21% of failures are unattributed.
- The audit error and how it was corrected (two sentences, with a pointer to the supplement).

## 9. Conclusion (about 250 words)

**Sources:** B §22, A conclusion. Three sentences on the findings, one on the honest mechanism, one on future work.

## Back matter

- Data and code availability: the repository plus an archived release (Zenodo DOI). Take this from A's specifications table.
- Author contributions (CRediT), funding, competing interests, acknowledgements: reuse A's and B's back matter, and confirm with the mentor whether they are a co-author.
- References: merge `refs_A.bib` and `refs_B.bib` and remove duplicates. Add Paper C's arXiv entry.

---

## Combined supplement (S1–S8)

| Section | Content | Source |
|---|---|---|
| S1 | Full derivations and proofs (Propositions 1–7, Lemmas 1–2, Invariant 1, decomposition, monotone non-restrictiveness) | A supplement §S1 |
| S2 | Exact joint posterior for the top-two comparison; numerical evaluation of C_w | A supplement §S3, A main text |
| S3 | Failure-study extras: full noise factorial, fault-injection detail, full attribution breakdown | B §§7–8 detail, B supplement §S1.5 |
| S4 | The audit error: before and after, the defect, three independent checks, coverage | B §11, B supplement §S1.1–S1.4 |
| S5 | Retracted module, carry check, full assumptions and parameters, protocol for reuse | A main text |
| S6 | Module D secondary outcomes and the constructed placement probe | A supplement §S4 |
| S7 | Calibration bin-level data; per-instance ablation; realised resources; exactness of the matched-resource controls | A supplement §S2, B supplement §S2 |
| S8 | Artifact map (script → result → section), determinism | A supplement §S5 (extend it to B's scripts) |

## Figures and tables for the main text (suggested 7 figures, 4 tables)

1. Fig. 1. Pipeline schematic with the three fixed constants and where the modules act (new or adapted from the graphical abstract).
2. Fig. 2. Joint starvation failure surface (B §6).
3. Fig. 3. Noise interaction heatmap (`Data/failure_study/plots/noise_combined_heatmap.png`).
4. Fig. 4. Failure attribution (`paper_fig2_failure_attribution.png`).
5. Fig. 5. Ablation and held-out success (B §§16–17).
6. Fig. 6. Retention coverage against k, adaptive vs fixed (A).
7. Fig. 7. C_w calibration / ROC (A).
- Table 1. Notation. Table 2. Theoretical properties summary. Table 3. Matched-resource controls. Table 4. Computational cost.

## Consistency checks before drafting

- [ ] **Module C:** A says "not recommended alone". Check how B's ablation treats it, and use one wording everywhere.
- [ ] **Two evidence bases:** every number must say whether it is trial level (144-cell grid) or window level (3,354 windows).
- [ ] **Statistics:** A and B define their test families separately. Make the p-value convention and multiple-comparison wording match.
- [ ] **Notation:** C_w, k, beam width and shot cap must use the same symbols in both sources.
- [ ] **Numbers:** compare every headline figure against B's verification audit (`PAPER_B_FINAL_VERIFICATION_AUDIT.md`) and A's freeze manifest.
- [ ] **Self-overlap:** if B's arXiv package was ever uploaded, the merged paper must replace it as a new version, not appear as a duplicate. Current evidence says it was never uploaded; confirm.
- [ ] **Journal template:** IEEE TQE uses the IEEE two-column class, so check its page and word limits before cutting.

## Suggested schedule (after Paper C is submitted, about 12–14 Nov)

| Days | Work |
|---|---|
| 1–2 | New branch, skeleton `AB_master.tex`, merge the bib files, notation table |
| 3–5 | §§3–4 (condense B Part I; port A's method) |
| 6–8 | §§5–6 (unify statistics; results) |
| 9–10 | §§1, 2, 7–9, abstract; build the supplement |
| 11–12 | Consistency checks, number audit, mentor read-through |
| 13 | arXiv post (about 24–27 Nov), then journal submission by 30 Nov |
