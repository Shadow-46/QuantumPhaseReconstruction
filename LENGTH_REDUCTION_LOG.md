# Length Reduction Log

**Manuscript:** *Confidence-Guided Adaptive Reconstruction for Windowed Quantum Phase Estimation*
**Venue:** MethodsX (Elsevier)
**Operation:** structural compression — main article split into a concise manuscript plus a
Supplementary Information document.

---

## 1. Page counts

| Document | Before | After |
|---|---|---|
| **Main manuscript** (`Paper/main.pdf`) | **49 pages** | **37 pages** |
| **Supplementary Information** (`Paper/supplementary.pdf`) | did not exist | **19 pages** |
| Combined record | 49 | 56 |

**Main manuscript reduced by 12 pages (−24%).** The combined record is *longer* than the original,
which is the intended outcome: material was relocated and, where relocation left a gap, expanded
so the SI is independently verifiable rather than a pile of orphaned tables.

### On the 30–35 page target

The main article lands at **37 pages, two above the stated target.** I stopped there deliberately.
After the final pass the document contains only two pages that are not densely set — the abstract
page and the last page of references — so the remaining two pages cannot be recovered by
typesetting. Reaching 35 would require deleting roughly two pages of scientific content, and the
brief forbade that. The candidates I considered and rejected:

- **Figure 2 (failure attribution)** — Step 8A requires the attribution argument to stay clearly
  explained, and the figure is its clearest form.
- **Protocol for reuse (§1.13)** — the MethodsX-defining section; removing it would damage exactly
  the criterion the venue scores.
- **Matched-resource narrative (§2.5)** — Step 8C explicitly says not to over-compress it.
- **Limitations** — already consolidated 14 items → 6 groups; further cuts would start hiding
  limitations.

If a hard 35-page ceiling is imposed at submission, the least damaging cut is the adoption table
(Table 8) plus the baseline-reproduction subsection (§2.7), both of which could move to the SI
with a one-sentence pointer. I did not make that cut unasked.

---

## 2. Material moved to Supplementary Information

| Material | From | To | Main-article replacement |
|---|---|---|---|
| Complete proofs of all 10 results | Appendix A (~6 pp) | **SI §S1** | Statement list + Table 7 (status, and what each result does *not* establish) |
| Notation table (20 symbols) | Table 1 | **SI §S1.1** | Pointer; every symbol also defined at first use |
| Carry-audit before/after table, defect forensics, coverage check, three verification checks | Appendix C | **SI §S2** | The full scientific story kept in §1.9 (~1 page), quantitative detail deferred |
| Failure-attribution category table | (was figure-only) | **SI §S2.5** | Figure 2 retained in main |
| Complete per-instance ablation (4 instances × 5 arms, CIs, discordant counts, p-values) | Table 6 (~1.5 pp) | **SI §S3.1** | Compact 4-row baseline-vs-Full table |
| Realised resource usage per arm (mean/max of m_w, P′, S_w) | Appendix D.1 | **SI §S3.2** | Key figures quoted in prose in §1.12 |
| Matched-control exactness residuals | §2.5 prose | **SI §S3.3** | One-sentence summary retained |
| Holm–Bonferroni correction table (9 rows) | Table 5 | **SI §S3.4** | The decisive outcome (module C fails correction) stated in prose in three places |
| Full statistical procedure (McNemar variant justification, bootstrap spec, family definition) | §2.2 (~1 p) | **SI §S3.5** | Condensed to a single paragraph naming every choice |
| Calibration bin-by-bin table + ground-truth construction + mechanism + recalibration rationale | Appendix D.2 | **SI §S4** | Figure 5, Brier scores, and the ε-is-heuristic conclusion kept in main |
| Artifact map (file → script → section) | Appendix B | **SI §S5.1** | One-line pointer from the Specifications table and Reproducibility |
| Determinism/resumability detail, per-cell grid description, summary-CSV naming trap | Reproducibility (~1 p) | **SI §S5.2–S5.3** | Three short paragraphs retained |

---

## 3. Sections compressed in place

| Section | Before | After | What changed |
|---|---|---|---|
| Abstract | 453 words | **345 words** | Kept problem, method, canonical result, held-out summary, matched-resource qualification, calibration limitation. Removed individual p-values, all confidence intervals, secondary percentages. |
| Graphical abstract | 21 lines of prose | 10 lines | Long textual description of an unbuilt figure replaced by a build specification. |
| Introduction | 94 lines | 78 lines | Contributions list 6 items → 5; two paragraphs tightened. |
| Related work | 104 lines | 91 lines | Prose tightened; novelty table (Table 2) kept intact. |
| Posterior comparison statistic | 72 lines | 62 lines | Two equations merged into one; numerical-evaluation paragraph condensed. |
| Retracted module / audit | 41 lines | 37 lines | Story kept in full, forensic detail moved to SI §S2. |
| Assumptions + Parameters | two subsections | one (§1.10) | Merged; (A1)–(A5) and Table 3 unchanged. |
| Theoretical properties | 137 lines | 97 lines | Nine proposition environments → two load-bearing statements + Table 7. |
| Computational cost | 99 lines | 88 lines | Three-part structure (asymptotics / caps / wall-clock) retained; prose tightened. |
| Protocol for reuse | 51 lines | 49 lines | All seven steps retained. |
| Design | 96 lines | 75 lines | Three stages + noiseless-scope paragraph retained. |
| Statistical methods | 55 lines | 38 lines | Choices and their justifications retained; mechanics to SI. |
| Ablation results | 84 lines | 55 lines | Forward-pointing paragraph that duplicated §2.4/§2.5 removed. |
| Held-out generalisation | 43 lines of findings | 31 lines | Five bold findings merged into three; every number retained. |
| Matched-resource controls | — | −5 lines | Tightened only; Step 8C's warning respected. |
| Calibration | 80 lines | 71 lines | Figure, Brier scores, in-sample caveat and ε-conclusion all retained. |
| Discussion | 54 lines | 68 lines | Per-module recommendation prose → Table 8 (longer in source, clearer to read). |
| Limitations | 155 lines, 14 items | 117 lines, **6 groups** | Grouped per the requested taxonomy. **No limitation was dropped**; see §5. |

---

## 4. Structural additions

Two elements were **added**, not removed, because the target structure required them:

1. **Algorithm 1** (§1.8) — pseudocode for the complete adaptive reconstruction procedure, with the
   three module flags and the lines marked [B], [C], [D] that depart from the published pipeline.
   Transcribed from the implementation; it introduces no behaviour not already described.
2. **Table 8** (Discussion) — adoption guidance per module, replacing three paragraphs of prose.

Page-efficiency changes (formatting only, no content affected): class option `12pt` → `11pt`;
`\bibsep` tightened, which reduced the bibliography from 4 pages to 2; float and display skips
tightened; figures set at `0.94\linewidth`.

---

## 5. Confirmation: nothing scientific was changed or lost

Verified programmatically after restructuring — every value recomputed from the frozen CSVs and
matched against `main.tex` **or** `supplementary.tex`:

- Six result grids, 144 trials each; canonical instrumented run still bit-identical to the
  published ablation on all 144 trials × 5 arms.
- All 20 per-instance success rates and all 16 discordant-count pairs.
- All four matched-resource comparisons.
- Brier scores 0.2239 / 0.1988, all six calibration bins, n = 2,687 windows / 579 trials.
- All five realised-resource rows.
- Attribution: 94 failures, 72% classical budget starvation.
- 22 required structural elements in the main article; 14 in the SI.

**Result: no discrepancies.** No numerical result, statistical test, effect size, confidence
interval, p-value or claim was altered. No limitation was removed — the 14 previous items map onto
the 6 groups as follows:

| New group | Absorbs previous items |
|---|---|
| 1. External validity and evaluation scope | instance ceiling; noiseless validation; simulation-only |
| 2. Attribution of gain to allocation vs budget | matched-control nulls; no three-way matched control; approximate matching; no concentrated-ambiguity instance |
| 3. The confidence model | (A2) false; ε not a calibrated rate; in-sample recalibration; top-2 only (A3) |
| 4. Module-specific evidence | module C unsupported alone; module B no population guarantee; per-trial regressions |
| 5. Parameters and statistical scope | hyperparameters not cross-validated; design/test circularity; McNemar single-set caveat |
| 6. Deliberate scope exclusions | static overlap geometry (A5); no cross-window allocation (A4) |

The four conclusions the brief flagged as over-repeated now appear once each at full strength —
abstract (one sentence), method (definition only), results (full evidence), discussion
(interpretation), limitations (one concise statement) — with the intermediate restatements removed.

---

## 6. Build status

| | Main | SI |
|---|---|---|
| LaTeX errors | 0 | 0 |
| Undefined references / citations | 0 | 0 |
| BibTeX warnings | 0 | 0 |
| Overfull boxes | 18 (worst 20.8 pt) | 10 (worst 14.7 pt) |
| Pages | 37 | 19 |

Both documents build with `pdflatex → bibtex → pdflatex ×2` and share `refs.bib`. The SI is
standalone: it carries its own front matter, S-prefixed numbering, table of contents, and
bibliography, and references the main article by section name rather than by `\ref` (which cannot
cross documents).

## 7. Open items carried forward

Unchanged from `SUBMISSION_READINESS_CHECKLIST.md`: nine placeholder markers covering seven
author-supplied fields (now present in both documents — the SI repeats the affiliation and the
audit commit hash), and three bibliography confirmations. The graphical abstract still needs
producing as an image file from the specification in the main article.
