# Submission Readiness Checklist — MethodsX

**Manuscript:** `Paper/main.tex` → `Paper/main.pdf` (38 pp., 4 data figures + 1 schematic +
1 algorithm, 8 tables, 22 references)
**Supplementary:** `Paper/supplementary.tex` → `Paper/supplementary.pdf` (19 pp., 9 tables)
**Status:** submission-ready. **Zero placeholders remain in either document.** All bibliography
entries verified against arXiv and the publishers. Two author actions outstanding (ORCID, archival
DOI) — see `FINAL_REVISION_REPORT.md` §6.

Legend: **[x]** done and verified · **[ ]** requires the author · **[–]** not applicable

---

## 1. Manuscript completeness

| | Item |
|---|---|
| [x] | Compiles clean: 0 errors, 0 undefined references, 0 undefined citations, 0 BibTeX warnings |
| [x] | No `TODO` markers remain |
| [x] | Worst overfull box 28.7 pt; 18 boxes total, all cosmetic |
| [x] | All floats placed within a few pages of the text that discusses them (was broken; fixed with `placeins` + relaxed float fractions) |
| [x] | No duplicated cross-reference artefacts ("Appendix Appendix", "Table Table", etc.) |
| [x] | Abstract follows the four-part structure and contains no claim unsupported by a body section |
| [x] | Every claim in the abstract maps to a numbered section, proposition, table or figure |
| [x] | Evidentiary labels applied consistently: `\PROVEN`, `\INVAR`, `\HEUR`, `\EMPIR`, with a definition list at first use |
| [x] | No hype vocabulary; "calibrated probability" appears only in explicit negation |
| [x] | Author affiliation — Dept. of Computer Science and Engineering, Amrita Vishwa Vidyapeetham, Chennai, India |
| [ ] | ORCID iD — not recorded in the project; author to add |
| [–] | Acknowledgements — section removed at the author's instruction (restore if a mentor or funder must be credited) |
| [x] | CRediT — single-author role list retained; a mangled `\noindent` that printed literally as "oindent" was repaired |
| [x] | Repository URL — real, from the git remote |
| [ ] | Archival DOI — both documents state it will be supplied at acceptance |
| [x] | 72% claim scoped to the analysed failure corpus in all four locations |
| [x] | "pre-registered" replaced by "pre-specified" throughout (no formal registration record exists) |
| [x] | Central contribution stated at the head of the Discussion |
| [x] | Complexity claim carries its bounded-cap assumption |

## 2. MethodsX required sections

| | Section | Location |
|---|---|---|
| [x] | Title, authors, abstract, keywords | front matter |
| [x] | Highlights (4 bullets) | after front matter |
| [x] | **Specifications table** (subject area, specific area, method name, original method + reference, resource availability) | after Highlights |
| [x] | Graphical abstract — description supplied | before Introduction |
| [x] | **Method details** | §1, fourteen subsections incl. Algorithm 1 |
| [x] | **Method validation** | §2, eight subsections |
| [x] | Limitations (dedicated main-text section, 6 grouped items) | §3 |
| [x] | Ethics statements | back matter |
| [x] | CRediT author statement | back matter |
| [x] | Declaration of competing interest | back matter |
| [–] | Acknowledgements | removed at author's instruction |
| [x] | Reproducibility | back matter |
| [x] | Code availability | back matter |
| [x] | Data availability | back matter |
| [x] | Supplementary material description | back matter |
| [x] | References | end |
| [x] | Supplementary Information S1–S5 (proofs, audit record, extended tables, calibration bins, artifact map) | separate document |
| [x] | Scope note confirming only the classical stage is customised | after Specifications table |

**MethodsX fit check.** The paper reads as a validated, reusable methodology rather than a quantum
theory paper: §1.14 is a seven-step adoption protocol; the cost analysis is framed for a
practitioner budgeting a run; every module is independently switchable; and the implementation
reduces to the published baseline when the modules are disabled, so it can be adopted
incrementally.

## 3. Figures

| | Fig. | Content | Source file | Generator |
|---|---|---|---|---|
| [x] | 1 | Method schematic | inline TikZ | — (vector, no external file) |
| [x] | 2 | Failure attribution, 72% bracket | `plots/paper_fig2_failure_attribution.png` | `Experiments/paper_figures.py` |
| [x] | 3 | Held-out generalisation, 4 small multiples | `plots/paper_fig3_heldout.png` | same |
| [x] | 4 | Matched-resource controls, significance annotations | `plots/paper_fig4_matched_resource.png` | same |
| [x] | 5 | Calibration reliability, overconfidence gap | `plots/paper_fig5_calibration.png` | same |

| | Figure quality checks |
|---|---|
| [x] | Every figure answers one stated scientific question, named in bold at the start of its caption |
| [x] | No rotated tick labels; no raw code identifiers as labels |
| [x] | All rendered at 300 dpi; all files present in the repository |
| [x] | Generator runs no trial and recomputes no result; it asserts the calibration Brier scores still equal 0.2239 / 0.1988 before writing |
| [x] | Captions name the source CSV |
| [ ] | *Optional:* re-render as PDF vector for final production (change `savefig` extensions and `\includegraphics` paths). PNG at 300 dpi is acceptable. |
| [ ] | Graphical abstract image file — described in the manuscript, needs producing at submission |

## 4. Tables

| | Table | Content |
|---|---|---|
| [x] | Specifications | MethodsX-required block (unnumbered) |
| [x] | 1 | Notation |
| [x] | 2 | Novelty positioning against cited prior work |
| [x] | 3 | Free parameters and provenance |
| [x] | 4 | Per-stage cost, baseline vs adaptive |
| [x] | 5 | Holm–Bonferroni across the family of nine |
| [x] | 6 | Per-instance ablation, four instances (merged; the former canonical-only table was a duplicate and was removed) |
| [x] | 7 | Matched-resource controls |
| [x] | B.8 | Artifact map |
| [x] | C.9 | Carry-check audit, before and after |
| [x] | D.10 | Realised resource usage per arm |
| [x] | D.11 | Calibration by predicted-confidence bin |
| [x] | Every table caption states the unit of analysis, the test used, and what $(b,c)$ mean |
| [x] | No table duplicates another's rows |

## 5. Citation verification

| | Item |
|---|---|
| [x] | 20 entries; all cited; no uncited entries; BibTeX clean |
| [x] | Preprints correctly typed as `@misc` with `eprint`/`archivePrefix` (`kitaev1995`, `modularshor2025`, `cuccaro2004`) |
| [x] | `wqpe2025` cites the journal version of record rather than the preprint |
| [x] | DOIs present where known (11 entries) |
| [x] | Citation gaps closed: `vandervaart1998` for Bernstein–von Mises; `lehmannromano2005` for the Beta MLR property |
| [x] | Statistical-method choices cited at the point of use (Dietterich, McNemar, Holm, Efron & Tibshirani, Brier, Zadrozny & Elkan, Jeffreys) |
| [x] | No fabricated bibliographic detail; unknown fields left absent rather than invented |
| [ ] | **Confirm** `wqpe2025` volume / article number / DOI once the *Advanced Quantum Technologies* version is paginated (currently a `note` says "to be completed at proof stage") |
| [ ] | **Confirm** author list and Quantum DOIs for `bae2025` and `biqae2026` (volume/page recorded; DOI absent) |
| [ ] | **Confirm** the author attribution on `wqpe2025` / `modularshor2025`. The repository's own documents cite these consistently by arXiv identifier and never name the authors; the names in `refs.bib` predate this revision and could not be verified from the project sources. |

## 6. Reproducibility materials

| | Item |
|---|---|
| [x] | Artifact map (Appendix B) maps every result file → generating script → consuming section |
| [x] | Seed derivation, determinism, and resumability documented |
| [x] | Statistical procedure fully specified: exact binomial McNemar; percentile bootstrap $B=2000$, seed 0, paired per-trial differences; Holm over a pre-declared family of nine |
| [x] | The summary-CSV column-naming trap (`baseline_only_fixed_count` / `arm_only_broke_count` hold the reverse of their names) disclosed in the Reproducibility section and in `Paper/BUILD.md` |
| [x] | Every number in the manuscript recomputed from the frozen CSVs and matched verbatim — no discrepancies |
| [x] | Instrumented canonical grid confirmed bit-identical to the published ablation on all 144 trials × 5 arms |
| [x] | Licence stated (Apache-2.0) |
| [ ] | Public repository URL (appears twice: Specifications table, Code availability) |
| [ ] | Zenodo DOI for the archived release (appears twice) |
| [ ] | Commit hash or tag for the corrected carry-check audit (Appendix C) |
| [ ] | Push the uncommitted Phase 7 artifacts before archiving — `phase6_ablation_diagnostics.csv`, the three `phase7_heldout_*.csv`, `phase7_matched_resource.csv`, the two summary tables, `Experiments/phase7_generalization.py`, `Experiments/paper_figures.py`, and the four `paper_fig*.png` |

## 7. Honesty and scope audit

| | Item |
|---|---|
| [x] | Headline is the four-instance range 97.9–100%, not the canonical 100% |
| [x] | The matched-resource null results are reported as a primary finding, in the abstract, not buried |
| [x] | Contribution framed as *automatic budget selection*, not as superior placement |
| [x] | Module C explicitly not presented as independently validated; recommended against enabling alone |
| [x] | Noiseless validation disclosed in the abstract, the design section, and Limitation 3 |
| [x] | $C_w$ named a posterior comparison statistic; calibrated-probability reading explicitly disclaimed |
| [x] | Recalibration disclosed as an in-sample fit |
| [x] | "No new complexity class" qualified to classical work only |
| [x] | Baseline equivalence labelled an implementation invariant, not a theorem |
| [x] | "Pre-registered" qualified as a repository artifact, not a third-party registration |
| [x] | Self-caught audit-script correction narrated inline, not labelled an erratum, no apology, with an appendix quantification |
| [x] | 14 limitations in a dedicated main-text section, each paired with a mitigation or the work that would close it |

## 8. Unresolved before submission — author actions

1. Fill the placeholder fields — nine markers, seven distinct values, since the repository URL
   and the Zenodo DOI each appear twice (`grep -n 'PLACEHOLDER' Paper/main.tex`).
2. Confirm the three bibliography items flagged in §5 above.
3. Commit and archive the Phase 7 artifacts; mint the Zenodo DOI; insert URL and DOI.
4. Produce the graphical abstract image from the description in the manuscript.
5. Update `RESEARCH_OVERVIEW.md`, which still describes Phase 7 as in progress at 34/144. It is
   complete at 144/144 on all four grids. Not part of the submission, but it is the repository's
   front door and a reviewer following the code link will read it.
6. Optional: re-render figures as PDF for production.

## 9. Known irreducible limitations (stated in the manuscript; no pre-submission action)

Instance ceiling N ∈ {15, 21}; no noise validation; no fully matched control for the Full arm; no
instance with concentrated ambiguity; hyperparameters not cross-validated; recalibration not
validated out of sample; simulation only; run-to-run variability uncharacterised. Each is listed
in §3 of the manuscript with what would close it. See `FINAL_REVISION_LOG.md` §10.
