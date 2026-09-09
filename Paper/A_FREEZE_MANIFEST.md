# Paper A — Final Scientific Freeze Manifest

## 1. Paper

**Confidence-Guided Adaptive Reconstruction for Windowed Quantum Phase Estimation**
Sanjay Lakshmanan R — Department of Computer Science and Engineering, Amrita Vishwa
Vidyapeetham, Chennai, India. Target venue: MethodsX.

## 2. Freeze date

2026-09-09

## 3. Git identification

| Item | Value |
|---|---|
| Branch | `paper-split` |
| Parent commit (state this freeze builds on) | `041a71c433df91a1e90cee11d44aa885665a0776` |
| Freeze commit | identified by the annotated tag below |
| Annotated tag | `paper-a-v1.0-frozen` |

Resolve the freeze commit with `git rev-parse paper-a-v1.0-frozen`. The manifest cannot
contain its own commit hash, since it is part of the commit it describes.

## 4. Main manuscript source

`Paper/A_method.tex`

## 5. Supplementary source

`Paper/A_supplementary.tex`

## 6. Final manuscript PDF names

- `Paper/A_method.pdf` — two-column submission layout, **17 pages**
- Preprint/single-column layout of the same source — **35 pages**

Both are build outputs and are gitignored (`.gitignore:13`, `Paper/*.pdf`). The layout is
selected by the documentclass toggle at `A_method.tex:27-28`; the frozen source is left in
the **two-column submission configuration** (`final,5p,times,longtitle`).

## 7. Final SI PDF name

`Paper/A_supplementary.pdf` — two-column **8 pages**, preprint **15 pages**. Its toggle is
at `A_supplementary.tex:12-13` and must be kept in step with the main document.

## 8. Required figures

| Path | Used in |
|---|---|
| `Paper/graphical_abstract.png` | graphical abstract |
| `Data/failure_study/plots/paper_figA_window_retention.png` | main, Fig. 2 |
| `Data/failure_study/plots/paper_fig5_calibration_raw.png` | main, Fig. 3 (raw reliability panel) |
| `Data/failure_study/plots/paper_fig5_calibration.png` | SI, Fig. S1 (raw + isotonic, in-sample) |

## 9. Required datasets and artifacts

| Path | Role |
|---|---|
| `Data/failure_study/raw_windows/*.json` (579 files) | primary per-window counts; 3,354 windows |
| `Data/failure_study/paperA_window_validation.csv` | Table 5, Fig. 2 |
| `Data/failure_study/phase6_calibration.csv` | calibration study, 2,687 observations |
| `Data/failure_study/phase6_calibration_dirichlet.csv` | Dirichlet diagnostic |
| `Data/failure_study/paperA_moduleD_window_validation.csv.gz` | module D replay, 201,240 rows, 20 replicates, caps {32,64,256} |
| `Data/failure_study/paperA_synthetic_grid_A.csv.gz` | Family A, 39 realisable configurations |
| `Data/failure_study/paperA_synthetic_grid_B.csv.gz` | Family B, constructed ensembles |
| `Data/failure_study/summary_tables/paperA_window_validation_summary.csv` | Table 5 summary |
| `Data/failure_study/summary_tables/paperA_dirichlet_comparison.csv` | Table S3 |
| `Data/failure_study/summary_tables/paperA_dirichlet_reliability.csv` | SI reliability bins |
| `Data/failure_study/summary_tables/paperA_moduleD_summary.csv` | Table 6 |
| `Data/failure_study/summary_tables/paperA_moduleD_tests.csv` | Table 6 confirmatory tests |
| `Data/failure_study/summary_tables/paperA_synthetic_A_summary.csv` | Family A summary |
| `Data/failure_study/summary_tables/paperA_synthetic_B_summary.csv` | Table S4 |
| `Data/failure_study/phase7_matched_resource.csv` | shot-cost accounting (9.3 to 53.0) |

## 10. Reproducibility scripts

| Path | Role |
|---|---|
| `Reconstruction/confidence.py` | shipped `C_w`; additive Dirichlet diagnostic |
| `Reconstruction/candidate_generation.py` | module B retention |
| `Algorithms/adaptive_reconstruction.py` | modules C and D, switchable path, Invariant 1 |
| `Experiments/recapture_window_counts.py` | produces `raw_windows/` |
| `Experiments/paperA_window_validation.py` | Table 5, Fig. 2 |
| `Experiments/phase6_calibration.py` | calibration corpus, exact marginals |
| `Experiments/paperA_dirichlet_recomputation.py` | Dirichlet diagnostic |
| `Experiments/replay.py` | circuit-free replay primitives, CRN streams |
| `Experiments/paperA_moduleD_validation.py` | module D window-level study, cap sweep |
| `Experiments/paperA_synthetic_grid.py` | Family A and Family B |
| `validate_environment.py`, `requirements.txt` | environment gate, pinned versions |

Generated but referenced by no manuscript, retained as analysis outputs:
`plots/paperA_figA_moduleD_window.png`, `plots/paperA_figS_dirichlet.png`,
`plots/paperA_figS_synthetic.png`.

## 11. Environment validation result

`python validate_environment.py` gives **22/22 modules pass**, `environment validation
passed`, with the seven pinned dependency versions verified (qiskit 2.4.1, qiskit-aer
0.17.2, numpy 2.4.6, scipy 1.17.1, matplotlib 3.11.0, pandas 3.0.3, scikit-learn 1.9.0).

## 12. Compilation result

| Build | Pages | Errors | Undefined refs | Undefined citations | BibTeX |
|---|---|---|---|---|---|
| `A_method` two-column | 17 | 0 | 0 | 0 | 0 |
| `A_supplementary` two-column | 8 | 0 | 0 | 0 | 0 |
| `A_method` preprint | 35 | 0 | 0 | 0 | 0 |
| `A_supplementary` preprint | 15 | 0 | 0 | 0 | 0 |

All figures resolve. Table 6 (`tab:dsweep`) and Table S3 (`tab:dirichlet`) both render as
full-width `widetable` floats. No table, figure, or page carries clipped or overprinted
content in either layout.

## 13. Scientific status

**Scientific content frozen; subsequent changes require explicit post-freeze approval.**

The frozen position, for the avoidance of doubt:

- Module B is the primary adaptive reconstruction mechanism.
- Module C is optional and documented, **not** a default or recommended module, and is not
  independently validated.
- Module D is an adaptive recapture/allocation mechanism. It is **not** claimed to
  outperform an equally-resourced uniform allocation; at the shipped cap it is not
  distinguishable from its matched control (+1.12 pp, Holm-adjusted p = 0.075).
- The shipped `C_w` is unchanged. The Dirichlet-joint treatment is a diagnostic and a
  pre-registered future change, **not** a replacement.
- Family A covers the 39 corpus-represented, realisable configurations. Family B is a
  constructed mechanism-characterisation probe, **not** physical or end-to-end validation.
- Attribution: 48/94 (51%) positively attributed to reconstruction-stage shot-budget
  starvation; 25/94 (27%) to noise; 20/94 (21%) unattributed; none to a defect in the
  reconstruction logic. No causal claim is made for the 20. The branching structure of the
  attribution procedure is described as implemented.
- Table 5: k=1 argmax-kept = 0.561; N=21, a=19 argmax-coverage = 0.932.

## 14. Known non-blocking issues

1. 168.0 pt overfull box on page 1 of each two-column document — elsarticle's title-block
   output routine; structural to the class, not content.
2. Notation table (SI): 27.7 pt margin overflow. Renders in a clean separate column; no
   overprint.
3. Main text: 42.8 pt and 17.3 pt prose margin overhangs (lines 1146-1151, 1168-1172).
4. SI: 6.58 pt and 1.343 pt margin overhangs.
5. Cosmetic short line-breaks at `A_method.tex:681` and `:1610`.
6. Three generated figures referenced by no manuscript (listed in section 10).
7. Stale, unreferenced `plots/calibration_reliability*.png` retained, not deleted.
8. Two uncited entries in `refs_A.bib` (`lehmannromano2005`, `vandervaart1998`); harmless
   under the numeric bibliography style.
9. The companion reference `@misc{companion}` carries no persistent identifier. It points
   at the public code archive. Once Paper B is posted, its `eprint`/`archivePrefix` must be
   added — the single intended post-freeze edit to Paper A.

## 15. Freeze-operation integrity

**No scientific dataset or experiment was modified during the freeze operation.** No
experiment was rerun, no dataset regenerated, no figure regenerated, and no analysis code
altered. The only manuscript change made during Phase 1 was the authorized two-line layout
correction converting Table S3 from a single-column `table` float to the project's existing
`widetable` environment, which fixed a rendering corruption in which the table's right-hand
column overprinted adjacent body text in the two-column build. No table value, caption, or
scientific statement changed.

The shipped `C_w` implementation was verified byte-identical to its pre-freeze state by AST
comparison (unparsed-function SHA-256 begins `761e45db0eb591440d9b57ec`).

Per-file SHA-256 digests for the frozen set are in `Paper/A_FILE_HASHES.sha256`; verify with
`sha256sum -c Paper/A_FILE_HASHES.sha256` from the repository root.

## 16. Scope note

This freeze covers **Paper A only**. Paper B sources carry uncommitted working-tree changes
from the authorized cross-paper corrections; they are deliberately excluded from the freeze
commit and left untouched. `Paper/build/` PDFs are untracked build outputs and are also
excluded.
