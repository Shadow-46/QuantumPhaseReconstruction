# Submission Readiness Checklist

**Paper A** — *Confidence-Guided Adaptive Reconstruction for Windowed Quantum
Phase Estimation*. `Paper/A_method.tex` (28 pp. preprint / 14 pp. two-column)
plus `Paper/A_supplementary.tex` (11 / 6 pp.). Target: **MethodsX**.

**Paper B** — *Where Windowed Shor's Algorithm Actually Fails, and What Fixing
It Buys: A Controlled Study*. `Paper/B_master.tex` + the `Paper/B/` tree
(30 pp.) plus `Paper/B_supplementary.tex` (11 pp.). Target: **not yet
selected**.

Last verified: 2026-08-27.

---

## 1. Build state

| | Paper A | Paper A SI | Paper B | Paper B SI |
|---|---|---|---|---|
| LaTeX errors | 0 | 0 | 0 | 0 |
| Undefined references | 0 | 0 | 0 | 0 |
| Undefined citations | 0 | 0 | 0 | 0 |
| Overfull boxes | 12 | 10 | 6 | 2 |

Anything other than 0/0/0 is a regression. Build instructions: `Paper/BUILD.md`.

---

## 2. MethodsX compliance — Paper A

| | Item | Requirement | Actual |
|---|---|---|---|
| [x] | Abstract | ≤250 words | 232 |
| [x] | Introduction | ≤500 words | 439 |
| [x] | References | ≤25 | 16 cited / 18 in file |
| [x] | Highlights | 3–5 bullets, ≤85 chars each | 5 bullets, 70–71 chars |
| [x] | Specifications table | required | present |
| [x] | Graphical abstract | ≥1328×531 px @ 300 dpi | 2413×809 px |
| [x] | Method details | required | present |
| [x] | Method validation | required | present |
| [x] | Limitations | required | present |

---

## 3. Required statements

| | Statement | Paper A | Paper B |
|---|---|---|---|
| [x] | Ethics | present | present |
| [x] | CRediT author statement | present | present |
| [x] | Declaration of competing interest | present | present |
| [x] | Funding | present (nil declaration) | present (nil declaration) |
| [x] | Data availability | present | present |
| [x] | Code availability | present | present |
| [x] | Reproducibility | present | present |
| [x] | Supplementary material | present | present |
| [x] | Relation to the companion article | present | present |
| [x] | ORCID iD | 0009-0006-9042-254X | 0009-0006-9042-254X |

---

## 4. Claim discipline

| | Check |
|---|---|
| [x] | No `\PLACE{}` or `\TODO{}` is invoked in any document |
| [x] | "pre-registered" does not appear — the term implies a third-party registration that does not exist; all sites read "pre-specified" |
| [x] | The matched-resource null is stated in Paper A's abstract and Limitations, not only in the companion |
| [x] | The headline is the 97.9–100% range across four instances, not a flat 100% |
| [x] | Baseline equivalence is stated as an implementation invariant, not a theorem |
| [x] | C_w is a posterior comparison statistic, explicitly not a calibrated probability; AUC 0.640 reported |
| [x] | Isotonic recalibration is labelled in-sample and is not wired into the shipped rule |
| [x] | Noiseless, simulation-only scope disclosed in both papers |
| [x] | The N ∈ {15,21} ceiling is attributed to a measured compilation cliff, not presented as a choice |
| [x] | Module C is not claimed as an independently validated intervention |
| [x] | Paper B's Conclusion restates only claims made verbatim elsewhere |

---

## 5. Contribution boundary

| | Check |
|---|---|
| [x] | No image file is shared between the papers |
| [x] | No figure or table is shared |
| [x] | Artifact-map rows are disjoint (A: 3 data files, B: 12) |
| [x] | Each paper has its own primary outcome variable — window-level vs. end-to-end |
| [x] | Both papers disclose that they analyse one released corpus at different units |
| [x] | Cross-references are citations, never duplicated content |

---

## 6. Reproducibility

| | Check |
|---|---|
| [x] | `requirements.txt` pins the exact versions that produced the results |
| [x] | `validate_environment.py` floors match those pins |
| [x] | `README.md` names every script needed to reproduce both papers |
| [x] | Working tree clean; all Phase 7 artifacts committed |
| [x] | `RESEARCH_OVERVIEW.md` numbers match the manuscripts |
| [x] | Every bibliography entry verified against Crossref or the publisher record — see §9 |
| [ ] | **Zenodo DOI** — not minted; both papers state it will be supplied at acceptance |
| [ ] | **Commit hash for the corrected carry-check audit** — not yet recorded |

---

## 7. Open before submission

| Item | Owner | Blocking? |
|---|---|---|
| Select Paper B's venue, then fork a profile in `Paper/B/config/` | author | yes for B |
| Mint the Zenodo DOI and insert it | author | at acceptance |
| Flatten relative figure paths for publisher upload | either | at upload |
| `widefigure`/`widetable` are defined but unused in Paper B — a two-column target would shrink wide tables rather than widen them | either | before any two-column conversion |

`IEEEtran.cls` and `ieeeaccess.cls` are **not installed** in this MiKTeX tree.
Install before attempting an IEEE target.

---

## 8. Known irreducible limitations

These are disclosed in both manuscripts and are not defects to fix before
submission.

- All validation is noiseless and simulation-only.
- The corpus is confined to N ∈ {15, 21} by a measured compilation cliff:
  0.14 s at 4 work qubits, 42.2 s at 7, timeout at 8.
- No control matches all three budgets simultaneously against the Full arm.
- C_w compares only the top two patterns, so a genuine three-way near-tie is
  not flagged.
- Module B carries no population coverage guarantee; the measured shortfall is
  reported.

---

## 9. Bibliography verification

The three entries previously flagged as unconfirmed were checked against
Crossref and the publisher records on 2026-08-27. All resolve, and all
author lists, volumes, article numbers, years and DOIs match what the `.bib`
files state.

| Entry | Verified against | Result |
|---|---|---|
| `wqpe2025` | Crossref `10.1002/qute.202500683` | Shukla & Vedula, *Adv. Quantum Technol.* **9**(3), e00683, 2026 — matches. **Title corrected** from "Towards" (the arXiv form) to "Toward" (the published form). |
| `modularshor2025` | Crossref `10.1140/epjp/s13360-026-07719-0` | Shukla & Vedula, *Eur. Phys. J. Plus* **141**, 474, 2026 — matches. An author page still lists it "to appear"; the DOI is registered and the record complete. |
| `bae2025` | quantum-journal.org, `10.22331/q-2025-09-11-1856` | Ramôa & Santos, *Quantum* **9**, 1856, 2025 — matches. The key name is unrelated to the authors, which is harmless. |
| `biqae2026` | quantum-journal.org, `10.22331/q-2026-01-14-1962` | Li, Vidwans, Wang & Soley, *Quantum* **10**, 1962, 2026 — matches. |

Attribution on both Shukla & Vedula entries is confirmed correct.
