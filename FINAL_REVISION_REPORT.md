# Final Precision Revision — Report

**Manuscript:** *Confidence-Guided Adaptive Reconstruction for Windowed Quantum Phase Estimation*
**Venue:** MethodsX (Elsevier)
**Scope of this pass:** targeted precision only. No method, algorithm, experiment, or numerical
result was touched.

> **Note on the input file.** The brief refers to `main2.pdf`. No such file exists in the
> repository; the current build is `Paper/main.pdf` (37 pp. before this pass). That is the version
> revised here.

---

## 1. Changes Made

### Priority 1 — placeholders (all ten resolved; none remain in either PDF)

Front matter affiliation inserted in both documents; specifications table and Code availability now
cite the real repository; Acknowledgements section removed; CRediT instruction text removed; SI
audit-provenance line rewritten to name artifacts rather than a hash; `refs.bib` proof-stage note
deleted. Detail in §3.

### Priority 2 — 72% attribution wording (4 sites)

Rewritten in Highlights, abstract, Introduction and the Figure 2 caption so the figure is scoped to
the analysed corpus. The underlying result (68 of 94 failures) is unchanged.

### Priority 3 — "pre-registered" → "pre-specified" (6 sites)

No formal preregistration record exists — the only record is a module docstring in
`Experiments/phase7_generalization.py`. The §2.1 paragraph that previously had to walk the term back
was rewritten around the accurate word and is now shorter.

### Priority 4 — central contribution

§2.8 Discussion now opens with the required statement, in bold, before any other material.

### Priority 5 — Module C positioning

The differing evidential status is now signalled at both places where the three modules are first
introduced as a set (abstract, contributions item 2), not only later. Existing positioning at the
module definition, in the adoption table and in Limitation group 4 is unchanged.

### Priority 6 — calibration section heading

Now *"Calibration analysis: $C_w$ provides useful ranking information but is not a calibrated
probability estimate."* The memorable phrasing survives in the contributions list echo and body.
No calibration result was weakened.

### Priority 7 — complexity claim

The bounded-cap assumption is now explicit, including what happens if the caps are allowed to scale.
Checked against §1.12(a)–(c) and `COMPLEXITY_ANALYSIS.md`; no new mathematical claim was introduced.

### Priority 8 — abstract

371 → **336 words**. Removed the Beta/Jeffreys construction detail (it belongs in §1.4). The
matched-resource qualification and the calibration limitation are retained in full.

### Priority 9 — carry-audit history

Trimmed by 2 lines; all six required story beats retained; forensic detail was already in SI §S2.

### Priority 10 — original-method citation metadata

**Verified against arXiv and Quantum, not assumed.** Four entries upgraded:

| Key | Now cites |
|---|---|
| `wqpe2025` | Shukla & Vedula, *Adv. Quantum Technol.* **9**(3), e00683 (2026), doi:10.1002/qute.202500683 |
| `modularshor2025` | Shukla & Vedula, *Eur. Phys. J. Plus* **141**, 474 (2026), doi:10.1140/epjp/s13360-026-07719-0 — was cited as a bare preprint; a journal version exists |
| `bae2025` | Ramôa & Santos, *Quantum* **9**, 1856 (2025), doi:10.22331/q-2025-09-11-1856 |
| `biqae2026` | Li, Vidwans, Wang & Soley, *Quantum* **10**, 1962 (2026), doi:10.22331/q-2026-01-14-1962 |

This also closes the "confirm author attribution" item that had been standing as unresolvable: the
Shukla–Vedula attribution is confirmed correct for both source papers.

### Priority 11 — length

See §5. Net compression happened, but page count did not fall; explained there.

### Priority 12 — formatting audit

**A live defect was found and fixed.** The CRediT statement contained a mangled `\noindent` — the
backslash-n had been eaten by a shell-escaping bug in an earlier session — so the shipped PDF
literally printed *"oindent Sanjay Lakshmanan R: Conceptualization…"*. Confirmed present in the
previous PDF via `pdftotext`, now corrected.

Otherwise: no orphaned headings, no broken cross-references (`??`) or citations (`[?]`), no raw
macro leakage, no table or figure overflow. Two sparse pages remain (main p38, SI p16), both the
tail of a reference/table list — normal.

### Priority 13 — MethodsX fit

All fourteen required elements verified present programmatically. The article still reads as a
reusable method: seven-step adoption protocol, practitioner-framed cost analysis, independently
switchable modules, verified reduction to the published baseline.

---

## 2. Claim Precision Changes

| Claim | Before | After |
|---|---|---|
| **72% attribution** (Highlights) | "Fixed classical reconstruction budgets, not noise, **cause 72% of windowed-QPE failures**." | "**In the analysed failure corpus**, 72% of observed failures were **attributed to** classical reconstruction-budget starvation." |
| **72%** (Introduction) | "…72% of observed failures traced to classical budget starvation" | "…of the 94 failures the campaign produced, 72% were attributed to classical budget starvation… **Both figures characterise the corpus we analysed, not windowed phase estimation in general.**" |
| **Central contribution** (Discussion) | opened with "Three findings sit uneasily together…" | opens with "**The central result of this study is not that confidence-guided allocation has been shown to outperform an equal uniform allocation of the same total resources.** The matched-resource controls do not support that claim. What the evidence shows is that the proposed rules can infer and deploy a workable reconstruction budget from the observed measurement data, whereas the original pipeline requires that budget to be selected manually in advance." |
| **Module C** (abstract) | "three independently switchable modules: …" | "…**B and D carry the stronger independent evidence; C is evaluated chiefly as a component of the full framework.**" |
| **Module C** (contributions) | listed flat with B and D | "**The three are *not* equally supported by the evidence:** B and D each carry a large independent effect, whereas C's isolated contribution is small and does not survive multiple-comparison correction…" |
| **Complexity** | "In this restricted sense… no module introduces a new complexity class." | "**Under the bounded-cap implementation analysed here, the adaptive modules do not alter the asymptotic complexity class…; they modify bounded constants and realised resource usage. If the caps $\delta_{cov}$, $\beam\_max$ and $S_{max}$ were instead allowed to scale with problem size, the resulting complexity would require separate analysis, which we do not attempt.**" |
| **Held-out instances** | "pre-registered" (×6) | "pre-specified" (×6), with an accurate description of what record exists |
| **Calibration heading** | "$C_w$ is a good ranking signal and a poor probability" | "Calibration analysis: $C_w$ provides useful ranking information but is not a calibrated probability estimate" |

Matched-resource interpretation is unchanged in substance and remains prominent: it is a primary
finding in the abstract, a contributions bullet, §2.5 with its own figure and table, the opening of
the Discussion, and Limitation group 2.

---

## 3. Placeholder Resolution

| Placeholder | Where | Resolution |
|---|---|---|
| `DEPARTMENT, INSTITUTION, CITY, COUNTRY` | main + SI front matter | **Resolved** — "Department of Computer Science and Engineering, Amrita Vishwa Vidyapeetham, Chennai, India" (author-supplied) |
| `ORCID iD` | main front matter | **Removed.** It was a LaTeX comment and never rendered. Author will add after confirming; see §6 |
| `public repository URL` ×2 | specifications table, Code availability | **Resolved** — `https://github.com/Shadow-46/QuantumPhaseReconstruction`, recovered from the repository's own git remote. Not invented |
| `Zenodo DOI` ×2 | specifications table, Code availability | **Replaced with truthful wording**: an archived snapshot will be deposited in a persistent repository and its DOI supplied at acceptance. No DOI invented |
| Acknowledgements | back matter | **Section removed**, per author instruction |
| CRediT "confirm the role list…" | back matter | **Removed.** Single-author paper; the existing role list is accurate and retained. Also fixed the mangled `\noindent` on the adjacent line |
| `commit hash or tag` | SI §S2 | **Replaced** with the artifact names (`phase3_adversarial.py`, `phase3_carry_truth_table.csv`, the instrumented-replay record). A hash would have misled: `HEAD` = `e9ff12a` predates 26 uncommitted files including all Phase 7 data |
| "Volume and article number to be completed at proof stage" | `refs.bib` | **Removed** — real volume, article number and DOI now present |

**Nothing was left unresolved and nothing was invented.** Verified: `grep` for
`PLACEHOLDER` / `TODO` / `TBD` / "to be completed" returns zero hits across `main.tex`,
`supplementary.tex` and `refs.bib`, and `pdftotext` finds no `PLACEHOLDER` in either PDF.

---

## 4. Scientific Integrity Check

Confirmed by automated re-verification against the frozen result CSVs:

- **No numerical result changed.** Every reported value was recomputed from
  `Data/failure_study/*.csv` and matched against the manuscript sources: six grids × 144 trials,
  the canonical instrumented run still bit-identical to the published ablation, all 20 per-instance
  success rates, all 16 discordant-count pairs, four matched-resource comparisons, Brier scores
  0.2239 → 0.1988, all six calibration bins, all five realised-resource rows, and the 94-failure
  attribution. **No discrepancies.**
- **No experiment was rerun.** No script under `Experiments/` was executed.
- **No method was changed.** No file under `Algorithms/`, `Circuits/`, `Reconstruction/`,
  `Simulation/` or `Evaluation/` was touched.
- **No unsupported claim was added.** Every edit either narrowed a claim, made an existing
  assumption explicit, or resolved a placeholder from verified data. The one addition of new
  factual content — the four bibliography entries — was verified against arXiv and the publisher.
- **No limitation was removed.** All six limitation groups are intact; Module C's weaker evidence
  is now stated *earlier* as well as in Limitations.

Automated check of the brief's must-remain list: 16 claims, 14 MethodsX sections, 4 DOIs and 2
affiliations — all present.

---

## 5. Page Count

| Document | Before this pass | After |
|---|---|---|
| Main manuscript | 37 pp. | **38 pp.** |
| Supplementary Information | 19 pp. | **19 pp.** |

**The main manuscript grew by one page, and I did not force it back down.** The reason is that
Priorities 4, 5 and 7 each *mandated additional text*: the Discussion opening statement (~11 lines),
the Module C evidential clauses (~7 lines), and the bounded-cap conditional (~7 lines). Against
that I removed ~23 lines of genuine redundancy — the duplicated central-contribution statement in
the Introduction, the third restatement of the calibration limitation, the repeated
matched-resource explanation, and the abstract's implementation detail.

The document is therefore denser in content per page than before, but the same length. Reaching
34–35 pp. was not achievable naturally: after this pass only two pages are not densely set (p38, six
lines of references; p2, the abstract tail), so the remaining three pages would have to come out of
evidence — figures, the reuse protocol, the matched-resource narrative, or limitations — all of
which the brief protects. Per the brief's own instruction that "scientific quality is more important
than an arbitrary page count", I stopped.

---

## 6. Remaining Manual Actions

Only two items genuinely require human input:

1. **ORCID iD** — not recorded anywhere in the project. Add to both `main.tex` and
   `supplementary.tex` front matter once confirmed. (Registration is free at orcid.org; Elsevier
   will prompt for it at submission in any case.)
2. **Archival DOI** — the snapshot has not been deposited yet. Both documents currently state that
   it will be supplied at acceptance, which is accurate; replace with the DOI once minted.

### One judgement call to review

The **Acknowledgements section was removed** as instructed. If the mentor referred to in that
instruction should be acknowledged, or if any institutional or funding support needs declaring, the
section must be restored — Elsevier will not add it, and an omitted funding declaration is a
correction-worthy defect after publication. Restoring it is a four-line change.

### Not blocking, but worth knowing

The repository has 26 uncommitted files, including all Phase 7 data, both manuscript sources and the
figures. Before the archival snapshot is cut, these need committing so the deposited state matches
what the article describes.
