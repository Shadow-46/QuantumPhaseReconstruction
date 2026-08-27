# Paper build notes

The work is written up as **two cross-cited manuscripts**, each with its own
supplementary. They share one released corpus and analyse it at different units
— Paper A at the window, Paper B at the trial — and no figure, table, or result
is common to both.

## Files

| File | What it is |
|---|---|
| `A_method.tex` | **Paper A**, the method article. Targets MethodsX. `elsarticle`. |
| `A_supplementary.tex` | Paper A's SI: derivations (S1), calibration bins (S2), artifact map (S3). |
| `refs_A.bib` | Paper A bibliography, 18 entries (MethodsX allows 25). |
| `graphical_abstract.tex` | Generates Paper A's graphical abstract from its Figure 1. |
| `B_master.tex` | **Paper B**, the study article. Journal-neutral driver. |
| `B/` | Paper B's content tree — see below. |
| `B_supplementary.tex` | Paper B's SI: carry audit (S1), statistical tables (S2), artifact map (S3). |
| `si_heldout_table.tex` | Full per-instance ablation table, `\input`-ed by Paper B's SI. |
| `refs_B.bib` | Paper B bibliography, 21 entries, 12 cited. |
| `BUILD.md` | This file. |

`A_method.tex` and the `B/` tree are the **source of truth** and are edited
directly. (Earlier drafts were generated from a since-retired `main.tex`; that
file, `supplementary.tex` and `refs.bib` were removed once superseded and are
recoverable at the `pre-paper-split` tag.)

## Building

Each document is a separate LaTeX root. From this directory:

```bash
pdflatex -interaction=nonstopmode A_method
bibtex A_method
pdflatex -interaction=nonstopmode A_method
pdflatex -interaction=nonstopmode A_method
```

Same four-step cycle for `A_supplementary`, `B_master`, and `B_supplementary`.
Two passes after BibTeX are required — the first resolves citations, the second
the page cross-references.

The graphical abstract is built separately and converted to PNG:

```bash
pdflatex -interaction=nonstopmode graphical_abstract
pdftocairo -png -r 300 -singlefile graphical_abstract.pdf graphical_abstract
```

Expected clean state: **0 errors and 0 undefined references or citations** in
all four documents. Anything else is a regression.

| Document | Pages | Overfull boxes |
|---|---|---|
| `A_method` | 28 (preprint) / 14 (two-column) | 12 / 3 |
| `A_supplementary` | 11 / 6 | 10 / 4 |
| `B_master` | 30 | 6 |
| `B_supplementary` | 11 | 2 |

## Paper A: the layout switch

`A_method.tex` and `A_supplementary.tex` each carry two `\documentclass` lines
at the top, exactly one active:

```latex
\documentclass[preprint,11pt]{elsarticle}              % active default
%\documentclass[final,5p,times,longtitle]{elsarticle}  % two-column
```

`preprint` is the format Elsevier asks for — *"Articles are prepared and
submitted in single column format even if the final printed article will come
in a double column format journal."* The two-column line shows the published
length and is permitted for LaTeX submissions.

Nothing in the body depends on the choice. Floats too wide for one column use
`widefigure` / `widetable` / `widealgorithm`, which resolve once at preamble
time to starred or unstarred floats. `placeins` is loaded in one-column mode
only, since a section-level float barrier cannot be satisfied by a starred
float.

## Paper B: the journal profile

`B_master.tex` is journal-neutral. Retargeting touches only two lines there
plus one file in `B/config/`:

```latex
\documentclass[preprint,11pt]{elsarticle}   % LaTeX requires this first
\newcommand{\JournalProfile}{neutral}       % selects B/config/neutral.tex
```

| Path | Contents |
|---|---|
| `B/config/neutral.tex` | The only file a new target needs forking: `\journal`, `\BibStyle`, float fractions, spacing, wide-float definitions. |
| `B/config/README.md` | How to add a profile, and what each likely target requires. |
| `B/preamble.tex` | Journal-independent packages, theorem environments, symbol macros. |
| `B/frontmatter.tex` | Title, author, abstract, keywords. |
| `B/sections/*.tex` | One file per section, 22 of them, in reading order. |
| `B/backmatter.tex` | Declarations, availability statements, companion relation. |

**No file under `B/sections/` contains a journal-specific command.**

No target journal is selected yet. `B/config/README.md` records the candidates
and their prerequisites — note that `IEEEtran.cls` and `ieeeaccess.cls` are not
currently installed in this MiKTeX tree.

## What lives where

| Paper A owns | Paper B owns |
|---|---|
| The C_w statistic, modules B/C/D, Algorithm 1 | The Phase 1–5 failure characterisation |
| Theoretical properties, cost analysis, reuse protocol | The five-arm ablation and held-out generalisation |
| Window-level validation and the calibration study | The matched-resource controls |

Cross-references between the papers are citations (`companion` in `refs_A.bib`,
`methodpaper` in `refs_B.bib`), never duplicated content. Both papers disclose
that they analyse one released corpus at different units.

## Figures

| Paper | Figure | Source |
|---|---|---|
| A | method schematic | inline TikZ in `A_method.tex` |
| A | window retention | `Experiments/paperA_window_validation.py` |
| A | calibration reliability | `Experiments/paper_figures.py` |
| B | logistic effect sizes, noise heatmap, adversarial conditions, scaling | phase scripts, rendered during the campaign |
| B | failure attribution, held-out, matched-resource | `Experiments/paper_figures.py` |

Figure paths are relative (`../Data/failure_study/plots/...`). **Flatten them
before uploading to Overleaf or a publisher system**, which will not resolve a
parent directory.

## Cross-document macros

The `\SIproofs` / `\SIcalib` / `\SIartifacts` (Paper A) and `\SIaudit` /
`\SItables` / `\SIartifacts` (Paper B) macros expand to **literal** section
numbers — LaTeX cannot resolve a `\ref` across documents. If either
supplementary is renumbered, update the macro definitions in the corresponding
preamble by hand.

## Remaining open items

| Item | Status |
|---|---|
| ORCID iD | Not recorded; required by MethodsX for the submitting author. Add to both front matters. |
| Archival DOI | Not yet minted; both documents state it will be supplied at acceptance. |
| Paper B venue | Not selected. |
| Bibliography verification | Three entries unconfirmed: `wqpe2025` volume/DOI, `bae2025`/`biqae2026` author lists, attribution on `wqpe2025`/`modularshor2025`. |

## Checking the numbers

Every numeric claim in either manuscript is a mechanical summary of a released
result file. `Experiments/paper_figures.py` renders figures from frozen CSVs
only — it runs no trial and recomputes no outcome, and asserts the reported
Brier scores before writing. `Experiments/paperA_window_validation.py` does the
same for Paper A's window-level table, asserting the calibration corpus size
and raw Brier score first.

The artifact map in each supplementary maps every result file to the script
that produced it and the section that consumes it.
