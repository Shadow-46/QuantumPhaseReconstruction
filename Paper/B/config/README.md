# Journal profiles

`B_master.tex` is journal-neutral. Retargeting touches only:

1. the `\documentclass` line in `B_master.tex` (LaTeX requires it first), and
2. `\JournalProfile`, which selects a file in this directory.

Nothing under `B/sections/` refers to a journal, so no section file changes.

## Adding a profile

Copy `neutral.tex`, adjust, and point `\JournalProfile` at it.

A profile owns: `\journal{}`, `\BibStyle`, float fractions and spacing, and
the `widefigure` / `widetable` definitions. Sections use those two
environments for anything too wide for a single column, so a two-column
profile needs no edits to the content.

## Known requirements of likely targets

- **IEEE Access** and **IEEE TQE** need `IEEEtran.cls`; IEEE Access also needs
  `ieeeaccess.cls`. Neither is currently installed in this MiKTeX tree --
  install before attempting either. Both are two-column, use
  `\bibliographystyle{IEEEtran}`, replace `\begin{keyword}` with index terms,
  and IEEE Access additionally requires an author biography with a photograph.
  Note IEEEtran has no `\part`, so the two Part dividers need replacing.
- **Springer QINP** uses `svjour3.cls`, is single-column, and wants
  `\keywords{}` rather than the elsarticle `keyword` environment.
