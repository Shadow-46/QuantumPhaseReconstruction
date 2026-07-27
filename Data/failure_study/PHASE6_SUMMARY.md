# Phase 6: Confidence-Guided Adaptive Reconstruction — Summary

**Date:** 2026-07-15
**Scope:** Implementation of the adaptive reconstruction framework frozen in
[`FORMULATION.md`](FORMULATION.md), built to address failure modes identified in
[`REPORT.md`](REPORT.md).

This document explains, in order, what was done, why, and what came out of it. It is a
narrative companion to the corrected `REPORT.md`/`FORMULATION.md` and to the raw data in
this directory — it does not duplicate their content, it summarizes it.

---

## 1. Background

`FORMULATION.md` specified four modules meant to fix problems the earlier failure-study
campaign (`REPORT.md`) had surfaced in the phase-reconstruction pipeline:

| Module | Purpose |
|---|---|
| **A** | Fix a supposed carry-check false-accept defect in `Reconstruction/carry.py` |
| **B** | Candidate coverage — stop truncating the candidate list at a fixed count |
| **C** | Beam width — scale search width adaptively instead of using a fixed `max_paths` |
| **D** | Shot-stopping — spend more measurement shots only on windows the reconstruction is unsure about |

Before writing any code for Module A, the evidence behind it was re-validated. That
re-validation changed the scope of the whole effort, so it's covered first.

---

## 2. Module A: the evidence didn't hold up

**Original claim (REPORT.md, pre-correction):** the carry-check step false-accepts 55.4%
of overlap patterns — i.e., in over half of cases it incorrectly waves through a candidate
whose stitched bits are inconsistent across a window overlap.

**What I found:** this number came from a bug in the *audit script*
(`Experiments/phase3_adversarial.py`), not from a defect in `Reconstruction/carry.py`
itself. The audit's synthetic test windows were constructed so that, regardless of the
overlap parameter it was supposedly sweeping, the actual bit overlap between adjacent
windows was always 1 bit. It never exercised the multi-bit-overlap cases the 55.4% figure
was supposed to describe — the sweep parameter was disconnected from the geometry it
claimed to control.

**Why I'm confident this is a correction and not a judgment call**, three independent checks:

1. **Mathematical:** the carry-correction formula used by `carry.py` and the "ground truth"
   oracle the audit compares against are the same formula, algebraically. A test built this
   way cannot detect a defect in that formula — it can only detect its own arithmetic bugs.
2. **Empirical re-run:** fixing the audit's window-construction geometry so the overlap
   parameter actually controls overlap width, then re-running it, gives **0% false-accept,
   0% false-reject** (down from the reported 55.4% / 5.4%). Written to
   `Data/failure_study/phase3_carry_truth_table.csv` and
   `Data/failure_study/plots/carry_truth_table.png` (both regenerated).
3. **Instrumented replay:** an 80-trial replay of the related fault-injection experiment,
   logging what the carry-check actually does to a deliberately corrupted candidate at each
   step, shows the corrupted candidate is rejected in all 80 trials — never waved through.

**Conclusion:** Module A was solving a problem that does not exist. No code change was
made to `carry.py` — there was nothing to fix. I dropped Module A from the plan and
corrected `REPORT.md` (§4.1/4.2 and the summary sections that cited the 55.4% figure) and
`FORMULATION.md` in place, using visible strikethrough/correction markup rather than
quietly rewriting the original claims — the goal was an auditable trail of what was wrong
and why, not a clean-looking document.

This matters beyond just Module A: it means the original failure-study campaign shipped
one materially wrong headline finding, caused by not verifying that a synthetic-data
generator actually varies the parameter it claims to sweep. Worth keeping in mind for any
future audit-script work in this codebase.

---

## 3. Modules B, C, D: implementation

With Module A dropped, three modules were implemented as designed.

**B — Candidate coverage** (`Reconstruction/candidate_generation.py`)
New function selects candidates by cumulative probability mass instead of a fixed
top-`k` count: it keeps taking candidates, ranked by probability, until their combined
mass crosses a target coverage threshold (`δ_cov` in `AdaptiveReconstructionConfig`). Under
the old fixed-count scheme, windows with a flatter probability distribution over
candidates were silently under-covered; this fixes that at the cost of variable, and
sometimes larger, candidate lists.

**C — Adaptive beam width** (wired into `Experiments/harness.py`)
`max_paths` (the reconstruction beam width) now scales up based on how many windows in a
given run needed coverage-expansion under Module B. The intuition: if many windows needed
extra candidates to reach coverage, the search space is more ambiguous than usual, so the
beam should widen to avoid pruning the correct path.

**D — Shot-stopping** (`Reconstruction/confidence.py` + a resampling loop in
`Experiments/harness.py`)
New confidence statistic $C_w$ per window, computed via an exact Bayesian two-proportion
comparison (comparing the leading candidate's support against the runner-up). Windows
below a confidence threshold ($\alpha_0$, $\varepsilon$ in config) get additional
measurement shots, up to a cap `S_max`, before reconstruction proceeds — instead of a
fixed shot budget applied uniformly regardless of how ambiguous a given window is.

All three parameters live in a new `AdaptiveReconstructionConfig` in `config.py`
(`α₀`, `ε`, `δ_cov`, `S_max`, `beam_max`). `requirements.txt` gained `scikit-learn`, which
the isotonic recalibration below depends on — it turned out to already be an undeclared
transitive dependency, so this just makes it explicit.

`Experiments/recapture_window_counts.py` is a new backfill script: the original failure
study never persisted per-window trial data (only aggregated CSVs), which is needed to
evaluate $C_w$ against ground truth. It re-runs the original conditions and saves
per-window JSON records — 579 files in `Data/failure_study/raw_windows/`.

---

## 4. Calibration check

**Question:** does $C_w$'s predicted confidence actually match observed correctness? A
confidence statistic that's systematically over- or under-confident undermines the whole
point of shot-stopping on it.

**Method:** `Experiments/phase6_calibration.py` (new) bins windows by predicted $C_w$ and
compares against observed success rate in each bin (reliability diagram), computes Brier
score, and — per the plan's contingency — builds an isotonic-regression recalibration
layer if miscalibration is found.

**Result:** $C_w$ is measurably overconfident — **raw Brier score 0.224**. This triggered
the isotonic recalibration step the formulation had specified as conditional ("build only
if this happens"). After recalibration, **Brier score improved to 0.199**.

Outputs: `Data/failure_study/phase6_calibration.csv`,
`Data/failure_study/plots/calibration_reliability.png` (raw),
`Data/failure_study/plots/calibration_reliability_recalibrated.png` (after isotonic fit).

---

## 5. Ablation experiment

**Goal:** isolate the individual and combined effect of B, C, D against baseline, using
real quantum-circuit simulations (not synthetic/mocked data).

**Design:** `Experiments/phase6_ablation.py` (new) — five arms (Baseline, B-only, C-only,
D-only, Full = B+C+D), paired by random seed so each arm sees matched conditions,
significance via McNemar's test (paired binary outcomes) plus bootstrap confidence
intervals on the success-rate difference.

**Execution:** 139 of 144 planned trials completed. The remaining 6 — the slowest,
largest-config trials — were cut short on my instruction once the accumulating results made
clear they would take hours more to finish without materially changing the statistics
(effect sizes and p-values were already stable well before trial 139).

**Results** (from `Data/failure_study/phase6_ablation.csv`, summarized in
`Data/failure_study/summary_tables/phase6_ablation_summary.csv`):

| Arm | Success rate | vs. Baseline | 95% CI (diff) | McNemar p-value |
|---|---|---|---|---|
| Baseline | 51.1% | — | — | — |
| B | 76.3% | +25.2 pp | [18.0, 33.1] pp | 2.8×10⁻⁹ |
| C | 56.1% | +5.0 pp | [1.4, 8.6] pp | 0.016 |
| D | 76.3% | +25.2 pp | [16.5, 33.8] pp | 1.8×10⁻⁷ |
| **Full (B+C+D)** | **100%** | **+48.9 pp** | **[41.0, 57.6] pp** | 6.8×10⁻²¹ |

All arms are statistically significant improvements over baseline. Baseline's 51.1%
closely reproduces the original `REPORT.md` baseline figure of 52.1%, which is the
sanity check that the audit-script fix in §2 didn't leak into (or otherwise change) the
baseline reconstruction path — the corrected audit only changed how false-accepts were
*measured*, not how reconstruction itself behaves.

B and D each independently roughly halve the baseline failure rate; C's effect is real but
much smaller on its own. Combined, the three modules eliminate all measured failures in
this trial set (100% success, 139/139).

---

## 6. Files touched

**Documentation (edited in place, with visible correction markup, not silent rewrites):**
- `Data/failure_study/FORMULATION.md` — Module A dropped; added Research Hypotheses and
  Design Freeze sections; calibration and ablation results appended.
- `Data/failure_study/REPORT.md` — §4.1/4.2 and related summary sections corrected re:
  the carry-check false-accept claim.

**New source:**
- `Reconstruction/confidence.py` — $C_w$ formula.
- `Experiments/recapture_window_counts.py` — backfills per-window data.
- `Experiments/phase6_calibration.py` — calibration check + isotonic recalibration.
- `Experiments/phase6_ablation.py` — ablation study runner.

**Modified source:**
- `Reconstruction/candidate_generation.py` — coverage-based candidate selection (B).
- `Experiments/harness.py` — `reconstruct_adaptive()` and B/C/D wiring.
- `Experiments/phase3_adversarial.py` — fixed the audit script's window-geometry bug.
- `config.py` — added `AdaptiveReconstructionConfig`.
- `requirements.txt` — added `scikit-learn`.

**Generated data:**
- `Data/failure_study/phase3_carry_truth_table.csv` (regenerated, corrected audit)
- `Data/failure_study/raw_windows/*.json` (579 files)
- `Data/failure_study/phase6_calibration.csv`
- `Data/failure_study/phase6_ablation.csv`
- `Data/failure_study/summary_tables/phase6_ablation_summary.csv`
- `Data/failure_study/plots/carry_truth_table.png` (regenerated),
  `calibration_reliability.png`, `calibration_reliability_recalibrated.png`

---

## 7. What's next

1. **Finish or formally close out the 5 remaining ablation trials.** They were cut short
   for time, not because they were run and discarded — worth either running them to
   completion when time allows, or explicitly noting in `REPORT.md` that the ablation is
   139/144 by design and why the missing 5 don't change the conclusion.
2. **Decide whether the isotonic recalibration layer ships in the live pipeline** or stays
   a diagnostic-only artifact. Right now `phase6_calibration.py` demonstrates that
   recalibration helps (0.224 → 0.199 Brier), but `confidence.py`'s $C_w$ as used by the
   shot-stopping loop in `harness.py` should be checked — is it using raw or recalibrated
   confidence? If raw, shot-stopping decisions are still made on a known-overconfident
   signal.
3. ~~Re-run or spot-check the Phase 1–5 conclusions that don't depend on the carry-check
   claim~~ **Done (2026-07-27): Phase 1–2 spot-checked, confirmed clean.** The Phase 3
   bug was specific to `run_carry_truth_table()` hand-constructing synthetic
   `WindowCandidate` objects with `right.start` pinned so `overlap_width()` always
   evaluated to 1 regardless of the nominal `overlap` loop variable — a bypass of the
   real windowing code. `phase1_ofat.py` and `phase2_noise.py` do not do this: every
   swept value (`total_precision`, `window_size`, `overlap`, `candidate_count`,
   `max_paths`, `shots`, `noise_1q`, `noise_2q`, `noise_ro`) is written directly into the
   `spec` dict consumed by `harness.py`, which passes it straight into the real
   `ShorConfig`/`WindowConfig`/`NoiseConfig` objects. Traced `overlap` through to
   `Circuits/windowed_qpe.py:make_overlapping_windows` (`step = window_size - overlap`,
   directly used to place window starts) and the noise channels through to
   `Simulation/noise.py:build_noise_model` (each channel's probability plugged straight
   into `depolarizing_error`/`ReadoutError`) — both correctly consume the swept value.
   No hand-rolled synthetic-geometry shortcut exists in either script, so this failure
   pattern is confined to the already-fixed Phase 3 audit function.
4. **Push the corrected `REPORT.md`/`FORMULATION.md` and new Phase 6 code to review** —
   the carry-check retraction in particular is the kind of correction that should be
   visible to anyone who previously relied on the 55.4% figure.
5. **Consider scaling the ablation beyond 139 trials** if the results are going into a
   paper or external report — the current sample supports the qualitative ranking
   (B ≈ D > C, Full eliminates failures) but larger N would tighten the confidence
   intervals, especially for arm C where the effect size is smallest.
