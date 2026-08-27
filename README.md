# QuantumPhaseReconstruction

Research-grade Python platform for reproducing and extending windowed quantum
phase reconstruction for Shor's algorithm.

The repository is organized around a hard boundary between quantum sampling and
classical reconstruction:

- `Circuits/`: modular multiplication, inverse QFT, standard QPE, and windowed
  measurement circuits.
- `Reconstruction/`: window candidates, carry/overlap checks, stitching, and
  continued-fraction order recovery.
- `Algorithms/`: standard Shor (`standard_shor.py`), the paper-style windowed
  algorithm (`paper_algorithm.py`), the confidence-guided adaptive method
  itself (`adaptive_reconstruction.py`), and a small window-geometry policy
  (`window_policy.py`).
- `Simulation/`: Aer backends and noise models.
- `Evaluation/`: experiment sweeps, metrics, plots, tables, and CSV logs.
- `Experiments/`: the failure-study phase campaign (see below) and its shared
  harness.

## Install

```bash
python3.12 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

`requirements.txt` pins exact versions rather than ranges. The recorded
per-trial outcomes reproduce bit-for-bit under those versions; a different
Aer build may resample differently.

Then verify the environment before running anything. This checks every pinned
dependency and runs the self-test of all 18 modules:

```bash
python validate_environment.py
```

## Smoke run

```bash
python main.py --N 15 --a 2 --phase-qubits 8 --shots 2048
```

## Module self-tests

Every executable module carries a self-test that runs without arguments;
`validate_environment.py` above runs all of them. To run one directly:

```bash
python Algorithms/adaptive_reconstruction.py   # asserts the baseline-equivalence invariant
python Reconstruction/confidence.py            # symmetry, monotonicity, sample-size sensitivity
python Circuits/windowed_qpe.py                # block-local vs full-register agreement
```

## Experiment sweep

```bash
python Evaluation/experiments.py --quick
```

The driver writes CSV logs, summary tables, and plots under `Data/`.

## Failure-study phase campaign

`Experiments/` holds a separate, larger campaign investigating and fixing
practical failure modes of the windowed reconstruction pipeline, sharing
`Experiments/harness.py` as a common driver. Each phase script writes its
CSVs and plots under `Data/failure_study/`:

```bash
python Experiments/phase1_ofat.py             # one-factor-at-a-time sweeps
python Experiments/phase1b_combined_stress.py # combined-parameter stress grid
python Experiments/phase2_noise.py            # single/combined Aer noise sweeps
python Experiments/phase3_adversarial.py       # adversarial fault injection, carry audit
python Experiments/phase4_scaling.py           # scaling behaviour by N/precision/window
python Experiments/phase5_attribution.py       # failure-mode attribution
python Experiments/recapture_window_counts.py  # backfills per-window raw data
python Experiments/phase6_calibration.py       # confidence-statistic calibration
python Experiments/phase6_ablation.py          # ablation of the adaptive modules (B/C/D)
python Experiments/phase7_generalization.py    # held-out instances + matched-resource controls
python Experiments/analyze.py                  # aggregate report + Pareto plot
```

Phase 7 is long-running and resumable. `Experiments/resume_campaign.py`
restarts it safely (finished stages are skipped) and
`Experiments/check_status.py` reports progress, since a single heavy trial can
take 20-25 minutes.

Analyses and figures for the manuscripts are then rendered from the frozen
CSVs, without re-running any trial:

```bash
python Experiments/paperA_window_validation.py  # window-level retention + C_w discrimination
python Experiments/paper_figures.py             # publication figures from frozen result files
```

Phase 6 introduces confidence-guided adaptive reconstruction
(`Reconstruction/confidence.py`, `AdaptiveReconstructionConfig` in
`config.py`): coverage-based candidate retention, ambiguity-scaled beam
width, and confidence-gated shot-stopping. See
`Data/failure_study/PHASE6_SUMMARY.md` for the narrative writeup and
`Data/failure_study/REPORT.md` / `COMPILED_REPORT.md` for the full findings,
and `Paper/` for the two manuscripts built on this data — `A_method.tex`
(the method and its window-level validation, targeting MethodsX) and
`B_master.tex` (the failure characterisation and end-to-end evaluation,
journal-neutral). See `Paper/BUILD.md` for how to build them.

## Licence

Apache-2.0. See `LICENSE`.
