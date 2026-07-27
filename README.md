# QuantumPhaseReconstruction

Research-grade Python platform for reproducing and extending windowed quantum
phase reconstruction for Shor's algorithm.

The repository is organized around a hard boundary between quantum sampling and
classical reconstruction:

- `Circuits/`: modular multiplication, inverse QFT, standard QPE, and windowed
  measurement circuits.
- `Reconstruction/`: window candidates, carry/overlap checks, stitching, and
  continued-fraction order recovery.
- `Algorithms/`: standard Shor, the paper-style windowed algorithm, and a clean
  adaptive extension point.
- `Simulation/`: Aer backends and noise models.
- `Evaluation/`: experiment sweeps, metrics, plots, tables, and CSV logs.
- `Experiments/`: the failure-study phase campaign (see below) and its shared
  harness.

## Install

```bash
python3.12 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Smoke run

```bash
python main.py --N 15 --a 2 --phase-qubits 8 --shots 2048
```

## Module self-tests

Each executable module has a lightweight self-test:

```bash
python Circuits/modular_multiplication.py
python Circuits/iqft.py
python Circuits/qpe.py
python Algorithms/standard_shor.py
python Circuits/windowed_qpe.py
python Algorithms/paper_algorithm.py
python Evaluation/experiments.py
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
python Experiments/analyze.py                  # aggregate report + Pareto plot
```

Phase 6 introduces confidence-guided adaptive reconstruction
(`Reconstruction/confidence.py`, `AdaptiveReconstructionConfig` in
`config.py`): coverage-based candidate retention, ambiguity-scaled beam
width, and confidence-gated shot-stopping. See
`Data/failure_study/PHASE6_SUMMARY.md` for the narrative writeup and
`Data/failure_study/REPORT.md` / `COMPILED_REPORT.md` for the full findings,
and `Paper/` for the MethodsX manuscript draft built on this data.
