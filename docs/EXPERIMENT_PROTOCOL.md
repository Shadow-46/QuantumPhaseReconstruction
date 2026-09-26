# Experiment protocol (living document; frozen per phase before test-split runs)

Status: **draft v0.1 (2026-09-27)**. The sections for Phases 5–8 must be frozen here, with a git commit and config hashes recorded in `DECISIONS.md`, before any `split: test` run of those phases.

## Units, pairing and splits
- **Unit of inference: the phase (phase_id).** Replicates within a phase are correlated, because they share the same kernel distributions. Confidence intervals therefore use a phase-cluster bootstrap (≥ 2000 resamples for final analyses). Wilson intervals are reported only as within-phase sampling precision.
- **Pairing.**
  - Competing methods are run on the same phases with common random numbers. A block geometry's shot i is identical across policies (`sim/oracle_sim.py`).
  - In fixed-schedule experiments the counts are literally shared.
- **Splits.**
  - `dev` is used for characterisation and tuning.
  - `test` is used only for frozen, predeclared comparisons.
  - The splits use disjoint seed domains (`phases.py`, D-009).
- **Strata.** S0–S5 are defined in `phases.py`. Every result is reported per stratum as well as pooled. The pooled weights are the design weights (equal per stratum) and are stated with the result. They are never re-weighted after seeing results.

## Endpoints
- **Mode A (fixed total budget B shots).**
  - Primary: P(|φ̂−φ|_circ ≤ τ).
  - Secondary: circular RMSE, median error, and exact-best-n-bit rate (excluding S1).
- **Mode B (fixed tolerance τ).**
  - Primary: total shots used to stop. Co-reported: U-queries.
  - Validity condition: realised coverage P(|φ̂−φ| ≤ τ) ≥ 1−α on the test split. A method that fails coverage fails, whatever its savings.
- **τ grid.** τ ∈ {2⁻³, …, 2⁻ⁿ}, keeping only τ ≥ 2⁻ⁿ, because τ < 2⁻⁽ⁿ⁺¹⁾ is unachievable with n bits.

## Tests (for final, test-split comparisons)
- **Binary endpoints.**
  - Paired difference, with a phase-cluster bootstrap CI.
  - Exact McNemar test on trial-level discordant pairs as a secondary check; it ignores clustering, and that limitation is stated.
- **Continuous endpoints.** Wilcoxon signed-rank on per-phase paired means, with the Hodges–Lehmann estimate.
- **Multiplicity.** Holm correction across the predeclared primary family of each phase.
- **Equivalence.** TOST with margins declared before the run: ±1 percentage point for success rates and ±5% for shot ratios. "No difference" is a reportable outcome.
- **Practical effect.** Always report the absolute shot saving (median and IQR) next to any p-value.

## Leakage rules
- Policies and stopping rules may use only the `MeasurementOracle`, declared configuration, budget and τ, and model quantities.
- `evaluation/metrics.py` and the simulator's truth accessors must never be imported by `allocation/`, `overlap/`, `stopping/` or `controller.py`. This is enforced by a test (added with the controller).
- Oracle or reference policies live in `oracle/` and are tagged `is_oracle=True` in every row.

## Tuning
- Policy constants are tuned on `dev`, then frozen: batch size ΔS, thresholds, stopping α-calibration, and the uniform baseline's per-block shots at each budget.
- The uniform baseline is tuned with the same effort as the adaptive policies.
