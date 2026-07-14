"""
Purpose
    Phase 5: counterfactual failure-ladder attribution for every failed
    trial recorded in Phases 1-4.
Theory
    A single failed trial is ambiguous: it could be quantum shot noise,
    injected noise, beam-search pruning, candidate-count starvation,
    stitching collapse, or continued-fraction failure. Replaying a fixed,
    ordered sequence of counterfactuals -- each relaxing exactly one
    candidate confound -- and recording the *first* one that flips the
    outcome to success gives an evidence-based, not guessed, root cause.
Inputs
    Data/failure_study/phase1_ofat.csv and phase2_noise_*.csv (failed rows).
Outputs
    Data/failure_study/phase5_attribution.csv with a `failure_category`
    column per failed trial.
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from concurrent.futures import ProcessPoolExecutor, TimeoutError as FutureTimeoutError

import pandas as pd

from Algorithms.paper_algorithm import run_windowed_shor, run_windowed_shor_exact
from Simulation.noise import build_noise_model
from Experiments import harness
from config import NoiseConfig, ShorConfig, WindowConfig

# Kept modest deliberately: combined with candidate_count=2**window_size in
# the "exhaustive" replay and up to 32 windows (phase1b's total_precision=32,
# overlap=3 configs), a proxy in the thousands caused multi-GB path-list
# growth in stitch_candidates and made attribution impractically slow (see
# REPORT.md Appendix A). 64 is still >>32x the largest max_paths actually
# used in any failing trial (phase1b's grid never exceeds 2), so it remains
# a meaningful "much wider beam" counterfactual without the blowup, and a
# per-item timeout below guards against any remaining pathological case.
MAX_PATHS_PROXY = 64
ATTRIBUTION_TIMEOUT_S = 30


def _row_to_configs(row: pd.Series) -> tuple[ShorConfig, WindowConfig, dict]:
    shor = ShorConfig(
        N=int(row["spec_N"]), a=int(row["spec_a"]), phase_qubits=int(row["spec_total_precision"]),
        shots=int(row["spec_shots"]), random_seed=int(row["spec_seed"]) if pd.notna(row["spec_seed"]) else None,
    )
    window = WindowConfig(
        total_precision=int(row["spec_total_precision"]), window_size=int(row["spec_window_size"]),
        overlap=int(row["spec_overlap"]), candidate_count=int(row["spec_candidate_count"]),
        max_paths=int(row["spec_max_paths"]),
    )
    noise = {"noise_1q": row.get("spec_noise_1q", 0.0), "noise_2q": row.get("spec_noise_2q", 0.0), "noise_ro": row.get("spec_noise_ro", 0.0)}
    return shor, window, noise


def attribute_failure(row: pd.Series) -> str:
    shor, window, noise = _row_to_configs(row)

    exact = run_windowed_shor_exact(shor, window)
    if not exact.success:
        # Classical reconstruction would fail even with perfect measurements.
        wide_beam = WindowConfig(window.total_precision, window.window_size, window.overlap, window.candidate_count, MAX_PATHS_PROXY)
        if run_windowed_shor_exact(shor, wide_beam).success:
            return "beam_search_limitation"
        exhaustive = WindowConfig(window.total_precision, window.window_size, window.overlap, 2**window.window_size, MAX_PATHS_PROXY)
        if run_windowed_shor_exact(shor, exhaustive).success:
            return "candidate_pruning"
        from Circuits.windowed_qpe import deterministic_window_samples, make_overlapping_windows
        from fractions import Fraction
        from Algorithms.standard_shor import multiplicative_order

        order = multiplicative_order(shor.a, shor.N)
        if order is None:
            return "unknown"
        specs = make_overlapping_windows(window.total_precision, window.window_size, window.overlap)
        counts = deterministic_window_samples(Fraction(1, order), specs)
        stitched, recovered = harness.reconstruct_from_counts(
            shor.N, shor.a, window.total_precision, window.window_size, window.overlap, 2**window.window_size, MAX_PATHS_PROXY, counts
        )
        if not stitched:
            return "stitching_failure"
        return "continued_fraction_failure"

    # Classical logic is fine given perfect bits: the fault is on the quantum-sampling side.
    noise_enabled = any(v > 0 for v in noise.values())
    if noise_enabled:
        clean_result = run_windowed_shor(shor, window, noise_model=None)
        if clean_result.success:
            return "noise_sensitivity"
    more_shots_config = ShorConfig(shor.N, shor.a, shor.phase_qubits, shor.shots * 10, shor.simulator_method, shor.random_seed)
    noise_model = build_noise_model(NoiseConfig(enabled=noise_enabled, depolarizing_1q=noise["noise_1q"], depolarizing_2q=noise["noise_2q"], measurement_error=noise["noise_ro"]))
    more_shots_result = run_windowed_shor(more_shots_config, window, noise_model=noise_model)
    if more_shots_result.success:
        return "quantum_estimation_shot_limited"
    return "unknown"


def run_attribution(csv_names: list[str], sample_cap: int = 150) -> pd.DataFrame:
    frames = []
    for name in csv_names:
        path = harness.FAILURE_STUDY_DIR / f"{name}.csv"
        if not path.exists():
            continue
        frame = pd.read_csv(path)
        frame["source_csv"] = name
        frames.append(frame)
    if not frames:
        raise FileNotFoundError("no phase CSVs found to attribute")
    combined = pd.concat(frames, ignore_index=True)
    failed = combined[combined["result_success"] == False]  # noqa: E712
    if len(failed) > sample_cap:
        failed = failed.sample(sample_cap, random_state=0)
    print(f"attributing {len(failed)} failed trials out of {len(combined)} total")

    categories = []
    for i, (_, row) in enumerate(failed.iterrows()):
        try:
            with ProcessPoolExecutor(max_workers=1) as pool:
                future = pool.submit(attribute_failure, row)
                try:
                    categories.append(future.result(timeout=ATTRIBUTION_TIMEOUT_S))
                except FutureTimeoutError:
                    categories.append("unknown_attribution_timeout")
                    pool.shutdown(wait=False, cancel_futures=True)
        except Exception as exc:  # noqa: BLE001
            categories.append(f"attribution_error: {exc}")
        if (i + 1) % 10 == 0 or (i + 1) == len(failed):
            print(f"  attributed {i + 1}/{len(failed)}", flush=True)
    failed = failed.copy()
    failed["failure_category"] = categories
    failed.to_csv(harness.FAILURE_STUDY_DIR / "phase5_attribution.csv", index=False)
    return failed


def main() -> None:
    result = run_attribution(["phase1_ofat", "phase2_noise_single", "phase2_noise_combined", "phase1b_combined_stress"])
    print(result["failure_category"].value_counts())


if __name__ == "__main__":
    main()
