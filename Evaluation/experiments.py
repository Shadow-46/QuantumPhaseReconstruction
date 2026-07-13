"""
Purpose
    Run reproducible parameter sweeps for windowed phase reconstruction.
Theory
    The experiment driver varies noise labels, window size, overlap, candidate
    count, precision, number of windows implied by geometry, N, and a while
    recording success, runtime, and phase reconstruction metrics.
Inputs
    Command-line sweep controls and configuration dataclasses.
Outputs
    CSV logs, aggregated tables, and plot images under Data/.
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import time
from dataclasses import asdict
from fractions import Fraction
from itertools import product

import pandas as pd

from Algorithms.paper_algorithm import run_windowed_shor
from Algorithms.standard_shor import multiplicative_order
from Evaluation.metrics import circular_phase_error, experiment_metrics
from Simulation.noise import build_noise_model
from config import DATA_DIR, NoiseConfig, ShorConfig, WindowConfig


def run_single_experiment(N: int, a: int, precision: int, window_size: int, overlap: int, candidate_count: int, noise: float = 0.0) -> dict[str, object]:
    """Run one Aer-sampled reconstruction experiment and return a log row.

    `noise` is applied as a depolarizing/readout error rate through an Aer
    noise model, so the noise sweep dimension actually perturbs the sampled
    circuits rather than only labeling the row.
    """
    started = time.perf_counter()
    shor = ShorConfig(N=N, a=a, phase_qubits=precision)
    window = WindowConfig(total_precision=precision, window_size=window_size, overlap=overlap, candidate_count=candidate_count)
    noise_model = build_noise_model(
        NoiseConfig(enabled=noise > 0.0, depolarizing_1q=noise, depolarizing_2q=noise, measurement_error=noise)
    )
    result = run_windowed_shor(shor, window, noise_model=noise_model)
    runtime = time.perf_counter() - started
    order = multiplicative_order(a, N)
    truth = Fraction(1, order) if order else None
    phase_error = circular_phase_error(result.best_phase, truth) if result.best_phase is not None and truth is not None else None
    row = {
        "N": N,
        "a": a,
        "precision": precision,
        "window_size": window_size,
        "overlap": overlap,
        "candidate_budget": candidate_count,
        "noise": noise,
        "num_windows": len(result.reconstructed_phases[0].bitstring) if result.reconstructed_phases else 0,
        "order": result.order,
        "factors": str(result.factors),
        **experiment_metrics(result.success, runtime, len(result.reconstructed_phases), phase_error),
    }
    row.update({f"shor_{key}": value for key, value in asdict(shor).items() if key in {"shots", "random_seed"}})
    return row


def run_sweep(quick: bool = False) -> pd.DataFrame:
    """Run the configured sweep and return a Pandas DataFrame."""
    moduli = [(15, 2), (21, 2)] if quick else [(15, 2), (21, 2), (33, 5)]
    precisions = [6, 8] if quick else [6, 8, 10]
    window_sizes = [3, 4] if quick else [3, 4, 5]
    overlaps = [1, 2]
    candidate_counts = [1, 2, 4]
    noise_levels = [0.0] if quick else [0.0, 0.001, 0.01]
    rows: list[dict[str, object]] = []
    for (N, a), precision, window_size, overlap, candidates, noise in product(moduli, precisions, window_sizes, overlaps, candidate_counts, noise_levels):
        if overlap < window_size and window_size <= precision:
            rows.append(run_single_experiment(N, a, precision, window_size, overlap, candidates, noise))
    return pd.DataFrame(rows)


def write_outputs(frame: pd.DataFrame, output_dir: Path = DATA_DIR) -> tuple[Path, Path]:
    """Write CSV logs and aggregate tables for an experiment DataFrame."""
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "experiment_log.csv"
    table_path = output_dir / "summary_table.csv"
    frame.to_csv(log_path, index=False)
    summary = frame.groupby(["N", "a", "precision", "window_size", "overlap", "candidate_budget"], as_index=False).agg(success_rate=("success", "mean"), mean_runtime=("runtime_seconds", "mean"), mean_phase_error=("phase_error", "mean"))
    summary.to_csv(table_path, index=False)
    return log_path, table_path


def build_parser() -> argparse.ArgumentParser:
    """Create the experiment command-line parser."""
    parser = argparse.ArgumentParser(prog="experiments")
    parser.add_argument("--quick", action="store_true")
    return parser


def main() -> int:
    """Run the sweep, write tables, and generate plots."""
    args = build_parser().parse_args()
    frame = run_sweep(quick=args.quick)
    log_path, table_path = write_outputs(frame)
    from Evaluation.plots import generate_standard_plots
    plot_paths = generate_standard_plots(frame, DATA_DIR)
    print(f"wrote {log_path}")
    print(f"wrote {table_path}")
    for path in plot_paths:
        print(f"wrote {path}")
    return 0


def self_test() -> None:
    """Verify the quick sweep produces successful rows."""
    frame = run_sweep(quick=True)
    assert not frame.empty
    assert frame["success"].any()


if __name__ == "__main__":
    raise SystemExit(main())
