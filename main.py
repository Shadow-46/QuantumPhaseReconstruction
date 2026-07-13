"""
Purpose
    Command-line entry point for a fast repository smoke run.
Theory
    The research workflow separates quantum sampling from classical recovery.
    This entry point exercises both paths on small Shor instances without
    launching long experiments.
Inputs
    Optional command-line flags for modulus, base, shots, and precision.
Outputs
    A compact report containing recovered order, factors, and reconstruction
    diagnostics.
Author
    Sanjay
"""

from __future__ import annotations

import argparse

from Algorithms.paper_algorithm import run_windowed_shor
from Algorithms.standard_shor import run_standard_shor
from Simulation.noise import build_noise_model
from config import DEFAULT_SHOR, DEFAULT_WINDOW, NoiseConfig, ShorConfig, WindowConfig


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line parser for the smoke run."""

    parser = argparse.ArgumentParser(prog="QuantumPhaseReconstruction")
    parser.add_argument("--N", type=int, default=DEFAULT_SHOR.N)
    parser.add_argument("--a", type=int, default=DEFAULT_SHOR.a)
    parser.add_argument("--phase-qubits", type=int, default=DEFAULT_SHOR.phase_qubits)
    parser.add_argument("--shots", type=int, default=DEFAULT_SHOR.shots)
    parser.add_argument("--window-size", type=int, default=DEFAULT_WINDOW.window_size)
    parser.add_argument("--overlap", type=int, default=DEFAULT_WINDOW.overlap)
    parser.add_argument("--candidate-count", type=int, default=DEFAULT_WINDOW.candidate_count)
    parser.add_argument("--noise", type=float, default=0.0, help="depolarizing/readout error rate for Aer sampling")
    return parser


def main() -> int:
    """Run a short standard and windowed reconstruction smoke test."""

    args = build_parser().parse_args()
    shor = ShorConfig(N=args.N, a=args.a, phase_qubits=args.phase_qubits, shots=args.shots)
    window = WindowConfig(
        total_precision=args.phase_qubits,
        window_size=args.window_size,
        overlap=args.overlap,
        candidate_count=args.candidate_count,
    )
    noise_model = build_noise_model(
        NoiseConfig(enabled=args.noise > 0.0, depolarizing_1q=args.noise, depolarizing_2q=args.noise, measurement_error=args.noise)
    )
    standard = run_standard_shor(shor.N, shor.a)
    windowed = run_windowed_shor(shor, window, noise_model=noise_model)
    print("QuantumPhaseReconstruction")
    print(f"N={shor.N}, a={shor.a}, phase_qubits={shor.phase_qubits}, shots={shor.shots}")
    print(f"standard_order={standard.order}, standard_factors={standard.factors}")
    print(f"windowed_order={windowed.order}, windowed_factors={windowed.factors}")
    print(f"best_phase={windowed.best_phase} candidates={len(windowed.reconstructed_phases)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
