"""
Purpose
    Reproduce the paper-style windowed QPE reconstruction pipeline.
Theory
    The algorithm samples local phase windows on an Aer backend (one shallow
    block-local circuit per window, per Algorithm 2 of arXiv:2509.05010),
    generates the top candidates per window, checks carry-aware overlaps,
    stitches consistent paths, and only then invokes continued fractions for
    Shor order recovery.
Inputs
    ShorConfig, WindowConfig, and an optional Aer noise model.
Outputs
    WindowedShorResult containing reconstructed phases, order, and factors.
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from dataclasses import dataclass
from fractions import Fraction
from math import gcd

from config import ShorConfig, WindowConfig
from Circuits.modular_multiplication import ModularMultiplicationOperator
from Circuits.windowed_qpe import WindowSpec, deterministic_window_samples, make_overlapping_windows, run_windowed_qpe_block
from Reconstruction.candidate_generation import generate_window_candidates
from Reconstruction.continued_fraction import recover_order_and_factors
from Reconstruction.stitching import StitchedPhase, stitch_candidates


@dataclass(frozen=True)
class WindowedShorResult:
    """Complete result for the windowed reconstruction pipeline."""
    N: int
    a: int
    reconstructed_phases: tuple[StitchedPhase, ...]
    order: int | None
    factors: tuple[int, int] | None

    @property
    def best_phase(self) -> Fraction | None:
        """Return the highest-scoring reconstructed phase, if available."""
        return self.reconstructed_phases[0].phase if self.reconstructed_phases else None

    @property
    def success(self) -> bool:
        """Return True when non-trivial factors were recovered."""
        return self.factors is not None


def reconstruct_from_window_counts(N: int, a: int, total_precision: int, window_size: int, overlap: int, candidate_count: int, max_paths: int, window_counts: list[dict[str, int]]) -> WindowedShorResult:
    """Run candidate generation, overlap stitching, and continued fractions."""
    specs = make_overlapping_windows(total_precision, window_size, overlap)
    if len(specs) != len(window_counts):
        raise ValueError("window_counts length does not match generated windows.")
    candidate_groups = [
        generate_window_candidates(counts, spec.start, spec.width, total_precision, candidate_count)
        for spec, counts in zip(specs, window_counts, strict=True)
    ]
    phases = tuple(stitch_candidates(candidate_groups, max_paths=max_paths))
    for phase in phases:
        recovered = recover_order_and_factors(phase.phase, a, N)
        if recovered.factors is not None:
            return WindowedShorResult(N, a, phases, recovered.order, recovered.factors)
    return WindowedShorResult(N, a, phases, None, None)


def sample_window_counts(shor: ShorConfig, window: WindowConfig, noise_model=None) -> tuple[list[dict[str, int]], list[WindowSpec]]:
    """Sample every window's block-local QPE circuit on Aer and return raw counts."""
    op = ModularMultiplicationOperator(shor.a, shor.N)
    specs = make_overlapping_windows(window.total_precision, window.window_size, window.overlap)
    counts = [
        run_windowed_qpe_block(
            op,
            eigenstate=1,
            spec=spec,
            shots=shor.shots,
            method=shor.simulator_method,
            seed=shor.random_seed,
            noise_model=noise_model,
        )
        for spec in specs
    ]
    return counts, specs


def run_windowed_shor(shor: ShorConfig, window: WindowConfig, noise_model=None) -> WindowedShorResult:
    """Run the windowed reconstruction pipeline using genuine Aer-sampled QPE
    blocks (no classical shortcut through the already-known order)."""
    if gcd(shor.a, shor.N) != 1:
        divisor = gcd(shor.a, shor.N)
        return WindowedShorResult(shor.N, shor.a, tuple(), None, tuple(sorted((divisor, shor.N // divisor))))
    counts, _specs = sample_window_counts(shor, window, noise_model=noise_model)
    return reconstruct_from_window_counts(shor.N, shor.a, window.total_precision, window.window_size, window.overlap, window.candidate_count, window.max_paths, counts)


def run_windowed_shor_exact(shor: ShorConfig, window: WindowConfig) -> WindowedShorResult:
    """Run the reconstruction pipeline against exact, noiseless window bits
    sliced from the already-known order.

    This bypasses quantum sampling entirely and is intended only as a fast
    unit-test oracle for the classical candidate/stitching/continued-fraction
    logic; it is not a physical simulation and must not be used to claim
    reproduction of the paper's quantum measurement process.
    """
    if gcd(shor.a, shor.N) != 1:
        divisor = gcd(shor.a, shor.N)
        return WindowedShorResult(shor.N, shor.a, tuple(), None, tuple(sorted((divisor, shor.N // divisor))))
    from Algorithms.standard_shor import multiplicative_order
    order = multiplicative_order(shor.a, shor.N)
    if order is None:
        return WindowedShorResult(shor.N, shor.a, tuple(), None, None)
    phase = Fraction(1, order)
    specs = make_overlapping_windows(window.total_precision, window.window_size, window.overlap)
    counts = deterministic_window_samples(phase, specs)
    return reconstruct_from_window_counts(shor.N, shor.a, window.total_precision, window.window_size, window.overlap, window.candidate_count, window.max_paths, counts)


def paper_example() -> WindowedShorResult:
    """Run the canonical N=15, a=2, phase=1/4 windowed reconstruction example."""
    return run_windowed_shor(ShorConfig(N=15, a=2, phase_qubits=8), WindowConfig(total_precision=8, window_size=4, overlap=2, candidate_count=2))


def self_test() -> None:
    """Verify the paper example reproduces the N=15 factors from genuine
    Aer-sampled circuits, and that the exact reference oracle agrees.

    The standard Shor initial state |1> is an equal superposition of the
    order-4 eigenstates with phases k/4 for k=0..3, so the top-scoring
    stitched phase is not deterministically 1/4; only membership in that
    eigenphase set and successful factor recovery are guaranteed.
    """
    result = paper_example()
    assert result.best_phase in {Fraction(0, 1), Fraction(1, 4), Fraction(1, 2), Fraction(3, 4)}
    assert result.factors == (3, 5)

    exact = run_windowed_shor_exact(ShorConfig(N=15, a=2, phase_qubits=8), WindowConfig(total_precision=8, window_size=4, overlap=2, candidate_count=2))
    assert exact.best_phase == Fraction(1, 4)
    assert exact.factors == (3, 5)


if __name__ == "__main__":
    self_test()
    print("paper_algorithm self-test passed")
