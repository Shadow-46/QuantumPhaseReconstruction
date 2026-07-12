"""
Purpose
    Compute accuracy, runtime, and success metrics for reconstruction experiments.
Theory
    Windowed reconstruction is evaluated by exact phase agreement, recovered
    order validity, factor success, candidate volume, and wall-clock runtime.
Inputs
    Experiment outcomes, expected phases, elapsed seconds, and recovered factors.
Outputs
    Metric dictionaries suitable for Pandas tables and CSV logs.
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from fractions import Fraction
from typing import Iterable


def circular_phase_error(estimate: Fraction, truth: Fraction) -> float:
    """Return wrap-around absolute phase error on the unit interval."""
    delta = abs(float(estimate - truth)) % 1.0
    return min(delta, 1.0 - delta)


def success_rate(values: Iterable[bool]) -> float:
    """Return the fraction of successful boolean outcomes."""
    outcomes = list(values)
    if not outcomes:
        return 0.0
    return sum(bool(value) for value in outcomes) / len(outcomes)


def experiment_metrics(success: bool, runtime_seconds: float, candidate_count: int, phase_error: float | None) -> dict[str, float | int | bool | None]:
    """Build a normalized metric row for experiment logging."""
    if runtime_seconds < 0:
        raise ValueError("runtime_seconds must be non-negative.")
    return {
        "success": bool(success),
        "runtime_seconds": float(runtime_seconds),
        "candidate_count": int(candidate_count),
        "phase_error": phase_error,
    }


def self_test() -> None:
    """Verify metric edge cases."""
    assert circular_phase_error(Fraction(7, 8), Fraction(1, 8)) == 0.25
    assert success_rate([True, False, True]) == 2 / 3


if __name__ == "__main__":
    self_test()
    print("metrics self-test passed")
