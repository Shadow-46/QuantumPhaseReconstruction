"""
Purpose
    Provide a clean adaptive extension point for reconstruction-stage research.
Theory
    The baseline paper reproduction is fixed. Adaptive experiments are isolated
    here so future changes can vary only window geometry and candidate budgets
    before delegating to the unchanged reconstruction pipeline.
Inputs
    ShorConfig plus optional minimum and maximum window controls.
Outputs
    WindowedShorResult produced by the selected adaptive configuration.
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

from config import ShorConfig, WindowConfig
from Algorithms.paper_algorithm import WindowedShorResult, run_windowed_shor


@dataclass(frozen=True)
class AdaptivePolicy:
    """Deterministic policy for choosing reconstruction window parameters."""
    min_window_size: int = 3
    max_window_size: int = 6
    overlap_fraction: float = 0.5
    candidate_count: int = 4

    def choose(self, precision: int) -> WindowConfig:
        """Choose a valid WindowConfig for the requested phase precision."""
        if precision < 1:
            raise ValueError("precision must be positive.")
        width = min(self.max_window_size, max(self.min_window_size, precision // 2 or 1), precision)
        overlap = min(width - 1, max(0, int(width * self.overlap_fraction)))
        return WindowConfig(total_precision=precision, window_size=width, overlap=overlap, candidate_count=self.candidate_count)


def run_adaptive_shor(shor: ShorConfig, policy: AdaptivePolicy | None = None) -> WindowedShorResult:
    """Run the unchanged windowed algorithm with policy-selected parameters."""
    active_policy = AdaptivePolicy() if policy is None else policy
    return run_windowed_shor(shor, active_policy.choose(shor.phase_qubits))


def self_test() -> None:
    """Verify the adaptive entry point delegates cleanly to reconstruction."""
    assert run_adaptive_shor(ShorConfig(N=15, a=2, phase_qubits=8)).factors == (3, 5)


if __name__ == "__main__":
    self_test()
    print("adaptive_algorithm self-test passed")
