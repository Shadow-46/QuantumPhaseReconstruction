"""
Purpose
    Central configuration for QuantumPhaseReconstruction experiments.
Theory
    Reproducible phase-estimation studies require every experimental degree of
    freedom to be recorded: modulus, base, precision, window geometry, backend,
    shots, and reconstruction search limits.
Inputs
    Imported by command-line entry points, algorithms, and experiment drivers.
Outputs
    Frozen dataclass instances and conservative default constants.
Author
    Sanjay
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent
DATA_DIR = PROJECT_ROOT / "Data"
DEFAULT_RANDOM_SEED = 314159


@dataclass(frozen=True)
class ShorConfig:
    """Configuration for a single Shor/QPE reconstruction run."""

    N: int = 15
    a: int = 2
    phase_qubits: int = 8
    shots: int = 2048
    simulator_method: str = "automatic"
    random_seed: int = DEFAULT_RANDOM_SEED


@dataclass(frozen=True)
class WindowConfig:
    """Windowed phase reconstruction parameters."""

    total_precision: int = 8
    window_size: int = 4
    overlap: int = 2
    candidate_count: int = 4
    max_paths: int = 32


@dataclass(frozen=True)
class NoiseConfig:
    """Noise model parameters used by Aer simulation experiments."""

    enabled: bool = False
    depolarizing_1q: float = 0.0
    depolarizing_2q: float = 0.0
    measurement_error: float = 0.0


@dataclass(frozen=True)
class AdaptiveReconstructionConfig:
    """Free parameters for confidence-guided adaptive reconstruction
    (Data/failure_study/FORMULATION.md). Placeholder defaults pending the
    cross-validation step in FORMULATION.md Section 7; not tuned guesses."""

    alpha0: float = 0.5
    epsilon: float = 0.05
    delta_cov: float = 0.05
    # Bounded to stay well above every shots value the combined-stress grid
    # actually starves (4-16) while keeping per-trial resampling cost
    # tractable: each doubling round is a real Aer transpile+run per window
    # (~3.3s measured), and 65536 (14 doublings) made a single D-arm trial
    # take ~120s; 64 (4 doublings from shots=4) keeps 16x headroom over the
    # starved grid's shots values while bounding worst-case per-window cost.
    s_max: int = 64
    beam_max: int = 4096


DEFAULT_SHOR = ShorConfig()
DEFAULT_WINDOW = WindowConfig()
DEFAULT_NOISE = NoiseConfig()
DEFAULT_ADAPTIVE = AdaptiveReconstructionConfig()
