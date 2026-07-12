"""
Purpose
    Build Qiskit Aer noise models from explicit scalar parameters.
Theory
    Reconstruction robustness experiments vary physical error channels while
    preserving the same circuit and classical recovery logic.
Inputs
    NoiseConfig containing depolarizing and measurement error probabilities.
Outputs
    Aer NoiseModel objects or None for noiseless simulation.
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from config import NoiseConfig


def build_noise_model(config: NoiseConfig):
    """Create an Aer NoiseModel from the supplied noise configuration."""
    if not config.enabled:
        return None
    from qiskit_aer.noise import NoiseModel, ReadoutError, depolarizing_error
    model = NoiseModel()
    if config.depolarizing_1q:
        error = depolarizing_error(config.depolarizing_1q, 1)
        model.add_all_qubit_quantum_error(error, ["h", "x", "rz", "sx"])
    if config.depolarizing_2q:
        error = depolarizing_error(config.depolarizing_2q, 2)
        model.add_all_qubit_quantum_error(error, ["cx", "cp", "swap"])
    if config.measurement_error:
        p = config.measurement_error
        model.add_all_qubit_readout_error(ReadoutError([[1 - p, p], [p, 1 - p]]))
    return model


def self_test() -> None:
    """Verify disabled noise produces no model and enabled noise is constructible."""
    assert build_noise_model(NoiseConfig()) is None
    assert build_noise_model(NoiseConfig(enabled=True, measurement_error=0.01)) is not None


if __name__ == "__main__":
    self_test()
    print("noise self-test passed")
