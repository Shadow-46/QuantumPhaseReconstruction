"""
Purpose
    Provide Aer backend construction and count execution helpers.
Theory
    Simulation is an infrastructure concern separate from circuit construction
    and classical reconstruction. This module owns transpilation, seeding, and
    backend execution.
Inputs
    QuantumCircuit objects, shot counts, simulator method, seed, and noise model.
Outputs
    Measurement count dictionaries.
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from typing import Mapping

from qiskit import QuantumCircuit, transpile


def aer_simulator(method: str = "automatic", seed: int | None = None, noise_model=None):
    """Return a configured Qiskit AerSimulator."""
    from qiskit_aer import AerSimulator
    return AerSimulator(method=method, seed_simulator=seed, noise_model=noise_model)


def run_counts(circuit: QuantumCircuit, shots: int, method: str = "automatic", seed: int | None = None, noise_model=None) -> dict[str, int]:
    """Execute a measured circuit on Aer and return counts."""
    if shots < 1:
        raise ValueError("shots must be positive.")
    backend = aer_simulator(method=method, seed=seed, noise_model=noise_model)
    compiled = transpile(circuit, backend)
    result = backend.run(compiled, shots=shots).result()
    counts: Mapping[str, int] = result.get_counts()
    return {key.replace(" ", ""): int(value) for key, value in counts.items()}


def self_test() -> None:
    """Verify backend execution on a one-qubit deterministic circuit."""
    qc = QuantumCircuit(1, 1)
    qc.x(0)
    qc.measure(0, 0)
    counts = run_counts(qc, shots=16, seed=7)
    assert counts == {"1": 16}


if __name__ == "__main__":
    self_test()
    print("backend self-test passed")
