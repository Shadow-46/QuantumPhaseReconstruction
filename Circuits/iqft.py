"""
Purpose
    Construct an inverse quantum Fourier transform without deprecated helpers.
Theory
    QPE converts phase kickback into a computational-basis estimate by applying
    the inverse QFT on the phase register. This module expands the controlled
    phase and Hadamard decomposition explicitly.
Inputs
    Number of phase-register qubits and an optional swap policy.
Outputs
    A Qiskit QuantumCircuit implementing IQFT and numerical verification data.
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator


def inverse_qft_circuit(num_qubits: int, do_swaps: bool = True) -> QuantumCircuit:
    """Build an inverse QFT circuit over num_qubits qubits."""
    if num_qubits < 1:
        raise ValueError("num_qubits must be positive.")
    circuit = QuantumCircuit(num_qubits, name=f"IQFT_{num_qubits}")
    if do_swaps:
        for qubit in range(num_qubits // 2):
            circuit.swap(qubit, num_qubits - qubit - 1)
    for target in range(num_qubits):
        for control in range(target):
            circuit.cp(-np.pi / (2 ** (target - control)), control, target)
        circuit.h(target)
    return circuit


def dense_inverse_qft(num_qubits: int) -> np.ndarray:
    """Return the analytical inverse QFT matrix in Qiskit's basis ordering."""
    dim = 2**num_qubits
    omega = np.exp(-2j * np.pi / dim)
    mat = np.empty((dim, dim), dtype=complex)
    for row in range(dim):
        for col in range(dim):
            mat[row, col] = omega ** (row * col) / np.sqrt(dim)
    return mat


def verify_decomposition(num_qubits: int, atol: float = 1e-9) -> bool:
    """Compare the circuit decomposition with the analytical IQFT matrix."""
    circuit_matrix = Operator(inverse_qft_circuit(num_qubits)).data
    return bool(np.allclose(circuit_matrix, dense_inverse_qft(num_qubits), atol=atol))


def self_test() -> None:
    """Verify IQFT decompositions for small phase registers."""
    for qubits in (1, 2, 3, 4):
        assert verify_decomposition(qubits)


if __name__ == "__main__":
    self_test()
    print("iqft self-test passed")
