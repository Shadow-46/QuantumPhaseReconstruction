"""
Purpose
    Build the modular multiplication unitary U_a |y> = |a y mod N>.
Theory
    Shor phase estimation uses eigenphases of the permutation induced by
    multiplication by a modulo N. For gcd(a, N) = 1 this map is reversible on
    residues modulo N; computational states outside the modulus range are fixed
    to embed the operation in a power-of-two Hilbert space.
Inputs
    Integers a and N, with gcd(a, N) = 1.
Outputs
    NumPy matrices, Qiskit operators, gates, and verification diagnostics.
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
from math import ceil, gcd, log2

import numpy as np
from qiskit.circuit import Gate
from qiskit.circuit.library import UnitaryGate
from qiskit.quantum_info import Operator


def work_qubits_for_modulus(N: int) -> int:
    """Return the minimum number of work qubits needed to represent 0..N-1."""

    if N < 2:
        raise ValueError("N must be at least 2.")
    return ceil(log2(N))


@dataclass(frozen=True)
class ModularMultiplicationOperator:
    """Permutation representation of modular multiplication by a modulo N."""

    a: int
    N: int
    num_work_qubits: int | None = None

    def __post_init__(self) -> None:
        """Validate the modular multiplication parameters."""

        if self.a <= 0 or self.N <= 1:
            raise ValueError("Require a > 0 and N > 1.")
        if gcd(self.a, self.N) != 1:
            raise ValueError("Modular multiplication is unitary only when gcd(a, N) = 1.")
        if self.num_work_qubits is not None and 2**self.num_work_qubits < self.N:
            raise ValueError("num_work_qubits is too small for N.")

    @property
    def qubits(self) -> int:
        """Return the number of work qubits."""

        return self.num_work_qubits if self.num_work_qubits is not None else work_qubits_for_modulus(self.N)

    @property
    def dimension(self) -> int:
        """Return the matrix dimension of the embedded operation."""

        return 2**self.qubits

    def mapping(self) -> dict[int, int]:
        """Return the computational-basis permutation implemented by U_a."""

        return {
            y: (self.a * y) % self.N if y < self.N else y
            for y in range(self.dimension)
        }

    def matrix(self) -> np.ndarray:
        """Return the dense permutation matrix for modular multiplication."""

        mat = np.zeros((self.dimension, self.dimension), dtype=complex)
        for src, dst in self.mapping().items():
            mat[dst, src] = 1.0
        return mat

    def power_matrix(self, exponent: int) -> np.ndarray:
        """Return the dense matrix for U_a raised to a non-negative exponent."""

        if exponent < 0:
            raise ValueError("exponent must be non-negative.")
        powered_a = pow(self.a, exponent, self.N)
        return ModularMultiplicationOperator(powered_a, self.N, self.qubits).matrix()

    def operator(self) -> Operator:
        """Return a Qiskit Operator for U_a."""

        return Operator(self.matrix())

    def gate(self, label: str | None = None) -> Gate:
        """Return U_a as a Qiskit gate."""

        gate_label = label if label is not None else f"M_{self.a} mod {self.N}"
        return UnitaryGate(self.matrix(), label=gate_label)

    def controlled_power_gate(self, exponent: int) -> Gate:
        """Return a controlled gate for U_a raised to exponent."""

        powered_a = pow(self.a, exponent, self.N)
        gate = ModularMultiplicationOperator(powered_a, self.N, self.qubits).gate(
            label=f"M_{powered_a} mod {self.N}"
        )
        return gate.control(1)

    def is_permutation(self) -> bool:
        """Check that every matrix row and column has exactly one nonzero entry."""

        mat = self.matrix()
        return bool(np.allclose(mat.sum(axis=0), 1) and np.allclose(mat.sum(axis=1), 1))

    def is_unitary(self) -> bool:
        """Check the unitary condition U dagger U = I."""

        mat = self.matrix()
        return bool(np.allclose(mat.conj().T @ mat, np.eye(self.dimension)))


def self_test() -> None:
    """Verify reversibility and unitarity for representative Shor moduli."""

    for N, a in [(15, 2), (21, 2)]:
        op = ModularMultiplicationOperator(a, N)
        assert op.is_permutation()
        assert op.is_unitary()
        assert op.matrix().shape == (2**op.qubits, 2**op.qubits)


if __name__ == "__main__":
    self_test()
    print("modular_multiplication self-test passed")
