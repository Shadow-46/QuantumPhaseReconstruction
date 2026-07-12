"""
Purpose
    Provide the standard Shor order-finding baseline.
Theory
    Shor succeeds by finding the multiplicative order r of a modulo N and then
    extracting factors from gcd(a^(r/2) +/- 1, N). The baseline keeps quantum
    circuit construction separate from exact classical validation utilities.
Inputs
    Composite modulus N and coprime base a.
Outputs
    StandardShorResult containing order, factors, and success status.
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
from math import gcd

from Circuits.modular_multiplication import ModularMultiplicationOperator
from Circuits.qpe import build_qpe_circuit
from Reconstruction.continued_fraction import factors_from_order


@dataclass(frozen=True)
class StandardShorResult:
    """Result of standard order finding for one (N, a) pair."""
    N: int
    a: int
    order: int | None
    factors: tuple[int, int] | None

    @property
    def success(self) -> bool:
        """Return True when non-trivial factors were recovered."""
        return self.factors is not None


def multiplicative_order(a: int, N: int) -> int | None:
    """Compute the multiplicative order of a modulo N by exact validation."""
    if gcd(a, N) != 1:
        return None
    value = 1
    for r in range(1, N + 1):
        value = (value * a) % N
        if value == 1:
            return r
    return None


def build_standard_shor_qpe_circuit(N: int, a: int, phase_qubits: int):
    """Build the standard QPE circuit for modular multiplication by a modulo N."""
    op = ModularMultiplicationOperator(a, N)
    return build_qpe_circuit(op.gate(), phase_qubits, op.qubits, eigenstate=1, measure=True)


def run_standard_shor(N: int, a: int) -> StandardShorResult:
    """Run exact standard Shor post-processing for a small reproducible instance."""
    divisor = gcd(a, N)
    if 1 < divisor < N:
        return StandardShorResult(N, a, None, tuple(sorted((divisor, N // divisor))))
    order = multiplicative_order(a, N)
    return StandardShorResult(N, a, order, factors_from_order(a, N, order))


def self_test() -> None:
    """Verify the required N=15 and N=21 standard Shor cases."""
    assert run_standard_shor(15, 2).factors == (3, 5)
    assert run_standard_shor(21, 2).factors == (3, 7)


if __name__ == "__main__":
    self_test()
    print("standard_shor self-test passed")
