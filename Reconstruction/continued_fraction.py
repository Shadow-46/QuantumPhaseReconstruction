"""
Purpose
    Recover Shor orders and factors from phase estimates.
Theory
    A measured phase approximates s/r, where r is the order of a modulo N.
    Continued fractions recover candidate denominators, which are validated by
    modular exponentiation before factors are extracted with gcd(a^(r/2) +/- 1, N).
Inputs
    Phase fractions or measured bit strings, modulus N, and base a.
Outputs
    Validated orders and non-trivial factors.
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


@dataclass(frozen=True)
class OrderRecoveryResult:
    """Validated order-recovery result for one phase estimate."""
    phase: Fraction
    order: int | None
    factors: tuple[int, int] | None


def order_from_phase(phase: Fraction, a: int, N: int, max_denominator: int | None = None) -> int | None:
    """Recover a valid multiplicative order from a phase estimate."""
    limit = N if max_denominator is None else max_denominator
    frac = phase.limit_denominator(limit)
    for multiplier in range(1, N + 1):
        candidate = frac.denominator * multiplier
        if candidate > 0 and pow(a, candidate, N) == 1:
            return candidate
    return None


def factors_from_order(a: int, N: int, order: int | None) -> tuple[int, int] | None:
    """Extract non-trivial factors of N from a validated even order."""
    if order is None or order % 2:
        return None
    root = pow(a, order // 2, N)
    if root in (1, N - 1):
        return None
    factors = tuple(sorted((gcd(root - 1, N), gcd(root + 1, N))))
    if factors[0] in (1, N) or factors[1] in (1, N):
        return None
    return factors


def recover_order_and_factors(phase: Fraction, a: int, N: int) -> OrderRecoveryResult:
    """Recover both the modular order and Shor factors from a phase estimate."""
    order = order_from_phase(phase, a, N)
    return OrderRecoveryResult(phase=phase, order=order, factors=factors_from_order(a, N, order))


def self_test() -> None:
    """Verify continued-fraction recovery on N=15 and N=21 examples."""
    assert recover_order_and_factors(Fraction(1, 4), 2, 15).factors == (3, 5)
    assert recover_order_and_factors(Fraction(1, 6), 2, 21).order == 6


if __name__ == "__main__":
    self_test()
    print("continued_fraction self-test passed")
