"""Evaluator-side metrics. These functions take the true phase and must never
be imported by allocation, overlap, stopping or controller code (enforced by
research/tests/test_leakage.py)."""

from __future__ import annotations

import numpy as np


def circular_error(estimate, phi) -> np.ndarray:
    """Distance on the unit circle between phase estimates and truth (in turns)."""
    d = np.mod(np.asarray(estimate, dtype=float) - np.asarray(phi, dtype=float), 1.0)
    return np.minimum(d, 1.0 - d)


def best_nbit_integer(phi, n: int) -> np.ndarray:
    """AWQPE's target: floor(2^n phi + 0.5) mod 2^n (Theorem 3.7)."""
    return np.mod(np.floor(np.ldexp(np.asarray(phi, dtype=float), n) + 0.5), 1 << n).astype(np.int64)


def final_residual(phi, n: int) -> np.ndarray:
    """frac(2^n phi): 0.5 means the best n-bit approximation itself is a tie."""
    return np.mod(np.ldexp(np.asarray(phi, dtype=float), n), 1.0)


def exact_success(estimate_int, phi, n: int, tie_tol: float = 1e-12) -> np.ndarray:
    """Estimate equals the best n-bit approximation.

    When frac(2^n phi) is within tie_tol of 0.5 both neighbouring grid points
    are equally best, and either counts as success.
    """
    est = np.asarray(estimate_int, dtype=np.int64)
    best = best_nbit_integer(phi, n)
    tie = np.abs(final_residual(phi, n) - 0.5) < tie_tol
    lower = np.mod(best - 1, 1 << n)
    return (est == best) | (tie & (est == lower))


def boundary_remainders(phi, widths) -> np.ndarray:
    """delta_j = frac(2^(K_j) phi) at every internal chunk boundary K_j.

    delta_j near 0.5 is the regime where chunk j's rounding is ambiguous and
    chunk j+1 reads close to 10..0 (the special-chunk case). Shape (..., B-1).
    """
    cuts = np.cumsum(widths)[:-1]
    phi = np.asarray(phi, dtype=float)
    if len(cuts) == 0:
        return np.zeros(phi.shape + (0,))
    return np.stack([np.mod(np.ldexp(phi, int(K)), 1.0) for K in cuts], axis=-1)


def boundary_hardness(phi, widths) -> np.ndarray:
    """min_j |delta_j - 0.5| over internal boundaries (0 = hardest); 0.5 if B = 1."""
    rem = boundary_remainders(phi, widths)
    if rem.shape[-1] == 0:
        return np.full(np.shape(phi), 0.5)
    return np.min(np.abs(rem - 0.5), axis=-1)
