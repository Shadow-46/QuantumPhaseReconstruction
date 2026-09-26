"""Classical Fisher information of AWQPE blocks with respect to phi.

    I_b(phi) = sum_y (dp_y/dphi)^2 / p_y .

For the ideal kernel, p_y has exact double zeros whenever 2^m theta_y is a
nonzero integer. Near such a zero p ~ c x^2 and dp ~ 2 c x, so the term tends
to 4c = 2 d2p/dphi2, a finite nonzero value. Terms whose p_y falls below
`zero_tol` use that limit instead of dividing by round-off. Fisher
information is a local quantity; whether it predicts where extra shots help is
an empirical question (Phase 4), not an assumption.
"""

from __future__ import annotations

import numpy as np

from research.awqpe.model.noise import NoiseSpec, noisy_block_derivatives, noisy_block_probabilities


def fisher_from_derivatives(p: np.ndarray, d1: np.ndarray, d2: np.ndarray, zero_tol: float = 1e-10) -> np.ndarray:
    """Sum over the last axis with the double-zero limit handled analytically."""
    small = p < zero_tol
    safe = np.where(small, 1.0, p)
    terms = np.where(small, 2.0 * d2, d1 * d1 / safe)
    # A zero that is not a double zero would make the limit negative; clamp so
    # information is never negative.
    return np.sum(np.clip(terms, 0.0, None), axis=-1)


def block_fisher(phi, offset: int, width: int, noise: NoiseSpec | None = None) -> np.ndarray:
    """Per-shot Fisher information of one block about phi."""
    return fisher_from_derivatives(*noisy_block_derivatives(phi, offset, width, noise))


def block_fisher_finite_difference(phi, offset: int, width: int, noise: NoiseSpec | None = None, h: float = 1e-7) -> np.ndarray:
    """Central-difference Fisher information, for validation only.

    Terms with p_y below 1e-9 are dropped because the difference quotient is
    dominated by round-off there; compare against block_fisher only at phases
    away from exact kernel zeros.
    """
    phi = np.asarray(phi, dtype=float)
    p = noisy_block_probabilities(phi, offset, width, noise)
    dp = (noisy_block_probabilities(phi + h, offset, width, noise) - noisy_block_probabilities(phi - h, offset, width, noise)) / (2 * h)
    keep = p > 1e-9
    return np.sum(np.where(keep, dp * dp / np.where(keep, p, 1.0), 0.0), axis=-1)


def quantum_fisher_bound(offset: int, width: int) -> float:
    """QFI of the pre-measurement control state about phi (ideal, pure).

    The control register holds M^-1/2 sum_t e^{2 pi i 2^k t phi}|t>, whose QFI
    is 4 Var(2 pi 2^k t) = 4 (2 pi 2^k)^2 (M^2 - 1)/12 for t uniform on 0..M-1.
    Every measurement's classical Fisher information is bounded by it.
    """
    M = 1 << width
    return 4.0 * (2.0 * np.pi * (1 << offset)) ** 2 * (M * M - 1) / 12.0
