"""Squared Dirichlet-kernel measurement model of an AWQPE block.

AWQPE (arXiv:2507.22460v3, Sec. 3.2) models the outcome of an m-control-qubit
block as

    P(j | delta) = (1 / 2^(2m)) * sin^2(2^m pi theta) / sin^2(pi theta),
    theta = delta - j / 2^m,

where delta is the fractional phase the block sees. A block whose controls
apply U^(2^(k+p)), p = 0..m-1, sees delta = frac(2^k phi) (Algorithm 1,
lines 8-9).

Implementation note. With M = 2^m the kernel equals the Fejer form

    K_M(theta) = M^-2 * sum_{|d|<M} (M - |d|) cos(2 pi d theta),

a trigonometric polynomial with no removable singularity at theta in Z. We
evaluate that form, which gives exact analytic first and second derivatives
(needed for Fisher information) and lets phase-jitter noise act as a
per-frequency damping factor g_d (see model/noise.py). The closed sin^2 ratio
is kept as `kernel_closed_form` purely for cross-validation.
"""

from __future__ import annotations

import numpy as np

TWO_PI = 2.0 * np.pi


def _fourier_terms(M: int, damping: np.ndarray | None) -> tuple[np.ndarray, np.ndarray]:
    """Return frequencies d = 1..M-1 and weights (M-d)/M^2 * g_d."""
    if M < 2:
        raise ValueError("M must be at least 2.")
    d = np.arange(1, M, dtype=float)
    w = (M - d) / float(M * M)
    if damping is not None:
        damping = np.asarray(damping, dtype=float)
        if damping.shape != (M - 1,):
            raise ValueError("damping must have shape (M-1,).")
        w = w * damping
    return d, w


def kernel(theta, M: int, damping: np.ndarray | None = None) -> np.ndarray:
    """K_M(theta) in Fejer form; broadcast over any theta shape."""
    theta = np.asarray(theta, dtype=float)
    d, w = _fourier_terms(M, damping)
    phase = TWO_PI * theta[..., None] * d
    return 1.0 / M + 2.0 * np.sum(w * np.cos(phase), axis=-1)


def kernel_d1(theta, M: int, damping: np.ndarray | None = None) -> np.ndarray:
    """dK_M/dtheta, exact."""
    theta = np.asarray(theta, dtype=float)
    d, w = _fourier_terms(M, damping)
    phase = TWO_PI * theta[..., None] * d
    return -2.0 * np.sum(w * TWO_PI * d * np.sin(phase), axis=-1)


def kernel_d2(theta, M: int, damping: np.ndarray | None = None) -> np.ndarray:
    """d^2K_M/dtheta^2, exact."""
    theta = np.asarray(theta, dtype=float)
    d, w = _fourier_terms(M, damping)
    phase = TWO_PI * theta[..., None] * d
    return -2.0 * np.sum(w * (TWO_PI * d) ** 2 * np.cos(phase), axis=-1)


def kernel_closed_form(theta, M: int) -> np.ndarray:
    """The paper's sin^2 ratio, with the removable singularity set to 1."""
    theta = np.asarray(theta, dtype=float)
    s = np.sin(np.pi * theta)
    num = np.sin(M * np.pi * theta) ** 2
    out = np.empty_like(theta)
    small = np.abs(s) < 1e-12
    out[~small] = num[~small] / (M * M * s[~small] ** 2)
    out[small] = 1.0
    return out


def block_delta(phi, offset: int) -> np.ndarray:
    """Fractional phase frac(2^offset * phi) seen by a block at `offset`."""
    return np.mod(np.ldexp(np.asarray(phi, dtype=float), offset), 1.0)


def block_thetas(phi, offset: int, width: int) -> np.ndarray:
    """theta_y = frac(2^k phi) - y/2^m for every outcome y (last axis)."""
    M = 1 << width
    delta = block_delta(phi, offset)
    return delta[..., None] - np.arange(M) / M


def block_probabilities(phi, offset: int, width: int, damping: np.ndarray | None = None) -> np.ndarray:
    """Ideal outcome distribution p(y | phi) of one block; last axis is y.

    Outcome y is the integer read from the m measured control bits with the
    same convention as Qiskit's counts keys (verified against statevector in
    research/tests/test_kernel.py).
    """
    M = 1 << width
    p = kernel(block_thetas(phi, offset, width), M, damping)
    # Clip float round-off (values of order -1e-17 at exact zeros) and renormalise.
    p = np.clip(p, 0.0, None)
    return p / p.sum(axis=-1, keepdims=True)


def block_probability_derivatives(phi, offset: int, width: int, damping: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (p, dp/dphi, d2p/dphi2) for one block, unclipped.

    Chain rule: theta = 2^k phi - y/M, so d/dphi = 2^k d/dtheta.
    """
    M = 1 << width
    th = block_thetas(phi, offset, width)
    scale = float(1 << offset)
    return (
        kernel(th, M, damping),
        scale * kernel_d1(th, M, damping),
        scale * scale * kernel_d2(th, M, damping),
    )
