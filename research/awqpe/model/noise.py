"""Analytic (Dirichlet-level) noise channels for an AWQPE block.

These are phenomenological channels applied to the block's outcome
distribution, used for large Monte Carlo campaigns. Their fidelity to
circuit-level noise is checked against Qiskit Aer in Phase 10; nothing here is
claimed to be the noise model of the published AWQPE paper (see
docs/PAPER_VERSION_NOTES.md).

Channels (applied in this order, all linear in the ideal kernel so the same
maps act on derivatives):

1. Phase jitter: each shot sees theta + xi, xi ~ N(0, s^2). In Fejer form this
   multiplies frequency d by g_d = exp(-2 pi^2 d^2 s^2), exactly. With
   jitter_mode="scaled", s = sigma * 2^k (jitter on phi itself, amplified by
   the block's power); with "fixed", s = sigma for every block.
2. Global depolarisation: p -> (1 - lam) p + lam / M, with
   lam = 1 - exp(-gamma * q), q = controlled-U applications per shot
   (depol_mode="per_query") or lam = gamma (depol_mode="fixed").
3. Readout: independent per-bit confusion, P(read 1 | 0) = p01,
   P(read 0 | 1) = p10.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from research.awqpe.model import kernel as K


@dataclass(frozen=True)
class NoiseSpec:
    """Declarative noise configuration; all-zero means ideal."""

    jitter_sigma: float = 0.0
    jitter_mode: str = "fixed"  # "fixed" | "scaled"
    depol_gamma: float = 0.0
    depol_mode: str = "per_query"  # "per_query" | "fixed"
    readout_p01: float = 0.0
    readout_p10: float = 0.0

    @property
    def is_ideal(self) -> bool:
        return self.jitter_sigma == 0 and self.depol_gamma == 0 and self.readout_p01 == 0 and self.readout_p10 == 0

    def validate(self) -> None:
        if self.jitter_mode not in {"fixed", "scaled"}:
            raise ValueError("jitter_mode must be 'fixed' or 'scaled'.")
        if self.depol_mode not in {"per_query", "fixed"}:
            raise ValueError("depol_mode must be 'per_query' or 'fixed'.")
        if self.jitter_sigma < 0 or self.depol_gamma < 0:
            raise ValueError("noise strengths must be non-negative.")
        if self.depol_mode == "fixed" and self.depol_gamma > 1:
            raise ValueError("fixed depolarising strength must be <= 1.")
        for p in (self.readout_p01, self.readout_p10):
            if not 0 <= p <= 0.5:
                raise ValueError("readout flip probabilities must lie in [0, 0.5].")


def u_queries_per_shot(offset: int, width: int) -> int:
    """Controlled-U applications in one execution of a block: 2^k (2^m - 1)."""
    return (1 << offset) * ((1 << width) - 1)


def jitter_damping(noise: NoiseSpec, offset: int, width: int) -> np.ndarray | None:
    if noise.jitter_sigma == 0:
        return None
    s = noise.jitter_sigma * (float(1 << offset) if noise.jitter_mode == "scaled" else 1.0)
    d = np.arange(1, 1 << width, dtype=float)
    return np.exp(-2.0 * np.pi**2 * d**2 * s**2)


def depolarising_lambda(noise: NoiseSpec, offset: int, width: int) -> float:
    if noise.depol_gamma == 0:
        return 0.0
    if noise.depol_mode == "fixed":
        return float(noise.depol_gamma)
    return float(-np.expm1(-noise.depol_gamma * u_queries_per_shot(offset, width)))


def apply_readout(p: np.ndarray, width: int, p01: float, p10: float) -> np.ndarray:
    """Apply independent per-bit readout confusion along the last axis."""
    if p01 == 0 and p10 == 0:
        return p
    c = np.array([[1.0 - p01, p10], [p01, 1.0 - p10]])  # c[read, true]
    lead = p.shape[:-1]
    t = p.reshape(lead + (2,) * width)
    for axis in range(width):
        ax = len(lead) + axis
        t = np.moveaxis(np.tensordot(t, c, axes=([ax], [1])), -1, ax)
    return t.reshape(lead + (1 << width,))


def _linear_noise(x: np.ndarray, noise: NoiseSpec, offset: int, width: int, constant: bool) -> np.ndarray:
    lam = depolarising_lambda(noise, offset, width)
    if lam:
        x = (1.0 - lam) * x + (lam / (1 << width) if constant else 0.0)
    return apply_readout(x, width, noise.readout_p01, noise.readout_p10)


def noisy_block_probabilities(phi, offset: int, width: int, noise: NoiseSpec | None = None) -> np.ndarray:
    """Outcome distribution of one block under `noise` (ideal if None)."""
    if noise is None or noise.is_ideal:
        return K.block_probabilities(phi, offset, width)
    noise.validate()
    p = K.block_probabilities(phi, offset, width, jitter_damping(noise, offset, width))
    p = _linear_noise(p, noise, offset, width, constant=True)
    p = np.clip(p, 0.0, None)
    return p / p.sum(axis=-1, keepdims=True)


def noisy_block_derivatives(phi, offset: int, width: int, noise: NoiseSpec | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(p, dp/dphi, d2p/dphi2) under noise; derivatives of the constant vanish."""
    damping = None if noise is None else jitter_damping(noise, offset, width)
    p, d1, d2 = K.block_probability_derivatives(phi, offset, width, damping)
    if noise is None or noise.is_ideal:
        return p, d1, d2
    return (
        _linear_noise(p, noise, offset, width, constant=True),
        _linear_noise(d1, noise, offset, width, constant=False),
        _linear_noise(d2, noise, offset, width, constant=False),
    )
