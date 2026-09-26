"""Decoder D2: grid likelihood / posterior over phi from all observed shots.

This is NOT part of AWQPE. It is included as a crossed experimental factor so
that an improvement from a better decoder is never attributed to a better
allocation policy, and because it is the principled way to fuse evidence from
overlapping blocks: every shot enters the log-likelihood exactly once, and
blocks are separately executed circuits, so their outcomes are conditionally
independent given phi. Only model quantities (no truth) are used.

Grid: phi_g = g / G, G = 2^(n + refine). A block at offset k sees
delta = frac(2^k phi_g) = ((g * 2^k) mod G) / G, so its per-outcome
log-probabilities are tabulated once on the G-point delta grid and gathered.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

from research.awqpe.blocks.geometry import BlockSpec
from research.awqpe.model.noise import NoiseSpec, noisy_block_probabilities

LOG_FLOOR = 1e-300
MAX_TABLE_ELEMENTS = 1 << 24  # above this (G x 2^w) the sparse observed-outcome path is used


@lru_cache(maxsize=128)
def _log_table(G: int, offset: int, width: int, noise: NoiseSpec | None) -> np.ndarray:
    """log p(y | phi_g) for all grid points g, shape (G, 2^width)."""
    delta_grid = np.arange(G) / G
    # frac(2^k * (delta / 2^k)) = delta exactly, and offset-dependent noise is
    # evaluated at the block's true offset.
    p = noisy_block_probabilities(np.ldexp(delta_grid, -offset), offset, width, noise)
    logp = np.log(np.maximum(p, LOG_FLOOR))
    idx = (np.arange(G, dtype=np.int64) << offset) % G
    return np.ascontiguousarray(logp[idx])


def _sparse_loglik(G: int, spec: BlockSpec, counts: np.ndarray, noise: NoiseSpec | None) -> np.ndarray:
    """Same log-likelihood as the dense table, evaluated only at observed outcomes.

    For wide blocks (2^w large) with few shots most outcomes have zero count;
    this path bounds memory by the number of distinct observed outcomes.
    Ideal kernel only (noise channels mix outcomes, so they need the dense path).
    """
    if noise is not None and not noise.is_ideal:
        raise NotImplementedError("sparse likelihood path supports the ideal kernel only.")
    from research.awqpe.model.kernel import kernel_closed_form

    M = spec.M
    delta = ((np.arange(G, dtype=np.int64) << spec.offset) % G) / G
    out = np.zeros((counts.shape[0], G))
    for y in np.flatnonzero(counts.sum(axis=0)):
        logk = np.log(np.maximum(kernel_closed_form(delta - y / M, M), LOG_FLOOR))
        nz = counts[:, y] > 0
        out[nz] += counts[nz, y][:, None] * logk[None, :]
    return out


class GridPosterior:
    """Running log-likelihood over the phi grid for T independent runs."""

    def __init__(self, n: int, T: int = 1, refine: int = 4, noise: NoiseSpec | None = None):
        self.n, self.refine, self.noise = int(n), int(refine), noise
        self.G = 1 << (self.n + self.refine)
        self.loglik = np.zeros((T, self.G))

    def add_counts(self, spec: BlockSpec, counts) -> None:
        """Add counts of shape (2^w,) or (T, 2^w) observed on block `spec`."""
        c = np.atleast_2d(np.asarray(counts, dtype=float))
        if self.G * spec.M <= MAX_TABLE_ELEMENTS:
            self.loglik += c @ _log_table(self.G, spec.offset, spec.width, self.noise).T
        else:
            self.loglik += _sparse_loglik(self.G, spec, c, self.noise)

    def posterior(self) -> np.ndarray:
        z = self.loglik - self.loglik.max(axis=1, keepdims=True)
        w = np.exp(z)
        return w / w.sum(axis=1, keepdims=True)

    def map_estimate(self) -> np.ndarray:
        return np.argmax(self.loglik, axis=1) / self.G

    def credible_mass(self, center, tau: float) -> np.ndarray:
        """Posterior probability that |phi - center|_circle <= tau."""
        post = self.posterior()
        grid = np.arange(self.G) / self.G
        d = np.abs(np.mod(grid[None, :] - np.atleast_1d(np.asarray(center, dtype=float))[:, None], 1.0))
        d = np.minimum(d, 1.0 - d)
        return np.sum(np.where(d <= tau + 1e-15, post, 0.0), axis=1)


def likelihood_decode(block_specs, block_counts, n: int, refine: int = 4, noise: NoiseSpec | None = None) -> np.ndarray:
    """MAP phase for T runs; block_counts[i] has shape (T, 2^w_i)."""
    T = np.asarray(block_counts[0]).shape[0]
    gp = GridPosterior(n, T, refine, noise)
    for spec, c in zip(block_specs, block_counts):
        gp.add_counts(spec, c)
    return gp.map_estimate()
