"""Posterior-credible stopping for a declared phase tolerance (Phase 8, D-028).

Decision code: sees only the running grid log-likelihood, never the true phase
(enforced by research/tests/test_leakage.py).

Rule: stop when P(|phi - phi_MAP|_circle <= tau | data) >= 1 - alpha, where the
posterior is the uniform-prior grid posterior of decode.likelihood.GridPosterior.
`credible_mass_at_map` is an O(T G) prefix-sum version of
GridPosterior.credible_mass(map_estimate(), tau); equality is tested.
"""

from __future__ import annotations

import numpy as np


def tau_bits(tau: float) -> int:
    """Leading bits K with 2^-(K+1) <= tau: rounding to K bits already meets tolerance tau."""
    return max(1, int(np.ceil(-np.log2(tau) - 1e-12)) - 1)


def truncated_widths(widths, tau: float) -> list[int]:
    """Shortest MSB-first block prefix carrying at least tau_bits(tau) bits (tau-matched truncation)."""
    need, out, K = tau_bits(tau), [], 0
    for m in widths:
        out.append(int(m))
        K += int(m)
        if K >= need:
            break
    return out


def credible_mass_at_map(loglik: np.ndarray, tau: float) -> tuple[np.ndarray, np.ndarray]:
    """(mass, map_index) for each row of a (T, G) log-likelihood on the grid i / G."""
    T, G = loglik.shape
    rows = np.arange(T)
    c = np.argmax(loglik, axis=1)
    post = np.exp(loglik - loglik[rows, c][:, None])
    post /= post.sum(axis=1, keepdims=True)
    h = int(np.floor(tau * G + 1e-9))  # grid points within tau of the MAP on each side
    if 2 * h + 1 >= G:
        return np.ones(T), c
    cs = np.concatenate([np.zeros((T, 1)), np.cumsum(post, axis=1)], axis=1)
    lo, hi = c - h, c + h + 1  # half-open [lo, hi) modulo G
    zero, full = np.zeros_like(lo), np.full_like(lo, G)

    def seg(a, b):  # sum of post[:, a:b] for 0 <= a <= b <= G
        return cs[rows, b] - cs[rows, a]

    mass = np.where(lo < 0, seg(zero, np.clip(hi, 0, G)) + seg(np.clip(lo + G, 0, G), full),
                    np.where(hi > G, seg(np.clip(lo, 0, G), full) + seg(zero, np.clip(hi - G, 0, G)),
                             seg(np.clip(lo, 0, G), np.clip(hi, 0, G))))
    return np.clip(mass, 0.0, 1.0), c


def should_stop(loglik: np.ndarray, tau: float, alpha: float) -> tuple[np.ndarray, np.ndarray]:
    """(stop mask, credible mass) for the rule mass >= 1 - alpha."""
    mass, _ = credible_mass_at_map(loglik, tau)
    return mass >= 1.0 - alpha, mass
