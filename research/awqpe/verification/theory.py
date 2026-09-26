"""Infinite-shot failure theory for AWQPE Algorithms 1-2 (arXiv-v3 reading).

PRELIMINARY / arXiv-v3 result: a candidate decoder limitation, not a claim
about the published algorithm until the published PDF is checked.

Mechanism. Take an internal boundary after chunk j (1-indexed upper chunk),
with K bits above the boundary and k = n - K bits below it, and let
delta = frac(2^K phi). In the infinite-shot limit the lower chunks decode to
the best k-bit approximation of delta. If that approximation is exactly
10..0 (b_k = 2^(k-1), i.e. |delta - 1/2| <= 2^-(k+1)), Lemma 3.2's excluded
case occurs: the lower chunks do not say whether the upper block rounded up
(delta > 1/2) or down (delta < 1/2).

* If the upper block is flagged ambiguous (C2/C1 > eps), Algorithm 1 takes
  min(t1, t2) = floor, which is correct either way, and no borrow is applied.
* If it is not flagged, its raw value is round(.), and chunk j+1 is the
  special chunk 10..0, so the borrow is suppressed. The result is wrong iff
  delta > 1/2. (With the special-chunk rule disabled the borrow is applied
  and the result is wrong iff delta < 1/2.)

The upper block's top-two ratio in this situation is
    r(delta) = K_M(max(delta, 1-delta)/M) / K_M(min(delta, 1-delta)/M),
which decreases as |delta - 1/2| grows, so every phase in the danger band is
flagged iff eps < eps*(k, m) = r(1/2 + 2^-(k+1)) (the band's edge).
"""

from __future__ import annotations

import numpy as np

from research.awqpe.model.kernel import kernel


def top_two_ratio(delta, width: int) -> np.ndarray:
    """Exact kernel ratio P(second)/P(first) for fractional offset delta."""
    M = 1 << width
    delta = np.asarray(delta, dtype=float)
    near = np.minimum(delta, 1.0 - delta)
    far = np.maximum(delta, 1.0 - delta)
    return kernel(far / M, M) / kernel(near / M, M)


def epsilon_star(k: int, m: int) -> float:
    """Largest eps for which the whole danger band of a (k lower bits, width-m upper block) boundary is flagged."""
    return float(top_two_ratio(0.5 + 2.0 ** -(k + 1), m))


def partition_epsilon_star(widths) -> float:
    """min over internal boundaries of eps*(k_j, m_j): below it, no infinite-shot failure is predicted."""
    n, K, out = sum(widths), 0, np.inf
    for m in widths[:-1]:
        K += m
        out = min(out, epsilon_star(n - K, m))
    return float(out)


def predicted_failure(phi, widths, eps: float, special_rule: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """Per-phase prediction of an infinite-shot failure and the 1-indexed upper chunk responsible (0 = none).

    Exact ties (delta exactly 1/2 at some level, measure zero) are outside
    the prediction; they are reported separately by the verifier.
    """
    phi = np.asarray(phi, dtype=float)
    n = int(sum(widths))
    fail = np.zeros(phi.shape, dtype=bool)
    block = np.zeros(phi.shape, dtype=np.int64)
    K = 0
    for j, m in enumerate(widths[:-1], start=1):
        K += m
        k = n - K
        delta = np.mod(np.ldexp(phi, K), 1.0)
        danger = np.floor(delta * 2.0**k + 0.5) == 2 ** (k - 1)
        unflagged = top_two_ratio(delta, m) <= eps
        wrong_side = delta > 0.5 if special_rule else delta < 0.5
        f = danger & unflagged & wrong_side & ~fail
        block = np.where(f, j, block)
        fail |= f
    return fail, block
