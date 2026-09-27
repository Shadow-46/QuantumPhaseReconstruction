"""Observable overlap triggers (decision side; no truth, covered by test_leakage.py).

lowerhalf_mask: for each overlap candidate at boundary j, whether the current
faithful (eps 0.9) chunk decode gives a lower part (chunks j+1..B) of exactly
10..0 -- the observable form of the situation in which Algorithm 2 cannot
decide the borrow into chunk j (RESULTS_LOG V1).
"""

from __future__ import annotations

import numpy as np

from research.awqpe.baseline.awqpe import awqpe_vectorised


def lowerhalf_mask(chunk_counts, widths, jitter, cands, epsilon: float = 0.9) -> np.ndarray:
    widths = [int(m) for m in widths]
    corr = awqpe_vectorised(chunk_counts, widths, epsilon, jitter=jitter)["corrected"]
    T, B = corr.shape
    out = np.zeros((T, len(cands)), dtype=bool)
    for i, (j, _, _) in enumerate(cands):
        m = corr[:, j + 1] == (1 << (widths[j + 1] - 1))
        for k in range(j + 2, B):
            m &= corr[:, k] == 0
        out[:, i] = m
    return out
