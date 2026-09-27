"""AWQPE decoding with O-ext overlap blocks: a LABELLED NON-PAPER variant ("awqpe_ext").

Faithful AWQPE (Algorithms 1-2, baseline/awqpe.py) has no input for extra
blocks, so an overlap block can only influence a D1-type result through a
modified decoder. This module defines the minimal such modification for the
O-ext mechanism (docs/DECISIONS.md D-020).

Why the widened block is consulted only in one situation. Let chunk j have
real value x = p + delta (p integer, delta = the remainder below the
boundary). Algorithm 2 decides the borrow into chunk j from the decoded lower
part (chunks j+1..B). When that lower part is exactly 10..0 (Lemma 3.2's
excluded case; the source of the eps = 0.9 floor, RESULTS_LOG V1), the lower
chunks cannot say whether delta was above or below 1/2. A block re-measuring
chunk j with v extra bits reads round(2^v x); because delta is within the
lower chunks' half-LSB of 1/2, 2^v delta rounds to 2^(v-1) without a carry,
so floor(reading / 2^v) = p exactly. Outside that situation the widened
reading's own rounding can carry into chunk j's bits (probability ~2^-(v+1)),
and substituting it there makes decoding WORSE (verified: 16-54% failures
for unconditional substitution). So:

  when processing boundary j in Algorithm 2 (LSB -> MSB), if a widened block
  for chunk j has shots AND the already-corrected lower part equals 10..0,
  set chunk j := floor(t1_ext / 2^v) mod 2^(m_j) (no borrow);
  otherwise apply Algorithm 2's rule unchanged.

With no widened blocks the result equals faithful AWQPE exactly (tested).
O-bridge blocks have no D1 interpretation without new stitching logic, so
they are used only by the likelihood decoder D2.
"""

from __future__ import annotations

import numpy as np

from research.awqpe.baseline.awqpe import awqpe_vectorised


def awqpe_ext_decode(chunk_counts, widths, epsilon, jitter, ext: dict) -> np.ndarray:
    """ext maps boundary index j (0-based upper chunk) -> (ext_counts (T, 2^(m_j+v)), v, has (T,) bool, ext_jitter)."""
    widths = [int(m) for m in widths]
    r = awqpe_vectorised(chunk_counts, widths, epsilon, jitter=jitter)
    raw, flags = r["raw"], r["flags"]
    T, B = raw.shape
    # Algorithm 2, step 1 (special chunk), unchanged.
    special = np.zeros(T, dtype=np.int64)
    decided = np.zeros(T, dtype=bool)
    for j in range(B, 0, -1):
        x = raw[:, j - 1]
        nz = (x != 0) & ~decided
        special = np.where(nz & (x == (1 << (widths[j - 1] - 1))), j, special)
        decided |= nz
    # Step 2 with the O-ext rule.
    corr = raw.copy()
    ext_sub = {}
    for j0, (counts, v, has, jit) in ext.items():
        t1 = np.argmax(np.asarray(counts, dtype=float) + jit, axis=1)
        ext_sub[j0] = ((t1 >> v) % (1 << widths[j0]), np.asarray(has, dtype=bool))
    for j in range(B - 1, 0, -1):  # 1-indexed upper chunk j, lower chunks j+1..B
        b_corr = (corr[:, j] >> (widths[j] - 1)) & 1
        b_corr = np.where(flags[:, j - 1] | (special == j + 1), 0, b_corr)
        new = (corr[:, j - 1] - b_corr) % (1 << widths[j - 1])
        if (j - 1) in ext_sub:
            sub, has = ext_sub[j - 1]
            lower_half = corr[:, j] == (1 << (widths[j] - 1))
            for i in range(j + 1, B):
                lower_half &= corr[:, i] == 0
            use = has & lower_half
            new = np.where(use, sub, new)
        corr[:, j - 1] = new
    est = np.zeros(T, dtype=np.int64)
    for i, m in enumerate(widths):
        est = (est << m) | corr[:, i]
    return est
