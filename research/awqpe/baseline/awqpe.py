"""Faithful AWQPE decoding: Algorithm 1 (outcome selection and ambiguity
flagging) and Algorithm 2 (AWQPEAmbiguityResolution, LSB-to-MSB correction)
of Shukla & Vedula, arXiv:2507.22460v3 (journal-ref Adv. Quantum Technol.
9(3) e00683, 2026).

Two implementations are kept on purpose:

* `awqpe_reference`: a line-by-line transcription of the pseudocode for one
  run, used as the specification in tests;
* `awqpe_vectorised`: the same logic over T independent runs at once, used by
  Monte Carlo campaigns. Tests assert the two agree exactly.

The quantum part of Algorithm 1 (lines 4-13) is performed by whoever supplies
the counts (Dirichlet simulator or Qiskit); this module is purely classical.

Interpretation choices where the paper is ambiguous (docs/DECISIONS.md D-004):
  (a) Ties in lines 14-15 are "picked randomly": we add U(0, 0.5) jitter from a
      caller-supplied generator to integer counts before ranking, which picks
      uniformly among tied outcomes and never reorders distinct counts. The
      ratio test on line 18 uses the unjittered counts.
  (b) Eq. 2.1 defines min(a, b) mod n for the wrap pair {0, n-1} (-> n-1) and
      for a < b with a, b not in {0, n-1}; it is silent on a non-wrap pair
      containing 0 or n-1 (e.g. {0, 2}). We return the ordinary minimum there.
  (c) Algorithm 2, line 17 tests A[j] with 1-indexed chunks; the text says the
      list A is 0-indexed, so chunk j's flag is A[j-1] in Python.
  (d) Line 16 reads the MSB of chunk j+1 as it stands when chunk j is
      processed, i.e. after chunk j+1's own correction (the loop runs
      j = B-1 .. 1 and overwrites in place). Implemented verbatim.
  (e) If the final block is ambiguous its flag is set but b_ml stays t1*
      (line 20 guard); Algorithm 2 never reads the final block's flag.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

DEFAULT_EPSILON = 0.9  # "typically around 0.9" (Sec. 2; Alg. 1 line 2)


def modular_min(a: int, b: int, n: int) -> int:
    """Eq. 2.1: minimum of two distinct residues on the cycle Z/nZ."""
    if a == b:
        raise ValueError("modular_min expects distinct residues.")
    if {a, b} == {0, n - 1}:
        return n - 1
    return min(a, b)


def _msb(x: int, width: int) -> int:
    return (x >> (width - 1)) & 1


def _top_two(counts: np.ndarray, jitter: np.ndarray) -> tuple[int, int]:
    key = np.asarray(counts, dtype=float) + jitter
    order = np.argsort(-key, kind="stable")
    return int(order[0]), int(order[1])


@dataclass(frozen=True)
class AWQPEResult:
    estimate: int  # n-bit integer; phi_est = estimate / 2^n
    raw: tuple[int, ...]  # b_ml per chunk (Algorithm 1 output)
    corrected: tuple[int, ...]  # chunks after Algorithm 2
    flags: tuple[bool, ...]
    special_index: int | None  # 1-indexed S_idx, or None
    top1: tuple[int, ...]
    top2: tuple[int, ...]
    ratios: tuple[float, ...]

    def phase(self, n: int) -> float:
        return self.estimate / float(1 << n)


def ambiguity_resolution(raw: list[int], widths: list[int], flags: list[bool]) -> tuple[list[int], int | None]:
    """Algorithm 2 (AWQPEAmbiguityResolution), transcribed line by line."""
    B = len(widths)
    chunks = list(raw)  # line 1 (already partitioned); chunk j is chunks[j-1]
    s_idx = None  # line 2
    for j in range(B, 0, -1):  # lines 3-14: rightmost non-zero chunk equal to 10..0
        x = chunks[j - 1]
        if x == 0:
            continue
        if x == 1 << (widths[j - 1] - 1):
            s_idx = j
        break
    for j in range(B - 1, 0, -1):  # lines 15-23
        b_corr = _msb(chunks[j], widths[j])  # MSB of chunk j+1
        if flags[j - 1] or s_idx == j + 1:
            b_corr = 0
        x = chunks[j - 1]
        chunks[j - 1] = (x - b_corr) % (1 << widths[j - 1])
    return chunks, s_idx


def awqpe_reference(block_counts, widths, epsilon: float = DEFAULT_EPSILON, jitter=None) -> AWQPEResult:
    """Algorithm 1 lines 14-26 followed by Algorithm 2, for one run.

    block_counts[i] holds counts (or probabilities, for the infinite-shot
    limit) over the 2^(m_i) outcomes of chunk i, MSB chunk first. `jitter`
    supplies the tie-breaking noise (zeros if None: lowest index wins ties).
    """
    widths = [int(m) for m in widths]
    n = sum(widths)
    k = 0
    raw, flags, t1s, t2s, ratios = [], [], [], [], []
    for i, m in enumerate(widths):
        c = np.asarray(block_counts[i], dtype=float)
        if c.shape != (1 << m,):
            raise ValueError(f"block {i} counts must have length 2^{m}.")
        jit = np.zeros(1 << m) if jitter is None else np.asarray(jitter[i], dtype=float)
        t1, t2 = _top_two(c, jit)  # lines 14-15
        b_ml = t1  # line 16
        flag = False  # line 17
        ratio = c[t2] / c[t1] if c[t1] > 0 else 0.0
        if ratio > epsilon:  # line 18
            flag = True  # line 19
            if k + m < n:  # line 20
                b_ml = modular_min(t1, t2, 1 << m)  # line 21
        flags.append(flag)  # line 24
        raw.append(b_ml)  # line 25
        t1s.append(t1)
        t2s.append(t2)
        ratios.append(float(ratio))
        k += m  # line 26
    corrected, s_idx = ambiguity_resolution(raw, widths, flags)
    est = 0
    for x, m in zip(corrected, widths):
        est = (est << m) | x
    return AWQPEResult(est, tuple(raw), tuple(corrected), tuple(flags), s_idx, tuple(t1s), tuple(t2s), tuple(ratios))


def awqpe_vectorised(block_counts, widths, epsilon: float = DEFAULT_EPSILON, rng: np.random.Generator | None = None, jitter=None, special_chunk_rule: bool = True) -> dict[str, np.ndarray]:
    """AWQPE over T runs. block_counts[i] has shape (T, 2^(m_i)).

    Tie-breaking jitter is drawn from `rng` (or taken from `jitter`, arrays of
    the same shapes, which is how tests pin it). Returns arrays: estimate (T,),
    raw/corrected/flags/top1/top2/ratio (T, B), special_index (T,) with
    1-indexed chunk numbers and 0 meaning None.

    special_chunk_rule=False is a diagnostic ABLATION, not the published
    algorithm: S_idx is still reported but no longer suppresses the borrow
    (Algorithm 2, line 17 reduced to the ambiguity-flag test).
    """
    widths = [int(m) for m in widths]
    B, n = len(widths), sum(widths)
    T = np.asarray(block_counts[0]).shape[0]
    raw = np.zeros((T, B), dtype=np.int64)
    flags = np.zeros((T, B), dtype=bool)
    top1 = np.zeros((T, B), dtype=np.int64)
    top2 = np.zeros((T, B), dtype=np.int64)
    ratio = np.zeros((T, B))
    rows = np.arange(T)
    k = 0
    for i, m in enumerate(widths):
        c = np.asarray(block_counts[i], dtype=float)
        if c.shape != (T, 1 << m):
            raise ValueError(f"block {i} counts must have shape (T, 2^{m}).")
        if jitter is not None:
            jit = np.asarray(jitter[i], dtype=float)
        elif rng is not None:
            jit = rng.random(c.shape) * 0.5
        else:
            jit = np.zeros(c.shape)
        order = np.argsort(-(c + jit), axis=1, kind="stable")
        t1, t2 = order[:, 0], order[:, 1]
        c1, c2 = c[rows, t1], c[rows, t2]
        r = np.divide(c2, c1, out=np.zeros(T), where=c1 > 0)
        amb = r > epsilon
        b = t1.copy()
        if k + m < n:
            M = 1 << m
            wrap = ((t1 == 0) & (t2 == M - 1)) | ((t1 == M - 1) & (t2 == 0))
            mm = np.where(wrap, M - 1, np.minimum(t1, t2))
            b = np.where(amb, mm, t1)
        raw[:, i], flags[:, i], top1[:, i], top2[:, i], ratio[:, i] = b, amb, t1, t2, r
        k += m
    # Algorithm 2, step 1: rightmost non-zero chunk equal to 10..0.
    special = np.zeros(T, dtype=np.int64)
    decided = np.zeros(T, dtype=bool)
    for j in range(B, 0, -1):
        x = raw[:, j - 1]
        nz = (x != 0) & ~decided
        special = np.where(nz & (x == (1 << (widths[j - 1] - 1))), j, special)
        decided |= nz
    # Step 2: LSB-to-MSB corrections.
    corr = raw.copy()
    for j in range(B - 1, 0, -1):
        b_corr = (corr[:, j] >> (widths[j] - 1)) & 1
        suppress = flags[:, j - 1] | ((special == j + 1) if special_chunk_rule else False)
        b_corr = np.where(suppress, 0, b_corr)
        corr[:, j - 1] = (corr[:, j - 1] - b_corr) % (1 << widths[j - 1])
    est = np.zeros(T, dtype=np.int64)
    for i, m in enumerate(widths):
        est = (est << m) | corr[:, i]
    return {"estimate": est, "raw": raw, "corrected": corr, "flags": flags, "special_index": special, "top1": top1, "top2": top2, "ratio": ratio}
