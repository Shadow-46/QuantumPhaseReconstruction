"""Independent re-implementation of AWQPE Algorithms 1-2 (arXiv:2507.22460v3)
for verification only.

Deliberately shares no code with research/awqpe/baseline/awqpe.py: it works
on binary strings exactly as the paper describes (phi_raw is a concatenated
bit string, Algorithm 2 partitions it into substrings), uses Python ints and
explicit loops, and was written from the pseudocode rather than from the
primary implementation.

Every point where the paper is ambiguous is a switch in `Variant`. The default
Variant reproduces the interpretation used by the faithful baseline
(docs/DECISIONS.md D-004); all other settings are LABELLED ALTERNATIVE
READINGS, not AWQPE as implemented for the main experiments.

    flag_index     "own"  : chunk j's correction is gated by chunk j's own flag
                            (paper text: A[j-1] with 0-indexed A)  [default]
                   "next" : literal A[j] on a 0-indexed list, i.e. the flag of
                            chunk j+1
    minmod         "min_wrap"      : ordinary min, wrap pair {0,n-1} -> n-1 [default]
                   "strict_eq21"   : only the two cases Eq. 2.1 defines; any
                                     other pair keeps t1*
                   "circular_lower": adjacent pair -> circular predecessor;
                                     non-adjacent pair keeps t1*
    borrow_source  "corrected": MSB of chunk j+1 after its own correction [default]
                   "raw"      : MSB of chunk j+1 as produced by Algorithm 1
    final_block    "keep_t1"  : final block keeps t1* even if flagged [default]
                   "min"      : final block also takes the modular min when flagged
    special        "paper"        : S_idx = rightmost non-zero chunk equal to 10..0 [default]
                   "off"          : no special chunk
                   "flagged_only" : S_idx only if that chunk was itself flagged
    ties           "random" (seeded) [default] | "low" (lowest outcome) | "high"
"""

from __future__ import annotations

import random
from dataclasses import dataclass, fields


@dataclass(frozen=True)
class Variant:
    flag_index: str = "own"
    minmod: str = "min_wrap"
    borrow_source: str = "corrected"
    final_block: str = "keep_t1"
    special: str = "paper"
    ties: str = "random"

    @property
    def label(self) -> str:
        default = Variant()
        diffs = [f"{f.name}={getattr(self, f.name)}" for f in fields(self) if getattr(self, f.name) != getattr(default, f.name)]
        return "default" if not diffs else ",".join(diffs)


OPTIONS = {
    "flag_index": ("own", "next"),
    "minmod": ("min_wrap", "strict_eq21", "circular_lower"),
    "borrow_source": ("corrected", "raw"),
    "final_block": ("keep_t1", "min"),
    "special": ("paper", "off", "flagged_only"),
    "ties": ("random", "low", "high"),
}


def _modmin(a: int, b: int, n: int, rule: str, t1: int) -> int:
    wrap = (a == 0 and b == n - 1) or (a == n - 1 and b == 0)
    if rule == "min_wrap":
        return n - 1 if wrap else min(a, b)
    if rule == "strict_eq21":
        if wrap:
            return n - 1
        if a not in (0, n - 1) and b not in (0, n - 1):
            return min(a, b)
        return t1
    if rule == "circular_lower":
        if (a - b) % n == 1:
            return b
        if (b - a) % n == 1:
            return a
        return t1
    raise ValueError(rule)


def _top_two(values, ties: str, rnd: random.Random) -> tuple[int, int]:
    idx = list(range(len(values)))
    if ties == "random":
        keys = {i: rnd.random() for i in idx}
        idx.sort(key=lambda i: (-values[i], keys[i]))
    elif ties == "low":
        idx.sort(key=lambda i: (-values[i], i))
    elif ties == "high":
        idx.sort(key=lambda i: (-values[i], -i))
    else:
        raise ValueError(ties)
    return idx[0], idx[1]


def algorithm1(block_values, widths, epsilon: float, variant: Variant = Variant(), rnd: random.Random | None = None):
    """Return (phi_raw bit string, flags, per-block (t1, t2, C1, C2))."""
    rnd = rnd or random.Random(0)
    n = sum(widths)
    k = 0
    phi_raw = ""
    flags = []
    detail = []
    for values, m in zip(block_values, widths):
        values = [float(v) for v in values]
        t1, t2 = _top_two(values, variant.ties, rnd)
        c1, c2 = values[t1], values[t2]
        b_ml = t1
        flag = False
        if c1 > 0 and c2 / c1 > epsilon:
            flag = True
            if k + m < n or variant.final_block == "min":
                b_ml = _modmin(t1, t2, 2 ** m, variant.minmod, t1)
        flags.append(flag)
        phi_raw += format(b_ml, "b").zfill(m)
        detail.append((t1, t2, c1, c2))
        k += m
    return phi_raw, flags, detail


def algorithm2(phi_raw: str, widths, flags, variant: Variant = Variant()):
    """Return (phi_est bit string, S_idx 1-indexed or None)."""
    chunks, pos = [], 0
    for m in widths:
        chunks.append(phi_raw[pos:pos + m])
        pos += m
    raw_chunks = list(chunks)
    B = len(chunks)
    s_idx = None
    if variant.special != "off":
        for j in range(B, 0, -1):
            x = int(chunks[j - 1], 2)
            if x == 0:
                continue
            if x == 2 ** (widths[j - 1] - 1) and (variant.special == "paper" or flags[j - 1]):
                s_idx = j
            break
    for j in range(B - 1, 0, -1):
        source = chunks if variant.borrow_source == "corrected" else raw_chunks
        b_corr = int(source[j][0])  # MSB of chunk j+1 (1-indexed)
        gate = flags[j - 1] if variant.flag_index == "own" else flags[j]
        if gate or s_idx == j + 1:
            b_corr = 0
        x_new = (int(chunks[j - 1], 2) - b_corr) % (2 ** widths[j - 1])
        chunks[j - 1] = format(x_new, "b").zfill(widths[j - 1])
    return "".join(chunks), s_idx


def decode(block_values, widths, epsilon: float, variant: Variant = Variant(), seed: int = 0) -> dict:
    rnd = random.Random(seed)
    phi_raw, flags, detail = algorithm1(block_values, widths, epsilon, variant, rnd)
    phi_est, s_idx = algorithm2(phi_raw, widths, flags, variant)
    return {"raw": phi_raw, "est": phi_est, "flags": flags, "special_index": s_idx, "detail": detail}
