"""Deterministic seed derivation.

Every random stream is a numpy SeedSequence built from a tuple of
non-negative integers. Python's built-in hash() is never used: it is salted
per process for str/bytes (PYTHONHASHSEED), which made some Paper A seeds
process-dependent (docs/HISTORICAL_SETUP_RECOVERY.md). Strings are mapped to
integers with SHA-256 instead.
"""

from __future__ import annotations

import hashlib

import numpy as np


def stable_int(text: str, bits: int = 63) -> int:
    """Process-independent non-negative integer for a string."""
    digest = hashlib.sha256(text.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") & ((1 << bits) - 1)


def seed_sequence(*parts: int) -> np.random.SeedSequence:
    ints = [int(p) for p in parts]
    if any(p < 0 for p in ints):
        raise ValueError("seed parts must be non-negative integers.")
    return np.random.SeedSequence(ints)


def generator(*parts: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(seed_sequence(*parts)))
