"""
Purpose
    Stitch ranked window candidates into full phase estimates.
Theory
    The paper reconstruction stage forms globally consistent paths through
    overlapping window candidates. This implementation performs deterministic
    beam stitching with exact overlap checks and probability-product scoring.
Inputs
    Candidate lists ordered by window position and a maximum path count.
Outputs
    Ranked reconstructed full-precision phase integers.
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from dataclasses import dataclass
from fractions import Fraction

from Reconstruction.candidate_generation import WindowCandidate
from Reconstruction.carry import candidates_compatible


@dataclass(frozen=True)
class StitchedPhase:
    """A globally stitched phase candidate and its score."""
    value: int
    total_precision: int
    bitstring: str
    score: float

    @property
    def phase(self) -> Fraction:
        """Return the phase candidate as a fraction."""
        return Fraction(self.value, 2**self.total_precision)


def stitch_candidates(windows: list[list[WindowCandidate]], max_paths: int = 32) -> list[StitchedPhase]:
    """Return ranked globally consistent phase candidates."""
    if not windows or any(not group for group in windows):
        raise ValueError("each window must contain at least one candidate.")
    total_precision = windows[0][0].total_precision
    paths: list[tuple[list[WindowCandidate], float]] = [([candidate], candidate.weight) for candidate in windows[0]]
    for group in windows[1:]:
        next_paths: list[tuple[list[WindowCandidate], float]] = []
        for path, score in paths:
            for candidate in group:
                if candidates_compatible(path[-1], candidate):
                    next_paths.append(([*path, candidate], score * candidate.weight))
        paths = sorted(next_paths, key=lambda item: item[1], reverse=True)[:max_paths]
        if not paths:
            return []
    stitched: list[StitchedPhase] = []
    for path, score in paths:
        bits = ["0"] * total_precision
        for candidate in path:
            for offset, bit in enumerate(candidate.local_bits):
                bits[candidate.start + offset] = bit
        bitstring = "".join(bits)
        stitched.append(StitchedPhase(int(bitstring, 2), total_precision, bitstring, score))
    return sorted(stitched, key=lambda item: item.score, reverse=True)


def self_test() -> None:
    """Verify stitching of three overlapping deterministic windows."""
    from Reconstruction.candidate_generation import exact_window_counts, generate_window_candidates
    groups = [generate_window_candidates(exact_window_counts(0b101100, 6, s, 3), s, 3, 6, 1) for s in (0, 2, 3)]
    stitched = stitch_candidates(groups)
    assert stitched[0].bitstring == "101100"


if __name__ == "__main__":
    self_test()
    print("stitching self-test passed")
