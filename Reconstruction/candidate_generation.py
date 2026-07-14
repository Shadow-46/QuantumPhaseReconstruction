"""
Purpose
    Generate phase candidates from window-local measurement counts.
Theory
    A windowed QPE experiment estimates contiguous phase bits. Candidate
    generation lifts each local bit pattern into full-precision bit positions;
    stitching later resolves the unknown positions through overlap agreement.
    generate_window_candidates_by_coverage (Data/failure_study/FORMULATION.md
    Section 4.1) is a coverage-based alternative to the fixed candidate_count
    truncation: it takes ranked candidates until cumulative probability mass
    reaches 1 - delta_cov, so ambiguous windows automatically receive more
    candidates and confident windows fewer, with no posterior smoothing.
Inputs
    Window counts, window start/width, total precision, and either a
    candidate limit or a coverage target.
Outputs
    Ranked WindowCandidate objects for downstream carry and stitching checks.
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
from typing import Mapping


@dataclass(frozen=True)
class WindowCandidate:
    """Ranked full-precision phase candidate induced by one window outcome."""
    value: int
    total_precision: int
    start: int
    width: int
    local_bits: str
    weight: float

    @property
    def phase(self) -> Fraction:
        """Return the candidate value as a phase fraction."""
        return Fraction(self.value, 2**self.total_precision)

    @property
    def end(self) -> int:
        """Return the exclusive end bit index of the represented window."""
        return self.start + self.width


def normalize_counts(counts: Mapping[str, int]) -> dict[str, float]:
    """Normalize a count dictionary into probabilities."""
    total = sum(counts.values())
    if total <= 0:
        raise ValueError("counts must contain positive total mass.")
    return {bits.replace(" ", ""): count / total for bits, count in counts.items()}


def generate_window_candidates(counts: Mapping[str, int], start: int, width: int, total_precision: int, candidate_count: int) -> list[WindowCandidate]:
    """Generate ranked full-precision candidates matching one observed window."""
    if start < 0 or width < 1 or start + width > total_precision:
        raise ValueError("invalid window geometry.")
    if candidate_count < 1:
        raise ValueError("candidate_count must be positive.")
    ranked = sorted(normalize_counts(counts).items(), key=lambda item: item[1], reverse=True)[:candidate_count]
    suffix_width = total_precision - start - width
    candidates: list[WindowCandidate] = []
    for bits, probability in ranked:
        window_bits = bits[-width:] if len(bits) != width else bits
        local = int(window_bits, 2)
        value = local << suffix_width
        candidates.append(WindowCandidate(value, total_precision, start, width, window_bits, probability))
    return candidates


def generate_window_candidates_by_coverage(counts: Mapping[str, int], start: int, width: int, total_precision: int, delta_cov: float) -> list[WindowCandidate]:
    """Generate ranked full-precision candidates covering cumulative
    probability mass >= 1 - delta_cov, in place of a fixed candidate_count
    (Data/failure_study/FORMULATION.md Section 4.1, module B). Ambiguous
    windows (flat count distributions) automatically receive more
    candidates; confident windows receive fewer, without any posterior
    smoothing of the ranking itself."""
    if start < 0 or width < 1 or start + width > total_precision:
        raise ValueError("invalid window geometry.")
    if not 0 < delta_cov <= 1:
        raise ValueError("delta_cov must satisfy 0 < delta_cov <= 1.")
    ranked = sorted(normalize_counts(counts).items(), key=lambda item: item[1], reverse=True)
    target = 1.0 - delta_cov
    suffix_width = total_precision - start - width
    candidates: list[WindowCandidate] = []
    cumulative = 0.0
    for bits, probability in ranked:
        if cumulative >= target and candidates:
            break
        window_bits = bits[-width:] if len(bits) != width else bits
        local = int(window_bits, 2)
        value = local << suffix_width
        candidates.append(WindowCandidate(value, total_precision, start, width, window_bits, probability))
        cumulative += probability
    return candidates


def exact_window_counts(phase_value: int, total_precision: int, start: int, width: int) -> dict[str, int]:
    """Return deterministic counts for a full-precision value restricted to one window."""
    if not 0 <= phase_value < 2**total_precision:
        raise ValueError("phase_value does not fit the requested precision.")
    bits = format(phase_value, f"0{total_precision}b")
    return {bits[start:start + width]: 1}


def self_test() -> None:
    """Verify deterministic candidate extraction for a known phase value."""
    counts = exact_window_counts(4, 4, 1, 2)
    candidates = generate_window_candidates(counts, 1, 2, 4, 1)
    assert candidates[0].local_bits == "10"
    assert candidates[0].phase == Fraction(4, 16)

    # A confident (near-deterministic) window needs only its top candidate
    # to cover 1 - delta_cov of the mass.
    confident_counts = {"00": 990, "01": 10}
    confident = generate_window_candidates_by_coverage(confident_counts, 0, 2, 4, delta_cov=0.05)
    assert len(confident) == 1 and confident[0].local_bits == "00"

    # An ambiguous (near-flat) window needs multiple candidates to reach the
    # same coverage target.
    ambiguous_counts = {"00": 260, "01": 250, "10": 245, "11": 245}
    ambiguous = generate_window_candidates_by_coverage(ambiguous_counts, 0, 2, 4, delta_cov=0.05)
    assert len(ambiguous) > 1


if __name__ == "__main__":
    self_test()
    print("candidate_generation self-test passed")
