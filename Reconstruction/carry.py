"""
Purpose
    Check carry-aware overlap consistency between adjacent window candidates.
Theory
    Exact overlap equality is always valid when both windows measure
    noiseless slices of the same true bits, but a naive equality-only check
    rejects legitimate stitches whenever a borrow/carry propagates across the
    block boundary from independent per-block rounding (see arXiv:2509.05010,
    Section 2.2 and Algorithm 3). Following that algorithm, the right
    candidate's bit immediately past the shared overlap acts as the carry
    indicator c, and a candidate pair is also accepted whenever
    (int(tail_left, 2) - c) mod 2**overlap == int(head_right, 2), strictly
    widening acceptance beyond plain equality rather than replacing it.

    An original exhaustive audit (Experiments/phase3_adversarial.py:
    run_carry_truth_table) reported a 55.4% false-accept rate for this
    carry-corrected branch, but that audit's synthetic window construction
    had a geometry bug (right.start=overlap instead of width-overlap): it
    always produced a 1-bit true overlap regardless of the loop's nominal
    `overlap` value, so tail/head were never actually compared at their
    intended width. With the bug fixed, the same 168-combination audit shows
    this function is an exact classifier against the "no borrow, or exactly
    one legitimate single-bit borrow" ground truth: 0% false-accept, 0%
    false-reject (see Data/failure_study/phase3_carry_truth_table.csv,
    regenerated; REPORT.md Section 4.2, corrected). No fix to this function
    was needed or made.
Inputs
    WindowCandidate objects or explicit bit strings and window coordinates.
Outputs
    Boolean compatibility decisions and merged bit strings.
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Reconstruction.candidate_generation import WindowCandidate


def overlap_width(left: WindowCandidate, right: WindowCandidate) -> int:
    """Return the number of shared bit positions between two candidates."""
    return max(0, min(left.end, right.end) - max(left.start, right.start))


def carry_bit(right: WindowCandidate, overlap: int) -> int:
    """Return the paper's carry indicator: the bit of `right` just past the
    shared overlap (0 when the overlap consumes the entire right window)."""
    if overlap < 0:
        raise ValueError("overlap must be non-negative.")
    if overlap >= right.width:
        return 0
    return int(right.local_bits[overlap])


def candidates_compatible(left: WindowCandidate, right: WindowCandidate) -> bool:
    """Return True when the overlapping bits of two candidates agree exactly,
    or agree modulo the paper's single-bit carry correction."""
    shared = overlap_width(left, right)
    if shared == 0:
        return True
    left_offset = max(left.start, right.start) - left.start
    right_offset = max(left.start, right.start) - right.start
    tail = left.local_bits[left_offset:left_offset + shared]
    head = right.local_bits[right_offset:right_offset + shared]
    if tail == head:
        return True
    c = carry_bit(right, right_offset + shared)
    corrected = (int(tail, 2) - c) % (2**shared)
    return corrected == int(head, 2)


def merge_window_bits(left: str, right: str, overlap: int) -> str:
    """Merge two bit strings after validating their overlap."""
    if overlap < 0 or overlap > min(len(left), len(right)):
        raise ValueError("invalid overlap.")
    if overlap and left[-overlap:] != right[:overlap]:
        raise ValueError("window bits are incompatible.")
    return left + right[overlap:]


def self_test() -> None:
    """Verify overlap acceptance/rejection and the carry-aware overlap check."""
    assert merge_window_bits("101", "111", 1) == "10111"
    try:
        merge_window_bits("101", "011", 1)
    except ValueError:
        pass
    else:
        raise AssertionError("incompatible windows were accepted")

    # A borrow of 1 across the block boundary: tail("01") disagrees with
    # head("00") under plain equality, but the paper's carry-aware check
    # accepts it because the bit past the overlap (right.local_bits[2]) is 1.
    left = WindowCandidate(value=0, total_precision=6, start=0, width=4, local_bits="1101", weight=1.0)
    right_with_carry = WindowCandidate(value=0, total_precision=6, start=2, width=4, local_bits="0011", weight=1.0)
    assert candidates_compatible(left, right_with_carry)

    # Without the carry bit set, the same tail/head pair must be rejected.
    right_without_carry = WindowCandidate(value=0, total_precision=6, start=2, width=4, local_bits="0001", weight=1.0)
    assert not candidates_compatible(left, right_without_carry)


if __name__ == "__main__":
    self_test()
    print("carry self-test passed")
