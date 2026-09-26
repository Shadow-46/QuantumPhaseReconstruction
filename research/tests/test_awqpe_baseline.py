"""Faithful AWQPE baseline: top-two, ambiguity, Algorithm 2, special chunk, golden cases (tests 6-9)."""

import numpy as np
import pytest

from research.awqpe.baseline.awqpe import ambiguity_resolution, awqpe_reference, awqpe_vectorised, modular_min
from research.awqpe.blocks.geometry import partition_blocks
from research.awqpe.model.kernel import block_probabilities


def _infinite_shot(phi, widths):
    return [block_probabilities(phi, s.offset, s.width) for s in partition_blocks(widths)]


def _bits(values, widths):
    return "".join(format(int(x), f"0{m}b") for x, m in zip(values, widths))


# Paper Sec. 6.1 walkthrough and Table 2. Case 5 is listed in the paper as
# m = [3]*9 but its 30-bit strings require ten 3-bit blocks (docs/PAPER_VERSION_NOTES.md).
GOLDEN = [
    (0.8203125, [3, 2, 3], "11110010", "11010010"),
    (0.3, [2, 2], "0101", "0101"),
    (np.pi / 6, [3, 2, 2, 3], "1000111000", "1000011000"),
    (0.671875, [4, 4], "10111100", "10101100"),
    (1 / np.sqrt(2), [3] * 10, "110101010000010100110011010101", "101101010000010011110011001101"),
    (np.sin(np.pi / 12), [5, 6, 7, 4], "0100001001000010001110", "0100001001000001111110"),
]


@pytest.mark.parametrize("phi,widths,raw,final", GOLDEN)
def test_paper_golden_cases_infinite_shot(phi, widths, raw, final):
    r = awqpe_reference(_infinite_shot(phi, widths), widths)
    assert _bits(r.raw, widths) == raw
    assert format(r.estimate, f"0{sum(widths)}b") == final


def test_walkthrough_counts_from_paper():
    """The paper's printed Sec. 6.1 counts (top entries) decode to the stated result."""
    b1 = np.zeros(8)
    for s, c in (("111", 5180), ("110", 3284), ("000", 534), ("101", 459), ("001", 237)):
        b1[int(s, 2)] = c
    b2 = np.zeros(4)
    for s, c in (("10", 8374), ("11", 1078), ("01", 467), ("00", 321)):
        b2[int(s, 2)] = c
    b3 = np.zeros(8)
    b3[2] = 10240
    r = awqpe_reference([b1, b2, b3], [3, 2, 3])
    assert _bits(r.raw, [3, 2, 3]) == "11110010" and r.estimate == 0b11010010
    assert r.flags == (False, False, False)


def test_modular_min_eq_2_1():
    assert modular_min(0, 7, 8) == 7 and modular_min(7, 0, 8) == 7
    assert modular_min(3, 4, 8) == 3 and modular_min(4, 3, 8) == 3
    assert modular_min(0, 1, 8) == 0  # non-wrap pair containing 0 (paper silent; D-004b)
    with pytest.raises(ValueError):
        modular_min(2, 2, 8)


def test_ambiguity_flag_and_min_selection():
    # Chunk 1 ambiguous between 5 and 6 (ratio 0.95 > 0.9): raw picks min = 5 and no correction.
    b1 = np.zeros(8)
    b1[6], b1[5] = 100, 95
    b2 = np.zeros(8)
    b2[4] = 100  # '100' would normally trigger a borrow
    r = awqpe_reference([b1, b2], [3, 3], epsilon=0.9)
    assert r.flags == (True, False) and r.raw[0] == 5 and r.corrected[0] == 5
    # Just below threshold: not ambiguous, t1 kept; chunk 2 = '100' is special -> no borrow either.
    b1[5] = 89
    r = awqpe_reference([b1, b2], [3, 3], epsilon=0.9)
    assert r.flags == (False, False) and r.special_index == 2 and r.corrected[0] == 6


def test_final_block_ambiguity_keeps_top1():
    b1 = np.zeros(4)
    b1[1] = 10
    b2 = np.zeros(4)
    b2[3], b2[2] = 10, 10  # tie; with zero jitter the lower index (2) wins the argsort
    r = awqpe_reference([b1, b2], [2, 2], epsilon=0.9)
    assert r.flags[1] is True and r.raw[1] == r.top1[1]


def test_borrow_and_wrap():
    # chunk2 MSB = 1 -> chunk1 decremented mod 2^m; chunk1 = 0 wraps to 2^m - 1.
    chunks, s = ambiguity_resolution([0, 3], [2, 2], [False, False])
    assert chunks == [3, 3] and s is None


def test_special_chunk_detection_rightmost_nonzero():
    chunks, s = ambiguity_resolution([3, 2, 0], [2, 2, 2], [False, False, False])
    assert s == 2 and chunks == [3, 2, 0]  # no borrow into chunk 1 because chunk 2 is special
    chunks, s = ambiguity_resolution([3, 3, 0], [2, 2, 2], [False, False, False])
    assert s is None and chunks == [2, 3, 0]


def test_vectorised_matches_reference_exactly():
    rng = np.random.default_rng(7)
    for widths in ([2, 2], [3, 2, 3], [4, 4], [2, 3, 2, 3]):
        T = 400
        phis = rng.random(T)
        counts = []
        for s in partition_blocks(widths):
            p = block_probabilities(phis, s.offset, s.width)
            counts.append(rng.multinomial(rng.integers(3, 40), p))  # small shots -> many ties
        jitter = [rng.random(c.shape) * 0.5 for c in counts]
        for eps in (0.5, 0.9):
            v = awqpe_vectorised(counts, widths, eps, jitter=jitter)
            for t in range(T):
                r = awqpe_reference([c[t] for c in counts], widths, eps, jitter=[j[t] for j in jitter])
                assert v["estimate"][t] == r.estimate
                assert tuple(v["raw"][t]) == r.raw and tuple(v["flags"][t]) == r.flags
                assert v["special_index"][t] == (r.special_index or 0)


def test_tie_breaking_is_seeded_and_uniform():
    counts = [np.tile([5.0, 5.0, 0.0, 0.0], (4000, 1)), np.tile([9.0, 0.0, 0.0, 0.0], (4000, 1))]
    a = awqpe_vectorised(counts, [2, 2], rng=np.random.default_rng(3))
    b = awqpe_vectorised(counts, [2, 2], rng=np.random.default_rng(3))
    assert np.array_equal(a["top1"], b["top1"])
    assert 0.45 < np.mean(a["top1"][:, 0] == 0) < 0.55
