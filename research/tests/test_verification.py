"""Independent AWQPE implementation and infinite-shot failure theory (PRELIMINARY, arXiv-v3)."""

import numpy as np
import pytest

from research.awqpe.baseline.awqpe import awqpe_vectorised
from research.awqpe.blocks.geometry import partition_blocks
from research.awqpe.evaluation.metrics import exact_success
from research.awqpe.model.kernel import block_probabilities
from research.awqpe.run_verify_epsilon import GOLDEN
from research.awqpe.verification import awqpe_strings as S
from research.awqpe.verification.theory import epsilon_star, partition_epsilon_star, predicted_failure


@pytest.mark.parametrize("phi,widths,raw,final", GOLDEN)
def test_independent_decoder_reproduces_paper_examples(phi, widths, raw, final):
    P = [block_probabilities(phi, s.offset, s.width).tolist() for s in partition_blocks(widths)]
    r = S.decode(P, widths, 0.9)
    assert r["raw"] == raw and r["est"] == final


def test_raw_borrow_source_breaks_paper_examples():
    """The worked examples select the in-place ('corrected') reading of Algorithm 2."""
    reproduced = []
    for phi, widths, raw, final in GOLDEN:
        P = [block_probabilities(phi, s.offset, s.width).tolist() for s in partition_blocks(widths)]
        reproduced.append(S.decode(P, widths, 0.9, S.Variant(borrow_source="raw"))["est"] == final)
    assert not all(reproduced)


def test_independent_matches_primary_on_grid():
    widths = [3, 2, 3]
    phis = (np.arange(1024) + 0.37) / 1024
    P = [block_probabilities(phis, s.offset, s.width) for s in partition_blocks(widths)]
    for eps in (0.5, 0.9):
        prim = awqpe_vectorised(P, widths, eps)["estimate"]
        ind = [int(S.decode([p[t] for p in P], widths, eps, S.Variant(ties="low"))["est"], 2) for t in range(len(phis))]
        assert np.array_equal(prim, np.array(ind))


def test_hand_traced_failure_case():
    """phi = 4.125/256, [3,2,3], eps = 0.9: block 2 unflagged (ratio 0.89), chunk 3 special -> 00001100."""
    widths = [3, 2, 3]
    phi = 4.125 / 256
    P = [block_probabilities(phi, s.offset, s.width).tolist() for s in partition_blocks(widths)]
    r = S.decode(P, widths, 0.9)
    assert r["est"] == "00001100" and r["special_index"] == 3 and r["flags"] == [False, False, False]
    assert S.decode(P, widths, 0.5)["est"] == "00000100"


def test_theory_predicts_failures_exactly_on_generic_grid():
    for widths in ([4, 4], [3, 2, 3], [2, 2, 2, 2]):
        n = sum(widths)
        phis = (np.arange(1 << (n + 3)) + 0.37) / (1 << (n + 3))
        P = [block_probabilities(phis, s.offset, s.width) for s in partition_blocks(widths)]
        for eps in (0.4, 0.62, 0.8, 0.9):
            obs = ~exact_success(awqpe_vectorised(P, widths, eps)["estimate"], phis, n)
            pred, _ = predicted_failure(phis, widths, eps)
            assert np.array_equal(obs, pred)


def test_epsilon_star_values():
    assert epsilon_star(2, 3) == pytest.approx(0.365, abs=2e-3)
    assert epsilon_star(4, 4) == pytest.approx(0.779, abs=2e-3)
    assert epsilon_star(6, 6) > 0.9
    assert partition_epsilon_star([4, 4]) == pytest.approx(epsilon_star(4, 4))
