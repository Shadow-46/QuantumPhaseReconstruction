"""Phase 7 (D-031): unified controller."""

import numpy as np

from research.awqpe.blocks.geometry import partition_blocks
from research.awqpe.controller import unified_run
from research.awqpe.overlap.candidates import overlap_candidates
from research.awqpe.phases import make_phase_table
from research.awqpe.run_adaptive_shots import make_streams
from research.awqpe.run_tolerance import run_stopping
from research.awqpe.run_unified import unified_shard


def _setup(widths=(4, 4, 4)):
    chunks = partition_blocks(list(widths))
    cands = [c[1] for c in overlap_candidates(list(widths), "bridge", shift="1")]
    phis = make_phase_table(list(widths), 2, "dev", 7, strata=("S2_boundary", "S4_uniform")).phi.to_numpy()
    return chunks, cands, phis, make_streams(phis, chunks + cands, 64, np.random.default_rng(2))


def test_no_overlap_equals_p8_stop_eig():
    chunks, _, _, streams = _setup()
    B = len(chunks)
    a = unified_run(chunks, [], streams[:B], 12, 3, 4, 64, 0, 2**-8, 0.05, np.random.default_rng(11))
    b = run_stopping(chunks, streams[:B], 12, 3, 4, 64, 2**-8, 0.05, "eig", np.random.default_rng(11))
    for x, y in zip(a[:3], b[:3]):
        np.testing.assert_array_equal(x, y)


def test_overlap_budget_and_stop_rule():
    chunks, cands, _, streams = _setup()
    est, sh, stp, nov = unified_run(chunks, cands, streams, 12, 3, 4, 64, 2, 2**-12, 0.05, np.random.default_rng(1))
    B = len(chunks)
    assert (nov <= 2).all() and (sh[:, B:].sum(axis=1) == 4 * nov).all()
    assert (sh[:, :B] <= 64).all() and (sh[:, :B] >= 4).all()


def test_shard_runs():
    cfg = {"replicates": 1, "master_seed": 4, "experiment_id": "t", "S0": 4, "A": 2, "cap_per_block": 32, "alpha": 0.05,
           "tau_exponents": ["half", "n"], "likelihood_refine": {}, "split": "dev"}
    t = make_phase_table([4, 4], 1, "dev", 4, strata=("S2_boundary", "S4_uniform"))
    out = unified_shard(cfg, {"widths": [4, 4], "chunk_index": 0, "phases": t[["phase_id", "stratum", "phi"]].to_dict("records")})
    assert set(out.arm) == {"stop_uniform", "stop_eig", "unified_bridge", "unified_ext"}
    assert set(out.tau_exp) == {4, 8}
    assert (out[out.arm.str.startswith("stop")].overlap_shots == 0).all()
