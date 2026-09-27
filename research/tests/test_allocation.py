"""Shot allocation (test area 11): policies, support-based EIG, fair budgets, CRN pairing."""

import numpy as np

from research.awqpe.allocation.policies import POLICIES, AllocationState, eig_cell_support
from research.awqpe.blocks.geometry import partition_blocks
from research.awqpe.decode.likelihood import GridPosterior
from research.awqpe.information.signals import global_eig_signals
from research.awqpe.oracle import greedy
from research.awqpe.run_adaptive_shots import adaptive_shard, batch_counts, make_streams
from research.awqpe.sim.oracle_sim import batch_block_counts


def _state(widths, S, T=40, seed=0):
    rng = np.random.default_rng(seed)
    specs = partition_blocks(widths)
    n = sum(widths)
    phis = rng.random(T)
    counts = [batch_block_counts(phis, s, S, rng) for s in specs]
    gp = GridPosterior(n, T, 4)
    for s, c in zip(specs, counts):
        gp.add_counts(s, c)
    return AllocationState(specs, n, counts, np.full((T, len(specs)), S), gp.loglik, gp.G, 10, rng), gp


def test_support_eig_matches_dense_eig():
    state, gp = _state([3, 2, 3], 4)
    sparse, dropped = eig_cell_support(state)
    dense = global_eig_signals(gp, state.specs, state.n)["eig_cell"]
    assert dropped.max() < 1e-5
    assert np.allclose(sparse, dense, atol=1e-4)


def test_policies_return_valid_blocks_and_fisher_picks_lsb():
    state, _ = _state([4, 4, 4], 4, T=20)
    for name, cls in POLICIES.items():
        chosen, sig = cls().choose(state)
        assert chosen.shape == (20,) and chosen.min() >= 0 and chosen.max() < 3, name
    assert (POLICIES["fisher"]().choose(state)[0] == 2).all()


def test_streams_give_policy_independent_counts():
    rng = np.random.default_rng(3)
    specs = partition_blocks([3, 3])
    phis = rng.random(10)
    s = make_streams(phis, specs, 20, rng)
    a = batch_counts(s[0], np.zeros(10, dtype=np.int64), 5, 8) + batch_counts(s[0], np.full(10, 5), 7, 8)
    b = batch_counts(s[0], np.zeros(10, dtype=np.int64), 12, 8)
    assert np.array_equal(a, b)


def test_oracle_is_tagged_and_fair_budget():
    assert greedy.IS_ORACLE
    cfg = {"experiment_id": "t", "master_seed": 1, "replicates": 2, "likelihood_refine": {"6": 3},
           "policies": ["uniform", "eig_cell", "fisher"], "oracle_decoders": ["likelihood"]}
    spec = {"widths": [3, 3], "S0": 2, "r": 2, "dS": "S0", "chunk_index": 0, "split": "dev",
            "phases": [{"phase_id": "p0", "stratum": "S4_uniform", "phi": 0.31, "final_residual": 0.0, "boundary_hardness": 0.1},
                       {"phase_id": "p1", "stratum": "S4_uniform", "phi": 0.77, "final_residual": 0.0, "boundary_hardness": 0.1}]}
    out = adaptive_shard(cfg, spec)
    fin = out[out.row_type == "final"]
    assert (fin.total_shots == 2 * 2 * 2).all()  # r * S0 * B for every policy
    uni = fin[fin.policy == "uniform"]
    assert (uni.shots_b1 == uni.shots_b2).all()
    assert fin[fin.policy.str.startswith("greedy_oracle")].is_oracle.all()
    assert not fin[~fin.policy.str.startswith("greedy_oracle")].is_oracle.any()
