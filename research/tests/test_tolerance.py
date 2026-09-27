"""Phase 8 (D-028): tolerance stopping rule and runner."""

import numpy as np
import pytest

from research.awqpe.blocks.geometry import partition_blocks
from research.awqpe.decode.likelihood import GridPosterior
from research.awqpe.phases import make_phase_table
from research.awqpe.run_adaptive_shots import make_streams
from research.awqpe.run_tolerance import run_stopping, tol_shard
from research.awqpe.stopping.rules import credible_mass_at_map, tau_bits, truncated_widths


@pytest.mark.parametrize("tau", [2**-3, 2**-6, 2**-8, 0.3, 0.49])
def test_credible_mass_matches_grid_posterior(tau):
    rng = np.random.default_rng(1)
    gp = GridPosterior(8, 7, 2)
    gp.loglik = rng.normal(0, 6, (7, gp.G))
    gp.loglik[0, 0] += 50  # MAP at the wrap point
    gp.loglik[1, -1] += 50
    mass, c = credible_mass_at_map(gp.loglik, tau)
    np.testing.assert_allclose(mass, gp.credible_mass(gp.map_estimate(), tau), atol=1e-12)
    np.testing.assert_array_equal(c / gp.G, gp.map_estimate())


def test_tau_bits_and_prefix():
    assert tau_bits(2**-3) == 2 and tau_bits(2**-8) == 7 and tau_bits(0.4) == 1
    assert truncated_widths([4, 4, 4, 4], 2**-3) == [4]
    assert truncated_widths([4, 4, 4, 4], 2**-5) == [4]
    assert truncated_widths([4, 4, 4, 4], 2**-8) == [4, 4]
    assert truncated_widths([3, 2, 3], 2**-8) == [3, 2, 3]
    assert truncated_widths([3, 2, 3], 2**-6) == [3, 2]


@pytest.mark.parametrize("policy", ["uniform", "eig"])
def test_stopping_respects_rule_and_cap(policy):
    widths = [4, 4]
    specs = partition_blocks(widths)
    phis = make_phase_table(widths, 2, "dev", 5, strata=("S1_final_half", "S4_uniform")).phi.to_numpy()
    streams = make_streams(phis, specs, 64, np.random.default_rng(3))
    est, shots, stopped, mass = run_stopping(specs, streams, 8, 3, 4, 64, 2**-6, 0.05, policy, np.random.default_rng(4))
    assert (shots <= 64).all() and (shots >= 4).all()
    assert (mass[stopped] >= 0.95).all()
    assert ((~stopped) <= (shots + 4 > 64).all(axis=1)).all()  # unstopped only when capped
    _, sh2, st2, _ = run_stopping(specs, streams, 8, 3, 4, 64, 2**-6, 1e-12, policy, np.random.default_rng(4))
    assert (sh2.sum(axis=1) >= shots.sum(axis=1)).all()  # a stricter alpha never stops earlier


def test_shard_rows_and_crn():
    cfg = {"S0": 4, "replicates": 2, "cap_per_block": 32, "alpha": 0.05, "master_seed": 9, "experiment_id": "t",
           "fixed_S": [4, 32], "stop_policies": ["uniform", "eig"], "tau_exponents": [3, 8], "likelihood_refine": {}, "split": "dev"}
    t = make_phase_table([4, 4], 1, "dev", 9, strata=("S0_dyadic", "S4_uniform"))
    spec = {"widths": [4, 4], "chunk_index": 0, "phases": t[["phase_id", "stratum", "phi"]].to_dict("records")}
    out = tol_shard(cfg, spec)
    arms = set(out.arm)
    assert {"fixed_S4", "fixed_S32", "stop_uniform", "stop_eig", "fixed_S4_trunc", "stop_eig_trunc"} <= arms
    assert "stop_eig_trunc" not in set(out[out.tau_exp == 8].arm)  # prefix == full partition at tau 2^-8
    tr = out[out.arm == "fixed_S32_trunc"]
    assert (tr.blocks_used == 1).all() and (tr.u_queries == 32 * 15).all()
    fx = out[(out.arm == "fixed_S32") & (out.decoder == "likelihood")]
    assert (fx.shots == 64).all() and fx.covered.mean() >= 0.75
    again = tol_shard(cfg, spec)
    assert out.drop(columns=["credible_mass"]).equals(again.drop(columns=["credible_mass"]))  # deterministic
