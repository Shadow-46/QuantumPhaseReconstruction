"""Phase 9 (D-030): noise runner."""

import numpy as np

from research.awqpe.blocks.geometry import partition_blocks
from research.awqpe.model.noise import NoiseSpec, noisy_block_probabilities
from research.awqpe.phases import make_phase_table
from research.awqpe.run_noise import noise_shard, noisy_streams


def test_noisy_streams_follow_noisy_distribution():
    spec = partition_blocks([4, 4])[1]
    nz = NoiseSpec(readout_p01=0.05, readout_p10=0.05, depol_gamma=1e-3)
    phis = np.array([0.3137])
    y = noisy_streams(phis, [spec], 40000, nz, np.random.default_rng(0))[0][0]
    emp = np.bincount(y, minlength=16) / len(y)
    p = noisy_block_probabilities(phis, spec.offset, spec.width, nz)[0]
    assert np.abs(emp - p).max() < 0.01


def test_shard_ideal_and_noisy():
    cfg = {"replicates": 2, "master_seed": 3, "experiment_id": "t", "cap_per_block": 32, "S0": 4, "alpha": 0.05,
           "fixed_S": [16], "stop_tau_exponents": ["half"], "likelihood_refine": {}, "split": "dev"}
    t = make_phase_table([4, 4], 1, "dev", 3, strata=("S0_dyadic", "S4_uniform"))
    base = {"widths": [4, 4], "chunk_index": 0, "phases": t[["phase_id", "stratum", "phi"]].to_dict("records")}
    ideal = noise_shard(cfg, {**base, "noise_id": "ideal", "channel": "none", "noise": {}})
    assert set(ideal.decoder) == {"awqpe", "awqpe_eps_safe", "likelihood_ideal"}
    noisy = noise_shard(cfg, {**base, "noise_id": "ro05", "channel": "readout", "noise": {"readout_p01": 0.05, "readout_p10": 0.05}})
    assert "likelihood_aware" in set(noisy.decoder)
    st = noisy[noisy.arm == "stop_tau4"]
    assert (st.shots <= 64).all() and st.tau.eq(2**-4).all()
    assert ideal.equals(noise_shard(cfg, {**base, "noise_id": "ideal", "channel": "none", "noise": {}}))
