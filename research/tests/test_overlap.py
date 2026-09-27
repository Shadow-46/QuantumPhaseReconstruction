"""Phase 6 overlap machinery (test area 12): O-ext decoder, geometry, fair budgets, CRN pairing."""

import numpy as np

from research.awqpe.baseline.awqpe import awqpe_vectorised
from research.awqpe.blocks.geometry import BlockSpec, partition_blocks
from research.awqpe.decode.awqpe_overlap import awqpe_ext_decode
from research.awqpe.evaluation.metrics import exact_success
from research.awqpe.model.kernel import block_probabilities
from research.awqpe.overlap.candidates import overlap_candidates
from research.awqpe.run_adaptive_overlap import overlap_shard


def test_candidate_geometry():
    ext = overlap_candidates([3, 2, 3], "ext", v=2)
    assert [(j, s, v) for j, s, v in ext] == [(0, BlockSpec(0, 5, "ext"), 2), (1, BlockSpec(3, 4, "ext"), 2)]
    br = overlap_candidates([4, 4], "bridge", shift="half")
    assert br == [(0, BlockSpec(2, 4, "bridge"), 0)]
    assert overlap_candidates([4, 4], "bridge", shift="1")[0][1] == BlockSpec(3, 4, "bridge")


def test_ext_decoder_identity_and_floor_removal():
    for w in ([3, 2, 3], [2, 2, 2, 2], [4, 4]):
        n = sum(w)
        phi = (np.arange(1 << (n + 4)) + 0.37) / (1 << (n + 4))
        T = len(phi)
        P = [block_probabilities(phi, s.offset, s.width) for s in partition_blocks(w)]
        J = [np.zeros_like(p) for p in P]
        base = awqpe_vectorised(P, w, 0.9, jitter=J)["estimate"]
        assert np.array_equal(base, awqpe_ext_decode(P, w, 0.9, J, {}))
        ext = {j: (block_probabilities(phi, s.offset, s.width), v, np.ones(T, bool), np.zeros((T, s.M))) for j, s, v in overlap_candidates(w, "ext", 1)}
        assert exact_success(awqpe_ext_decode(P, w, 0.9, J, ext), phi, n).all()


def test_shard_fairness_and_pairing():
    cfg = {"experiment_id": "t6", "master_seed": 3, "replicates": 2, "likelihood_refine": {"8": 3},
           "rescue_oracle_decoders": ["likelihood"],
           "variants": [{"label": "ext_v1", "mechanism": "ext", "v": 1, "A_max": 1, "trigger": "eig"},
                        {"label": "bridge_half", "mechanism": "bridge", "shift": "half", "A_max": 2, "trigger": "eig"}]}
    spec = {"widths": [3, 2, 3], "S0": 2, "r": 2, "dS": "S0", "chunk_index": 0, "split": "dev",
            "phases": [{"phase_id": "p0", "stratum": "S2_boundary", "phi": 0.1234, "final_residual": 0.0, "boundary_hardness": 0.0},
                       {"phase_id": "p1", "stratum": "S4_uniform", "phi": 0.6789, "final_residual": 0.0, "boundary_hardness": 0.1}]}
    out = overlap_shard(cfg, spec)
    f = out[out.row_type == "final"]
    main = f[f.arm.isin(["B1_uniform", "B2_p5_eig", "B3_overlap_alone", "B4_eig_plus_overlap"])]
    assert main.total_shots.nunique() == 1  # equal shots for every shot-budget arm
    for v in ("ext_v1", "bridge_half"):
        g = f[(f.variant == v) & (f.decoder == "likelihood")].set_index(["arm", "trial"])
        assert (g.loc["B1u_uniform_Umatched_B3"].u_queries.values >= g.loc["B3_overlap_alone"].u_queries.values).all()
        assert (g.loc["B2u_p5_eig_Umatched_B4"].u_queries.values >= g.loc["B4_eig_plus_overlap"].u_queries.values).all()
        assert (g.loc["B3_overlap_alone"].n_overlap_actions <= (1 if v == "ext_v1" else 2)).all()
    # chunk streams identical across variants -> the no-overlap baselines are identical
    b = f[(f.arm == "B2_p5_eig") & (f.decoder == "likelihood")].pivot(index="trial", columns="variant", values="error")
    assert np.array_equal(b["ext_v1"].values, b["bridge_half"].values)
