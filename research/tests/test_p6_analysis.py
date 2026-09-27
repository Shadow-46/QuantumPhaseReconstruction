"""P6 analysis conventions (D-022, D-026) and the lowerhalf_gated trigger."""

import numpy as np
import pandas as pd
import pytest

from research.analysis.scripts import summarize_p6 as S
from research.awqpe.run_adaptive_overlap import match_to_target, overlap_shard


def test_bracket_selection():
    # U path: 3 trials x 5 points
    U = np.array([[10, 10, 10], [20, 12, 11], [30, 14, 12], [40, 16, 13], [50, 18, 14]], dtype=float)
    lo, hi = match_to_target(U, np.array([25.0, 16.0, 99.0]))
    assert list(lo) == [1, 3, 4] and list(hi) == [2, 3, -1]  # exact hit -> lower == upper; beyond cap -> invalid upper
    assert U[lo, np.arange(3)].tolist() == [20, 16, 14]


def test_upper_bracket_is_primary_and_no_legacy_arms():
    x, y, acct, role = S.COMPARISONS[S.PRIMARY_COMPARISON]
    assert y.startswith("B2m_upper_") and acct == "equal_U_upper" and role == "PRIMARY"
    assert all(S.LEGACY_MARKER not in c[1] for c in S.COMPARISONS.values())
    with pytest.raises(ValueError):
        S.assert_no_legacy_arms(pd.DataFrame({"arm": ["B4_eig_plus_overlap", "B2u_p5_eig_Umatched_B4"]}))


def _slice():
    rows = []
    for c in range(6):
        for rep in range(2):
            base = {"cluster": f"w|p{c}", "widths": "w", "S0": 2, "r": 2, "dS": 2, "replicate_id": rep}
            rows.append({**base, "arm": S.B4, "tol": 1.0, "u_queries": 100.0, "total_shots": 10})
            invalid = c == 0 and rep == 0
            rows.append({**base, "arm": f"B2m_upper_{S.B4}", "tol": np.nan if invalid else 0.0,
                         "u_queries": np.nan if invalid else 110.0, "total_shots": np.nan if invalid else 11})
    return pd.DataFrame(rows)


def test_invalid_upper_brackets_are_excluded_and_counted():
    r = S.paired_comparison(_slice(), S.B4, f"B2m_upper_{S.B4}", boot=200, perm=200)
    assert r["invalid_matches_excluded"] == 1 and r["pairs"] == 11
    assert r["diff_p_tol"] == pytest.approx(1.0) and r["mean_u_diff"] == pytest.approx(10.0)


def test_selection_rule_uses_primary_and_u_tiebreak():
    common = {"ci_lo": 0, "ci_hi": 0.04, "p_holm_across_variants": 0.1}
    comp = pd.DataFrame([
        {"variant": "a", "mechanism": "ext", "decoder": "likelihood", "comparison": S.PRIMARY_COMPARISON, "diff_p_tol": 0.02, "mean_u_overlap_arm": 500, **common},
        {"variant": "b", "mechanism": "ext", "decoder": "likelihood", "comparison": S.PRIMARY_COMPARISON, "diff_p_tol": 0.02, "mean_u_overlap_arm": 400, **common},
        {"variant": "c", "mechanism": "ext", "decoder": "likelihood", "comparison": "B4_vs_B2_equal_shots", "diff_p_tol": 0.9, "mean_u_overlap_arm": 1, **common},
        {"variant": "d", "mechanism": "ext", "decoder": "awqpe", "comparison": S.PRIMARY_COMPARISON, "diff_p_tol": 0.5, "mean_u_overlap_arm": 1, **common},
    ])
    sel = S.select_variants(comp)
    assert len(sel) == 1  # faithful decoder is a control, never selected; equal-shot rows never used for selection
    assert sel.iloc[0].selected_variant == "b" and bool(sel.iloc[0].tie_broken_by_u)


def test_gated_trigger_acts_only_where_rule_applies_and_keeps_budget():
    cfg = {"experiment_id": "tg", "master_seed": 5, "replicates": 3, "likelihood_refine": {"8": 3},
           "rescue_oracle_decoders": ["awqpe_ext"],
           "variants": [{"label": "ext_gated_A2", "mechanism": "ext", "v": 1, "A_max": 2, "trigger": "lowerhalf_gated"}]}
    phases = [{"phase_id": f"p{i}", "stratum": "S2_boundary", "phi": float(p), "final_residual": 0.0, "boundary_hardness": 0.0}
              for i, p in enumerate([0.0161, 0.0479, 0.3, 0.5156, 0.8203125, 0.123])]
    out = overlap_shard(cfg, {"widths": [3, 2, 3], "S0": 2, "r": 2, "dS": "S0", "chunk_index": 0, "split": "dev", "phases": phases})
    act = out[(out.row_type == "action") & out.arm.isin(["B3_overlap_alone", "B4_eig_plus_overlap"])]
    assert len(act) > 0 and act.lowerhalf09.astype(bool).all()
    f = out[out.row_type == "final"]
    assert f[f.arm.isin(["B1_uniform", "B2_p5_eig", "B3_overlap_alone", "B4_eig_plus_overlap"])].total_shots.nunique() == 1


def test_cached_eig_is_bit_identical_to_p5_signal():
    from research.awqpe.allocation import eig_cached
    from research.awqpe.allocation.policies import AllocationState, eig_cell_support
    from research.awqpe.blocks.geometry import partition_blocks
    from research.awqpe.decode.likelihood import GridPosterior
    from research.awqpe.overlap.candidates import overlap_candidates
    from research.awqpe.sim.oracle_sim import batch_block_counts

    for widths, refine in (([3, 2, 3], 3), ([4, 4, 4, 4], 2)):
        rng = np.random.default_rng(4)
        n = sum(widths)
        chunks = partition_blocks(widths)
        specs = chunks + [c[1] for c in overlap_candidates(widths, "ext", 1)] + [c[1] for c in overlap_candidates(widths, "bridge", shift="1")]
        phis = rng.random(12)
        counts = [batch_block_counts(phis, s, 3, rng) for s in chunks] + [np.zeros((12, s.M), dtype=np.int64) for s in specs[len(chunks):]]
        gp = GridPosterior(n, 12, refine)
        for s, c in zip(chunks, counts):
            gp.add_counts(s, c)
        st = AllocationState(specs, n, counts, np.zeros((12, len(specs))), gp.loglik, gp.G, 0, rng)
        a, da = eig_cell_support(st)
        b, db = eig_cached.eig_cell_support(st)
        assert np.array_equal(a, b) and np.array_equal(da, db)
