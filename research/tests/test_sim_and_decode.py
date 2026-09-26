"""Simulator reproducibility, common random numbers, truth firewall, phases, likelihood decoder (tests 14, 28)."""

import numpy as np
import pytest

from research.awqpe.blocks.geometry import BlockSpec, bridge_block, extension_block, partition_blocks
from research.awqpe.decode.likelihood import GridPosterior, likelihood_decode
from research.awqpe.evaluation.metrics import best_nbit_integer, circular_error, exact_success
from research.awqpe.model.kernel import block_probabilities
from research.awqpe.phases import make_phase_table
from research.awqpe.sim.oracle_sim import PhaseSimulator, batch_block_counts
from research.awqpe.sim.seeds import stable_int


def test_seeded_reproducibility():
    spec = BlockSpec(2, 3)
    a = PhaseSimulator(0.3141, seed_parts=(1, 2, 3)).oracle().measure(spec, 500)
    b = PhaseSimulator(0.3141, seed_parts=(1, 2, 3)).oracle().measure(spec, 500)
    c = PhaseSimulator(0.3141, seed_parts=(1, 2, 4)).oracle().measure(spec, 500)
    assert np.array_equal(a, b) and not np.array_equal(a, c)


def test_common_random_numbers_are_order_and_batch_independent():
    s1, s2 = BlockSpec(0, 3), BlockSpec(3, 3)
    x = PhaseSimulator(0.71, seed_parts=(9,)).oracle()
    first = x.measure(s1, 30) + x.measure(s1, 70)
    y = PhaseSimulator(0.71, seed_parts=(9,)).oracle()
    y.measure(s2, 50)  # another block requested first must not disturb s1's stream
    assert np.array_equal(first, y.measure(s1, 100))


def test_sampling_matches_distribution():
    sim = PhaseSimulator(0.4321, seed_parts=(5,))
    spec = BlockSpec(1, 3)
    counts = sim.oracle().measure(spec, 400000)
    assert np.allclose(counts / counts.sum(), sim.distribution(spec), atol=4e-3)


def test_cost_ledger():
    o = PhaseSimulator(0.2, seed_parts=(1,)).oracle()
    o.measure(BlockSpec(0, 3), 10)
    o.measure(BlockSpec(3, 2), 5)
    assert o.ledger.shots == 15 and o.ledger.u_queries == 10 * 7 + 5 * 8 * 3
    assert o.ledger.max_control_qubits == 3 and o.ledger.max_power == 16


def test_oracle_exposes_no_truth():
    o = PhaseSimulator(0.123456, seed_parts=(1,)).oracle()
    public = {a for a in dir(o) if not a.startswith("__")}
    assert public == {"_measure_fn", "ledger", "measure"}
    for value in vars(o).values():
        assert not isinstance(value, (float, np.floating)) or value != 0.123456


def test_mixture_state_is_componentwise_mixture():
    sim = PhaseSimulator([0.25, 0.5], weights=[0.5, 0.5], seed_parts=(2,))
    spec = BlockSpec(0, 2)
    assert np.allclose(sim.distribution(spec), [0, 0.5, 0.5, 0])


def test_batch_counts_shapes_and_mixtures():
    rng = np.random.default_rng(0)
    c = batch_block_counts(np.array([0.1, 0.2]), BlockSpec(0, 3), 64, rng)
    assert c.shape == (2, 8) and (c.sum(1) == 64).all()
    cm = batch_block_counts(np.array([[0.25, 0.5]]), BlockSpec(0, 2), 1000, rng, weights=np.array([[1.0, 0.0]]))
    assert cm[0, 1] == 1000


def test_geometry_and_costs():
    specs = partition_blocks([3, 2, 3])
    assert [(s.offset, s.width) for s in specs] == [(0, 3), (3, 2), (5, 3)]
    assert specs[2].u_queries_per_shot == 32 * 7 and specs[2].max_power == 128
    with pytest.raises(ValueError):
        partition_blocks([1, 3])
    assert extension_block(specs[0], 2) == BlockSpec(0, 5, "ext")
    assert bridge_block(specs[0], specs[1], 3, 1) == BlockSpec(2, 3, "bridge")


def test_metrics():
    assert circular_error(0.99, 0.01) == pytest.approx(0.02)
    assert best_nbit_integer(0.8203125, 8) == 0b11010010
    assert best_nbit_integer(0.999, 4) == 0  # rounds up and wraps
    tie = (5 + 0.5) / 16
    assert exact_success(6, tie, 4) and exact_success(5, tie, 4) and not exact_success(7, tie, 4)


def test_phase_table_is_deterministic_and_splits_are_disjoint():
    a = make_phase_table([3, 3], 20, "dev", 11)
    b = make_phase_table([3, 3], 20, "dev", 11)
    t = make_phase_table([3, 3], 20, "test", 11)
    assert a.equals(b)
    # S0 draws from only 2^n dyadic values, so dev/test may share phases there by construction.
    cont_a, cont_t = a[a.stratum != "S0_dyadic"], t[t.stratum != "S0_dyadic"]
    assert not set(np.round(cont_a.phi, 15)) & set(np.round(cont_t.phi, 15))
    s2 = a[a.stratum == "S2_boundary"]
    assert (s2.boundary_hardness <= 0.05 + 1e-12).all()
    assert (np.abs(a[a.stratum == "S1_final_half"].final_residual - 0.5) <= 0.05 + 1e-12).all()
    assert "S5_paper" not in set(t.stratum)
    assert stable_int("x") == stable_int("x")


def test_likelihood_decoder_recovers_phase():
    rng = np.random.default_rng(1)
    widths = [3, 3]
    phis = rng.random(200)
    specs = partition_blocks(widths)
    counts = [batch_block_counts(phis, s, 400, rng) for s in specs]
    est = likelihood_decode(specs, counts, 6)
    assert np.quantile(circular_error(est, phis), 0.95) < 2 ** -6


def test_grid_posterior_credible_mass():
    gp = GridPosterior(4, 1)
    p = block_probabilities(0.3, 0, 4)
    gp.add_counts(BlockSpec(0, 4), p * 3)
    est = gp.map_estimate()
    assert gp.credible_mass(est, 0.5)[0] == pytest.approx(1.0)
    assert 0 < gp.credible_mass(est, 2 ** -6)[0] < 1


def test_unitaries_eigenstates_and_paper_a_mixture():
    """Eigenstate verification (test area: state preparation) and the Paper A |1> decomposition."""
    from research.awqpe.circuits import unitaries as UN
    from research.awqpe.circuits.qiskit_blocks import block_statevector_probabilities, general_unitary_block

    e = UN.modmul_eigenstate_input(2, 21, 1)
    assert e.is_eigenstate and e.eigen_residual < 1e-12 and np.allclose(e.phases, [1 / 6])
    one = UN.modmul_one_input(2, 21)
    assert not one.is_eigenstate
    assert np.allclose(one.phases, np.arange(6) / 6) and np.allclose(one.weights, 1 / 6)
    c = UN.conjugated_input([0.25, 0.1, 0.3, 0.7], 2, np.random.default_rng(0))
    assert c.is_eigenstate and np.allclose(c.phases, [0.3])
    sv = block_statevector_probabilities(general_unitary_block(one.unitary, one.state, 0, 3, measure=False), 3)
    mix = one.weights @ block_probabilities(one.phases, 0, 3)
    assert np.allclose(sv, mix, atol=1e-12)
    padded = UN.controlled_mixture_input(0.3, 0.5, [0.1, 0.2, 0.4, 0.6, 0.8])
    assert padded.unitary.shape == (8, 8) and np.isclose(padded.weights.sum(), 1.0)


def test_sparse_likelihood_path_matches_dense():
    import research.awqpe.decode.likelihood as L

    rng = np.random.default_rng(0)
    spec = BlockSpec(2, 5)
    c = rng.multinomial(12, np.ones(32) / 32, size=7).astype(float)
    dense = c @ L._log_table(1 << 10, 2, 5, None).T
    sparse = L._sparse_loglik(1 << 10, spec, c, None)
    finite = dense > -600
    assert np.allclose(dense[finite], sparse[finite], atol=1e-9)
    assert (sparse[~finite] < -600).all()
