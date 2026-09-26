"""Dirichlet-kernel model, noise channels and Fisher information (tests 1-5, 10, 15)."""

import numpy as np
import pytest

from research.awqpe.circuits.qiskit_blocks import block_statevector_probabilities, phase_gate_block
from research.awqpe.model import kernel as K
from research.awqpe.model.fisher import block_fisher, block_fisher_finite_difference, quantum_fisher_bound
from research.awqpe.model.noise import NoiseSpec, apply_readout, depolarising_lambda, noisy_block_probabilities

RNG = np.random.default_rng(20260927)


@pytest.mark.parametrize("m", [1, 2, 3, 4, 6, 8])
def test_normalisation(m):
    phi = RNG.random(500)
    for k in (0, 3, 7):
        p = K.block_probabilities(phi, k, m)
        assert np.allclose(p.sum(axis=-1), 1.0, atol=1e-12)
        assert (p >= 0).all()


def test_fejer_equals_closed_form_sin_ratio():
    theta = RNG.random(2000) - 0.5
    for M in (2, 4, 8, 16, 64):
        assert np.allclose(K.kernel(theta, M), K.kernel_closed_form(theta, M), atol=1e-13)


def test_symmetry_and_periodicity():
    theta = RNG.random(300)
    for M in (4, 8, 32):
        assert np.allclose(K.kernel(theta, M), K.kernel(-theta, M), atol=1e-13)
        assert np.allclose(K.kernel(theta, M), K.kernel(theta + 1.0, M), atol=1e-12)
    # Reflecting phi -> -phi reflects outcomes y -> -y mod M.
    phi = RNG.random(50)
    p, q = K.block_probabilities(phi, 2, 3), K.block_probabilities(-phi, 2, 3)
    assert np.allclose(p, q[:, (-np.arange(8)) % 8], atol=1e-12)


def test_exact_grid_phase_is_deterministic():
    for m in (2, 3, 5):
        M = 1 << m
        for j in range(M):
            p = K.block_probabilities(j / M, 0, m)
            assert p[j] == pytest.approx(1.0, abs=1e-12)


def test_halfway_phase_splits_top_two_equally_and_lemma_bounds():
    for m in (2, 3, 4, 6):
        M = 1 << m
        j = M // 2 - 1
        p = K.block_probabilities((j + 0.5) / M, 0, m)
        assert p[j] == pytest.approx(p[j + 1], rel=1e-12)
        assert p[j] + p[j + 1] >= 8 / np.pi**2 - 1e-12
    # Lemma 3.3(a): peak probability >= 4/pi^2 for every phase.
    for m in (2, 3, 4, 6):
        assert K.block_probabilities(RNG.random(2000), 0, m).max(axis=-1).min() >= 4 / np.pi**2 - 1e-12


def test_large_M_limit_peak():
    # As M grows the halfway peak decreases towards 4/pi^2.
    peaks = [K.block_probabilities(0.5 / (1 << m), 0, m)[0] for m in (4, 8, 10)]
    assert all(a >= b for a, b in zip(peaks, peaks[1:]))
    assert peaks[-1] == pytest.approx(4 / np.pi**2, abs=1e-4)


def test_offset_sees_fractional_power():
    phi = 0.8203125
    assert K.block_delta(phi, 3) == pytest.approx(0.5625)
    assert K.block_delta(phi, 5) == pytest.approx(0.25)


@pytest.mark.parametrize("m", [2, 3, 4])
@pytest.mark.parametrize("k", [0, 2, 5])
def test_kernel_matches_qiskit_statevector(m, k):
    for phi in list(RNG.random(3)) + [0.8203125, 0.5, 0.0]:
        sv = block_statevector_probabilities(phase_gate_block(phi, k, m, measure=False), m)
        assert np.allclose(sv, K.block_probabilities(phi, k, m), atol=1e-12)


def test_fisher_analytic_vs_finite_difference_and_qfi():
    phi = RNG.random(200)
    for m, k in ((2, 0), (3, 1), (4, 2)):
        f = block_fisher(phi, k, m)
        assert np.allclose(f, block_fisher_finite_difference(phi, k, m), rtol=1e-6)
        assert np.all(f <= quantum_fisher_bound(k, m) * (1 + 1e-8))


def test_ideal_fisher_is_phase_independent_and_equals_qfi():
    """Recorded finding: for the ideal kernel the block CFI equals the QFI for every phi,
    including exact grid points where the double-zero limit is required."""
    phi = np.concatenate([np.linspace(0, 1, 4001), [0.0, 0.25, 0.5]])
    for m, k in ((2, 0), (3, 0), (5, 3)):
        f = block_fisher(phi, k, m)
        assert np.allclose(f, quantum_fisher_bound(k, m), rtol=1e-8)


def test_fisher_scales_as_four_to_the_offset():
    phi = RNG.random(20)
    assert np.allclose(block_fisher(phi, 3, 3), 64 * block_fisher(np.mod(phi * 8, 1.0), 0, 3), rtol=1e-8)


def test_noisy_fisher_vanishes_at_grid_points_with_depolarisation():
    nz = NoiseSpec(depol_gamma=0.1, depol_mode="fixed")
    assert block_fisher(np.array([0.25]), 0, 2, nz)[0] == pytest.approx(0.0, abs=1e-9)
    assert block_fisher(np.array([0.125]), 0, 2, nz)[0] > 0


def test_noise_channels():
    phi = RNG.random(40)
    ideal = K.block_probabilities(phi, 1, 3)
    # Readout: stochastic, and a 50% flip on every bit gives the uniform distribution.
    ro = noisy_block_probabilities(phi, 1, 3, NoiseSpec(readout_p01=0.02, readout_p10=0.05))
    assert np.allclose(ro.sum(-1), 1.0)
    e0 = np.zeros(8)
    e0[0] = 1.0
    assert np.allclose(apply_readout(e0, 3, 0.5, 0.5), 1 / 8)
    one_bit = apply_readout(e0, 3, 0.1, 0.0)
    assert one_bit[0] == pytest.approx(0.9**3) and one_bit[1] == pytest.approx(0.9**2 * 0.1)
    # Depolarising: exact convex mixture with uniform; lam grows with U-queries.
    nz = NoiseSpec(depol_gamma=0.01)
    lam = depolarising_lambda(nz, 1, 3)
    assert lam == pytest.approx(1 - np.exp(-0.01 * 2 * 7))
    assert np.allclose(noisy_block_probabilities(phi, 1, 3, nz), (1 - lam) * ideal + lam / 8)
    assert depolarising_lambda(nz, 4, 3) > depolarising_lambda(nz, 1, 3)
    # Jitter equals Monte Carlo averaging of the kernel over Gaussian phase noise.
    s = 0.03
    xi = np.random.default_rng(1).normal(0, s, 200000)
    mc = K.kernel(0.37 + xi[:, None] - np.arange(8)[None, :] / 8, 8).mean(axis=0)
    an = noisy_block_probabilities(0.37, 0, 3, NoiseSpec(jitter_sigma=s))
    assert np.allclose(an, mc, atol=3e-3)
    with pytest.raises(ValueError):
        NoiseSpec(readout_p01=0.7).validate()
