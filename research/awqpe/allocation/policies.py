"""Sequential global shot-allocation policies (Phase 5).

A policy sees an `AllocationState` (observed counts per block, shots per
block, block geometry, the running grid log-likelihood over phi, remaining
budget, and a private tie-breaking generator) and returns, for each of T
trials, the index of the block that receives the next batch.

Nothing here may see the true phase, evaluator outputs, or future shots
(enforced by research/tests/test_leakage.py). The greedy oracle lives in
research/awqpe/oracle/ and is NOT one of these policies.

Signal orientation (higher = more need) follows docs/DECISIONS.md D-016.
The information-gain signal is the exact mutual information of ONE further
shot on the block with the n-bit cell of phi, evaluated on the smallest set
of grid points holding >= 1 - SUPPORT_TAIL posterior mass (capped at
SUPPORT_MAX points); the discarded mass is recorded.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from research.awqpe.decode.likelihood import MAX_TABLE_ELEMENTS, _log_table
from research.awqpe.information.signals import count_signals, legacy_cw
from research.awqpe.model.kernel import kernel_closed_form

SUPPORT_TAIL = 1e-6
SUPPORT_MAX = 8192


@dataclass
class AllocationState:
    specs: list
    n: int
    counts: list  # per block, (T, 2^w) int
    shots: np.ndarray  # (T, B)
    loglik: np.ndarray  # (T, G) running grid log-likelihood (uniform prior)
    G: int
    remaining: int
    rng: np.random.Generator
    step: int = 0
    diagnostics: dict = field(default_factory=dict)


def _block_probs_at(G: int, spec, idx: np.ndarray) -> np.ndarray:
    """p(y | phi_g) for grid indices idx (any shape), returns idx.shape + (M,)."""
    if G * spec.M <= MAX_TABLE_ELEMENTS:
        return np.exp(_log_table(G, spec.offset, spec.width, None))[idx]
    delta = ((idx.astype(np.int64) << spec.offset) % G) / G
    y = np.arange(spec.M) / spec.M
    return kernel_closed_form(delta[..., None] - y, spec.M)


def _entropy_bits(p, axis=-1):
    with np.errstate(divide="ignore", invalid="ignore"):
        return -np.sum(np.where(p > 0, p * np.log2(p), 0.0), axis=axis)


def eig_cell_support(state: AllocationState) -> tuple[np.ndarray, np.ndarray]:
    """(T, B) mutual information of one extra shot per block with the n-bit cell; plus dropped mass (T,)."""
    ll = state.loglik
    T, G = ll.shape
    post = np.exp(ll - ll.max(axis=1, keepdims=True))
    post /= post.sum(axis=1, keepdims=True)
    order = np.argsort(-post, axis=1)
    csum = np.cumsum(np.take_along_axis(post, order, axis=1), axis=1)
    K = int(min(SUPPORT_MAX, G, max(1, np.max(np.argmax(csum >= 1 - SUPPORT_TAIL, axis=1)) + 1)))
    idx = order[:, :K]
    w = np.take_along_axis(post, idx, axis=1)
    dropped = 1.0 - w.sum(axis=1)
    w = w / w.sum(axis=1, keepdims=True)
    sub = G >> state.n
    cells = np.mod((idx + sub // 2) // sub, 1 << state.n)
    key = (np.arange(T)[:, None] * (1 << state.n) + cells).ravel()
    uniq, inv = np.unique(key, return_inverse=True)
    mass = np.zeros(len(uniq))
    np.add.at(mass, inv, w.ravel())
    out = np.zeros((T, len(state.specs)))
    for b, spec in enumerate(state.specs):
        P = _block_probs_at(G, spec, idx)  # (T, K, M)
        pred = np.einsum("tk,tkm->tm", w, P)
        pc = np.zeros((len(uniq), spec.M))
        np.add.at(pc, inv, (w[..., None] * P).reshape(-1, spec.M))
        cond = pc / np.maximum(mass[:, None], 1e-300)
        h_cond = np.zeros(T)
        np.add.at(h_cond, uniq // (1 << state.n), mass * _entropy_bits(cond))
        out[:, b] = _entropy_bits(pred) - h_cond
    return out, dropped


def _per_block(state: AllocationState, fn) -> np.ndarray:
    return np.stack([fn(c) for c in state.counts], axis=1)


class Policy:
    name = "policy"
    is_oracle = False

    def signal(self, state: AllocationState) -> np.ndarray:
        raise NotImplementedError

    def choose(self, state: AllocationState) -> tuple[np.ndarray, np.ndarray]:
        s = self.signal(state)
        jitter = state.rng.random(s.shape) * 1e-9 * (np.abs(s).max() + 1.0)  # seeded tie-breaking only
        return np.argmax(s + jitter, axis=1), s


class UniformPolicy(Policy):
    """Round-robin over blocks; with the configured budgets every block ends with equal shots."""

    name = "uniform"

    def choose(self, state):
        T, B = state.shots.shape
        return np.full(T, state.step % B), np.zeros((T, B))


class EIGPolicy(Policy):
    name = "eig_cell"

    def signal(self, state):
        s, dropped = eig_cell_support(state)
        state.diagnostics["eig_dropped_mass_max"] = float(dropped.max())
        return s


class EntropyPolicy(Policy):
    name = "entropy"

    def signal(self, state):
        return _per_block(state, lambda c: count_signals(c)["entropy"])


class RatioPolicy(Policy):
    name = "ratio"

    def signal(self, state):
        return _per_block(state, lambda c: count_signals(c)["ratio"])


class LegacyCwPolicy(Policy):
    name = "legacy_cw"

    def signal(self, state):
        return -_per_block(state, legacy_cw)  # low confidence = more need


class FisherPolicy(Policy):
    """Negative control: ideal per-shot Fisher information about phi, 4^k (4 pi^2/3)(M^2-1)."""

    name = "fisher"

    def signal(self, state):
        T = state.shots.shape[0]
        f = np.array([4.0**s.offset * (4 * np.pi**2 / 3) * (s.M**2 - 1) for s in state.specs])
        return np.tile(f, (T, 1))


class RandomPolicy(Policy):
    name = "random"

    def signal(self, state):
        return state.rng.random(state.shots.shape)


POLICIES = {p.name: p for p in (UniformPolicy, EIGPolicy, EntropyPolicy, RatioPolicy, LegacyCwPolicy, FisherPolicy, RandomPolicy)}
