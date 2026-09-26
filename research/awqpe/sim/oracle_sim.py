"""Dirichlet-kernel measurement simulator with a truth firewall.

`PhaseSimulator` holds the ground truth (eigenphase, or an eigen-decomposition
mixture for non-eigenstate inputs) and the noise model. Adaptive algorithms
never receive it. They receive a `MeasurementOracle`, whose only capability is
measure(spec, shots) -> counts; truth is captured in a closure, not stored as
an attribute, and the oracle records the quantum cost it incurred.

Non-eigenstate inputs. If the target register is prepared in
sum_j c_j |u_j>, every shot of every block collapses independently onto one
eigencomponent, so each block's outcome distribution is the mixture
sum_j |c_j|^2 p(y | phi_j). Different blocks do not share a component: this is
exactly what happened to Paper A's |1> input to U_a (docs/HISTORICAL_SETUP_RECOVERY.md).

Common random numbers. In "crn" mode each block geometry (offset, width) in a
trial owns a private uniform stream; shot i of that block is always
searchsorted(cdf, u_i). Two policies that take the first s shots of the same
block therefore observe identical outcomes, whatever order they request
blocks in. "multinomial" mode draws counts per request from the same private
stream and is faster for large fixed allocations, but only pairs requests of
identical batch sizes.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from research.awqpe.blocks.geometry import BlockSpec
from research.awqpe.model.noise import NoiseSpec, noisy_block_probabilities
from research.awqpe.sim.seeds import generator


@dataclass
class CostLedger:
    """Quantum resources consumed through one oracle."""

    shots: int = 0
    u_queries: int = 0
    circuits: int = 0
    max_control_qubits: int = 0
    max_power: int = 0
    per_block_shots: dict = field(default_factory=dict)

    def charge(self, spec: BlockSpec, shots: int) -> None:
        self.shots += shots
        self.u_queries += shots * spec.u_queries_per_shot
        self.circuits += 1
        self.max_control_qubits = max(self.max_control_qubits, spec.width)
        self.max_power = max(self.max_power, spec.max_power)
        self.per_block_shots[spec.key] = self.per_block_shots.get(spec.key, 0) + shots


class MeasurementOracle:
    """The only interface an (adaptive) AWQPE algorithm may use."""

    def __init__(self, measure_fn, ledger: CostLedger):
        self._measure_fn = measure_fn
        self.ledger = ledger

    def measure(self, spec: BlockSpec, shots: int) -> np.ndarray:
        if shots < 1:
            raise ValueError("shots must be positive.")
        counts = self._measure_fn(spec, int(shots))
        self.ledger.charge(spec, int(shots))
        return counts


class PhaseSimulator:
    """Ground-truth holder. Only simulators and evaluators may touch it."""

    def __init__(self, phases, weights=None, noise: NoiseSpec | None = None, seed_parts: tuple[int, ...] = (0,), mode: str = "crn"):
        self._phases = np.atleast_1d(np.asarray(phases, dtype=float))
        w = np.ones_like(self._phases) if weights is None else np.asarray(weights, dtype=float)
        if w.shape != self._phases.shape or np.any(w < 0) or w.sum() <= 0:
            raise ValueError("weights must be non-negative and match phases.")
        self._weights = w / w.sum()
        self._noise = noise
        self._seed_parts = tuple(int(s) for s in seed_parts)
        if mode not in {"crn", "multinomial"}:
            raise ValueError("mode must be 'crn' or 'multinomial'.")
        self._mode = mode
        self._streams: dict[tuple[int, int], np.random.Generator] = {}
        self._cdfs: dict[tuple[int, int], np.ndarray] = {}

    # --- truth side (evaluator only) ---
    @property
    def true_phase(self) -> float:
        """Dominant eigenphase (the estimation target)."""
        return float(self._phases[int(np.argmax(self._weights))])

    def distribution(self, spec: BlockSpec) -> np.ndarray:
        p = noisy_block_probabilities(self._phases, spec.offset, spec.width, self._noise)
        return np.clip(self._weights @ p, 0.0, None)

    # --- measurement side ---
    def _stream(self, spec: BlockSpec) -> np.random.Generator:
        key = spec.key
        if key not in self._streams:
            self._streams[key] = generator(*self._seed_parts, *key)
        return self._streams[key]

    def _measure(self, spec: BlockSpec, shots: int) -> np.ndarray:
        rng = self._stream(spec)
        if self._mode == "multinomial":
            p = self.distribution(spec)
            return rng.multinomial(shots, p / p.sum()).astype(np.int64)
        key = spec.key
        if key not in self._cdfs:
            cdf = np.cumsum(self.distribution(spec))
            self._cdfs[key] = cdf / cdf[-1]
        y = np.searchsorted(self._cdfs[key], rng.random(shots), side="right")
        return np.bincount(np.minimum(y, spec.M - 1), minlength=spec.M).astype(np.int64)

    def oracle(self) -> MeasurementOracle:
        """A fresh oracle sharing this simulator's shot streams."""
        return MeasurementOracle(self._measure, CostLedger())


def batch_block_counts(phis: np.ndarray, spec: BlockSpec, shots: int, rng: np.random.Generator, noise: NoiseSpec | None = None, weights: np.ndarray | None = None) -> np.ndarray:
    """Vectorised multinomial counts for many independent trials.

    phis has shape (T,) for eigenstate inputs or (T, J) with `weights` (T, J)
    for mixtures. Returns int64 counts of shape (T, 2^width). Intended for
    non-adaptive Monte Carlo where every trial uses the same fixed schedule.
    """
    phis = np.asarray(phis, dtype=float)
    p = noisy_block_probabilities(phis, spec.offset, spec.width, noise)
    if p.ndim == 3:
        if weights is None:
            raise ValueError("mixture phases need weights.")
        w = np.asarray(weights, dtype=float)
        p = np.einsum("tj,tjy->ty", w / w.sum(axis=1, keepdims=True), p)
    p = np.clip(p, 0.0, None)
    p = p / p.sum(axis=-1, keepdims=True)
    return rng.multinomial(shots, p).astype(np.int64)
