"""Unified adaptive controller (Phase 7, docs/DECISIONS.md D-031).

One decision per step, per trial, among: a batch of S0 shots on chunk b, a batch of S0 shots
on an overlap block (O-ext or O-bridge candidate; at most A overlap batches per trial), or
stop. The action is the eligible block with the largest expected information gain about the
tau-resolution cell (P5 signal, cell bits = tau_bits(tau)); the stop is the P8 posterior-credible
rule. Decoder: D2 grid posterior over all shots of all blocks.

Decision code: sees counts and the posterior only (test_leakage.py). With no overlap candidates
it performs exactly the arithmetic of run_tolerance.run_stopping(policy="eig") (tested).
"""

from __future__ import annotations

import numpy as np

from research.awqpe.allocation.eig_cached import eig_cell_support
from research.awqpe.allocation.policies import AllocationState
from research.awqpe.decode.likelihood import GridPosterior
from research.awqpe.stopping.rules import should_stop, tau_bits


def _batch(stream, start, dS, M):
    T = stream.shape[0]
    ys = np.take_along_axis(stream, start[:, None] + np.arange(dS)[None, :], axis=1)
    c = np.zeros((T, M), dtype=np.int64)
    np.add.at(c, (np.repeat(np.arange(T), dS), ys.ravel()), 1)
    return c


def unified_run(chunks, cands, streams, n, refine, S0, cap, A, tau, alpha, rng):
    """Returns MAP (T,), shots (T, B + C), stopped (T,), overlap batches used (T,)."""
    specs = list(chunks) + list(cands)
    B, K = len(chunks), len(chunks) + len(cands)
    T = streams[0].shape[0]
    shots = np.zeros((T, K), dtype=np.int64)
    shots[:, :B] = S0
    counts = [_batch(streams[b], np.zeros(T, dtype=np.int64), S0, specs[b].M) if b < B else np.zeros((T, specs[b].M), dtype=np.int64)
              for b in range(K)]
    gp = GridPosterior(n, T, refine)
    for s, c in zip(specs[:B], counts[:B]):
        gp.add_counts(s, c)
    stopped, _ = should_stop(gp.loglik, tau, alpha)
    cell_bits = min(n, tau_bits(tau))
    n_ov = np.zeros(T, dtype=np.int64)
    step = 0
    while True:
        E = np.empty((T, K), dtype=bool)
        E[:, :B] = shots[:, :B] + S0 <= cap
        E[:, B:] = (n_ov < A)[:, None]
        active = ~stopped & E.any(axis=1)
        if not active.any():
            break
        ia = np.flatnonzero(active)
        st = AllocationState(specs, cell_bits, [c[ia] for c in counts], shots[ia], gp.loglik[ia], gp.G, 0, rng, step)
        sig, _ = eig_cell_support(st)
        sig = np.where(E[ia], sig + rng.random(sig.shape) * 1e-12, -np.inf)
        chosen = np.full(T, -1)
        chosen[ia] = np.argmax(sig, axis=1)
        for b in range(K):
            m = active & (chosen == b)
            if not m.any():
                continue
            d = _batch(streams[b], np.minimum(shots[:, b], streams[b].shape[1] - S0), S0, specs[b].M)
            d[~m] = 0
            counts[b] += d
            shots[m, b] += S0
            if b >= B:
                n_ov[m] += 1
            gp.add_counts(specs[b], d)
        step += 1
        live = ~stopped
        s_new, _ = should_stop(gp.loglik[live], tau, alpha)
        stopped[np.flatnonzero(live)[s_new]] = True
    return gp.map_estimate(), shots, stopped, n_ov
