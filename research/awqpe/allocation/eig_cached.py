"""P6-local, memoised copy of the frozen P5 signal `policies.eig_cell_support`.

The frozen P5 function recomputes np.exp(log-table) over the whole (G x M)
grid table on every call; at n = 16 that dominates Phase 6 runtime. This copy
performs exactly the same arithmetic in the same order, but computes each
exponentiated table once per (G, offset, width). It is used only by the Phase 6
runner; P5 code is not modified. Bit-identity with the P5 function is asserted
in research/tests/test_p6_analysis.py.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

from research.awqpe.allocation.policies import SUPPORT_MAX, SUPPORT_TAIL, AllocationState, _entropy_bits
from research.awqpe.decode.likelihood import MAX_TABLE_ELEMENTS, _log_table
from research.awqpe.model.kernel import kernel_closed_form


@lru_cache(maxsize=64)
def _prob_table(G: int, offset: int, width: int) -> np.ndarray:
    return np.exp(_log_table(G, offset, width, None))


def _block_probs_at(G: int, spec, idx: np.ndarray) -> np.ndarray:
    if G * spec.M <= MAX_TABLE_ELEMENTS:
        return _prob_table(G, spec.offset, spec.width)[idx]
    delta = ((idx.astype(np.int64) << spec.offset) % G) / G
    y = np.arange(spec.M) / spec.M
    return kernel_closed_form(delta[..., None] - y, spec.M)


def eig_cell_support(state: AllocationState) -> tuple[np.ndarray, np.ndarray]:
    """Identical to allocation.policies.eig_cell_support (see module docstring)."""
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
