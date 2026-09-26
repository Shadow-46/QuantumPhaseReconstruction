"""Observable per-window information signals.

Every function here takes only observed counts, the declared block geometry
and model quantities. None may import evaluation/metrics or read the true
phase (enforced by research/tests/test_leakage.py).

Signal families (all vectorised over T trials):
  count      c1, c2 (as fractions of the block's shots), ratio c2/c1,
             margin (c1-c2), empirical entropy (bits)
  local      block-local posterior over the fractional phase delta_b that the
             block sees (uniform prior, kernel likelihood of that block only):
             circular posterior sd (in outcome units), P(chunk correct) =
             P(round(M delta) = t1 | data), boundary risk =
             P(|frac(M delta) - 1/2| < BOUNDARY_HALF_WIDTH | data)
  fisher     per-shot Fisher information about phi at the local MLE
             (ideal kernel: constant 4^k (4 pi^2/3)(M^2-1), see RESULTS_LOG F-1)
  global     from the joint grid posterior over phi (decoder D2's model):
             eig_phi  = I(Y_b ; phi | data)   (bits per extra shot on block b)
             eig_cell = I(Y_b ; n-bit cell of phi | data)
  legacy     Paper A's C_w (independent Beta posteriors, Jeffreys prior),
             imported read-only from the frozen Reconstruction/confidence.py
"""

from __future__ import annotations

import sys
from functools import lru_cache
from pathlib import Path

import numpy as np

from research.awqpe.blocks.geometry import BlockSpec
from research.awqpe.decode.likelihood import GridPosterior, _log_table
from research.awqpe.model.fisher import block_fisher
from research.awqpe.model.noise import NoiseSpec, noisy_block_probabilities

_REPO = Path(__file__).resolve().parents[3]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

BOUNDARY_HALF_WIDTH = 0.15  # in outcome units; fixed a priori (docs/DECISIONS.md D-015)
LOCAL_REFINE = 6  # local delta grid has 2^(m + LOCAL_REFINE) points


def _entropy_bits(p: np.ndarray, axis: int = -1) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        return -np.sum(np.where(p > 0, p * np.log2(p), 0.0), axis=axis)


def count_signals(counts: np.ndarray) -> dict[str, np.ndarray]:
    c = np.asarray(counts, dtype=float)
    S = c.sum(axis=1)
    top = -np.sort(-c, axis=1)
    c1, c2 = top[:, 0] / S, top[:, 1] / S
    return {
        "c1": c1, "c2": c2, "ratio": np.divide(c2, c1, out=np.zeros_like(c1), where=c1 > 0),
        "margin": c1 - c2, "entropy": _entropy_bits(c / S[:, None]),
    }


@lru_cache(maxsize=64)
def _local_log_table(width: int, noise: NoiseSpec | None, offset: int) -> tuple[np.ndarray, np.ndarray]:
    M = 1 << width
    grid = np.arange(M << LOCAL_REFINE) / (M << LOCAL_REFINE)
    p = noisy_block_probabilities(np.ldexp(grid, -offset), offset, width, noise)
    return grid, np.log(np.maximum(p, 1e-300))


def local_posterior_signals(spec: BlockSpec, counts: np.ndarray, noise: NoiseSpec | None = None) -> dict[str, np.ndarray]:
    counts = np.asarray(counts, dtype=float)
    M = spec.M
    grid, logp = _local_log_table(spec.width, noise, spec.offset)
    ll = counts @ logp.T
    ll -= ll.max(axis=1, keepdims=True)
    post = np.exp(ll)
    post /= post.sum(axis=1, keepdims=True)
    x = grid * M  # delta in outcome units
    R = np.abs(post @ np.exp(2j * np.pi * grid))
    circ_sd = np.sqrt(np.maximum(-2 * np.log(np.clip(R, 1e-300, 1.0)), 0.0)) / (2 * np.pi) * M
    t1 = np.argmax(counts, axis=1)  # ties: lowest index (posterior mass is what matters)
    nearest = np.mod(np.floor(x + 0.5), M).astype(np.int64)
    p_correct = np.sum(post * (nearest[None, :] == t1[:, None]), axis=1)
    frac = x - np.floor(x)
    boundary = np.sum(post * (np.abs(frac - 0.5) < BOUNDARY_HALF_WIDTH)[None, :], axis=1)
    mle = grid[np.argmax(ll, axis=1)]
    fisher = block_fisher(np.ldexp(mle, -spec.offset), spec.offset, spec.width, noise)
    return {"local_post_sd": circ_sd, "p_chunk_correct": p_correct, "chunk_risk": 1.0 - p_correct,
            "boundary_risk": boundary, "fisher_phi": fisher}


def global_eig_signals(gp: GridPosterior, specs, n: int, noise: NoiseSpec | None = None) -> dict[str, np.ndarray]:
    """Per-block expected information of ONE extra shot, from the current grid posterior.

    Returns arrays of shape (T, B): eig_phi (bits about phi) and eig_cell
    (bits about the n-bit cell round(2^n phi), which decides tolerance-2^-n success).
    """
    post = gp.posterior()  # (T, G)
    T, G = post.shape
    sub = G >> n
    cell_of = np.mod((np.arange(G) + sub // 2) // sub, 1 << n)  # cell centred on j/2^n
    order = np.argsort(cell_of, kind="stable")
    eig_phi = np.zeros((T, len(specs)))
    eig_cell = np.zeros((T, len(specs)))
    post_c = post[:, order].reshape(T, 1 << n, sub)
    mass_c = post_c.sum(axis=2)
    for b, s in enumerate(specs):
        P = np.exp(_log_table(G, s.offset, s.width, noise))  # (G, M)
        pred = post @ P
        H_pred = _entropy_bits(pred)
        eig_phi[:, b] = H_pred - post @ _entropy_bits(P)
        P_c = P[order].reshape(1 << n, sub, s.M)
        pred_c = np.einsum("tcs,csm->tcm", post_c, P_c)
        cond = np.divide(pred_c, mass_c[:, :, None], out=np.zeros_like(pred_c), where=mass_c[:, :, None] > 0)
        eig_cell[:, b] = H_pred - np.sum(mass_c * _entropy_bits(cond), axis=1)
    return {"eig_phi": eig_phi, "eig_cell": eig_cell}


def global_risk(gp: GridPosterior, n: int, tau: float | None = None) -> np.ndarray:
    """1 - posterior mass within tau (default 2^-n) of the MAP: a trial-level signal."""
    tau = 2.0**-n if tau is None else tau
    return 1.0 - gp.credible_mass(gp.map_estimate(), tau)


@lru_cache(maxsize=200000)
def _cw(n1: int, n2: int, s: int) -> float:
    from Reconstruction.confidence import top1_vs_top2_confidence  # frozen Paper A code, read-only

    return float(top1_vs_top2_confidence(n1, n2, s))


def legacy_cw(counts: np.ndarray) -> np.ndarray:
    c = np.asarray(counts, dtype=np.int64)
    top = -np.sort(-c, axis=1)
    S = c.sum(axis=1)
    return np.array([_cw(int(a), int(b), int(s)) for a, b, s in zip(top[:, 0], top[:, 1], S)])


def all_block_signals(specs, counts_list, n: int, gp: GridPosterior, noise: NoiseSpec | None = None, with_cw: bool = True) -> dict[str, np.ndarray]:
    """Every signal for every block: dict of arrays shaped (T, B)."""
    out: dict[str, list] = {}
    for s, c in zip(specs, counts_list):
        sig = {**count_signals(c), **local_posterior_signals(s, c, noise)}
        if with_cw:
            sig["legacy_cw"] = legacy_cw(c)
        for k, v in sig.items():
            out.setdefault(k, []).append(v)
    res = {k: np.stack(v, axis=1) for k, v in out.items()}
    res.update(global_eig_signals(gp, specs, n, noise))
    return res
