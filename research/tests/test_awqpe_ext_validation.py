"""Validation of the EXPERIMENTAL decoder extension `awqpe_ext` (D-020). Not a P6 performance result.

A. no overlap block                      -> bit-identical to faithful AWQPE
B. overlap block present, trigger false  -> bit-identical to faithful AWQPE
C. trigger true                          -> only chunks at/above the triggered boundary can change
D. infinite shots                        -> every known eps-0.9 failure phase is repaired; no other phase changes
E. no true phase enters the trigger or the decoder
"""

import ast
import inspect
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from research.awqpe.baseline.awqpe import awqpe_vectorised
from research.awqpe.blocks.geometry import partition_blocks
from research.awqpe.decode import awqpe_overlap
from research.awqpe.decode.awqpe_overlap import awqpe_ext_decode
from research.awqpe.evaluation.metrics import best_nbit_integer, exact_success
from research.awqpe.model.kernel import block_probabilities
from research.awqpe.overlap.candidates import overlap_candidates
from research.awqpe.sim.oracle_sim import batch_block_counts

PARTITIONS = ([3, 2, 3], [2, 2, 2, 2], [4, 4], [3, 3, 3, 3])
TABLES = Path(__file__).resolve().parents[1] / "analysis" / "tables"


def _finite(widths, shots, T=3000, seed=0):
    """Finite-shot chunk and widened-block counts for random phases, plus tie jitter."""
    rng = np.random.default_rng(seed)
    phis = rng.random(T)
    chunks = [batch_block_counts(phis, s, shots, rng) for s in partition_blocks(widths)]
    jit = [rng.random(c.shape) * 0.5 for c in chunks]
    ext = {}
    for j, s, v in overlap_candidates(widths, "ext", 1):
        c = batch_block_counts(phis, s, shots, rng)
        ext[j] = (c, v, np.ones(T, bool), rng.random(c.shape) * 0.5)
    return phis, chunks, jit, ext


def _faithful_trigger(chunks, widths, jit, eps):
    """(T, B-1): faithful corrected lower part below boundary j is exactly 10..0."""
    corr = awqpe_vectorised(chunks, widths, eps, jitter=jit)["corrected"]
    B = len(widths)
    out = np.zeros((corr.shape[0], B - 1), dtype=bool)
    for j in range(B - 1):
        m = corr[:, j + 1] == (1 << (widths[j + 1] - 1))
        for k in range(j + 2, B):
            m &= corr[:, k] == 0
        out[:, j] = m
    return out


@pytest.mark.parametrize("widths", PARTITIONS)
@pytest.mark.parametrize("eps", [0.9, 0.6])
def test_A_no_overlap_is_faithful(widths, eps):
    _, chunks, jit, ext = _finite(widths, 4)
    base = awqpe_vectorised(chunks, widths, eps, jitter=jit)["estimate"]
    assert np.array_equal(base, awqpe_ext_decode(chunks, widths, eps, jit, {}))
    none_avail = {j: (c, v, np.zeros_like(h), jj) for j, (c, v, h, jj) in ext.items()}
    assert np.array_equal(base, awqpe_ext_decode(chunks, widths, eps, jit, none_avail))


@pytest.mark.parametrize("widths", PARTITIONS)
@pytest.mark.parametrize("shots", [2, 4, 16])
def test_B_and_C_trigger_semantics(widths, shots):
    eps = 0.9
    phis, chunks, jit, ext = _finite(widths, shots, seed=shots)
    base = awqpe_vectorised(chunks, widths, eps, jitter=jit)["estimate"]
    new = awqpe_ext_decode(chunks, widths, eps, jit, ext)
    trig = _faithful_trigger(chunks, widths, jit, eps)
    none = ~trig.any(axis=1)  # B: no triggered boundary -> identical
    assert none.sum() > 0 and np.array_equal(base[none], new[none])
    n = sum(widths)
    cuts = np.cumsum(widths)[:-1]
    changed = base != new  # C: a change implies a triggered boundary; bits below it untouched
    assert trig[changed].any(axis=1).all()
    for t in np.flatnonzero(changed):
        j = int(np.flatnonzero(trig[t]).max())  # least significant triggered boundary
        mask = (1 << (n - int(cuts[j]))) - 1
        assert (int(base[t]) & mask) == (int(new[t]) & mask)


def test_D_infinite_shot_repairs_known_failures_only():
    rec = pd.read_csv(TABLES / "verify_failure_records.csv.gz")
    rec = rec[(rec.epsilon == 0.9) & (rec.grid_offset == 0.37)]
    for w in PARTITIONS:
        n = sum(w)
        fail_phis = rec[rec.widths == "-".join(map(str, w))].phi.to_numpy()
        grid = (np.arange(1 << (n + 3)) + 0.37) / (1 << (n + 3))
        phis = np.concatenate([fail_phis, grid])
        T = len(phis)
        P = [block_probabilities(phis, s.offset, s.width) for s in partition_blocks(w)]
        J = [np.zeros_like(p) for p in P]
        ext = {j: (block_probabilities(phis, s.offset, s.width), v, np.ones(T, bool), np.zeros((T, s.M))) for j, s, v in overlap_candidates(w, "ext", 1)}
        base = awqpe_vectorised(P, w, 0.9, jitter=J)["estimate"]
        new = awqpe_ext_decode(P, w, 0.9, J, ext)
        assert len(fail_phis) > 0
        assert not exact_success(base[: len(fail_phis)], fail_phis, n).any()  # they really are failures
        assert exact_success(new, phis, n).all()  # every failure repaired, nothing broken
        ok_before = exact_success(base, phis, n)
        assert np.array_equal(base[ok_before], new[ok_before])  # correct phases untouched
        assert np.array_equal(new[: len(fail_phis)], best_nbit_integer(fail_phis, n))


def test_E_no_truth_in_trigger_or_decoder():
    forbidden = ("research.awqpe.evaluation", "research.awqpe.sim.oracle_sim", "research.awqpe.oracle")
    files = [Path(awqpe_overlap.__file__)] + list((Path(awqpe_overlap.__file__).parents[1] / "overlap").glob("*.py"))
    for f in files:
        for node in ast.walk(ast.parse(f.read_text(encoding="utf-8"))):
            mods = [a.name for a in node.names] if isinstance(node, ast.Import) else ([node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            assert not any(m.startswith(x) for m in mods for x in forbidden), (f.name, mods)
    assert set(inspect.signature(awqpe_ext_decode).parameters) == {"chunk_counts", "widths", "epsilon", "jitter", "ext"}
