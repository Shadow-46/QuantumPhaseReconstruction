"""Leakage firewall (test 28): decision-making code must not see the true phase.

Static check: modules under information/, allocation/, overlap/, stopping/
and controller.py may not import the evaluator (evaluation.metrics), the
simulator (sim.oracle_sim, whose PhaseSimulator holds truth), or the oracle
package, and may not reference truth attributes.
"""

import ast
from pathlib import Path

PKG = Path(__file__).resolve().parents[1] / "awqpe"
DECISION_DIRS = ("information", "allocation", "overlap", "stopping")
FORBIDDEN_MODULES = ("research.awqpe.evaluation", "research.awqpe.sim.oracle_sim", "research.awqpe.oracle")
FORBIDDEN_NAMES = {"true_phase", "_phases", "phi_true"}


def _decision_files():
    files = [p for d in DECISION_DIRS for p in (PKG / d).glob("*.py")]
    ctrl = PKG / "controller.py"
    return files + ([ctrl] if ctrl.exists() else [])


def test_decision_code_does_not_import_truth():
    offenders = []
    for path in _decision_files():
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            mods = []
            if isinstance(node, ast.Import):
                mods = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                mods = [node.module]
            for m in mods:
                if any(m == f or m.startswith(f + ".") for f in FORBIDDEN_MODULES):
                    offenders.append(f"{path.name}: imports {m}")
            if isinstance(node, ast.Attribute) and node.attr in FORBIDDEN_NAMES:
                offenders.append(f"{path.name}: uses .{node.attr}")
    assert not offenders, offenders


def test_signals_depend_only_on_counts():
    """Signals are a deterministic function of counts (no hidden state or truth)."""
    import numpy as np

    from research.awqpe.blocks.geometry import partition_blocks
    from research.awqpe.decode.likelihood import GridPosterior
    from research.awqpe.information.signals import all_block_signals

    specs = partition_blocks([3, 3])
    counts = [np.array([[5, 3, 0, 0, 0, 0, 0, 0]]), np.array([[0, 0, 7, 1, 0, 0, 0, 0]])]
    results = []
    for _ in range(2):
        gp = GridPosterior(6, 1)
        for s, c in zip(specs, counts):
            gp.add_counts(s, c)
        results.append(all_block_signals(specs, counts, 6, gp))
    assert all(np.array_equal(results[0][k], results[1][k]) for k in results[0])
