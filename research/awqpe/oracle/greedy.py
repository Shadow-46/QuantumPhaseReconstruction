"""Greedy one-step ORACLE allocation: an analysis reference, never a practical policy.

At each decision it inspects the simulator's realized counterfactuals (the
next dS shots that each block WOULD receive, scored against the true phase
under a given decoder) and picks the block with the largest realized
improvement: first in tolerance success, then in circular error. It is
myopic (one step), so it is an achievable-by-hindsight reference, not a
global optimum. Every row it produces is tagged is_oracle=True.
"""

from __future__ import annotations

import numpy as np

IS_ORACLE = True


def greedy_oracle_choice(gain_tol: np.ndarray, gain_err_lsb: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """gain_tol in {-1, 0, 1} and gain_err_lsb (error reduction in LSB units), both (T, B)."""
    score = gain_tol * 1e6 + gain_err_lsb + rng.random(gain_tol.shape) * 1e-9
    return np.argmax(score, axis=1)
