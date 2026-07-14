"""
Purpose
    Quantify how confidently a window's top-1 observed bit pattern beats its
    top-2 runner-up, given the shot count actually collected.
Theory
    Raw probability/margin/entropy over observed counts are sample-size-blind
    (a 2/4 split and a 500/1000 split score identically), which misses the
    dominant measured failure mode: shot-starvation, not genuine ambiguity
    (see Data/failure_study/REPORT.md Sections 2, 6 and FORMULATION.md
    Section 2). Modeling each candidate's true sampling probability with a
    Beta(alpha0+n_i, alpha0+S_w-n_i) posterior (Jeffreys prior, alpha0=0.5)
    and computing Pr[X_1 > X_2] for two such independent variables gives a
    single, exact, sample-size-aware confidence score with a closed form via
    the regularized incomplete Beta function (FORMULATION.md Section 2). Only
    the top-2 observed patterns are compared, a documented simplifying
    assumption, not a validated approximation (FORMULATION.md Section 2).
Inputs
    Observed counts n1, n2 for a window's top-1/top-2 bit patterns and the
    window's total shot count S_w.
Outputs
    C_w in [0, 1]: the posterior probability that pattern 1's true sampling
    probability exceeds pattern 2's.
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scipy.integrate import quad
from scipy.stats import beta as beta_dist


def top1_vs_top2_confidence(n1: int, n2: int, s_w: int, alpha0: float = 0.5) -> float:
    """Return Pr[X1 > X2] for independent Beta(alpha0+n1, alpha0+s_w-n1) and
    Beta(alpha0+n2, alpha0+s_w-n2) posteriors: Pr[X1 > X2] =
    integral_0^1 f1(x) F2(x) dx, where F2 is X2's CDF (the regularized
    incomplete Beta function, scipy's `betainc` under `beta.cdf`) and f1 is
    X1's density. Evaluated by deterministic numerical quadrature over both
    closed-form scipy functions — no Monte Carlo, no regime branching."""
    if n1 < 0 or n2 < 0 or s_w < 0 or n1 > s_w or n2 > s_w:
        raise ValueError("counts must satisfy 0 <= n1, n2 <= s_w.")
    a1, b1 = alpha0 + n1, alpha0 + s_w - n1
    a2, b2 = alpha0 + n2, alpha0 + s_w - n2

    def integrand(x: float) -> float:
        return beta_dist.pdf(x, a1, b1) * beta_dist.cdf(x, a2, b2)

    value, _ = quad(integrand, 0.0, 1.0, limit=200)
    return float(min(1.0, max(0.0, value)))


def self_test() -> None:
    """Verify symmetry, monotonicity, and sample-size sensitivity of C_w."""
    assert abs(top1_vs_top2_confidence(10, 10, 20) - 0.5) < 1e-6
    assert top1_vs_top2_confidence(3, 1, 4) > 0.5
    assert top1_vs_top2_confidence(1, 3, 4) < 0.5
    # Symmetry: swapping n1/n2 should reflect C_w around 0.5.
    c_fwd = top1_vs_top2_confidence(6, 4, 10)
    c_rev = top1_vs_top2_confidence(4, 6, 10)
    assert abs(c_fwd + c_rev - 1.0) < 1e-6
    # Sample-size sensitivity: the same 60/40 split is far more confident at
    # large S_w than at small S_w (the property raw margin/entropy lack).
    c_small = top1_vs_top2_confidence(3, 2, 5)
    c_large = top1_vs_top2_confidence(600, 400, 1000)
    assert c_small < c_large
    assert c_large > 0.999


if __name__ == "__main__":
    self_test()
    print("confidence self-test passed")
