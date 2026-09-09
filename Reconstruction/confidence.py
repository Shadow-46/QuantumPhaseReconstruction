"""
Purpose
    Quantify how confidently a window's top-1 observed bit pattern beats its
    top-2 runner-up, given the shot count actually collected.
Theory
    Raw probability/margin/entropy over observed counts are sample-size-blind
    (a 3-vs-1 split on 4 shots and a 750-vs-250 split on 1000 score identically
    under all three), which misses the dominant measured failure mode:
    shot-starvation, not genuine ambiguity (see Data/failure_study/REPORT.md
    Sections 2, 6 and FORMULATION.md Section 2). Modeling each candidate's true
    sampling probability with a Beta(alpha0+n_i, alpha0+S_w-n_i) posterior
    (Jeffreys prior, alpha0=0.5) and computing Pr[X_1 > X_2] for two such
    independent variables gives a single, sample-size-aware score. Only the
    top-2 observed patterns are compared, a documented simplifying assumption,
    not a validated approximation (FORMULATION.md Section 2).

    Two variants are provided. `top1_vs_top2_confidence` is the SHIPPED
    statistic: it treats the two marginals as independent Binomials, which is
    false -- they are negatively correlated coordinates of one K-way
    multinomial -- and it has no elementary closed form for the half-integer
    shape parameters that arise here, so it is evaluated by quadrature.
    `top1_vs_top2_confidence_dirichlet` is the EXACT posterior under the same
    sampling model, obtained by dropping that independence assumption; it does
    have a closed form. The Dirichlet version is diagnostic only and is
    deliberately not wired into any call site -- see its docstring.
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
from scipy.special import betainc
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


def top1_vs_top2_confidence_dirichlet(n1: int, n2: int, alpha0: float = 0.5) -> float:
    """Return Pr[p1 > p2] under the EXACT joint posterior, in closed form.

    DIAGNOSTIC ONLY. This is deliberately not called by any reconstruction
    path. The shipped method thresholds `top1_vs_top2_confidence`, and every
    frozen result in both articles was produced with that statistic; swapping
    it would change module D's behaviour on ~3% of windows and invalidate the
    released ablation and generalisation corpora. The swap is pre-registered
    for the next campaign, on the same reasoning that keeps the isotonic
    recalibration out of the shipped rule (see adaptive_reconstruction.py).

    Derivation. Under counts n ~ Multinomial(S_w, p) over K patterns with a
    symmetric Dirichlet(alpha0) prior, the posterior is p | n ~ Dir(alpha0+n).
    Aggregating the K-2 non-top-2 categories gives
        (p1, p2, rest) ~ Dir(alpha0+n1, alpha0+n2, (K-2)*alpha0 + S_w-n1-n2).
    By Dirichlet neutrality, V := p1/(p1+p2) ~ Beta(alpha0+n1, alpha0+n2) and
    is independent of p1+p2. Since p1+p2 > 0 a.s., {p1 > p2} = {V > 1/2}, so
        Pr[p1 > p2] = 1 - I_{1/2}(alpha0+n1, alpha0+n2),
    with I the regularized incomplete Beta function. One `betainc` call: no
    quadrature, no tolerance parameter, no clamping.

    Note the result depends only on (n1, n2) -- not on S_w and not on K. That
    is correct rather than a defect: conditional on n1+n2, the residual counts
    are ancillary for the ratio p1/(p1+p2), whose sign determines the event.
    Sample-size sensitivity is retained in the currency that matters for a
    pairwise comparison, n1+n2: the same 3:2 proportion scores 0.670 on five
    shots and 1.000 on a thousand. What is lost relative to the shipped
    statistic is an incidental penalty on windows carrying large mass outside
    the top two -- an accidental proxy for the top-2-only restriction (A3),
    not a modelled quantity."""
    if n1 < 0 or n2 < 0:
        raise ValueError("counts must satisfy 0 <= n1, n2.")
    # betainc(a, b, x) is the regularized incomplete Beta function I_x(a, b).
    return float(1.0 - betainc(alpha0 + n1, alpha0 + n2, 0.5))


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

    # --- Dirichlet-joint variant (diagnostic; not on any shipped path) ---
    import numpy as np

    # Symmetry and the tie value, matching the shipped statistic's identities.
    assert abs(top1_vs_top2_confidence_dirichlet(10, 10) - 0.5) < 1e-12
    d_fwd = top1_vs_top2_confidence_dirichlet(6, 4)
    d_rev = top1_vs_top2_confidence_dirichlet(4, 6)
    assert abs(d_fwd + d_rev - 1.0) < 1e-12
    assert top1_vs_top2_confidence_dirichlet(3, 1) > 0.5
    assert top1_vs_top2_confidence_dirichlet(1, 3) < 0.5

    # Sample-size sensitivity is retained in n1+n2, the informative count.
    assert top1_vs_top2_confidence_dirichlet(3, 2) < top1_vs_top2_confidence_dirichlet(600, 400)
    assert top1_vs_top2_confidence_dirichlet(600, 400) > 0.999

    # The closed form matches direct Dirichlet Monte Carlo, and is invariant
    # to the number of categories and to mass outside the top two -- the two
    # quantities the shipped statistic incorrectly responds to.
    rng = np.random.default_rng(0)
    for counts in ([3, 1], [3, 1, 0, 0], [5, 3, 1, 1] + [0] * 8, [10, 7, 2, 1] + [0] * 12):
        draws = rng.dirichlet(0.5 + np.asarray(counts, dtype=float), size=400_000)
        mc = float((draws[:, 0] > draws[:, 1]).mean())
        closed = top1_vs_top2_confidence_dirichlet(counts[0], counts[1])
        assert abs(mc - closed) < 5e-3, f"{counts}: MC {mc:.4f} vs closed {closed:.4f}"

    # The shipped statistic converges to the exact one as the mass outside the
    # top two grows, which is why the two agree on well-sampled windows and
    # differ most where the top two take all the shots.
    exact = top1_vs_top2_confidence_dirichlet(3, 1)
    gap_small = abs(top1_vs_top2_confidence(3, 1, 4) - exact)
    gap_large = abs(top1_vs_top2_confidence(3, 1, 2048) - exact)
    assert gap_large < gap_small
    assert gap_large < 1e-3


if __name__ == "__main__":
    self_test()
    print("confidence self-test passed")
