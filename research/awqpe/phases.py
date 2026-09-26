"""Stratified test-phase generation with a fixed development/test split.

Ambiguity in AWQPE is governed by where phi sits relative to each block's
measurement grid, so phases are not drawn only uniformly. Every stratum has
its own seed domain, and the development ("dev") and held-out ("test") pools
use disjoint seed domains, so no test phase can be seen while tuning.

Strata (n = sum of widths, K_j = cumulative bits through chunk j):
  S0_dyadic      phi = I / 2^n exactly.
  S1_final_half  phi = (I + t) / 2^n, t ~ U[0.45, 0.55]: best n-bit value near a tie.
  S2_boundary    phi = (I + d) / 2^(K_j), j uniform over internal boundaries,
                 d ~ U[0.45, 0.55]: chunk-j rounding ambiguity / special chunk.
  S3_near_grid   phi = (I + t) / 2^n, t ~ U[0.01, 0.10] or U[0.90, 0.99].
  S4_uniform     phi ~ U[0, 1).
  S5_paper       the phases used in the paper's Sec. 6 examples (dev only).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from research.awqpe.evaluation.metrics import boundary_hardness, final_residual
from research.awqpe.sim.seeds import generator, stable_int

STRATA = ("S0_dyadic", "S1_final_half", "S2_boundary", "S3_near_grid", "S4_uniform", "S5_paper")
SPLITS = {"dev": 0, "test": 1}
PAPER_PHASES = {
    "0.8203125": 0.8203125,
    "0.3": 0.3,
    "pi/6": float(np.pi / 6),
    "0.671875": 0.671875,
    "1/sqrt2": float(1 / np.sqrt(2)),
    "sin(pi/12)": float(np.sin(np.pi / 12)),
}


def _draw(stratum: str, count: int, widths, rng: np.random.Generator) -> np.ndarray:
    n = int(sum(widths))
    N = 1 << n
    if stratum == "S0_dyadic":
        return rng.integers(0, N, count) / N
    if stratum == "S1_final_half":
        return (rng.integers(0, N, count) + rng.uniform(0.45, 0.55, count)) / N
    if stratum == "S2_boundary":
        cuts = np.cumsum(widths)[:-1]
        if len(cuts) == 0:
            raise ValueError("S2_boundary needs at least two chunks.")
        K = cuts[rng.integers(0, len(cuts), count)]
        scale = np.ldexp(1.0, K)
        return (np.floor(rng.random(count) * scale) + rng.uniform(0.45, 0.55, count)) / scale
    if stratum == "S3_near_grid":
        t = np.where(rng.random(count) < 0.5, rng.uniform(0.01, 0.10, count), rng.uniform(0.90, 0.99, count))
        return (rng.integers(0, N, count) + t) / N
    if stratum == "S4_uniform":
        return rng.random(count)
    raise ValueError(f"unknown stratum {stratum!r}.")


def make_phase_table(widths, per_stratum: int, split: str, master_seed: int, strata=STRATA) -> pd.DataFrame:
    """Deterministic table of test phases with stratum features.

    Columns: phase_id, split, stratum, phi, final_residual, boundary_hardness.
    phase_id is stable across runs given (widths, split, master_seed).
    """
    if split not in SPLITS:
        raise ValueError("split must be 'dev' or 'test'.")
    wkey = stable_int("widths:" + ",".join(str(int(m)) for m in widths))
    rows = []
    for stratum in strata:
        if stratum == "S5_paper":
            if split == "test":
                continue  # the paper's own examples are reference cases, never held-out data
            phis, labels = np.array(list(PAPER_PHASES.values())), list(PAPER_PHASES)
        else:
            rng = generator(master_seed, SPLITS[split], stable_int(stratum), wkey)
            phis = np.mod(_draw(stratum, per_stratum, widths, rng), 1.0)
            labels = [str(i) for i in range(len(phis))]
        for label, phi in zip(labels, phis):
            rows.append({"phase_id": f"{split}:{stratum}:{label}", "split": split, "stratum": stratum, "phi": float(phi)})
    df = pd.DataFrame(rows)
    n = int(sum(widths))
    df["final_residual"] = final_residual(df["phi"].to_numpy(), n)
    df["boundary_hardness"] = boundary_hardness(df["phi"].to_numpy(), widths)
    return df
