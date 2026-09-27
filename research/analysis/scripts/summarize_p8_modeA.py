"""Phase 8 Mode A (D-028): tolerance success vs budget, a re-analysis of the frozen P5 held-out test.

No new data. For tau = 2^-j, j in {3, n/2, n}, success = circular error <= tau. For each
P5 cell and decoder, the paired difference eig_cell - uniform (points) is computed per phase
cluster, then pooled over cells with equal cell weight per (n, j, decoder); the CI is a
phase-cluster bootstrap of the within-cell-centred cluster differences.
Secondary, no multiplicity claims (the P5 primary family used tau = 2^-n only).

    python -m research.analysis.scripts.summarize_p8_modeA --run research/results/p5_adaptive_shots_test/<run>
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from research.analysis.scripts.summarize_p5 import boot_ci
from research.awqpe.runner.core import REPO_ROOT

OUT = REPO_ROOT / "research" / "analysis" / "tables"
CELL = ["widths", "n", "S0", "r", "dS"]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", required=True)
    ap.add_argument("--boot", type=int, default=2000)
    args = ap.parse_args()
    run = Path(args.run).resolve()
    cols = [*CELL, "decoder", "policy", "phase_id", "replicate_id", "error"]
    f = pd.concat([pd.read_parquet(x, columns=cols, filters=[("row_type", "==", "final"), ("policy", "in", ["uniform", "eig_cell"])])
                   for x in sorted((run / "shards").glob("*.parquet"))], ignore_index=True)
    f["cluster"] = f.widths + "|" + f.phase_id
    rows = []
    for n, g in f.groupby("n"):
        for j in sorted({3, int(n) // 2, int(n)}):
            g = g.assign(ok=(g.error <= 2.0**-j + 1e-15).astype(float))
            p = g.pivot_table(index=[*CELL, "decoder", "cluster", "replicate_id"], columns="policy", values="ok").dropna()
            p["d"] = p.eig_cell - p.uniform
            per_cell = p.groupby([*CELL, "decoder", "cluster"]).d.mean().reset_index()
            for dec, h in per_cell.groupby("decoder"):
                cm = h.groupby(CELL).d.transform("mean")
                point = float(h.groupby(CELL).d.mean().mean())
                lo, hi = boot_ci((h.d - cm).to_numpy() + point, args.boot)
                base = p.xs(dec, level="decoder")
                rows.append({"n": int(n), "tau_exp": j, "decoder": dec, "cells": h.groupby(CELL).ngroups,
                             "p_uniform": float(base.uniform.mean()), "p_eig_cell": float(base.eig_cell.mean()),
                             "diff_pts": 100 * point, "lo_pts": 100 * lo, "hi_pts": 100 * hi})
    out = pd.DataFrame(rows).assign(status="HELD-OUT TEST (re-analysis of frozen P5 test; secondary)", run=run.relative_to(REPO_ROOT).as_posix())
    out.to_csv(OUT / "p8_modeA_p5test_tolerance.csv", index=False)
    print(out.drop(columns=["status", "run"]).round(3).to_string())


if __name__ == "__main__":
    main()
