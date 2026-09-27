"""Post-hoc characterisation of the n = 16 decline in extended-decoder overlap rescue (P6 held-out test).

Descriptive only, on existing held-out data; no claim enters any primary family.
For each shot-fixable failure of the P5-eig baseline (B2_p5_eig, limit decode correct) under
awqpe_ext with the frozen variant ext_v1_gated_A2, the responsible boundary j* is the internal
boundary whose chunk-j LSB scale 2^(n - K_j) is nearest (in log2) to the error in n-bit LSB units.
Reported by n: how often the practical overlap arm placed overlap at all, how often at j*, and
the rescue rate conditional on that.

    python -m research.analysis.scripts.p6_ext_decline --run research/results/p6_adaptive_overlap_test/<run>
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from research.awqpe.runner.core import REPO_ROOT

OUT = REPO_ROOT / "research" / "analysis" / "tables"
IDX = ["widths", "S0", "r", "dS", "phase_id", "replicate_id", "trial"]


def responsible_boundary(widths: str, err_lsb: float) -> int:
    w = [int(m) for m in widths.split("-")]
    n, K = sum(w), np.cumsum(w)[:-1]
    return int(np.argmin(np.abs(np.log2(max(err_lsb, 1e-9)) - (n - K))))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", required=True)
    ap.add_argument("--variant", default="ext_v1_gated_A2")
    ap.add_argument("--decoder", default="awqpe_ext")
    args = ap.parse_args()
    run = Path(args.run).resolve()
    cols = [*IDX, "arm", "n", "error", "tol", "limit_correct", "overlap_boundary_mask", "n_overlap_actions"]
    filt = [("row_type", "==", "final"), ("variant", "==", args.variant), ("decoder", "==", args.decoder),
            ("arm", "in", ["B2_p5_eig", "R_trigger_overlap"])]
    f = pd.concat([pd.read_parquet(x, columns=cols, filters=filt) for x in sorted((run / "shards").glob("*.parquet"))], ignore_index=True)
    b = f[f.arm == "B2_p5_eig"].set_index(IDX)
    o = f[f.arm == "R_trigger_overlap"].set_index(IDX).loc[b.index]
    d = b[["n", "error", "tol", "limit_correct"]].assign(ov_tol=o.tol, mask=o.overlap_boundary_mask, nact=o.n_overlap_actions)
    d = d[(d.tol < 0.5) & d.limit_correct.astype(bool)].reset_index()
    d["err_lsb"] = d.error * np.ldexp(1.0, d.n.to_numpy().astype(int))
    d["jstar"] = [responsible_boundary(w, e) for w, e in zip(d.widths, d.err_lsb)]
    d["overlap_at_jstar"] = [pd.notna(m) and (int(m) >> j) & 1 == 1 for m, j in zip(d["mask"], d.jstar)]
    d["any_overlap"] = d.nact > 0
    rows = []
    for (n, j), g in d.groupby(["n", "jstar"]):
        rows.append({"n": n, "jstar": j, "failures": len(g), "any_overlap": g.any_overlap.mean(), "overlap_at_jstar": g.overlap_at_jstar.mean(),
                     "rescue": g.ov_tol.mean(), "rescue_if_at_jstar": g[g.overlap_at_jstar].ov_tol.mean() if g.overlap_at_jstar.any() else np.nan,
                     "rescue_if_not": g[~g.overlap_at_jstar].ov_tol.mean() if (~g.overlap_at_jstar).any() else np.nan})
    out = pd.DataFrame(rows).assign(variant=args.variant, decoder=args.decoder, status="HELD-OUT TEST, post hoc descriptive",
                                    run=run.relative_to(REPO_ROOT).as_posix())
    out.to_csv(OUT / "p6_test_ext_decline.csv", index=False)
    print(out.drop(columns=["run", "status"]).round(3).to_string())


if __name__ == "__main__":
    main()
