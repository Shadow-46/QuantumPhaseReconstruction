"""Summarise the epsilon-floor verification run (PRELIMINARY, arXiv-v3 reading).

    python -m research.analysis.scripts.summarize_verify [--run DIR]

Writes research/analysis/tables/verify_*.csv, including every infinite-shot
failure record and the a-priori safe-epsilon table (docs/DECISIONS.md D-011).
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd

from research.awqpe.runner.core import REPO_ROOT, latest_run
from research.awqpe.verification.theory import partition_epsilon_star

OUT = REPO_ROOT / "research" / "analysis" / "tables"


def safe_epsilon(widths) -> float:
    """D-011: eps_safe = floor_{0.01}(min_j eps*(k_j, m_j)) - 0.02, fixed a priori (no tuning)."""
    return round(np.floor(partition_epsilon_star(widths) * 100) / 100 - 0.02, 2)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", default=None)
    args = ap.parse_args()
    run = Path(args.run) if args.run else latest_run("verify_epsilon_floor", include_pilot=False)
    load = lambda pat: pd.concat([pd.read_parquet(f) for f in glob.glob(str(run / "shards" / pat))], ignore_index=True)
    OUT.mkdir(parents=True, exist_ok=True)
    tag = run.relative_to(REPO_ROOT).as_posix()

    c = load("cross_*")
    cs = c.groupby(["widths", "grid_offset", "epsilon"]).agg(decodes=("agree", "size"), agree=("agree", "sum")).reset_index()
    cs["run"] = tag
    cs.to_csv(OUT / "verify_crosscheck.csv", index=False)

    s = load("sweep_*")
    s["fail_rate"] = s.observed_fail / s.n_phases
    s["run"] = tag
    s.to_csv(OUT / "verify_epsilon_sweep.csv", index=False)

    f = load("fail_*")
    f = f[f["phi"].notna()].copy()
    f["error_lsb"] = f.circular_error * np.ldexp(1.0, f.n.to_numpy().astype(int))
    f["run"] = tag
    f.to_csv(OUT / "verify_failure_records.csv.gz", index=False, compression="gzip")

    v = load("var_*")
    vv = v.groupby(["variant", "widths", "epsilon"]).agg(failures=("failures", "sum"), beyond_1lsb=("failures_beyond_1lsb", "sum"), n=("n_phases", "sum")).reset_index()
    vv["fail_rate"] = vv.failures / vv.n
    gold = load("golden*")[["variant", "golden_all_reproduce"]]
    vv = vv.merge(gold, on="variant")
    vv["run"] = tag
    vv.to_csv(OUT / "verify_variants.csv", index=False)

    rows = []
    gen = s[(s.decoder == "awqpe") & (s.grid_offset == 0.37)]
    for w, g in gen.groupby("widths"):
        widths = [int(x) for x in w.split("-")]
        es, star = safe_epsilon(widths), partition_epsilon_star(widths)
        at = g[np.isclose(g.epsilon, es)]
        onset = g[g.observed_fail > 0].epsilon.min()
        rows.append({"widths": w, "eps_star_min": round(star, 4), "eps_safe": es, "empirical_onset_generic_grid": onset,
                     "infinite_shot_failures_at_eps_safe": int(at.observed_fail.sum()) if len(at) else None,
                     "fail_rate_at_0.9": float(g[np.isclose(g.epsilon, 0.9)].fail_rate.iloc[0])})
    pd.DataFrame(rows).assign(rule="floor_0.01(min eps*) - 0.02", run=tag).to_csv(OUT / "verify_safe_epsilon.csv", index=False)
    print(f"wrote verify_*.csv to {OUT}")


if __name__ == "__main__":
    main()
