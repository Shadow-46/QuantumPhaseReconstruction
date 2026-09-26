"""Summarise Phase 2 runs into committed CSV tables (PRELIMINARY: dev split).

    python -m research.analysis.scripts.summarize_p2 [--limit-run DIR] [--mc-run DIR]

Writes research/analysis/tables/p2a_limit_failure.csv and
p2b_mc_success.csv with Wilson 95% intervals. Phases, not replicates, are the
independent units in P2b, so the table also reports a phase-cluster bootstrap
interval.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from research.awqpe.runner.core import REPO_ROOT, latest_run, load_results

OUT = REPO_ROOT / "research" / "analysis" / "tables"


def wilson(k, n, z: float = 1.959964):
    k, n = np.asarray(k, dtype=float), np.asarray(n, dtype=float)
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return centre - half, centre + half


def cluster_bootstrap(df: pd.DataFrame, col: str, B: int = 400, seed: int = 0) -> tuple[float, float]:
    per_phase = df.groupby("phase_id")[col].mean().to_numpy()
    rng = np.random.default_rng(seed)
    means = per_phase[rng.integers(0, len(per_phase), (B, len(per_phase)))].mean(axis=1)
    return float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def decoder_label(d: pd.DataFrame) -> pd.Series:
    return d["decoder"] + np.where(d["epsilon"].notna(), "@" + d["epsilon"].astype(str), "")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--limit-run", default=None)
    ap.add_argument("--mc-run", default=None)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    lim_dir = Path(args.limit_run).resolve() if args.limit_run else latest_run("p2a_infinite_shot_limit", include_pilot=False)
    d = load_results(lim_dir)
    d["fail"] = ~d["exact"]
    t = d.groupby(["widths", "decoder", "epsilon"]).agg(n_phases=("fail", "size"), failures=("fail", "sum")).reset_index()
    t["fail_rate"] = t.failures / t.n_phases
    t["run"] = lim_dir.relative_to(REPO_ROOT).as_posix()
    t.to_csv(OUT / "p2a_limit_failure.csv", index=False)

    mc_dir = Path(args.mc_run).resolve() if args.mc_run else latest_run("p2b_awqpe_baseline_mc", include_pilot=False)
    m = load_results(mc_dir)
    m["tol_lsb"] = m["error"] <= np.ldexp(1.0, -m["n"].to_numpy()) + 1e-15
    m["decoder_label"] = decoder_label(m)
    m["strata_set"] = np.where(m["stratum"] == "S1_final_half", "S1", "all_but_S1")
    rows = []
    for (w, s, lab, sset), g in m.groupby(["widths", "shots_per_block", "decoder_label", "strata_set"]):
        k, n = int(g.tol_lsb.sum()), len(g)
        lo, hi = wilson(k, n)
        blo, bhi = cluster_bootstrap(g, "tol_lsb")
        rows.append({"widths": w, "shots_per_block": s, "decoder": lab, "strata": sset, "trials": n, "phases": g.phase_id.nunique(),
                     "p_tol_lsb": k / n, "wilson_lo": float(lo), "wilson_hi": float(hi), "cluster_boot_lo": blo, "cluster_boot_hi": bhi,
                     "median_error_lsb": float(np.median(g.error * np.ldexp(1.0, g.n.to_numpy())))})
    out = pd.DataFrame(rows)
    out["run"] = mc_dir.relative_to(REPO_ROOT).as_posix()
    out["status"] = "PRELIMINARY (dev split)"
    out.to_csv(OUT / "p2b_mc_success.csv", index=False)
    print(f"wrote {OUT / 'p2a_limit_failure.csv'} and {OUT / 'p2b_mc_success.csv'}")


if __name__ == "__main__":
    main()
