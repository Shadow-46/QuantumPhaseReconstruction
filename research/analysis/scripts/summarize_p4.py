"""Phase 4 analysis, implementing the plan declared in docs/DECISIONS.md D-016. PRELIMINARY (dev phases).

    python -m research.analysis.scripts.summarize_p4 [--run DIR] [--boot 300]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from research.awqpe.runner.core import REPO_ROOT, latest_run, load_results

OUT = REPO_ROOT / "research" / "analysis" / "tables"
# D-016 orientation: +1 means "higher value = more need for shots".
ORIENT = {"ratio": 1, "entropy": 1, "local_post_sd": 1, "chunk_risk": 1, "boundary_risk": 1, "eig_phi": 1, "eig_cell": 1,
          "fisher_phi": 1, "c1": -1, "margin": -1, "p_chunk_correct": -1, "legacy_cw": -1}
CELL = ["widths", "n", "B", "S0", "dS", "decoder"]


def auroc(score: np.ndarray, label: np.ndarray) -> float:
    label = label.astype(bool)
    npos, nneg = label.sum(), (~label).sum()
    if npos == 0 or nneg == 0:
        return np.nan
    r = rankdata(score)
    return float((r[label].sum() - npos * (npos + 1) / 2) / (npos * nneg))


def greedy_values(g: pd.DataFrame, sig: str) -> pd.DataFrame:
    """Per trial: P(tol1) for the argmax-signal block (ties averaged), a random block, and the oracle block."""
    x = g[["trial", "phase_id", "tol0", "tol1"]].copy()
    x["s"] = ORIENT[sig] * g[sig].to_numpy()
    mx = x.groupby("trial")["s"].transform("max")
    x["chosen"] = np.isclose(x["s"], mx, rtol=0, atol=1e-12)
    per = x.groupby("trial").agg(phase_id=("phase_id", "first"), tol0=("tol0", "first"), random=("tol1", "mean"), oracle=("tol1", "max"))
    per["greedy"] = x[x.chosen].groupby("trial")["tol1"].mean()
    return per


def cluster_boot_diff(per: pd.DataFrame, a: str, b: str, B: int, seed: int = 0) -> tuple[float, float]:
    ph = per.groupby("phase_id")[[a, b]].mean()
    diff = (ph[a] - ph[b]).to_numpy()
    rng = np.random.default_rng(seed)
    m = diff[rng.integers(0, len(diff), (B, len(diff)))].mean(axis=1)
    return float(np.quantile(m, 0.025)), float(np.quantile(m, 0.975))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", default=None)
    ap.add_argument("--boot", type=int, default=300)
    args = ap.parse_args()
    run = Path(args.run).resolve() if args.run else latest_run("p4_information_study", include_pilot=False)
    tag = run.relative_to(REPO_ROOT).as_posix()
    d = load_results(run)
    d["trial"] = d["phase_id"] + "#" + d["replicate_id"].astype(str)
    d["d_err"] = d["err0"] - d["err1"]
    d["fix"] = (~d["tol0"]) & d["tol1"]
    OUT.mkdir(parents=True, exist_ok=True)
    signals = [s for s in ORIENT if s in d.columns]

    rows, strat_rows, block_rows, head_rows = [], [], [], []
    for key, g in d.groupby(CELL):
        base = dict(zip(CELL, key))
        failing = g[~g.tol0]
        head = greedy_values(g, signals[0])
        fail_trials = ~head.tol0.astype(bool)
        head_rows.append({**base, "trials": len(head), "p_tol0": head.tol0.mean(), "p_random": head.random.mean(), "p_oracle": head.oracle.mean(),
                          "frac_failing_fixable_by_one_batch": float((head.oracle[fail_trials] > 0).mean()) if fail_trials.any() else np.nan})
        for sig in signals:
            s = ORIENT[sig] * g[sig].to_numpy()
            rho = spearmanr(s, g.d_err.to_numpy()).statistic if np.ptp(s) > 0 else np.nan
            per = greedy_values(g, sig)
            lo, hi = cluster_boot_diff(per, "greedy", "random", args.boot)
            headroom = per.oracle.mean() - per.random.mean()
            rows.append({**base, "signal": sig, "spearman_d_err": rho,
                         "auroc_fix_all": auroc(s, g.fix.to_numpy()),
                         "auroc_fix_within_failing": auroc(ORIENT[sig] * failing[sig].to_numpy(), failing.fix.to_numpy()) if len(failing) else np.nan,
                         "p_greedy": per.greedy.mean(), "p_random": per.random.mean(), "p_oracle": per.oracle.mean(),
                         "gain_vs_random": per.greedy.mean() - per.random.mean(), "gain_ci_lo": lo, "gain_ci_hi": hi,
                         "headroom_captured": (per.greedy.mean() - per.random.mean()) / headroom if headroom > 1e-12 else np.nan})
            for st, gs in g.groupby("stratum"):
                ps = greedy_values(gs, sig)
                fs = gs[~gs.tol0]
                strat_rows.append({**base, "signal": sig, "stratum": st, "trials": len(ps),
                                   "auroc_fix_within_failing": auroc(ORIENT[sig] * fs[sig].to_numpy(), fs.fix.to_numpy()) if len(fs) else np.nan,
                                   "gain_vs_random": ps.greedy.mean() - ps.random.mean(), "oracle_minus_random": ps.oracle.mean() - ps.random.mean()})
        for (blk, off, wid), gb in g.groupby(["block", "offset", "width"]):
            block_rows.append({**base, "block": blk, "offset": off, "width": wid, "p_fix": gb.fix.mean(),
                               "p_break": float((gb.tol0 & ~gb.tol1).mean()), "mean_d_err_lsb": float(gb.d_err.mean() * 2.0 ** gb.n.iloc[0])})
    pd.DataFrame(rows).assign(run=tag, status="PRELIMINARY").to_csv(OUT / "p4_signal_metrics.csv", index=False)
    pd.DataFrame(strat_rows).assign(run=tag, status="PRELIMINARY").to_csv(OUT / "p4_signal_by_stratum.csv", index=False)
    pd.DataFrame(block_rows).assign(run=tag, status="PRELIMINARY").to_csv(OUT / "p4_block_value.csv", index=False)
    pd.DataFrame(head_rows).assign(run=tag, status="PRELIMINARY").to_csv(OUT / "p4_headroom.csv", index=False)

    one = d[d.decoder == "likelihood"]  # signals do not depend on the decoder; avoid triple counting
    cal = []
    for (w, S0), g in one.groupby(["widths", "S0"]):
        g = g.drop_duplicates(["trial", "block"])
        bins = np.clip((g.p_chunk_correct * 10).astype(int), 0, 9)
        for bidx, gb in g.groupby(bins):
            cal.append({"widths": w, "S0": S0, "bin": int(bidx), "mean_predicted": gb.p_chunk_correct.mean(), "observed_t1_correct": gb.t1_correct.mean(), "count": len(gb)})
    pd.DataFrame(cal).assign(run=tag, status="PRELIMINARY").to_csv(OUT / "p4_calibration.csv", index=False)
    print(f"wrote p4_*.csv to {OUT}")


if __name__ == "__main__":
    main()
