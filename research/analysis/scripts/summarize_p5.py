"""Phase 5 analysis, implementing docs/DECISIONS.md D-018.

    python -m research.analysis.scripts.summarize_p5 --run DIR --tag dev|test [--boot 2000] [--perm 20000]

Unit of inference: the phase (cluster key = partition + phase_id). Paired
differences vs uniform are averaged per phase over cells and replicates,
then bootstrapped over phases (CI) and sign-flip permuted (p-value). Holm is
applied across the three decoders for the primary comparison eig_cell vs
uniform. The greedy oracle is reported only as a reference.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from research.awqpe.runner.core import REPO_ROOT

OUT = REPO_ROOT / "research" / "analysis" / "tables"
CELL = ["widths", "n", "B", "S0", "r", "dS"]
DECODERS = ["awqpe", "awqpe_eps_safe", "likelihood"]
PRIMARY = "eig_cell"


def holm(p: dict) -> dict:
    keys = sorted(p, key=lambda k: p[k])
    m, running, out = len(keys), 0.0, {}
    for i, k in enumerate(keys):
        running = max(running, min(1.0, (m - i) * p[k]))
        out[k] = running
    return out


def paired_phase_diffs(f: pd.DataFrame, policy: str, ref: str, decoder: str, col: str = "tol") -> pd.Series:
    x = f[(f.decoder == decoder) & (f.policy.isin([policy, ref]))]
    piv = x.pivot_table(index=["cluster", *CELL, "replicate_id"], columns="policy", values=col, aggfunc="first").dropna()
    return (piv[policy].astype(float) - piv[ref].astype(float)).groupby(level="cluster").mean()


def boot_ci(d: np.ndarray, B: int, seed: int = 0) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    m = d[rng.integers(0, len(d), (B, len(d)))].mean(axis=1)
    return float(np.quantile(m, 0.025)), float(np.quantile(m, 0.975))


def signflip_p(d: np.ndarray, N: int, seed: int = 1) -> float:
    rng = np.random.default_rng(seed)
    obs = abs(d.mean())
    signs = rng.choice([-1.0, 1.0], size=(N, len(d)))
    return float((np.abs((signs * d).mean(axis=1)) >= obs - 1e-15).mean())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--perm", type=int, default=20000)
    args = ap.parse_args()
    run = Path(args.run).resolve()
    run_tag = run.relative_to(REPO_ROOT).as_posix()
    status = "HELD-OUT TEST" if args.tag == "test" else "PRELIMINARY / DEV"
    files = sorted((run / "shards").glob("*.parquet"))
    f = pd.concat([pd.read_parquet(x, filters=[("row_type", "==", "final")]) for x in files], ignore_index=True)
    f = f.dropna(axis=1, how="all")
    f["cluster"] = f["widths"] + "|" + f["phase_id"]
    f["err_lsb"] = f.error * np.ldexp(1.0, f.n.to_numpy().astype(int))
    OUT.mkdir(parents=True, exist_ok=True)

    # 1. final summary per cell
    shot_cols = [c for c in f.columns if c.startswith("shots_b")]
    summ = f.groupby([*CELL, "decoder", "policy"]).agg(
        trials=("tol", "size"), p_tol=("tol", "mean"), rmse_lsb=("err_lsb", lambda e: float(np.sqrt(np.mean(e**2)))),
        median_err_lsb=("err_lsb", "median"), p_exact=("exact", "mean"), u_queries=("u_queries", "mean"),
        total_shots=("total_shots", "mean"), **{f"mean_{c}": (c, "mean") for c in shot_cols}).reset_index()
    summ.assign(run=run_tag, status=status).to_csv(OUT / f"p5_{args.tag}_final_summary.csv", index=False)

    # 2. primary + secondary paired comparisons vs uniform
    rows, pvals = [], {}
    policies = sorted(p for p in f.policy.unique() if p != "uniform" and not p.startswith("greedy_oracle"))
    for dec in DECODERS:
        oracle = f"greedy_oracle[{dec}]"
        d_or = paired_phase_diffs(f, oracle, "uniform", dec) if (f.policy == oracle).any() else None
        for pol in policies:
            arr = paired_phase_diffs(f, pol, "uniform", dec).to_numpy()
            lo, hi = boot_ci(arr, args.boot)
            p = signflip_p(arr, args.perm)
            or_mean = float(d_or.mean()) if d_or is not None else np.nan
            rows.append({"decoder": dec, "policy": pol, "phases": len(arr), "diff_p_tol": float(arr.mean()), "ci_lo": lo, "ci_hi": hi,
                         "p_signflip": p, "tost_equivalent_1pt": bool(lo > -0.01 and hi < 0.01),
                         "oracle_minus_uniform": or_mean, "efficiency_vs_oracle": float(arr.mean() / or_mean) if abs(or_mean) > 1e-12 else np.nan,
                         "family": "primary" if pol == PRIMARY else "secondary"})
            if pol == PRIMARY:
                pvals[dec] = p
    adj = holm(pvals)
    prim = pd.DataFrame(rows)
    prim["p_holm_primary"] = [adj.get(r.decoder) if r.policy == PRIMARY else np.nan for r in prim.itertuples()]
    prim.assign(run=run_tag, status=status).to_csv(OUT / f"p5_{args.tag}_primary.csv", index=False)

    # 3. regimes: mean paired difference vs uniform by factor
    reg = []
    x = f[~f.policy.str.startswith("greedy_oracle")]
    piv = x.pivot_table(index=["cluster", "stratum", *CELL, "replicate_id", "decoder", "limit_correct"], columns="policy", values="tol", aggfunc="first")
    for pol in policies:
        diff = (piv[pol].astype(float) - piv["uniform"].astype(float)).rename("diff").reset_index()
        for factor in ["stratum", "n", "widths", "S0", "r", "dS", "limit_correct"]:
            g = diff.groupby(["decoder", factor])["diff"].agg(["mean", "size"]).reset_index().rename(columns={factor: "level"})
            g["level"] = g["level"].astype(str)
            g["factor"], g["policy"] = factor, pol
            reg.append(g)
    pd.concat(reg).assign(run=run_tag, status=status).to_csv(OUT / f"p5_{args.tag}_regimes.csv", index=False)

    # 4. decision analysis on step rows
    dec_rows = []
    key = ["widths", "S0", "r", "dS", "trial"]
    step_cols = ["policy", "B", "widths", "S0", "r", "dS", "trial"] + [f"signal_b{b}" for b in range(1, 5)] + [
        f"{dec}__{k}" for dec in DECODERS for k in ("gain_tol_chosen", "gain_err_chosen", "gain_tol_best", "gain_err_best", "chose_best", "changed", "gain_err_b1", "gain_err_b2", "gain_err_b3", "gain_err_b4")]
    import pyarrow.parquet as pq

    schemas = {x: set(pq.read_schema(x).names) for x in files}
    for pol in sorted(f.policy.unique()):
        parts = [pd.read_parquet(x, columns=[c for c in step_cols if c in schemas[x]], filters=[("row_type", "==", "step"), ("policy", "==", pol)]) for x in files]
        g = pd.concat(parts, ignore_index=True)
        if g.empty:
            continue
        for dec in DECODERS:
            # signal-vs-realized rank correlation, centred within each decision, per block count B (no NaN padding)
            sc_all, gc_all = [], []
            for Bv, gb in g.groupby("B"):
                Bv = int(Bv)
                sig = gb[[f"signal_b{b}" for b in range(1, Bv + 1)]].to_numpy(float)
                gains = gb[[f"{dec}__gain_err_b{b}" for b in range(1, Bv + 1)]].to_numpy(float)
                sc_all.append((sig - sig.mean(axis=1, keepdims=True)).ravel())
                gc_all.append((gains - gains.mean(axis=1, keepdims=True)).ravel())
            sc, gc = np.concatenate(sc_all), np.concatenate(gc_all)
            rho = spearmanr(sc, gc).statistic if np.ptp(sc) > 0 else np.nan
            # a decision is informative if some block's next batch would change tolerance success or error
            informative = (g[f"{dec}__gain_tol_best"] > 0) | (g[f"{dec}__gain_err_best"] > 1e-9)
            gsum = g.groupby(key)[f"{dec}__gain_err_chosen"].sum()
            budget = g.groupby(key).size() * g.groupby(key).dS.first()
            dec_rows.append({"policy": pol, "decoder": dec, "decisions": int(len(g)),
                             "p_chose_oracle_best": float(g[f"{dec}__chose_best"].mean()),
                             "informative_decisions": int(informative.sum()),
                             "p_chose_best_when_informative": float(g.loc[informative, f"{dec}__chose_best"].mean()) if informative.any() else np.nan,
                             "mean_realized_gain_err_chosen": float(g[f"{dec}__gain_err_chosen"].mean()),
                             "mean_regret_tol": float((g[f"{dec}__gain_tol_best"] - g[f"{dec}__gain_tol_chosen"]).mean()),
                             "mean_regret_err_lsb": float((g[f"{dec}__gain_err_best"] - g[f"{dec}__gain_err_chosen"]).mean()),
                             "spearman_signal_vs_realized_gain": rho,
                             "err_reduction_lsb_per_extra_shot": float((gsum / budget).mean()),
                             "p_decoded_result_changed": float(g[f"{dec}__changed"].mean())})
    pd.DataFrame(dec_rows).assign(run=run_tag, status=status).to_csv(OUT / f"p5_{args.tag}_decisions.csv", index=False)
    print(f"wrote p5_{args.tag}_*.csv to {OUT}")


if __name__ == "__main__":
    main()
