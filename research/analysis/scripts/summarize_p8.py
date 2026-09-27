"""Phase 8 (D-028) analysis: coverage, resources, dev calibration of S*, paired comparisons.

    python -m research.analysis.scripts.summarize_p8 --run <run_dir> --tag dev
    python -m research.analysis.scripts.summarize_p8 --run <run_dir> --tag test --sstar research/configs/p8_frozen_sstar.csv

Dev mode writes the S* table (smallest fixed S per block reaching pooled coverage >= 1 - alpha).
Test mode must be given the frozen dev S* table and never re-calibrates.
Savings are fixed minus stop (positive = stopping uses fewer resources); per-trial
differences are averaged within phase clusters before the bootstrap and sign-flip test.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from research.analysis.scripts.summarize_p5 import boot_ci, holm, signflip_p
from research.awqpe.runner.core import REPO_ROOT

OUT = REPO_ROOT / "research" / "analysis" / "tables"
KEY = ["widths", "phase_id", "replicate_id", "tau_exp"]
PRIMARY_WIDTHS = {"4-4": 8, "4-4-4": 12, "4-4-4-4": 16}


def primary_cells():
    return {(w, j) for w, n in PRIMARY_WIDTHS.items() for j in (3, n // 2, n)}


def load(run: Path) -> pd.DataFrame:
    cols = [*KEY, "stratum", "n", "arm", "decoder", "fixed_S", "truncated", "stopped", "shots", "u_queries", "covered", "credible_mass", "alpha"]
    f = pd.concat([pd.read_parquet(x, columns=cols) for x in sorted((run / "shards").glob("*.parquet"))], ignore_index=True)
    f["cluster"] = f.widths + "|" + f.phase_id
    return f


def cluster_mean_ci(g: pd.DataFrame, col: str, B: int) -> tuple[float, float, float]:
    v = g.groupby("cluster")[col].mean().to_numpy(dtype=float)
    lo, hi = boot_ci(v, B)
    return float(v.mean()), lo, hi


def calibrate(f: pd.DataFrame, alpha: float) -> pd.DataFrame:
    fx = f[f.arm.str.startswith("fixed_S")]
    cov = fx.groupby(["widths", "tau_exp", "truncated", "decoder", "fixed_S"]).covered.mean().reset_index()
    rows = []
    for key, g in cov.groupby(["widths", "tau_exp", "truncated", "decoder"]):
        ok = g[g.covered >= 1 - alpha].sort_values("fixed_S")
        rows.append(dict(zip(["widths", "tau_exp", "truncated", "decoder"], key),
                         S_star=float(ok.fixed_S.iloc[0]) if len(ok) else np.nan,
                         dev_coverage=float(ok.covered.iloc[0]) if len(ok) else float(g.covered.max())))
    return pd.DataFrame(rows)


def paired(f, w, j, arm_a, dec_a, arm_b, dec_b, B, N):
    """Resources of arm_b minus arm_a (positive = arm_a cheaper), phase-clustered."""
    g = f[(f.widths == w) & (f.tau_exp == j)]
    a = g[(g.arm == arm_a) & (g.decoder == dec_a)].set_index(KEY)
    b = g[(g.arm == arm_b) & (g.decoder == dec_b)].set_index(KEY)
    if a.empty or b.empty:
        return None
    b = b.loc[a.index]
    d = pd.DataFrame({"cluster": a.cluster, "d_shots": b.shots - a.shots, "d_u": b.u_queries - a.u_queries,
                      "d_cov": a.covered.astype(float) - b.covered.astype(float)})
    c = d.groupby("cluster").mean()
    out = {"widths": w, "tau_exp": j, "arm": arm_a, "decoder": dec_a, "ref_arm": arm_b, "ref_decoder": dec_b,
           "clusters": len(c), "trials": len(d),
           "arm_shots": float(a.shots.mean()), "ref_shots": float(b.shots.mean()),
           "arm_u": float(a.u_queries.mean()), "ref_u": float(b.u_queries.mean()),
           "arm_coverage": float(a.covered.mean()), "ref_coverage": float(b.covered.mean())}
    for k in ("d_shots", "d_u", "d_cov"):
        v = c[k].to_numpy()
        lo, hi = boot_ci(v, B)
        out |= {k: float(v.mean()), f"{k}_lo": lo, f"{k}_hi": hi}
    out["p_shots"] = signflip_p(c.d_shots.to_numpy(), N)
    ca = a.groupby("cluster").covered.mean().to_numpy(dtype=float)
    out["arm_coverage_lo"], out["arm_coverage_hi"] = boot_ci(ca, B)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--sstar", default=None)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--perm", type=int, default=20000)
    args = ap.parse_args()
    run = Path(args.run).resolve()
    status = "HELD-OUT TEST" if args.tag == "test" else "PRELIMINARY / DEV"
    f = load(run)
    alpha = float(f.alpha.iloc[0])
    OUT.mkdir(parents=True, exist_ok=True)
    tag = f"p8_{args.tag}"

    # 1. per-cell coverage and resources
    rows = []
    for key, g in f.groupby(["widths", "tau_exp", "arm", "decoder"]):
        m, lo, hi = cluster_mean_ci(g.assign(cov=g.covered.astype(float)), "cov", args.boot)
        rows.append(dict(zip(["widths", "tau_exp", "arm", "decoder"], key), trials=len(g), coverage=m, coverage_lo=lo, coverage_hi=hi,
                         shots=g.shots.mean(), shots_median=g.shots.median(), u_queries=g.u_queries.mean(),
                         stopped=g.stopped.mean(), mean_credible_mass=g.credible_mass.mean()))
    pd.DataFrame(rows).assign(status=status).to_csv(OUT / f"{tag}_cells.csv", index=False)

    # 2. coverage by stratum for the stopping arms (under-coverage audit)
    st = f[f.arm.str.startswith("stop_")].groupby(["widths", "tau_exp", "arm", "stratum"]).agg(
        coverage=("covered", "mean"), shots=("shots", "mean"), stopped=("stopped", "mean")).reset_index()
    st.assign(status=status).to_csv(OUT / f"{tag}_strata.csv", index=False)

    # 3. calibration: credible mass (at stop) vs realised coverage, stopping arms
    s = f[f.arm.str.startswith("stop_") & f.stopped].copy()
    s["mass_bin"] = pd.cut(s.credible_mass, [0.95, 0.99, 0.999, 1.0 + 1e-12], include_lowest=True)
    cal = s.groupby(["arm", "mass_bin"], observed=True).agg(trials=("covered", "size"), mean_mass=("credible_mass", "mean"),
                                                            coverage=("covered", "mean")).reset_index()
    cal["mass_bin"] = cal.mass_bin.astype(str)
    cal.assign(status=status).to_csv(OUT / f"{tag}_calibration.csv", index=False)

    # 4. S*: calibrate on dev; test must use the frozen dev table
    if args.tag == "test":
        if not args.sstar:
            raise SystemExit("test analysis needs --sstar (frozen dev calibration); refusing to re-calibrate on test")
        sstar = pd.read_csv(args.sstar)
    else:
        sstar = calibrate(f, alpha)
        sstar.assign(status=status).to_csv(OUT / f"{tag}_sstar.csv", index=False)
    sstar = sstar.set_index(["widths", "tau_exp", "truncated", "decoder"]).S_star

    # 5. comparisons
    comps, pvals = [], {}
    prim = primary_cells()
    for (w, j), g in f.groupby(["widths", "tau_exp"]):
        has_trunc = bool(g.truncated.any())
        sfx = "_trunc" if has_trunc else ""
        stop = f"stop_eig{sfx}"
        todo = []
        S = sstar.get((w, j, has_trunc, "likelihood"), np.nan)
        if np.isfinite(S):
            todo.append(("primary" if (w, j) in prim else "secondary_vs_fixed", stop, "likelihood", f"fixed_S{int(S)}{sfx}", "likelihood"))
        S1 = sstar.get((w, j, has_trunc, "awqpe_eps_safe"), np.nan)
        if np.isfinite(S1):
            todo.append(("secondary_vs_awqpe_fixed", stop, "likelihood", f"fixed_S{int(S1)}{sfx}", "awqpe_eps_safe"))
        todo.append(("secondary_allocation", stop, "likelihood", f"stop_uniform{sfx}", "likelihood"))
        if has_trunc:
            todo.append(("secondary_truncation", "stop_eig_trunc", "likelihood", "stop_eig", "likelihood"))
        for fam, a, da, b, db in todo:
            r = paired(f, w, j, a, da, b, db, args.boot, args.perm)
            if r is None:
                continue
            r["family"] = fam
            comps.append(r)
            if fam == "primary":
                pvals[(w, j)] = r["p_shots"]
    comp = pd.DataFrame(comps)
    adj = holm(pvals)
    comp["p_holm"] = [adj.get((w, j)) if fam == "primary" else np.nan for w, j, fam in zip(comp.widths, comp.tau_exp, comp.family)]
    comp["positive"] = (comp.family == "primary") & (comp.p_holm < 0.05) & (comp.d_shots > 0) & (comp.arm_coverage >= 1 - alpha)
    missing = sorted(prim - set(pvals))
    comp.assign(status=status, missing_primary=str(missing)).to_csv(OUT / f"{tag}_comparisons.csv", index=False)
    print(f"wrote {tag}_*.csv; primary cells {len(pvals)}/{len(prim)}, missing baseline: {missing}")


if __name__ == "__main__":
    main()
