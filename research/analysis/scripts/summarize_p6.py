"""Phase 6 analysis (adaptive overlap). Same inference conventions as P5 (D-018):
phase clusters (partition + phase_id), paired per-phase mean differences,
phase-cluster bootstrap CIs, sign-flip permutation p-values, Holm within a
declared family.

    python -m research.analysis.scripts.summarize_p6 --run DIR --tag pilot|dev|test [--boot 2000] [--perm 20000]

Tables (research/analysis/tables/p6_<tag>_*.csv):
  arms         P(tol), shots, U-queries, overlap actions per variant x decoder x arm
  comparisons  paired differences: B4-B2 (equal shots), B4-B2u (equal U-queries),
               B3-B1 (equal shots), B3-B1u (equal U-queries); Holm across variants
               within each (decoder, comparison) family
  rescue       P5-policy (B2) failures classified by what rescues them with A extra
               batches: shots only / overlap only / both / neither (greedy-oracle
               capability and practical policies), split by decoder-limited status
  actions      overlap actions: count, boundary location, immediate improvement rate
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from research.analysis.scripts.summarize_p5 import boot_ci, holm, signflip_p
from research.awqpe.runner.core import REPO_ROOT

OUT = REPO_ROOT / "research" / "analysis" / "tables"
CELL = ["widths", "S0", "r", "dS"]
COMPARISONS = {"B4_vs_B2_equal_shots": ("B4_eig_plus_overlap", "B2_p5_eig"),
               "B4_vs_B2u_equal_U": ("B4_eig_plus_overlap", "B2u_p5_eig_Umatched_B4"),
               "B3_vs_B1_equal_shots": ("B3_overlap_alone", "B1_uniform"),
               "B3_vs_B1u_equal_U": ("B3_overlap_alone", "B1u_uniform_Umatched_B3")}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--perm", type=int, default=20000)
    args = ap.parse_args()
    run = Path(args.run).resolve()
    tag = run.relative_to(REPO_ROOT).as_posix()
    status = {"test": "HELD-OUT TEST", "dev": "PRELIMINARY / DEV"}.get(args.tag, "PILOT / DEV")
    files = sorted((run / "shards").glob("*.parquet"))
    f = pd.concat([pd.read_parquet(x, filters=[("row_type", "==", "final")]) for x in files], ignore_index=True).dropna(axis=1, how="all")
    a = pd.concat([pd.read_parquet(x, filters=[("row_type", "==", "action")]) for x in files], ignore_index=True).dropna(axis=1, how="all")
    f["cluster"] = f["widths"] + "|" + f["phase_id"]
    OUT.mkdir(parents=True, exist_ok=True)

    arms = f.groupby(["variant", "decoder", "arm"]).agg(trials=("tol", "size"), p_tol=("tol", "mean"), mean_total_shots=("total_shots", "mean"),
                                                        mean_overlap_shots=("overlap_shots", "mean"), mean_u_queries=("u_queries", "mean"),
                                                        mean_overlap_actions=("n_overlap_actions", "mean")).reset_index()
    arms.assign(run=tag, status=status).to_csv(OUT / f"p6_{args.tag}_arms.csv", index=False)

    rows = []
    for (variant, dec), g in f.groupby(["variant", "decoder"]):
        piv = g.pivot_table(index=["cluster", *CELL, "replicate_id"], columns="arm", values="tol", aggfunc="first")
        upiv = g.pivot_table(index=["cluster", *CELL, "replicate_id"], columns="arm", values="u_queries", aggfunc="first")
        for comp, (x, y) in COMPARISONS.items():
            if x not in piv or y not in piv:
                continue
            dd = (piv[x].astype(float) - piv[y].astype(float)).dropna()
            per = dd.groupby(level="cluster").mean().to_numpy()
            lo, hi = boot_ci(per, args.boot)
            rows.append({"variant": variant, "decoder": dec, "comparison": comp, "phases": len(per), "diff_p_tol": float(per.mean()),
                         "ci_lo": lo, "ci_hi": hi, "p_signflip": signflip_p(per, args.perm),
                         "u_ratio_x_over_y": float(upiv[x].mean() / upiv[y].mean())})
    comp = pd.DataFrame(rows)
    comp["p_holm_across_variants"] = np.nan
    for (dec, c), g in comp.groupby(["decoder", "comparison"]):
        for i, p in holm(dict(zip(g.index, g.p_signflip))).items():
            comp.loc[i, "p_holm_across_variants"] = p
    comp.assign(run=tag, status=status).to_csv(OUT / f"p6_{args.tag}_comparisons.csv", index=False)

    res = []
    for (variant, dec), g in f.groupby(["variant", "decoder"]):
        idx = ["cluster", *CELL, "replicate_id"]
        piv = g.pivot_table(index=idx, columns="arm", values="tol", aggfunc="first")
        lim = g[g.arm == "B2_p5_eig"].set_index(idx)["limit_correct"].reindex(piv.index).astype(bool)
        base_fail = ~piv["B2_p5_eig"].astype(bool)
        for kind, (sa, oa) in {"oracle_capability": (f"R_oracle_shots[{dec}]", f"R_oracle_overlap[{dec}]"),
                               "practical": ("R_eig_shots", "R_trigger_overlap")}.items():
            if sa not in piv or oa not in piv:
                continue
            s_ok, o_ok = piv[sa].astype(bool), piv[oa].astype(bool)
            for lim_state in (True, False):
                m = base_fail & (lim == lim_state)
                nfail = int(m.sum())
                if nfail == 0:
                    continue
                res.append({"variant": variant, "decoder": dec, "diagnostic": kind, "decoder_limited": not lim_state,
                            "p5_failures": nfail, "shots_only": int((m & s_ok & ~o_ok).sum()), "overlap_only": int((m & ~s_ok & o_ok).sum()),
                            "both": int((m & s_ok & o_ok).sum()), "neither": int((m & ~s_ok & ~o_ok).sum())})
    rs = pd.DataFrame(res)
    if len(rs):
        for c in ("shots_only", "overlap_only", "both", "neither"):
            rs[f"frac_{c}"] = rs[c] / rs["p5_failures"]
    rs.assign(run=tag, status=status).to_csv(OUT / f"p6_{args.tag}_rescue.csv", index=False)

    act = []
    gain_cols = [c for c in a.columns if c.endswith("__gain_err_lsb")]
    for (variant, arm), g in a.groupby(["variant", "arm"]):
        for c in gain_cols:
            v = g[c].dropna()
            if not len(v):
                continue
            act.append({"variant": variant, "arm": arm, "decoder": c.split("__")[0], "actions": len(v),
                        "frac_improved": float((v > 1e-9).mean()), "frac_worsened": float((v < -1e-9).mean()),
                        "mean_gain_err_lsb": float(v.mean()),
                        "boundary_distribution": ";".join(f"{int(k)}:{int(n)}" for k, n in g.boundary.value_counts().sort_index().items())})
    pd.DataFrame(act).assign(run=tag, status=status).to_csv(OUT / f"p6_{args.tag}_actions.csv", index=False)
    print(f"wrote p6_{args.tag}_*.csv to {OUT}")


if __name__ == "__main__":
    main()
