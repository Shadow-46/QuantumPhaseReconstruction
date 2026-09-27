"""Phase 6 analysis (adaptive overlap). Same inference conventions as P5 (D-018):
phase clusters (partition + phase_id), paired per-phase mean differences,
phase-cluster bootstrap CIs, sign-flip permutation p-values, Holm within a
declared family.

Equal-U-query comparisons use the bracketed matching of D-022 ONLY:
  *_upper_*  closest baseline point with U >= the overlap arm  (PRIMARY for claims
             that overlap improves performance: the baseline gets at least as much work)
  *_lower_*  closest baseline point with U <= the overlap arm  (sensitivity / primary
             for claims that overlap is worse)
Invalid upper brackets (baseline path cap reached; tol is NaN) are EXCLUDED from the
comparison and COUNTED, never substituted. The pilot's overshooting `*Umatched*` arms
are refused outright.

    python -m research.analysis.scripts.summarize_p6 --run DIR --tag pilot|dev|test [--boot 2000] [--perm 20000]

Tables (research/analysis/tables/p6_<tag>_*.csv):
  arms         P(tol), shots, U-queries, overlap actions per variant x decoder x arm
  comparisons  paired differences (equal shots; equal U upper = primary; equal U lower),
               with mean U of both arms, U difference, excluded invalid matches, Holm
  selection    predeclared dev selection per (mechanism, decoder) (D-026)
  rescue       P5-policy (B2) failures: shots only / overlap only / both / neither, with
               raw counts, by decoder x failure type x n x stratum x mechanism
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
PAIR_INDEX = ["cluster", *CELL, "replicate_id"]
B3, B4 = "B3_overlap_alone", "B4_eig_plus_overlap"
# comparison name -> (overlap arm, baseline arm, accounting, role)
COMPARISONS = {
    "B4_vs_B2_equal_shots": (B4, "B2_p5_eig", "equal_shots", "secondary"),
    "B4_vs_B2m_upper_equal_U": (B4, f"B2m_upper_{B4}", "equal_U_upper", "PRIMARY"),
    "B4_vs_B2m_lower_equal_U": (B4, f"B2m_lower_{B4}", "equal_U_lower", "sensitivity"),
    "B3_vs_B1_equal_shots": (B3, "B1_uniform", "equal_shots", "secondary"),
    "B3_vs_B1m_upper_equal_U": (B3, f"B1m_upper_{B3}", "equal_U_upper", "secondary"),
    "B3_vs_B1m_lower_equal_U": (B3, f"B1m_lower_{B3}", "equal_U_lower", "sensitivity"),
}
PRIMARY_COMPARISON = "B4_vs_B2m_upper_equal_U"
LEGACY_MARKER = "Umatched"
# decoders at which each mechanism's overlap information can be used (faithful decoders are controls)
USABLE_DECODERS = {"ext": ("awqpe_ext", "awqpe_eps_safe_ext", "likelihood"), "bridge": ("likelihood",)}


def assert_no_legacy_arms(f: pd.DataFrame) -> None:
    legacy = sorted(a for a in f["arm"].unique() if LEGACY_MARKER in a)
    if legacy:
        raise ValueError(f"refusing to analyse overshooting-U legacy arms {legacy}; use the D-022 bracketed arms")


def paired_comparison(g: pd.DataFrame, x: str, y: str, boot: int, perm: int) -> dict | None:
    """Paired per-phase difference x - y in P(tol) on one (variant, decoder) slice.

    Pairs whose baseline tol is NaN (invalid bracket) are excluded and counted.
    """
    piv = g.pivot_table(index=PAIR_INDEX, columns="arm", values="tol", aggfunc="first", dropna=False)
    if x not in piv or y not in piv:
        return None
    upiv = g.pivot_table(index=PAIR_INDEX, columns="arm", values="u_queries", aggfunc="first", dropna=False)
    spiv = g.pivot_table(index=PAIR_INDEX, columns="arm", values="total_shots", aggfunc="first", dropna=False)
    valid = piv[x].notna() & piv[y].notna()
    dd = (piv[x] - piv[y])[valid].astype(float)
    per = dd.groupby(level="cluster").mean().to_numpy()
    if len(per) == 0:
        return None
    lo, hi = boot_ci(per, boot)
    ux, uy = upiv[x][valid].astype(float), upiv[y][valid].astype(float)
    return {"pairs": int(valid.sum()), "invalid_matches_excluded": int((~valid & piv[x].notna()).sum()), "phases": len(per),
            "diff_p_tol": float(per.mean()), "ci_lo": lo, "ci_hi": hi, "p_signflip": signflip_p(per, perm),
            "p_tol_overlap_arm": float(piv[x][valid].mean()), "p_tol_baseline": float(piv[y][valid].mean()),
            "mean_u_overlap_arm": float(ux.mean()), "mean_u_baseline": float(uy.mean()), "mean_u_diff": float((uy - ux).mean()),
            "mean_shots_overlap_arm": float(spiv[x][valid].astype(float).mean()), "mean_shots_baseline": float(spiv[y][valid].astype(float).mean())}


def select_variants(comp: pd.DataFrame) -> pd.DataFrame:
    """D-026: per (mechanism, decoder) among usable decoders, the variant with the largest PRIMARY
    difference (B4 vs upper-bracket P5 eig at equal U); ties (to 1e-12) broken by lower mean U of B4."""
    prim = comp[comp.comparison == PRIMARY_COMPARISON].copy()
    prim = prim[[d in USABLE_DECODERS.get(m, ()) for m, d in zip(prim.mechanism, prim.decoder)]]
    out = []
    for (mech, dec), g in prim.groupby(["mechanism", "decoder"]):
        g = g.sort_values(["diff_p_tol", "mean_u_overlap_arm"], ascending=[False, True])
        best = g.iloc[0]
        tie = g[np.isclose(g.diff_p_tol, best.diff_p_tol, atol=1e-12)]
        out.append({"mechanism": mech, "decoder": dec, "selected_variant": best.variant, "diff_p_tol": best.diff_p_tol,
                    "ci_lo": best.ci_lo, "ci_hi": best.ci_hi, "p_holm_across_variants": best.p_holm_across_variants,
                    "mean_u_overlap_arm": best.mean_u_overlap_arm, "tie_broken_by_u": len(tie) > 1,
                    "runner_up": g.iloc[1].variant if len(g) > 1 else "", "runner_up_diff": g.iloc[1].diff_p_tol if len(g) > 1 else np.nan})
    return pd.DataFrame(out)


def rescue_table(f: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (variant, mech, dec), g in f.groupby(["variant", "mechanism", "decoder"]):
        piv = g.pivot_table(index=[*PAIR_INDEX, "n", "stratum"], columns="arm", values="tol", aggfunc="first")
        lim = g[g.arm == "B2_p5_eig"].set_index([*PAIR_INDEX, "n", "stratum"])["limit_correct"].reindex(piv.index).astype(bool)
        base_fail = piv["B2_p5_eig"] < 0.5
        for kind, (sa, oa) in {"oracle_capability": (f"R_oracle_shots[{dec}]", f"R_oracle_overlap[{dec}]"),
                               "practical": ("R_eig_shots", "R_trigger_overlap")}.items():
            if sa not in piv or oa not in piv:
                continue
            s_ok, o_ok = piv[sa] > 0.5, piv[oa] > 0.5
            d = pd.DataFrame({"fail": base_fail, "s": s_ok, "o": o_ok, "limited": ~lim}).reset_index()
            d = d[d.fail]
            for keys, gg in d.groupby(["limited", "n", "stratum"]):
                rows.append({"variant": variant, "mechanism": mech, "decoder": dec, "diagnostic": kind,
                             "failure_type": "decoder_limited" if keys[0] else "shot_fixable_class", "n": keys[1], "stratum": keys[2],
                             "p5_failures": len(gg), "shots_only": int((gg.s & ~gg.o).sum()), "overlap_only": int((~gg.s & gg.o).sum()),
                             "both": int((gg.s & gg.o).sum()), "neither": int((~gg.s & ~gg.o).sum())})
    return pd.DataFrame(rows)


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
    assert_no_legacy_arms(f)
    f["cluster"] = f["widths"] + "|" + f["phase_id"]
    OUT.mkdir(parents=True, exist_ok=True)
    meta = {"run": tag, "status": status}

    f.groupby(["variant", "mechanism", "decoder", "arm"]).agg(
        trials=("tol", "size"), p_tol=("tol", "mean"), invalid_rows=("tol", lambda s: int(s.isna().sum())),
        mean_total_shots=("total_shots", "mean"), mean_overlap_shots=("overlap_shots", "mean"),
        mean_u_queries=("u_queries", "mean"), mean_overlap_actions=("n_overlap_actions", "mean")).reset_index().assign(**meta).to_csv(OUT / f"p6_{args.tag}_arms.csv", index=False)

    rows = []
    for (variant, mech, dec), g in f.groupby(["variant", "mechanism", "decoder"]):
        acts = g[g.arm == B4].n_overlap_actions.mean()
        for comp, (x, y, acct, role) in COMPARISONS.items():
            r = paired_comparison(g, x, y, args.boot, args.perm)
            if r:
                rows.append({"variant": variant, "mechanism": mech, "decoder": dec, "comparison": comp, "accounting": acct, "role": role,
                             "overlap_usable_by_decoder": dec in USABLE_DECODERS.get(mech, ()), "mean_overlap_actions_B4": acts, **r})
    comp = pd.DataFrame(rows)
    comp["p_holm_across_variants"] = np.nan
    for (mech, dec, c), g in comp.groupby(["mechanism", "decoder", "comparison"]):
        for i, p in holm(dict(zip(g.index, g.p_signflip))).items():
            comp.loc[i, "p_holm_across_variants"] = p
    comp.assign(**meta).to_csv(OUT / f"p6_{args.tag}_comparisons.csv", index=False)
    select_variants(comp).assign(**meta).to_csv(OUT / f"p6_{args.tag}_selection.csv", index=False)

    rescue_table(f).assign(**meta).to_csv(OUT / f"p6_{args.tag}_rescue.csv", index=False)

    act = []
    gain_cols = [c for c in a.columns if c.endswith("__gain_err_lsb")]
    for (variant, arm), g in a.groupby(["variant", "arm"]):
        for c in gain_cols:
            v = g[c].dropna()
            if len(v):
                act.append({"variant": variant, "arm": arm, "decoder": c.split("__")[0], "actions": len(v),
                            "frac_improved": float((v > 1e-9).mean()), "frac_worsened": float((v < -1e-9).mean()),
                            "mean_gain_err_lsb": float(v.mean()),
                            "boundary_distribution": ";".join(f"{int(k)}:{int(n)}" for k, n in g.boundary.value_counts().sort_index().items())})
    pd.DataFrame(act).assign(**meta).to_csv(OUT / f"p6_{args.tag}_actions.csv", index=False)
    print(f"wrote p6_{args.tag}_*.csv to {OUT}")


if __name__ == "__main__":
    main()
