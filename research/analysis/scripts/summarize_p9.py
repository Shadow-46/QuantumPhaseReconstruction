"""Phase 9 (D-030) analysis: decoder robustness and stopping coverage under analytic noise.

    python -m research.analysis.scripts.summarize_p9 --run <run_dir> --tag dev|test
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from research.analysis.scripts.summarize_p5 import boot_ci, holm, signflip_p
from research.awqpe.runner.core import REPO_ROOT

OUT = REPO_ROOT / "research" / "analysis" / "tables"
KEY = ["widths", "noise_id", "phase_id", "replicate_id"]
PRIMARY_NOISE = ("ro03", "dp5", "jt02", "comb")
PRIMARY_ARM = "fixed_S16"  # D-032: S = 64 is at ceiling on dev
FAMILIES = {"F1_aware_vs_awqpe_safe": ("likelihood_aware", "awqpe_eps_safe"),
            "F2_aware_vs_ideal_model": ("likelihood_aware", "likelihood_ideal")}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--perm", type=int, default=20000)
    args = ap.parse_args()
    run = Path(args.run).resolve()
    status = "HELD-OUT TEST" if args.tag == "test" else "PRELIMINARY / DEV"
    cols = [*KEY, "stratum", "channel", "arm", "decoder", "tol", "tol_half", "covered", "stopped", "shots", "u_queries", "tau"]
    f = pd.concat([pd.read_parquet(x, columns=cols) for x in sorted((run / "shards").glob("*.parquet"))], ignore_index=True)
    f["cluster"] = f.widths + "|" + f.phase_id
    tag = f"p9_{args.tag}"
    OUT.mkdir(parents=True, exist_ok=True)

    cells = f.groupby(["widths", "noise_id", "channel", "arm", "decoder"]).agg(
        trials=("tol", "size"), p_tol=("tol", "mean"), p_tol_half=("tol_half", "mean"), coverage=("covered", "mean"),
        shots=("shots", "mean"), u_queries=("u_queries", "mean"), stopped=("stopped", "mean")).reset_index()
    cells.assign(status=status).to_csv(OUT / f"{tag}_cells.csv", index=False)

    # coverage audit of the stopping arms, with phase-cluster CI
    rows = []
    for key, g in f[f.arm.str.startswith("stop_")].groupby(["widths", "noise_id", "channel", "arm", "decoder"]):
        v = g.groupby("cluster").covered.mean().to_numpy(dtype=float)
        lo, hi = boot_ci(v, args.boot)
        rows.append(dict(zip(["widths", "noise_id", "channel", "arm", "decoder"], key), coverage=float(g.covered.mean()),
                         coverage_lo=lo, coverage_hi=hi, shots=float(g.shots.mean()), stopped=float(g.stopped.mean()),
                         meets_095=bool(g.covered.mean() >= 0.95)))
    pd.DataFrame(rows).assign(status=status).to_csv(OUT / f"{tag}_coverage.csv", index=False)

    # primary family: paired success differences at fixed S = 64
    # F3 (D-032): stopping coverage at tau = 2^-n, noise-aware minus ideal-likelihood posterior
    n_of = {w: sum(map(int, w.split("-"))) for w in f.widths.unique()}
    fx = pd.concat([f[f.arm == PRIMARY_ARM].assign(family_arm="fixed", y=lambda x: x.tol),
                    f[f.arm == f.widths.map(lambda w: f"stop_tau{n_of[w]}")].assign(family_arm="stop", y=lambda x: x.covered)])
    families = {**{k: (a, b, "fixed") for k, (a, b) in FAMILIES.items()},
                "F3_stop_coverage_aware_vs_ideal": ("likelihood_aware", "likelihood_ideal", "stop")}
    comps, pvals = [], {}
    for fam, (a, b, kind) in families.items():
        sub = fx[(fx.family_arm == kind) & fx.noise_id.isin(PRIMARY_NOISE)]
        for (w, nid), g in sub.groupby(["widths", "noise_id"]):
            A = g[g.decoder == a].set_index(KEY)
            Bf = g[g.decoder == b].set_index(KEY).loc[A.index]
            d = pd.DataFrame({"cluster": A.cluster, "d": A.y.astype(float) - Bf.y.astype(float)}).groupby("cluster").d.mean().to_numpy()
            lo, hi = boot_ci(d, args.boot)
            p = signflip_p(d, args.perm)
            pvals[(fam, w, nid)] = p
            comps.append({"family": fam, "widths": w, "noise_id": nid, "decoder": a, "ref_decoder": b, "clusters": len(d),
                          "p_tol": float(A.y.mean()), "ref_p_tol": float(Bf.y.mean()),
                          "diff_pts": 100 * float(d.mean()), "lo_pts": 100 * lo, "hi_pts": 100 * hi, "p": p})
    comp = pd.DataFrame(comps)
    adj = holm(pvals)
    comp["p_holm"] = [adj[(r.family, r.widths, r.noise_id)] for r in comp.itertuples()]
    comp.assign(status=status).to_csv(OUT / f"{tag}_primary.csv", index=False)
    print(f"wrote {tag}_*.csv ({len(comp)} primary comparisons)")


if __name__ == "__main__":
    main()
