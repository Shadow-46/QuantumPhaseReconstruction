"""Phase 7 (D-031) analysis: unified controller vs its ablations (shots to stop, U-queries, coverage).

    python -m research.analysis.scripts.summarize_p7 --run <run_dir> --tag dev|test
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from research.analysis.scripts.summarize_p5 import boot_ci, holm, signflip_p
from research.awqpe.runner.core import REPO_ROOT

OUT = REPO_ROOT / "research" / "analysis" / "tables"
KEY = ["widths", "phase_id", "replicate_id", "tau_exp"]
PRIMARY_WIDTHS = ("4-4", "4-4-4", "4-4-4-4")
COMPARISONS = [("unified_bridge", "stop_eig", "primary"), ("unified_ext", "stop_eig", "primary"),
               ("unified_bridge", "stop_uniform", "secondary"), ("unified_ext", "stop_uniform", "secondary"),
               ("stop_eig", "stop_uniform", "secondary")]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--perm", type=int, default=20000)
    args = ap.parse_args()
    run = Path(args.run).resolve()
    status = "HELD-OUT TEST" if args.tag == "test" else "PRELIMINARY / DEV"
    cols = [*KEY, "stratum", "arm", "shots", "u_queries", "covered", "stopped", "overlap_batches"]
    f = pd.concat([pd.read_parquet(x, columns=cols) for x in sorted((run / "shards").glob("*.parquet"))], ignore_index=True)
    f["cluster"] = f.widths + "|" + f.phase_id
    tag = f"p7_{args.tag}"
    OUT.mkdir(parents=True, exist_ok=True)

    f.groupby(["widths", "tau_exp", "arm"]).agg(
        trials=("shots", "size"), coverage=("covered", "mean"), shots=("shots", "mean"), u_queries=("u_queries", "mean"),
        overlap_batches=("overlap_batches", "mean"), used_overlap=("overlap_batches", lambda x: (x > 0).mean()),
        stopped=("stopped", "mean")).reset_index().assign(status=status).to_csv(OUT / f"{tag}_cells.csv", index=False)
    f.groupby(["arm", "stratum"]).covered.mean().reset_index().assign(status=status).to_csv(OUT / f"{tag}_strata.csv", index=False)

    rows, pvals = [], {}
    for (w, j), g in f.groupby(["widths", "tau_exp"]):
        for arm, ref, fam in COMPARISONS:
            A, R = g[g.arm == arm].set_index(KEY), g[g.arm == ref].set_index(KEY)
            R = R.loc[A.index]
            d = pd.DataFrame({"cluster": A.cluster, "d_shots": R.shots - A.shots, "d_u": R.u_queries - A.u_queries}).groupby("cluster").mean()
            fam = fam if w in PRIMARY_WIDTHS or fam != "primary" else "secondary_other_partition"
            r = {"family": fam, "widths": w, "tau_exp": j, "arm": arm, "ref_arm": ref, "clusters": len(d),
                 "arm_shots": float(A.shots.mean()), "ref_shots": float(R.shots.mean()), "arm_u": float(A.u_queries.mean()),
                 "ref_u": float(R.u_queries.mean()), "arm_coverage": float(A.covered.mean()), "ref_coverage": float(R.covered.mean())}
            for k in ("d_shots", "d_u"):
                lo, hi = boot_ci(d[k].to_numpy(), args.boot)
                r |= {k: float(d[k].mean()), f"{k}_lo": lo, f"{k}_hi": hi}
            r["p_shots"] = signflip_p(d.d_shots.to_numpy(), args.perm)
            if fam == "primary":
                pvals[(w, j, arm)] = r["p_shots"]
            rows.append(r)
    comp = pd.DataFrame(rows)
    adj = holm(pvals)
    comp["p_holm"] = [adj.get((r.widths, r.tau_exp, r.arm)) if r.family == "primary" else None for r in comp.itertuples()]
    comp["positive"] = (comp.family == "primary") & (comp.p_holm.fillna(1) < 0.05) & (comp.d_shots > 0) & (comp.arm_coverage >= 0.95)
    comp.assign(status=status).to_csv(OUT / f"{tag}_comparisons.csv", index=False)
    print(f"wrote {tag}_*.csv; primary cells {len(pvals)}")


if __name__ == "__main__":
    main()
