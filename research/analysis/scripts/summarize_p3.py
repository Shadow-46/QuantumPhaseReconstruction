"""Summarise Phase 3 (initial-state sensitivity) into committed CSV tables. PRELIMINARY (dev phases).

    python -m research.analysis.scripts.summarize_p3 [--run DIR]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from research.awqpe.evaluation.metrics import circular_error
from research.awqpe.runner.core import REPO_ROOT, latest_run, load_results

OUT = REPO_ROOT / "research" / "analysis" / "tables"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", default=None)
    args = ap.parse_args()
    run = Path(args.run).resolve() if args.run else latest_run("p3_state_sensitivity", include_pilot=False)
    tag = run.relative_to(REPO_ROOT).as_posix()
    d = load_results(run)
    OUT.mkdir(parents=True, exist_ok=True)

    cv = pd.DataFrame(json.loads((run / "manifest.json").read_text(encoding="utf-8"))["case_verification"])
    cv.assign(run=tag).to_csv(OUT / "p3_case_verification.csv", index=False)

    rows = []
    e = d[d.family == "E_modmul_one"].copy()
    e["r"] = e.instance.str.extract(r"r=(\d+)")[0].astype(int)
    for (inst, w, s, dec), g in e.groupby(["instance", "widths", "shots_per_block", "decoder"]):
        r, n = int(g.r.iloc[0]), int(g.n.iloc[0])
        comps = np.arange(1, r) / r
        err = np.min(circular_error(g.phi_hat.to_numpy()[:, None], comps[None, :]), axis=1)
        rows.append({"instance": inst, "widths": w, "shots_per_block": s, "decoder": dec, "state": "|1> (mixture)", "trials": len(g),
                     "p_hit_nonzero_s_over_r": float((err <= 2.0**-n).mean()), "p_output_zero": float((g.phi_hat == 0).mean()),
                     "flag_rate": float((g.flags_mask > 0).mean()) if dec != "likelihood" else np.nan,
                     "true_max_top2_ratio": float(g.true_max_top2_ratio.mean())})
    dm = d[d.family == "D_modmul_eigen"]
    for (inst, w, s, dec), g in dm.groupby(["instance", "widths", "shots_per_block", "decoder"]):
        rows.append({"instance": inst, "widths": w, "shots_per_block": s, "decoder": dec, "state": "eigenstates |u_s>", "trials": len(g),
                     "p_tol_target": float(g.tol_dominant.mean()), "flag_rate": float((g.flags_mask > 0).mean()) if dec != "likelihood" else np.nan,
                     "true_max_top2_ratio": float(g.true_max_top2_ratio.mean())})
    pd.DataFrame(rows).assign(run=tag, status="PRELIMINARY").to_csv(OUT / "p3_modmul_eigen_vs_one.csv", index=False)

    a = d[d.family == "A_phase_gate"].copy()
    a["grid_location"] = np.select([a.is_dyadic_n, a.boundary_hardness < 0.05, (a.final_residual - 0.5).abs() < 0.05],
                                   ["dyadic_n", "boundary_hard", "final_half"], "generic")
    ga = a.groupby(["widths", "grid_location", "shots_per_block", "decoder"]).agg(
        trials=("tol_dominant", "size"), phases=("case_id", "nunique"), p_tol=("tol_dominant", "mean"),
        true_max_top2_ratio=("true_max_top2_ratio", "mean")).reset_index()
    ga.assign(run=tag, status="PRELIMINARY").to_csv(OUT / "p3_eigenstate_grid_location.csv", index=False)

    fam = d.groupby(["family", "state_prep", "widths", "shots_per_block", "decoder"]).agg(
        trials=("tol_dominant", "size"), cases=("case_id", "nunique"), p_tol_dominant=("tol_dominant", "mean"),
        p_tol_any=("tol_any", "mean"), max_eigen_residual=("eigen_residual", "max"), max_qiskit_diff=("qiskit_max_abs_diff", "max")).reset_index()
    fam.assign(run=tag, status="PRELIMINARY").to_csv(OUT / "p3_family_summary.csv", index=False)
    print(f"wrote p3_*.csv to {OUT}")


if __name__ == "__main__":
    main()
