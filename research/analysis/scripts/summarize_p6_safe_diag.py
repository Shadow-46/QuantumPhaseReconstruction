"""P6 safe-epsilon MECHANISM DIAGNOSTIC (dev split; not a performance claim, no tuning).

    python -m research.analysis.scripts.summarize_p6_safe_diag --run DIR

Boundary state at the moment of an overlap action (pre-batch counts, observable):
  rule_inapplicable    decoded lower part below the boundary is not 10..0 ->
                       the awqpe_ext rule cannot use the widened block at all
  resolved_by_flag     lower part is 10..0 and Algorithm 1 already flagged the upper
                       chunk (min rule applies)
  ambiguous_unflagged  lower part is 10..0 and the upper chunk is NOT flagged
                       (the eps-floor situation)
evaluated under eps_safe ("safe") and eps 0.9 ("09").

Tables (research/analysis/tables/p6diag_*.csv):
  action_states        share of overlap actions in each state, per arm and variant
  gain_by_state        immediate error gain / improve / worsen rates by state
  displacement         B4 - B2 success difference, displaced chunk shots and added
                       U-queries, grouped by where B4's overlap actions landed
                       (post-treatment grouping: descriptive, not causal)
  rescue_by_state      P5-policy (B2) failures rescued by one practical trigger
                       overlap batch vs one extra eig shot batch, by acted state
  v_cost               B4 U-queries per trial and immediate gains, v = 1 vs v = 2
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from research.awqpe.runner.core import REPO_ROOT

OUT = REPO_ROOT / "research" / "analysis" / "tables"
KEY = ["widths", "S0", "r", "dS", "phase_id", "replicate_id", "variant"]
DEC_STATE = {"awqpe_eps_safe_ext": "safe", "awqpe_ext": "09", "likelihood": "safe"}


def state(df: pd.DataFrame, tag: str) -> np.ndarray:
    lh, fl = df[f"lowerhalf{tag}"].astype(bool), df[f"flag{tag}"].astype(bool)
    return np.where(~lh, "rule_inapplicable", np.where(fl, "resolved_by_flag", "ambiguous_unflagged"))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", required=True)
    args = ap.parse_args()
    run = Path(args.run).resolve()
    tag = run.relative_to(REPO_ROOT).as_posix()
    files = sorted((run / "shards").glob("*.parquet"))
    f = pd.concat([pd.read_parquet(x, filters=[("row_type", "==", "final")]) for x in files], ignore_index=True)
    a = pd.concat([pd.read_parquet(x, filters=[("row_type", "==", "action")]) for x in files], ignore_index=True)
    for t in ("safe", "09"):
        a[f"state_{t}"] = state(a, t)
    OUT.mkdir(parents=True, exist_ok=True)
    meta = {"run": tag, "status": "DIAGNOSTIC / DEV"}

    rows = []
    for t in ("safe", "09"):
        c = pd.crosstab([a.variant, a.arm], a[f"state_{t}"], normalize="index").reset_index()
        rows.append(c.assign(epsilon_state=t))
    pd.concat(rows).assign(**meta).to_csv(OUT / "p6diag_action_states.csv", index=False)

    rows = []
    for dec, t in DEC_STATE.items():
        g = a.groupby(["variant", "arm", f"state_{t}"])[f"{dec}__gain_err_lsb"].agg(
            actions="size", mean_gain_err_lsb="mean", frac_improved=lambda x: (x > 1e-9).mean(), frac_worsened=lambda x: (x < -1e-9).mean()).reset_index()
        rows.append(g.rename(columns={f"state_{t}": "state"}).assign(decoder=dec, epsilon_state=t))
    pd.concat(rows).assign(**meta).to_csv(OUT / "p6diag_gain_by_state.csv", index=False)

    rows = []
    for dec, t in DEC_STATE.items():
        x = f[f.decoder == dec]
        piv = x.pivot_table(index=KEY, columns="arm", values=["tol", "chunk_shots", "u_queries"], aggfunc="first")
        a4 = a[a.arm == "B4_eig_plus_overlap"].assign(app=lambda d: d[f"lowerhalf{t}"].astype(bool))
        act = a4.groupby(KEY)["app"].max().rename("any_applicable")
        d = pd.DataFrame({"dtol": piv["tol"]["B4_eig_plus_overlap"] - piv["tol"]["B2_p5_eig"],
                          "shots_displaced": piv["chunk_shots"]["B2_p5_eig"] - piv["chunk_shots"]["B4_eig_plus_overlap"],
                          "u_added": piv["u_queries"]["B4_eig_plus_overlap"] - piv["u_queries"]["B2_p5_eig"]}).join(act)
        d["group"] = np.select([d.any_applicable.isna(), d.any_applicable.fillna(False).astype(bool)],
                               ["no_overlap_action", "action_on_applicable"], "actions_all_inapplicable")
        agg = {"trials": ("dtol", "size"), "mean_dtol": ("dtol", "mean"), "shots_displaced": ("shots_displaced", "mean"), "u_added": ("u_added", "mean")}
        g = d.reset_index().groupby(["variant", "group"]).agg(**agg).reset_index()
        tot = d.reset_index().groupby("variant").agg(**agg).reset_index().assign(group="ALL")
        rows.append(pd.concat([g, tot]).assign(decoder=dec, epsilon_state=t))
    pd.concat(rows).assign(**meta).to_csv(OUT / "p6diag_displacement.csv", index=False)

    rows = []
    for dec, t in DEC_STATE.items():
        x = f[f.decoder == dec]
        piv = x.pivot_table(index=KEY, columns="arm", values="tol", aggfunc="first")
        r2 = a[a.arm == "R_trigger_overlap"].groupby(KEY)[f"state_{t}"].first()
        d = pd.DataFrame({"fail": piv["B2_p5_eig"] < 0.5, "rescued_by_overlap": piv["R_trigger_overlap"] > 0.5,
                          "rescued_by_extra_shots": piv["R_eig_shots"] > 0.5}).join(r2)
        d = d[d.fail]
        d[f"state_{t}"] = d[f"state_{t}"].fillna("no_overlap_action")
        g = d.groupby(f"state_{t}").agg(p5_failures=("fail", "size"), rescued_by_overlap=("rescued_by_overlap", "mean"),
                                        rescued_by_extra_shots=("rescued_by_extra_shots", "mean")).reset_index()
        rows.append(g.rename(columns={f"state_{t}": "acted_state"}).assign(decoder=dec, epsilon_state=t))
    pd.concat(rows).assign(**meta).to_csv(OUT / "p6diag_rescue_by_state.csv", index=False)

    b4 = a[a.arm == "B4_eig_plus_overlap"]
    x = f[(f.decoder == "likelihood") & (f.arm.isin(["B4_eig_plus_overlap", "B2_p5_eig"]))]
    u = x.pivot_table(index=KEY, columns="arm", values="u_queries", aggfunc="first").reset_index()
    u["u_added_per_trial"] = u["B4_eig_plus_overlap"] - u["B2_p5_eig"]
    v = u.groupby("variant").u_added_per_trial.mean().reset_index()
    gains = b4.groupby("variant").agg(overlap_actions=("step", "size"), **{f"mean_gain_{d}": (f"{d}__gain_err_lsb", "mean") for d in ("awqpe_ext", "awqpe_eps_safe_ext", "likelihood")}).reset_index()
    v.merge(gains, on="variant").assign(**meta).to_csv(OUT / "p6diag_v_cost.csv", index=False)
    print(f"wrote p6diag_*.csv to {OUT}")


if __name__ == "__main__":
    main()
