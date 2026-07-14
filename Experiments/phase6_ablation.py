"""
Purpose
    Ablation study for confidence-guided adaptive reconstruction
    (Data/failure_study/FORMULATION.md Section 5): Baseline -> B -> C -> D ->
    Full, evaluated on the existing combined-stress grid.
Theory
    Reuses phase1b_combined_stress's exact 72-combination x 2-repeat grid and
    seeds (same instance N=21, a=19, order 6) so every arm is paired by seed
    with Baseline -- the textbook precondition for McNemar's test. B (Section
    4.1), C (Section 4.2), and D (Section 4.3) are each run in isolation and
    combined (Full), letting B-vs-C isolate whether REPORT.md Section 6's
    "both budgets starved together" coupling is real.
Inputs
    None (grid defined below, matching phase1b_combined_stress.py exactly).
Outputs
    Data/failure_study/phase6_ablation.csv (per-trial success by arm) and
    console-printed McNemar's/bootstrap summary per arm vs. Baseline.
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import itertools
import random

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from Experiments import harness
from config import DEFAULT_ADAPTIVE, ShorConfig, WindowConfig

N, A = 21, 19
REPEATS = 2
GRID = {
    "candidate_count": [1, 2],
    "max_paths": [1, 2],
    "shots": [4, 8, 16],
    "total_precision": [16, 24, 32],
    "overlap": [0, 3],
}
WINDOW_SIZE = 4

ARMS = {
    "baseline": dict(use_coverage=False, use_beam_rule=False, use_shot_stopping=False),
    "B": dict(use_coverage=True, use_beam_rule=False, use_shot_stopping=False),
    "C": dict(use_coverage=False, use_beam_rule=True, use_shot_stopping=False),
    "D": dict(use_coverage=False, use_beam_rule=False, use_shot_stopping=True),
    "full": dict(use_coverage=True, use_beam_rule=True, use_shot_stopping=True),
}


def all_combos() -> list[tuple]:
    return list(
        itertools.product(GRID["candidate_count"], GRID["max_paths"], GRID["shots"], GRID["total_precision"], GRID["overlap"])
    )


def run_one_trial(combo: tuple, repeat: int) -> dict:
    candidate_count, max_paths, shots, total_precision, overlap = combo
    seed = hash((*combo, repeat)) % (2**31)
    shor = ShorConfig(N=N, a=A, phase_qubits=total_precision, shots=shots, random_seed=seed)
    window = WindowConfig(
        total_precision=total_precision, window_size=WINDOW_SIZE, overlap=overlap,
        candidate_count=candidate_count, max_paths=max_paths,
    )
    row = {
        "candidate_count": candidate_count, "max_paths": max_paths, "shots": shots,
        "total_precision": total_precision, "overlap": overlap, "repeat": repeat, "seed": seed,
    }
    for arm_name, flags in ARMS.items():
        _stitched, recovered = harness.reconstruct_adaptive(
            shor, window,
            delta_cov=DEFAULT_ADAPTIVE.delta_cov, beam_max=DEFAULT_ADAPTIVE.beam_max,
            epsilon=DEFAULT_ADAPTIVE.epsilon, s_max=DEFAULT_ADAPTIVE.s_max,
            **flags,
        )
        row[f"success_{arm_name}"] = recovered is not None and recovered.factors is not None
    return row


def run_ablation_grid(resume: bool = True) -> pd.DataFrame:
    """Run every (combo, repeat) trial (all 5 arms each) sequentially.

    Deliberately not parallelized: concurrent ProcessPoolExecutor workers
    each running real Aer/qiskit circuits proved unreliable in this
    environment (hangs with no output, not clean crashes -- consistent with
    the native-worker instability REPORT.md Appendix A already documents for
    this exact combination). Sequential is slower per trial but the only
    mode observed to make reliable forward progress; resumable so a session
    that's interrupted partway doesn't lose completed trials."""
    csv_path = harness.FAILURE_STUDY_DIR / "phase6_ablation.csv"
    done_seeds: set[int] = set()
    rows: list[dict] = []
    if resume and csv_path.exists():
        existing = pd.read_csv(csv_path)
        rows = existing.to_dict("records")
        done_seeds = set(existing["seed"])

    work_items = [(combo, repeat) for combo in all_combos() for repeat in range(REPEATS)]
    total = len(work_items)
    print(f"  {len(done_seeds)}/{total} already done, resuming", flush=True)
    for i, (combo, repeat) in enumerate(work_items):
        seed = hash((*combo, repeat)) % (2**31)
        if seed in done_seeds:
            continue
        rows.append(run_one_trial(combo, repeat))
        pd.DataFrame(rows).to_csv(csv_path, index=False)
        print(f"  progress: {i + 1}/{total}", flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(csv_path, index=False)
    return frame


def mcnemar_test(baseline: pd.Series, arm: pd.Series) -> dict:
    """Exact McNemar's test on paired binary outcomes, plus a bootstrap 95%
    CI on the paired success-rate difference."""
    baseline_only = int(((baseline) & (~arm)).sum())
    arm_only = int(((~baseline) & (arm)).sum())
    discordant = baseline_only + arm_only
    p_value = 1.0 if discordant == 0 else binomtest(min(baseline_only, arm_only), discordant, 0.5).pvalue

    diffs = (arm.astype(int) - baseline.astype(int)).to_numpy()
    rng = random.Random(0)
    n = len(diffs)
    boot_means = []
    for _ in range(2000):
        idx = [rng.randrange(n) for _ in range(n)]
        boot_means.append(np.mean(diffs[idx]))
    ci_low, ci_high = np.percentile(boot_means, [2.5, 97.5])

    return {
        "baseline_success_rate": float(baseline.mean()),
        "arm_success_rate": float(arm.mean()),
        "paired_diff": float(arm.mean() - baseline.mean()),
        "diff_ci_low": float(ci_low),
        "diff_ci_high": float(ci_high),
        "mcnemar_p_value": float(p_value),
        "baseline_only_fixed_count": baseline_only,
        "arm_only_broke_count": arm_only,
    }


def summarize(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for arm_name in ARMS:
        if arm_name == "baseline":
            continue
        stats = mcnemar_test(frame["success_baseline"], frame[f"success_{arm_name}"])
        rows.append({"arm": arm_name, **stats})
    b_vs_c = mcnemar_test(frame["success_B"], frame["success_C"])
    rows.append({"arm": "B_vs_C", **b_vs_c})
    return pd.DataFrame(rows)


def main() -> None:
    print(f"phase6_ablation: {len(all_combos()) * REPEATS} trials x {len(ARMS)} arms")
    frame = run_ablation_grid()
    summary = summarize(frame)
    summary.to_csv(harness.FAILURE_STUDY_DIR / "summary_tables" / "phase6_ablation_summary.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
