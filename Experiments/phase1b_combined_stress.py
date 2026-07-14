"""
Purpose
    Phase 1b: a combined (non-OFAT) worst-case stress test on the hardest
    small instance available under the work_qubits<=5 compilation-cliff
    ceiling, simultaneously starving candidate_count, max_paths, and shots
    while maximizing total_precision/window count.
Theory
    Phase 1's one-factor-at-a-time sweep found zero failures (435/435
    success) because varying a single parameter at a time never pushed the
    classical pipeline hard enough. This combines every "hard" regime at
    once to test whether starvation on multiple axes together induces
    failure even at N=21 (order 6), independent of any noise or adversarial
    corruption.
Inputs
    None (fixed instance N=21, a=19; grid defined below).
Outputs
    Data/failure_study/phase1b_combined_stress.csv and per-trial JSON under
    raw/. Resumable: scans existing raw records first and only runs the
    combinations not already covered, so a prior partial/interrupted run is
    not repeated.
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
import json

import pandas as pd

from Experiments import harness

N, A = 21, 19
REPEATS = 2
GRID = {
    "candidate_count": [1, 2],
    "max_paths": [1, 2],
    "shots": [4, 8, 16],
    "total_precision": [16, 24, 32],
    "overlap": [0, 3],
}


def _combo_key(spec: dict) -> tuple:
    return (spec["candidate_count"], spec["max_paths"], spec["shots"], spec["total_precision"], spec["overlap"])


def _existing_records() -> list[dict]:
    records = []
    for path in harness.RAW_DIR.glob("phase1b_combined_stress_*.json"):
        records.append(json.loads(path.read_text()))
    return records


def all_combos() -> list[tuple]:
    return list(
        itertools.product(GRID["candidate_count"], GRID["max_paths"], GRID["shots"], GRID["total_precision"], GRID["overlap"])
    )


def build_missing_specs(existing: list[dict]) -> list[dict]:
    """For each of the 72 parameter combinations, top up to REPEATS trials
    using however many are already present in `existing` (the original
    interrupted run had no stored repeat index, so completeness is tracked
    per-combo by count, not by matching an exact repeat number)."""
    existing_counts: dict[tuple, int] = {}
    for record in existing:
        key = _combo_key(record["spec"])
        existing_counts[key] = existing_counts.get(key, 0) + 1

    missing: list[dict] = []
    for combo in all_combos():
        candidate_count, max_paths, shots, total_precision, overlap = combo
        have = existing_counts.get(combo, 0)
        need = max(0, REPEATS - have)
        for extra in range(need):
            spec = {
                **harness.DEFAULT_SPEC,
                "N": N, "a": A,
                "total_precision": total_precision, "window_size": 4, "overlap": overlap,
                "candidate_count": candidate_count, "max_paths": max_paths, "shots": shots,
                "phase": "phase1b_combined_stress", "param_swept": "combined",
                "seed": hash((*combo, have + extra)) % (2**31),
            }
            missing.append(spec)
    return missing


def main() -> None:
    existing = _existing_records()
    missing = build_missing_specs(existing)
    total_target = len(all_combos()) * REPEATS
    print(f"phase1b_combined_stress: {total_target} target, {len(existing)} already done, {len(missing)} remaining")

    if missing:
        frame = harness.run_batch(missing, "phase1b_combined_stress_new")
        print(frame["result_success"].value_counts())

    # Aggregate ALL records (pre-existing + newly run) into the canonical CSV.
    all_records = _existing_records()
    rows = [harness._flatten(r) for r in all_records]  # noqa: SLF001 - reusing harness's own flattening for consistency
    frame = pd.DataFrame(rows).drop_duplicates(subset=["trial_id"])
    frame.to_csv(harness.FAILURE_STUDY_DIR / "phase1b_combined_stress.csv", index=False)
    print(f"wrote {len(frame)} total rows to phase1b_combined_stress.csv")
    print("overall success rate:", frame["result_success"].mean())


if __name__ == "__main__":
    main()
