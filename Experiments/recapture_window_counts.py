"""
Purpose
    Recapture per-window measurement counts for the existing phase1_ofat and
    phase1b_combined_stress trials, which were never persisted originally.
Theory
    The calibration check (Data/failure_study/FORMULATION.md Section 6)
    needs per-window bit counts paired with ground-truth window correctness,
    but the original raw JSON under Data/failure_study/raw/ only stored
    trial-level summary outcomes. Aer sampling is deterministic given a seed,
    so replaying each existing trial's exact stored (N, a, seed, shots,
    window geometry, noise) spec through sample_window_counts reproduces the
    same trial with the missing per-window data added -- this is not a new
    experimental design, only an instrumentation gap fix.
Inputs
    Existing Data/failure_study/raw/phase1_ofat_*.json and
    phase1b_combined_stress_*.json records.
Outputs
    Data/failure_study/raw_windows/{phase}_{trial_id}.json, each containing
    the per-window counts and window specs (start/width/total_precision) for
    one existing trial. Resumable: skips trial_ids already recaptured.
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import json

from Algorithms.paper_algorithm import sample_window_counts
from Experiments import harness

RAW_WINDOWS_DIR = harness.FAILURE_STUDY_DIR / "raw_windows"
RAW_WINDOWS_DIR.mkdir(parents=True, exist_ok=True)

SOURCE_PATTERNS = ("phase1_ofat_*.json", "phase1b_combined_stress_*.json")


def _source_records() -> list[dict]:
    records = []
    for pattern in SOURCE_PATTERNS:
        for path in harness.RAW_DIR.glob(pattern):
            records.append(json.loads(path.read_text()))
    return records


def recapture_one(record: dict) -> dict:
    """Replay one existing trial's exact spec and return its per-window counts."""
    spec = {**harness.DEFAULT_SPEC, **record["spec"]}
    shor, window = harness._to_shor_window(spec)  # noqa: SLF001 - reusing harness's own spec parsing for exact reproduction
    noise_model = harness._noise_model_from_spec(spec)  # noqa: SLF001
    counts, specs = sample_window_counts(shor, window, noise_model=noise_model)
    return {
        "trial_id": record["trial_id"],
        "phase": record["phase"],
        "spec": spec,
        "windows": [{"start": s.start, "width": s.width, "total_precision": s.total_precision} for s in specs],
        "window_counts": counts,
    }


def main() -> None:
    existing_ids = {p.stem[-32:] for p in RAW_WINDOWS_DIR.glob("*.json")}
    records = _source_records()
    todo = [r for r in records if r["trial_id"] not in existing_ids]
    print(f"recapture_window_counts: {len(records)} source trials, {len(todo)} remaining")

    for i, record in enumerate(todo):
        recaptured = recapture_one(record)
        path = RAW_WINDOWS_DIR / f"{record['phase']}_{record['trial_id']}.json"
        path.write_text(json.dumps(recaptured, indent=2))
        if (i + 1) % 50 == 0:
            print(f"  progress: {i + 1}/{len(todo)}")

    print(f"done. {len(list(RAW_WINDOWS_DIR.glob('*.json')))} total recaptured trials.")


if __name__ == "__main__":
    main()
