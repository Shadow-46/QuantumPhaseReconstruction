"""
Purpose
    Self-service status check for the W1/W2/W3 generalization campaign
    (Experiments/phase7_generalization.py): prints how many of each grid's
    144 trials are done and flags whether a grid looks stalled.
Theory
    Every grid is a resumable CSV (row count - 1 header = trials done);
    "stalled" is a heuristic (no file write in a long time while the grid
    isn't finished) since a single trial can legitimately take up to ~20-25
    minutes in the slowest ablation cells (large window count + shots=4).
Inputs
    None -- reads the fixed set of CSVs this campaign writes to.
Outputs
    Console report only.
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import time

PROJECT_ROOT = Path(__file__).resolve().parents[1]
FAILURE_STUDY_DIR = PROJECT_ROOT / "Data" / "failure_study"

# (label, csv filename, expected total trials)
STAGES = [
    ("Stage 1: canonical (21,19) instrumented re-run", "phase6_ablation_diagnostics.csv", 144),
    ("Stage 2: held-out (15, a=2, r=4)", "phase7_heldout_N15_a2.csv", 144),
    ("Stage 3: held-out (21, a=8, r=2)", "phase7_heldout_N21_a8.csv", 144),
    ("Stage 4: held-out (21, a=2, r=6)", "phase7_heldout_N21_a2.csv", 144),
    ("Stage 5: matched-resource control arms", "phase7_matched_resource.csv", 144),
]

STALL_WARNING_SECONDS = 60 * 60  # flag if no write in >1h and not finished


def count_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with open(path, "r", encoding="utf-8") as f:
        return max(0, sum(1 for _ in f) - 1)  # minus header


def main() -> None:
    now = time.time()
    print(f"Status as of {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    any_incomplete = False
    for label, filename, total in STAGES:
        path = FAILURE_STUDY_DIR / filename
        done = count_rows(path)
        if not path.exists():
            print(f"  {label}: not started")
            any_incomplete = True
            continue
        mtime = path.stat().st_mtime
        age_min = (now - mtime) / 60
        status = "DONE" if done >= total else f"{done}/{total}"
        line = f"  {label}: {status}"
        if done < total:
            any_incomplete = True
            line += f"  (last write {age_min:.0f} min ago)"
            if (now - mtime) > STALL_WARNING_SECONDS:
                line += "  <-- no progress in over an hour, check if the process is still running"
        print(line)

    print()
    if any_incomplete:
        print("Not finished yet. If nothing is running (check Task Manager for python.exe,")
        print("or run: powershell -Command \"Get-Process python -ErrorAction SilentlyContinue\"),")
        print("relaunch with: python Experiments/resume_campaign.py")
        print("It is always safe to run resume_campaign.py -- finished stages are skipped instantly.")
    else:
        print("All stages complete.")


if __name__ == "__main__":
    main()
