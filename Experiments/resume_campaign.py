"""
Purpose
    Single entry point to run or resume the full W1/W2/W3 generalization
    campaign (plan: C:\\Users\\sanja\\.claude\\plans\\handoff-prompt-goofy-hollerith.md).
    Safe to run any number of times, including after the machine slept,
    the process was killed, or a previous run finished some stages fully:
    every stage checks its own CSV first and only computes missing trials,
    so a completed stage costs a few seconds to confirm, not a re-run.
Theory
    Stage order matters once (canonical re-run must be verified against the
    published phase6_ablation.csv before anything downstream is trusted),
    but is otherwise just the pre-registered plan (Part 6): instrumented
    canonical re-run -> verification gate -> 3 held-out grids -> matched-
    resource control arms built from the canonical grid's diagnostics.
Inputs
    None.
Outputs
    Data/failure_study/phase6_ablation_diagnostics.csv
    Data/failure_study/phase7_heldout_N{15,21}_a{2,8,2}.csv
    Data/failure_study/phase7_matched_resource.csv
    Data/failure_study/summary_tables/phase7_heldout_summary.csv
    Data/failure_study/summary_tables/phase7_matched_resource_summary.csv
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import pandas as pd

from Experiments import harness
from Experiments import phase7_generalization as p7


def main() -> None:
    print("=== Stage 1: instrumented canonical (21,19) re-run ===", flush=True)
    canonical_diag_path = harness.FAILURE_STUDY_DIR / "phase6_ablation_diagnostics.csv"
    canonical = p7.run_instance_grid(21, 19, canonical_diag_path)
    print("=== Stage 1: verifying against published phase6_ablation.csv ===", flush=True)
    p7.verify_diagnostics_match_canonical(canonical)
    print("=== Stage 1: COMPLETE AND VERIFIED ===\n", flush=True)

    print("=== Stages 2-4: held-out instances ===", flush=True)
    heldout_summary = p7.run_held_out_instances()
    print(heldout_summary.to_string(index=False), flush=True)
    print("=== Stages 2-4: COMPLETE ===\n", flush=True)

    print("=== Stage 5: matched-resource control arms ===", flush=True)
    matched_path = harness.FAILURE_STUDY_DIR / "phase7_matched_resource.csv"
    matched = p7.run_matched_resource_grid(canonical, matched_path)
    matched_summary = p7.summarize_matched_resource(matched)
    matched_summary.to_csv(harness.FAILURE_STUDY_DIR / "summary_tables" / "phase7_matched_resource_summary.csv", index=False)
    print(matched_summary.to_string(index=False), flush=True)
    print("=== Stage 5: COMPLETE ===\n", flush=True)

    print("=== ALL STAGES COMPLETE ===", flush=True)


if __name__ == "__main__":
    main()
