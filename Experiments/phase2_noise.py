"""
Purpose
    Phase 2: single-source and combined Aer noise-model sweeps.
Theory
    Depolarizing 1q/2q error and readout error are varied independently
    (to find each channel's individual breakdown threshold) and then
    combined factorially at a few representative levels (to test for
    super-/sub-additive interaction between noise sources).
Inputs
    The cached ground-truth corpus (Experiments.harness.build_corpus).
Outputs
    Data/failure_study/phase2_noise_single.csv,
    Data/failure_study/phase2_noise_combined.csv, and per-trial JSON.
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from itertools import product

from Experiments import harness

# See the note in phase1_ofat.py: throughput under this environment's
# process-pool + Aer overhead (~0.25-1 trial/s) requires a reduced-resolution
# grid to fit a practical session budget; scope reduction is reported
# explicitly rather than silently assumed.
REPEATS_SINGLE = 4
REPEATS_COMBINED = 3
PAIR_COUNT = 3

SINGLE_LEVELS = [0.0, 0.003, 0.01, 0.03, 0.1, 0.2, 0.3]
COMBINED_LEVELS = {"none": 0.0, "med": 0.05, "high": 0.2}


def build_single_source_specs() -> list[dict]:
    corpus = harness.build_corpus()
    pairs = harness.representative_subset(corpus, k=PAIR_COUNT, max_work_qubits=5)
    channels = ["noise_1q", "noise_2q", "noise_ro"]
    specs: list[dict] = []
    for channel in channels:
        for level in SINGLE_LEVELS:
            for pair in pairs:
                for repeat in range(REPEATS_SINGLE):
                    spec = {
                        **harness.DEFAULT_SPEC,
                        "N": pair["N"],
                        "a": pair["a"],
                        "phase": "phase2_noise_single",
                        "param_swept": channel,
                        "sweep_value": level,
                        "pair_work_qubits": pair["work_qubits"],
                        "true_order": pair["order"],
                    }
                    spec[channel] = level
                    spec["seed"] = hash((channel, level, pair["N"], pair["a"], repeat)) % (2**31)
                    specs.append(spec)
    return specs


def build_combined_specs() -> list[dict]:
    corpus = harness.build_corpus()
    pairs = harness.representative_subset(corpus, k=PAIR_COUNT, max_work_qubits=5)
    specs: list[dict] = []
    for l1q, l2q, lro in product(COMBINED_LEVELS.items(), COMBINED_LEVELS.items(), COMBINED_LEVELS.items()):
        (n1, v1), (n2, v2), (n3, v3) = l1q, l2q, lro
        for pair in pairs:
            for repeat in range(REPEATS_COMBINED):
                spec = {
                    **harness.DEFAULT_SPEC,
                    "N": pair["N"],
                    "a": pair["a"],
                    "phase": "phase2_noise_combined",
                    "param_swept": "combined",
                    "noise_1q": v1,
                    "noise_2q": v2,
                    "noise_ro": v3,
                    "level_1q": n1,
                    "level_2q": n2,
                    "level_ro": n3,
                    "pair_work_qubits": pair["work_qubits"],
                    "true_order": pair["order"],
                }
                spec["seed"] = hash((n1, n2, n3, pair["N"], pair["a"], repeat)) % (2**31)
                specs.append(spec)
    return specs


def main() -> None:
    single_specs = build_single_source_specs()
    print(f"phase2_noise_single: {len(single_specs)} trials")
    single_frame = harness.run_batch(single_specs, "phase2_noise_single")
    print(single_frame["result_success"].mean(), "single-source overall success rate")

    combined_specs = build_combined_specs()
    print(f"phase2_noise_combined: {len(combined_specs)} trials")
    combined_frame = harness.run_batch(combined_specs, "phase2_noise_combined")
    print(combined_frame["result_success"].mean(), "combined-noise overall success rate")


if __name__ == "__main__":
    main()
