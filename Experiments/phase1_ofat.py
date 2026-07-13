"""
Purpose
    Phase 1: one-factor-at-a-time sensitivity sweeps around the reference
    operating point R0, across a representative subset of (N, a) pairs.
Theory
    Holding every other parameter fixed while sweeping one isolates that
    parameter's marginal effect on success rate, producing the first
    evidence for which knob the algorithm is most sensitive to.
Inputs
    The cached ground-truth corpus (Experiments.harness.build_corpus).
Outputs
    Data/failure_study/phase1_ofat.csv and per-trial JSON under raw/.
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Experiments import harness

REPEATS = 5
PAIR_COUNT = 3

# NOTE: the initial design (10 repeats x 5 pairs x full-resolution grids,
# ~2400 trials) was found to run at only ~0.25 trials/s under this
# environment's process-pool + Aer overhead (~2.7h wall clock), which is
# impractical for this campaign. This reduced-resolution grid (~350 trials,
# ~20-25 min) keeps every qualitatively important region (below/at/above the
# QPE precision threshold, minimal/exhaustive candidate budgets, narrow/wide
# beams, and a wide shots range spanning the paper's own Lemma 3.3 regime)
# while fitting a practical session budget; this trade-off is reported
# explicitly as an experimental-scope limitation.
SWEEPS: dict[str, list[int]] = {
    "total_precision": [4, 6, 8, 10, 12],
    "window_size": [2, 3, 4, 6],
    "overlap": [0, 1, 2, 3],
    "candidate_count": [1, 2, 4, 8, 16],
    "max_paths": [1, 4, 16, 64, 256],
    "shots": [16, 64, 256, 1024, 4096, 16384],
}


def spec_for(param: str, value: int, pair: dict) -> dict:
    spec = {**harness.DEFAULT_SPEC, "N": pair["N"], "a": pair["a"], "phase": "phase1_ofat", "param_swept": param}
    if param == "total_precision":
        spec["total_precision"] = value
    elif param == "window_size":
        spec["window_size"] = value
        spec["overlap"] = min(spec["overlap"], value - 1)
        spec["total_precision"] = max(spec["total_precision"], value)
    elif param == "overlap":
        spec["overlap"] = value
    elif param == "candidate_count":
        spec["candidate_count"] = value
    elif param == "max_paths":
        spec["max_paths"] = value
    elif param == "shots":
        spec["shots"] = value
    return spec


def build_specs() -> list[dict]:
    corpus = harness.build_corpus()
    pairs = harness.representative_subset(corpus, k=PAIR_COUNT, max_work_qubits=5)
    specs: list[dict] = []
    for param, values in SWEEPS.items():
        for pair in pairs:
            for value in values:
                for repeat in range(REPEATS):
                    spec = spec_for(param, value, pair)
                    spec["seed"] = hash((param, value, pair["N"], pair["a"], repeat)) % (2**31)
                    spec["sweep_value"] = value
                    spec["pair_work_qubits"] = pair["work_qubits"]
                    spec["true_order"] = pair["order"]
                    specs.append(spec)
    return specs


def main() -> None:
    specs = build_specs()
    print(f"phase1_ofat: {len(specs)} trials")
    frame = harness.run_batch(specs, "phase1_ofat")
    print(frame["result_success"].mean(), "overall success rate")


if __name__ == "__main__":
    main()
