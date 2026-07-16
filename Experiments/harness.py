"""
Purpose
    Shared black-box experiment harness for the failure-mode study.
Theory
    Every experiment drives the repository exclusively through its public
    APIs (ShorConfig/WindowConfig/NoiseConfig, run_windowed_shor,
    run_windowed_shor_exact, sample_window_counts, and the Reconstruction
    stage functions). No repository source is modified. This module owns
    trial execution, timing/memory instrumentation, ground-truth lookup,
    Wilson confidence intervals, and JSON/CSV persistence so every phase
    script only has to describe *which* configurations to run.
Inputs
    Trial specification dicts (plain, JSON-serializable).
Outputs
    Per-trial JSON records under Data/failure_study/raw/, aggregated CSVs,
    and summary helpers used by analyze.py.
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
import math
import time
import tracemalloc
import uuid
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any, Callable

import pandas as pd

from Algorithms.paper_algorithm import run_windowed_shor, run_windowed_shor_exact, sample_window_counts
from Algorithms.standard_shor import multiplicative_order, run_standard_shor
from Circuits.modular_multiplication import ModularMultiplicationOperator
from Circuits.windowed_qpe import make_overlapping_windows, run_windowed_qpe_block
from Reconstruction.candidate_generation import generate_window_candidates, generate_window_candidates_by_coverage, normalize_counts
from Reconstruction.confidence import top1_vs_top2_confidence
from Reconstruction.continued_fraction import recover_order_and_factors
from Reconstruction.stitching import stitch_candidates
from Simulation.noise import build_noise_model
from config import DEFAULT_ADAPTIVE, NoiseConfig, ShorConfig, WindowConfig

FAILURE_STUDY_DIR = PROJECT_ROOT / "Data" / "failure_study"
RAW_DIR = FAILURE_STUDY_DIR / "raw"
PLOTS_DIR = FAILURE_STUDY_DIR / "plots"
TABLES_DIR = FAILURE_STUDY_DIR / "summary_tables"
CORPUS_PATH = FAILURE_STUDY_DIR / "corpus.json"

for _dir in (RAW_DIR, PLOTS_DIR, TABLES_DIR):
    _dir.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------
# Statistics helpers
# --------------------------------------------------------------------------

def wilson_ci(successes: int, n: int, z: float = 1.96) -> tuple[float, float, float]:
    """Return (point_estimate, low, high) Wilson score interval for a binomial proportion."""
    if n == 0:
        return (float("nan"), float("nan"), float("nan"))
    p = successes / n
    denom = 1 + z**2 / n
    center = (p + z**2 / (2 * n)) / denom
    half = (z * math.sqrt((p * (1 - p) + z**2 / (4 * n)) / n)) / denom
    return (p, max(0.0, center - half), min(1.0, center + half))


# --------------------------------------------------------------------------
# Ground-truth corpus
# --------------------------------------------------------------------------

def _candidate_bases(N: int) -> list[int]:
    """Return up to two coprime bases for N: one small-order, one large-order."""
    from math import gcd

    orders: list[tuple[int, int]] = []
    for a in range(2, N):
        if gcd(a, N) != 1:
            continue
        order = multiplicative_order(a, N)
        if order is not None and order > 1:
            orders.append((order, a))
    if not orders:
        return []
    orders.sort()
    small = orders[0][1]
    large = orders[-1][1]
    return sorted({small, large})


def build_corpus(force: bool = False) -> list[dict[str, Any]]:
    """Build (or load cached) ground-truth corpus of (N, a, order, factors)."""
    if CORPUS_PATH.exists() and not force:
        return json.loads(CORPUS_PATH.read_text())

    semiprimes = [15, 21, 33, 35, 51, 55, 65, 77, 85, 91, 95, 119, 143, 187, 221, 247, 299, 323, 377, 437, 493, 589]
    corpus: list[dict[str, Any]] = []
    for N in semiprimes:
        for a in _candidate_bases(N):
            result = run_standard_shor(N, a)
            if result.order is None or result.factors is None:
                continue
            work_qubits = max(1, math.ceil(math.log2(N)))
            corpus.append(
                {
                    "N": N,
                    "a": a,
                    "order": result.order,
                    "factors": list(result.factors),
                    "work_qubits": work_qubits,
                }
            )
    CORPUS_PATH.write_text(json.dumps(corpus, indent=2))
    return corpus


def representative_subset(corpus: list[dict[str, Any]], k: int = 5, max_work_qubits: int | None = None) -> list[dict[str, Any]]:
    """Return k (N, a) pairs spanning small/medium/large work-register sizes.

    `max_work_qubits` caps the pool before sampling: transpiling the dense
    UnitaryGate-based controlled modular multiplier scales very poorly with
    work-register size (Qiskit's generic unitary synthesis, not the Aer
    simulation itself, dominates -- see Phase 4). Parameter/noise sweeps
    that are not specifically studying scaling should stay below that cliff
    so their per-trial cost reflects the swept parameter, not N.
    """
    pool = [row for row in corpus if max_work_qubits is None or row["work_qubits"] <= max_work_qubits]
    ordered = sorted(pool, key=lambda row: (row["work_qubits"], row["N"]))
    if len(ordered) <= k:
        return ordered
    step = (len(ordered) - 1) / (k - 1)
    indices = sorted({round(i * step) for i in range(k)})
    return [ordered[i] for i in indices]


# --------------------------------------------------------------------------
# Trial specification / execution
# --------------------------------------------------------------------------

DEFAULT_SPEC: dict[str, Any] = {
    "total_precision": 8,
    "window_size": 4,
    "overlap": 2,
    "candidate_count": 2,
    "max_paths": 32,
    "shots": 2048,
    "noise_1q": 0.0,
    "noise_2q": 0.0,
    "noise_ro": 0.0,
    "simulator_method": "automatic",
    "seed": None,
}


def _noise_model_from_spec(spec: dict[str, Any]):
    enabled = spec["noise_1q"] > 0.0 or spec["noise_2q"] > 0.0 or spec["noise_ro"] > 0.0
    return build_noise_model(
        NoiseConfig(
            enabled=enabled,
            depolarizing_1q=spec["noise_1q"],
            depolarizing_2q=spec["noise_2q"],
            measurement_error=spec["noise_ro"],
        )
    )


def _to_shor_window(spec: dict[str, Any]) -> tuple[ShorConfig, WindowConfig]:
    shor = ShorConfig(
        N=spec["N"],
        a=spec["a"],
        phase_qubits=spec["total_precision"],
        shots=spec["shots"],
        simulator_method=spec.get("simulator_method", "automatic"),
        random_seed=spec["seed"] if spec.get("seed") is not None else ShorConfig.random_seed,
    )
    window = WindowConfig(
        total_precision=spec["total_precision"],
        window_size=spec["window_size"],
        overlap=spec["overlap"],
        candidate_count=spec["candidate_count"],
        max_paths=spec["max_paths"],
    )
    return shor, window


def run_trial(spec: dict[str, Any]) -> dict[str, Any]:
    """Execute one Aer-sampled trial exactly as configured and record everything."""
    full_spec = {**DEFAULT_SPEC, **spec}
    trial_id = full_spec.get("trial_id") or uuid.uuid4().hex
    shor, window = _to_shor_window(full_spec)
    noise_model = _noise_model_from_spec(full_spec)

    tracemalloc.start()
    started = time.perf_counter()
    try:
        result = run_windowed_shor(shor, window, noise_model=noise_model)
        error = None
    except Exception as exc:  # noqa: BLE001 - deliberately broad: this is a black-box probe
        result = None
        error = f"{type(exc).__name__}: {exc}"
    runtime = time.perf_counter() - started
    _current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    record: dict[str, Any] = {
        "trial_id": trial_id,
        "phase": full_spec.get("phase", "unknown"),
        "param_swept": full_spec.get("param_swept"),
        "spec": {k: v for k, v in full_spec.items() if k not in {"phase", "param_swept"}},
        "runtime_seconds": runtime,
        "peak_memory_bytes": peak,
        "error": error,
    }
    if result is not None:
        record["result"] = {
            "order": result.order,
            "factors": list(result.factors) if result.factors else None,
            "success": result.success,
            "best_phase": str(result.best_phase) if result.best_phase is not None else None,
            "num_candidates": len(result.reconstructed_phases),
        }
    else:
        record["result"] = {"order": None, "factors": None, "success": False, "best_phase": None, "num_candidates": 0}
    return record


def write_raw(record: dict[str, Any]) -> Path:
    """Persist one trial record as JSON under Data/failure_study/raw/."""
    path = RAW_DIR / f"{record['phase']}_{record['trial_id']}.json"
    path.write_text(json.dumps(record, indent=2))
    return path


def _flatten(record: dict[str, Any]) -> dict[str, Any]:
    row = {
        "trial_id": record["trial_id"],
        "phase": record["phase"],
        "param_swept": record["param_swept"],
        "runtime_seconds": record["runtime_seconds"],
        "peak_memory_bytes": record["peak_memory_bytes"],
        "error": record["error"],
    }
    row.update({f"spec_{k}": v for k, v in record["spec"].items()})
    row.update({f"result_{k}": v for k, v in record["result"].items()})
    return row


DEFAULT_MAX_WORKERS = min(3, (__import__("os").cpu_count() or 2))
_CHUNK_SIZE = 60


def run_batch(
    specs: list[dict[str, Any]],
    csv_name: str,
    trial_fn: Callable[[dict[str, Any]], dict[str, Any]] = run_trial,
    max_workers: int | None = DEFAULT_MAX_WORKERS,
    persist_raw: bool = True,
) -> pd.DataFrame:
    """Run a batch of trial specs (optionally in parallel), write raw JSON + a CSV, return a DataFrame.

    Runs in bounded-size chunks with a conservative worker count and falls
    back to sequential execution for any chunk whose pool dies (e.g. a
    native Aer/BLAS crash under heavy concurrent process load), so a single
    bad chunk cannot lose an entire multi-thousand-trial batch.
    """
    records: list[dict[str, Any]] = []
    if max_workers == 1:
        for spec in specs:
            record = trial_fn(spec)
            if persist_raw:
                write_raw(record)
            records.append(record)
    else:
        for start in range(0, len(specs), _CHUNK_SIZE):
            chunk = specs[start : start + _CHUNK_SIZE]
            try:
                with ProcessPoolExecutor(max_workers=max_workers) as pool:
                    futures = [pool.submit(trial_fn, spec) for spec in chunk]
                    for future in as_completed(futures):
                        record = future.result()
                        if persist_raw:
                            write_raw(record)
                        records.append(record)
            except Exception as exc:  # noqa: BLE001 - BrokenProcessPool or similar; degrade gracefully
                print(f"chunk starting at {start} failed in parallel ({exc}); retrying sequentially")
                for spec in chunk:
                    record = trial_fn(spec)
                    if persist_raw:
                        write_raw(record)
                    records.append(record)
            print(f"  progress: {min(start + _CHUNK_SIZE, len(specs))}/{len(specs)}")
    frame = pd.DataFrame([_flatten(r) for r in records])
    frame.to_csv(TABLES_DIR.parent / f"{csv_name}.csv", index=False)
    return frame


# --------------------------------------------------------------------------
# Classical-pipeline helpers (used by Phase 3 adversarial tests and Phase 5)
# --------------------------------------------------------------------------

def reconstruct_from_counts(
    N: int,
    a: int,
    total_precision: int,
    window_size: int,
    overlap: int,
    candidate_count: int,
    max_paths: int,
    window_counts: list[dict[str, int]],
):
    """Run only the classical stages (candidates -> stitching -> continued fraction).

    Returns (stitched_phases, recovered_or_None): `recovered_or_None` is the
    first stitched candidate that yields non-trivial factors, or None if
    either stitching collapsed to zero candidates or none of them recovered
    valid factors (continued-fraction failure).
    """
    specs = make_overlapping_windows(total_precision, window_size, overlap)
    groups = [
        generate_window_candidates(counts, spec.start, spec.width, total_precision, candidate_count)
        for spec, counts in zip(specs, window_counts, strict=True)
    ]
    stitched = stitch_candidates(groups, max_paths=max_paths)
    for phase in stitched:
        recovered = recover_order_and_factors(phase.phase, a, N)
        if recovered.factors is not None:
            return stitched, recovered
    return stitched, None


# --------------------------------------------------------------------------
# Confidence-guided adaptive reconstruction (Data/failure_study/FORMULATION.md)
# --------------------------------------------------------------------------

def _windows_needing_expansion(window_counts: list[dict[str, int]], candidate_count: int, delta_cov: float) -> int:
    """Return how many windows' fixed top-`candidate_count` truncation covers
    less than `1 - delta_cov` of the observed probability mass -- the signal
    module C's beam-width rule scales on, computed independently of whether
    module B's coverage-based candidate generation is actually enabled (so B
    and C remain separately ablatable, FORMULATION.md Section 5)."""
    target = 1.0 - delta_cov
    expanded = 0
    for counts in window_counts:
        ranked = sorted(normalize_counts(counts).values(), reverse=True)
        if sum(ranked[:candidate_count]) < target:
            expanded += 1
    return expanded


def _resample_shot_stopping(
    shor: ShorConfig,
    window: WindowConfig,
    window_counts: list[dict[str, int]],
    specs,
    epsilon: float,
    s_max: int,
    noise_model=None,
) -> list[dict[str, int]]:
    """Resample any window below the confidence target `1 - epsilon` in
    doubling increments up to a hard cap `s_max` (module D, FORMULATION.md
    Section 4.3). Uses the same operator/eigenstate as the original sample;
    each resample round uses a distinct seed so repeated calls are
    deterministic but not degenerate repeats of the same shots.

    Thresholds on the RAW top1_vs_top2_confidence (C_w), not the isotonic
    recalibration fit in Experiments/phase6_calibration.py. That
    recalibration is known to correct real overconfidence in C_w (Brier
    0.224 -> 0.199, FORMULATION.md Section 6), but is not applied here: this
    is a deliberate scope decision (documented in FORMULATION.md Section 6),
    not an oversight -- swapping in the recalibrated value would change
    every ablation arm that uses D (D, Full), invalidating the already-run
    phase6_ablation.py comparison, which was executed entirely against raw
    C_w. Wiring recalibration into this function is a valid follow-up but
    requires re-running the ablation to attribute any resulting change
    correctly."""
    op = ModularMultiplicationOperator(shor.a, shor.N)
    updated = [dict(c) for c in window_counts]
    for i, spec in enumerate(specs):
        round_index = 0
        while True:
            ranked = sorted(updated[i].values(), reverse=True)
            n1 = ranked[0] if ranked else 0
            n2 = ranked[1] if len(ranked) > 1 else 0
            s_w = sum(updated[i].values())
            if top1_vs_top2_confidence(n1, n2, s_w) >= 1.0 - epsilon:
                break
            extra_shots = min(s_w, s_max - s_w)
            if extra_shots <= 0:
                break
            round_index += 1
            extra_counts = run_windowed_qpe_block(
                op, eigenstate=1, spec=spec, shots=extra_shots,
                method=shor.simulator_method, seed=shor.random_seed + round_index,
                noise_model=noise_model,
            )
            for bits, count in extra_counts.items():
                updated[i][bits] = updated[i].get(bits, 0) + count
    return updated


def reconstruct_adaptive(
    shor: ShorConfig,
    window: WindowConfig,
    noise_model=None,
    use_coverage: bool = False,
    use_beam_rule: bool = False,
    use_shot_stopping: bool = False,
    delta_cov: float = DEFAULT_ADAPTIVE.delta_cov,
    beam_max: int = DEFAULT_ADAPTIVE.beam_max,
    epsilon: float = DEFAULT_ADAPTIVE.epsilon,
    s_max: int = DEFAULT_ADAPTIVE.s_max,
):
    """Run the classical reconstruction pipeline with modules B (candidate
    coverage), C (beam width), and D (shot stopping) independently
    switchable, for the ablation study (FORMULATION.md Section 5). With all
    three flags False this reproduces Baseline exactly."""
    window_counts, specs = sample_window_counts(shor, window, noise_model=noise_model)

    if use_shot_stopping:
        window_counts = _resample_shot_stopping(shor, window, window_counts, specs, epsilon, s_max, noise_model)

    if use_coverage:
        groups = [
            generate_window_candidates_by_coverage(counts, spec.start, spec.width, window.total_precision, delta_cov)
            for counts, spec in zip(window_counts, specs)
        ]
    else:
        groups = [
            generate_window_candidates(counts, spec.start, spec.width, window.total_precision, window.candidate_count)
            for counts, spec in zip(window_counts, specs)
        ]

    if use_beam_rule:
        expanded = _windows_needing_expansion(window_counts, window.candidate_count, delta_cov)
        max_paths = min(beam_max, window.max_paths * (1 + expanded))
    else:
        max_paths = window.max_paths

    stitched = stitch_candidates(groups, max_paths=max_paths)
    for phase in stitched:
        recovered = recover_order_and_factors(phase.phase, shor.a, shor.N)
        if recovered.factors is not None:
            return stitched, recovered
    return stitched, None
