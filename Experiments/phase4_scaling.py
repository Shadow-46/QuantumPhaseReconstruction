"""
Purpose
    Phase 4: scaling study of runtime, memory, and circuit qubit count as N
    (work-register size), total_precision (number of windows), and
    window_size (per-block width) each grow.
Theory
    The paper's central resource claim is that per-block qubit count and
    depth stay decoupled from total precision. This phase measures whether
    that decoupling actually holds in wall-clock/memory practice, and
    separately measures the classical-side cost driven by growing the
    number of windows and beam-search combinatorics.
Inputs
    The cached ground-truth corpus.
Outputs
    Data/failure_study/phase4_scaling_by_N.csv,
    phase4_scaling_by_precision.csv, phase4_scaling_by_window_size.csv.
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from concurrent.futures import ProcessPoolExecutor, TimeoutError as FutureTimeoutError

from Circuits.modular_multiplication import ModularMultiplicationOperator
from Circuits.windowed_qpe import build_windowed_qpe_circuit, make_overlapping_windows
from Experiments import harness

REPEATS = 3
QUBIT_CAP = 5  # work-qubit ceiling for *executed* trials (see transpile_cliff finding below)
TRANSPILE_TIMEOUT_S = 45


def circuit_stats(N: int, a: int, window_size: int, overlap: int, total_precision: int) -> dict:
    op = ModularMultiplicationOperator(a, N)
    specs = make_overlapping_windows(total_precision, window_size, overlap)
    depths, qubit_counts = [], []
    for spec in specs:
        circuit = build_windowed_qpe_circuit(op, eigenstate=1, spec=spec)
        depths.append(circuit.depth())
        qubit_counts.append(circuit.num_qubits)
    return {
        "num_windows": len(specs),
        "max_circuit_qubits": max(qubit_counts),
        "max_circuit_depth": max(depths),
        "work_qubits": op.qubits,
    }


def scaling_by_N() -> None:
    corpus = harness.build_corpus()
    specs = []
    for pair in corpus:
        if pair["work_qubits"] > QUBIT_CAP:
            continue
        for repeat in range(REPEATS):
            spec = {
                **harness.DEFAULT_SPEC,
                "N": pair["N"],
                "a": pair["a"],
                "phase": "phase4_scaling_by_N",
                "param_swept": "N",
                "sweep_value": pair["N"],
                "pair_work_qubits": pair["work_qubits"],
                "true_order": pair["order"],
                "seed": hash(("scaleN", pair["N"], pair["a"], repeat)) % (2**31),
            }
            specs.append(spec)
    print(f"phase4_scaling_by_N: {len(specs)} trials")
    frame = harness.run_batch(specs, "phase4_scaling_by_N")
    stats_rows = []
    for pair in corpus:
        if pair["work_qubits"] > QUBIT_CAP:
            continue
        stats_rows.append({"N": pair["N"], **circuit_stats(pair["N"], pair["a"], 4, 2, 8)})
    import pandas as pd

    pd.DataFrame(stats_rows).to_csv(harness.FAILURE_STUDY_DIR / "phase4_circuit_stats_by_N.csv", index=False)
    print(frame.groupby("spec_N")["runtime_seconds"].mean())


def scaling_by_precision() -> None:
    """Fix window_size=4; vary total_precision (more windows, same block width)."""
    corpus = harness.build_corpus()
    pair = harness.representative_subset(corpus, k=3, max_work_qubits=QUBIT_CAP)[1]  # safe mid-size representative
    specs = []
    for precision in [4, 8, 12, 16, 20, 24, 28, 32]:
        for repeat in range(REPEATS):
            spec = {
                **harness.DEFAULT_SPEC,
                "N": pair["N"],
                "a": pair["a"],
                "total_precision": precision,
                "phase": "phase4_scaling_by_precision",
                "param_swept": "total_precision",
                "sweep_value": precision,
                "seed": hash(("scaleP", precision, repeat)) % (2**31),
            }
            specs.append(spec)
    print(f"phase4_scaling_by_precision: {len(specs)} trials")
    frame = harness.run_batch(specs, "phase4_scaling_by_precision")
    print(frame.groupby("spec_total_precision")["runtime_seconds"].mean())


def scaling_by_window_size() -> None:
    """Fix total_precision=16; vary window_size (fewer/larger windows)."""
    corpus = harness.build_corpus()
    pair = harness.representative_subset(corpus, k=3, max_work_qubits=QUBIT_CAP)[1]
    specs = []
    for window_size in [2, 4, 8, 16]:
        overlap = max(0, window_size // 2 - 1) if window_size > 1 else 0
        for repeat in range(REPEATS):
            spec = {
                **harness.DEFAULT_SPEC,
                "N": pair["N"],
                "a": pair["a"],
                "total_precision": 16,
                "window_size": window_size,
                "overlap": overlap,
                "phase": "phase4_scaling_by_window_size",
                "param_swept": "window_size",
                "sweep_value": window_size,
                "seed": hash(("scaleW", window_size, repeat)) % (2**31),
            }
            specs.append(spec)
    print(f"phase4_scaling_by_window_size: {len(specs)} trials")
    frame = harness.run_batch(specs, "phase4_scaling_by_window_size")
    print(frame.groupby("spec_window_size")["runtime_seconds"].mean())


def structural_stats_full_corpus() -> None:
    """Parameter-only (no gate/circuit construction at all) qubit-count
    accounting across the FULL corpus (work_qubits 4-10). Deliberately
    avoids calling build_windowed_qpe_circuit/controlled_power_gate here:
    those go through UnitaryGate(..., check_input=True), whose O(dim^3)
    unitarity validation is itself a nontrivial cost at dim=2^10, separate
    from (and in addition to) the transpilation cliff measured below. This
    function isolates the pure claim "per-block qubit count = window_size +
    work_qubits, independent of total_precision" from either of those costs.
    """
    corpus = harness.build_corpus()
    rows = []
    for pair in corpus:
        specs = make_overlapping_windows(8, 4, 2)
        rows.append(
            {
                "N": pair["N"],
                "work_qubits": pair["work_qubits"],
                "num_windows": len(specs),
                "max_circuit_qubits": 4 + pair["work_qubits"],
            }
        )
    import pandas as pd

    pd.DataFrame(rows).to_csv(harness.FAILURE_STUDY_DIR / "phase4_structural_stats_full_corpus.csv", index=False)
    print("wrote structural (parameter-only) qubit-count accounting for the full corpus")


def _transpile_one(N: int, a: int) -> dict:
    """Time gate construction (UnitaryGate(check_input=True) validation,
    O(dim^3)) plus transpile() together, as one combined wall-clock cost:
    this is the total one-time preparation cost a user actually pays before
    a single window circuit can run, whichever sub-step dominates."""
    import time

    from qiskit import transpile
    from qiskit_aer import AerSimulator

    from Circuits.windowed_qpe import WindowSpec

    spec = WindowSpec(start=0, width=4, total_precision=4)
    backend = AerSimulator()
    t0 = time.perf_counter()
    op = ModularMultiplicationOperator(a, N)
    circuit = build_windowed_qpe_circuit(op, eigenstate=1, spec=spec)
    transpile(circuit, backend)
    return {"N": N, "work_qubits": op.qubits, "total_qubits": circuit.num_qubits, "transpile_seconds": time.perf_counter() - t0}


def transpile_cliff() -> None:
    """Measure combined gate-construction + transpile() wall-clock time (no
    shots, no stitching) as work-register size grows, with a hard
    per-attempt timeout so a runaway synthesis/validation pass can never
    hang the campaign. This directly tests whether the paper's phase-side
    qubit-count decoupling claim survives contact with the arithmetic
    layer's own preparation cost.
    """
    corpus = harness.build_corpus()
    seen_work_qubits: set[int] = set()
    probe_pairs = []
    for pair in sorted(corpus, key=lambda r: r["work_qubits"]):
        if pair["work_qubits"] not in seen_work_qubits:
            seen_work_qubits.add(pair["work_qubits"])
            probe_pairs.append(pair)

    rows = []
    for pair in probe_pairs:
        with ProcessPoolExecutor(max_workers=1) as pool:
            future = pool.submit(_transpile_one, pair["N"], pair["a"])
            try:
                result = future.result(timeout=TRANSPILE_TIMEOUT_S)
                result["timed_out"] = False
            except FutureTimeoutError:
                result = {"N": pair["N"], "work_qubits": pair["work_qubits"], "total_qubits": pair["work_qubits"] + 5, "transpile_seconds": TRANSPILE_TIMEOUT_S, "timed_out": True}
                pool.shutdown(wait=False, cancel_futures=True)
        rows.append(result)
        print(result, flush=True)
        if result["timed_out"]:
            print(f"transpile_cliff: hit the {TRANSPILE_TIMEOUT_S}s timeout at work_qubits={pair['work_qubits']}; stopping the probe here")
            break
    import pandas as pd

    pd.DataFrame(rows).to_csv(harness.FAILURE_STUDY_DIR / "phase4_transpile_cliff.csv", index=False)
    print("wrote transpile-cliff probe")


def main() -> None:
    structural_stats_full_corpus()
    transpile_cliff()
    scaling_by_N()
    scaling_by_precision()
    scaling_by_window_size()


if __name__ == "__main__":
    main()
