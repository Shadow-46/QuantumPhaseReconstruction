"""
Purpose
    Phase 3: black-box adversarial fault injection on the classical
    reconstruction pipeline, plus an exhaustive carry-mechanism truth table
    and a measurement-ambiguity boundary sweep.
Theory
    These experiments corrupt real sampled outputs (plain dict/candidate
    manipulation) before feeding them into the unmodified public
    Reconstruction functions, isolating classical-pipeline robustness from
    quantum-sampling noise (already covered in Phase 2).
Inputs
    Real Aer-sampled window counts (via Algorithms.paper_algorithm) and the
    ground-truth corpus.
Outputs
    Data/failure_study/phase3_*.csv and per-trial JSON under raw/.
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
import random
import time
import tracemalloc
import uuid
from fractions import Fraction

import pandas as pd

from Algorithms.paper_algorithm import sample_window_counts
from Circuits.windowed_qpe import make_overlapping_windows
from Reconstruction.candidate_generation import WindowCandidate, generate_window_candidates
from Reconstruction.carry import candidates_compatible, carry_bit, overlap_width
from Reconstruction.continued_fraction import recover_order_and_factors
from Reconstruction.stitching import stitch_candidates
from Experiments import harness
from config import ShorConfig, WindowConfig

REPEATS = 8
R0 = {"total_precision": 8, "window_size": 4, "overlap": 2, "candidate_count": 2, "max_paths": 32, "shots": 4096}


def _base_counts(pair: dict, seed: int) -> tuple[list[dict[str, int]], list]:
    shor = ShorConfig(N=pair["N"], a=pair["a"], phase_qubits=R0["total_precision"], shots=R0["shots"], random_seed=seed)
    window = WindowConfig(
        total_precision=R0["total_precision"], window_size=R0["window_size"], overlap=R0["overlap"],
        candidate_count=R0["candidate_count"], max_paths=R0["max_paths"],
    )
    counts, specs = sample_window_counts(shor, window)
    return counts, specs


def _flip_bit(bitstring: str, index: int) -> str:
    chars = list(bitstring)
    chars[index] = "1" if chars[index] == "0" else "0"
    return "".join(chars)


def _reconstruct(N: int, a: int, counts: list[dict[str, int]]) -> dict:
    stitched, recovered = harness.reconstruct_from_counts(
        N, a, R0["total_precision"], R0["window_size"], R0["overlap"], R0["candidate_count"], R0["max_paths"], counts
    )
    success = recovered is not None and recovered.factors is not None
    return {
        "num_stitched": len(stitched),
        "success": success,
        "order": recovered.order if recovered else None,
        "factors": list(recovered.factors) if recovered and recovered.factors else None,
    }


def run_corruption_condition(condition: str, pair: dict, seed: int) -> dict:
    started = time.perf_counter()
    tracemalloc.start()
    counts, specs = _base_counts(pair, seed)
    corrupted = [dict(c) for c in counts]
    rng = random.Random(seed)

    if condition == "clean":
        pass
    elif condition == "overlap_corruption":
        idx = rng.randrange(len(corrupted) - 1) if len(corrupted) > 1 else 0
        top_bits = max(corrupted[idx], key=corrupted[idx].get)
        width = specs[idx].width
        flip_pos = rng.randrange(width)
        new_bits = _flip_bit(top_bits, flip_pos)
        corrupted[idx] = {new_bits: corrupted[idx][top_bits], **{k: v for k, v in corrupted[idx].items() if k != top_bits}}
    elif condition == "missing_true_candidate":
        idx = rng.randrange(len(corrupted))
        top_bits = max(corrupted[idx], key=corrupted[idx].get)
        del corrupted[idx][top_bits]
        if not corrupted[idx]:
            corrupted[idx] = {"0" * specs[idx].width: 1}
    elif condition == "incorrect_high_weight_candidate":
        idx = rng.randrange(len(corrupted))
        width = specs[idx].width
        fake_bits = format(rng.randrange(2**width), f"0{width}b")
        max_count = max(corrupted[idx].values())
        corrupted[idx][fake_bits] = max_count * 10
    elif condition == "bit_flip_post_hoc":
        flip_prob = 0.1
        for i, c in enumerate(corrupted):
            new_c: dict[str, int] = {}
            for bits, cnt in c.items():
                for _ in range(cnt):
                    b = bits
                    if rng.random() < flip_prob:
                        b = _flip_bit(b, rng.randrange(len(b)))
                    new_c[b] = new_c.get(b, 0) + 1
            corrupted[i] = new_c
    else:
        raise ValueError(condition)

    outcome = _reconstruct(pair["N"], pair["a"], corrupted)
    runtime = time.perf_counter() - started
    _current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return {
        "trial_id": uuid.uuid4().hex,
        "phase": "phase3_adversarial",
        "condition": condition,
        "N": pair["N"],
        "a": pair["a"],
        "seed": seed,
        "work_qubits": pair["work_qubits"],
        "runtime_seconds": runtime,
        "peak_memory_bytes": peak,
        **outcome,
    }


def run_adversarial_battery() -> pd.DataFrame:
    corpus = harness.build_corpus()
    pairs = harness.representative_subset(corpus, k=8, max_work_qubits=5)
    conditions = ["clean", "overlap_corruption", "missing_true_candidate", "incorrect_high_weight_candidate", "bit_flip_post_hoc"]
    rows = []
    for condition, pair in itertools.product(conditions, pairs):
        for repeat in range(REPEATS):
            seed = hash((condition, pair["N"], pair["a"], repeat)) % (2**31)
            record = run_corruption_condition(condition, pair, seed)
            (harness.RAW_DIR / f"phase3_adversarial_{record['trial_id']}.json").write_text(json.dumps(record, indent=2))
            rows.append(record)
    frame = pd.DataFrame(rows)
    frame.to_csv(harness.FAILURE_STUDY_DIR / "phase3_adversarial.csv", index=False)
    return frame


def run_carry_truth_table() -> pd.DataFrame:
    """Exhaustively enumerate every (tail, head, carry-bit) combination for
    overlap widths 1-3 and audit candidates_compatible against a
    ground-truth 'are these truly consistent slices of one integer, allowing
    at most a single legitimate borrow' oracle.

    `right.start` must be `width - overlap` (not `overlap`) so that
    `overlap_width(left, right)` actually equals the loop's nominal
    `overlap`; the previous `right.start=overlap` collapsed every case to a
    1-bit true overlap regardless of the loop variable, so `tail`/`head`
    were never fully compared (only their outermost bit was), which is what
    the previously reported 55.4% false-accept rate was actually measuring."""
    rows = []
    for overlap in (1, 2, 3):
        width = overlap + 1  # smallest window wider than the overlap itself
        right_start = width - overlap
        for tail_val in range(2**overlap):
            for head_val in range(2**overlap):
                for extra_bit in range(2):  # the bit of `right` just past the head (carry_bit input)
                    tail = format(tail_val, f"0{overlap}b")
                    head = format(head_val, f"0{overlap}b")
                    left = WindowCandidate(value=0, total_precision=20, start=0, width=width, local_bits=("0" * (width - overlap)) + tail, weight=1.0)
                    right = WindowCandidate(value=0, total_precision=20, start=right_start, width=width, local_bits=head + str(extra_bit) + "0" * (width - overlap - 1), weight=1.0)
                    assert overlap_width(left, right) == overlap
                    accepted = candidates_compatible(left, right)
                    # Ground truth: true iff tail == head (no borrow) OR tail-1 == head mod 2^overlap AND extra_bit==1
                    # (a legitimate single-unit borrow is only physically meaningful when the
                    # bit that would have carried into it, i.e. extra_bit, is set).
                    truly_consistent_no_borrow = tail_val == head_val
                    truly_consistent_with_borrow = ((tail_val - 1) % (2**overlap) == head_val) and extra_bit == 1
                    ground_truth = truly_consistent_no_borrow or truly_consistent_with_borrow
                    rows.append(
                        {
                            "overlap": overlap,
                            "tail": tail,
                            "head": head,
                            "carry_bit_input": extra_bit,
                            "accepted": accepted,
                            "ground_truth_consistent": ground_truth,
                            "outcome": (
                                "true_accept" if accepted and ground_truth else
                                "true_reject" if not accepted and not ground_truth else
                                "false_accept" if accepted and not ground_truth else
                                "false_reject"
                            ),
                        }
                    )
    frame = pd.DataFrame(rows)
    frame.to_csv(harness.FAILURE_STUDY_DIR / "phase3_carry_truth_table.csv", index=False)
    return frame


def run_ambiguity_boundary_sweep() -> pd.DataFrame:
    """Sweep how close a single window's true fractional phase sits to the
    0.5 bucket boundary and measure (a) the top1/top2 count ratio the paper
    uses for ambiguity flagging, and (b) how often independent shot-noise
    realizations disagree on the winning bucket -- the AWQPE 'special chunk'
    phenomenon (Remark 2.1). Uses a genuine single-qubit phase eigenstate
    (PhaseGate, eigenstate=1) via Circuits.qpe.build_qpe_circuit so the true
    phase is exactly controllable; a single window at offset 0 is
    numerically identical to this standard-QPE circuit (verified in
    Circuits/windowed_qpe.py's self-test), so this is a faithful probe of
    the windowed block's own IQFT+measurement stage.
    """
    import numpy as np
    from qiskit.circuit.library import PhaseGate

    from Circuits.qpe import build_qpe_circuit
    from Simulation.backend import run_counts

    rows = []
    width = 4
    bucket = 7  # arbitrary interior bucket, away from the 0/2^width wraparound edge
    deltas = [0.25, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0]
    shots = 4096
    repeats = 20
    for delta in deltas:
        true_phase = (bucket + 0.5 - delta) / (2**width)
        gate = PhaseGate(2 * np.pi * true_phase)
        circuit = build_qpe_circuit(gate, width, 1, eigenstate=1, measure=True)
        winners = []
        for repeat in range(repeats):
            seed = hash((delta, repeat)) % (2**31)
            counts = run_counts(circuit, shots=shots, seed=seed)
            total = sum(counts.values())
            ranked = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)
            top1_frac = ranked[0][1] / total
            second_frac = ranked[1][1] / total if len(ranked) > 1 else 0.0
            ambiguous = (second_frac / top1_frac) > 0.9 if top1_frac > 0 else True
            winners.append(ranked[0][0])
            rows.append(
                {
                    "delta_from_boundary": delta,
                    "true_phase": true_phase,
                    "seed": seed,
                    "top1_share": top1_frac,
                    "second_share": second_frac,
                    "ambiguous_by_0.9_threshold": ambiguous,
                    "winner": ranked[0][0],
                }
            )
        winner_flip_rate = len(set(winners)) > 1
        for row in rows[-repeats:]:
            row["winner_disagreement_across_repeats"] = winner_flip_rate
    frame = pd.DataFrame(rows)
    frame.to_csv(harness.FAILURE_STUDY_DIR / "phase3_ambiguity_boundary.csv", index=False)
    return frame


def main() -> None:
    print("running adversarial corruption battery...")
    adv = run_adversarial_battery()
    print(adv.groupby("condition")["success"].mean())

    print("running exhaustive carry truth table...")
    carry = run_carry_truth_table()
    print(carry["outcome"].value_counts())

    print("running measurement-ambiguity boundary sweep...")
    ambig = run_ambiguity_boundary_sweep()
    print(ambig.groupby("delta_from_boundary")[["ambiguous_by_0.9_threshold", "winner_disagreement_across_repeats"]].mean())


if __name__ == "__main__":
    main()
