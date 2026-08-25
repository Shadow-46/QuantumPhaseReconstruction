"""
Purpose
    Close W1/W3 (single-instance validation, circularity) with three
    pre-registered held-out (N,a,r) instances run on independent grids, and
    close W2 (no matched-resource control) by instrumenting
    reconstruct_adaptive's realized resource usage (m_w, P', s_w) on the
    canonical grid and comparing each adaptive module against a
    matched-budget uniform control that spends the same realized resources
    without confidence guidance.
Theory
    Under the work_qubits<=5 ceiling (Circuits/modular_multiplication.py
    work_qubits_for_modulus), only N in {15, 21} admit an odd two-prime-
    factor modulus with an even, non-trivial order (required by
    Reconstruction/continued_fraction.py factors_from_order). r=6 at N=21
    (the canonical instance) is already the largest reachable even order,
    so held-out instances broaden coverage but cannot be "harder" than
    canonical -- this closes W3 (circularity) but only narrows W1
    (external validity); see the plan handoff for the full argument.
    Held-out instances were chosen by a fixed pre-registered rule (lowest
    valid base achieving the target order, excluding a=19) before any
    trial in this script ran:
        (N=15, a=2,  r=4)  -- new order, new modulus
        (N=21, a=8,  r=2)  -- new order, canonical modulus
        (N=21, a=2,  r=6)  -- canonical order, new base
    Matched-resource arms reuse reconstruct_adaptive's existing Baseline
    code path (all three flags False) fed post-hoc uniform parameters
    derived from an adaptive arm's own realized diagnostics -- no new
    reconstruction module, only new orchestration:
        shots-matched (vs D, vs Full): total extra shots the adaptive arm
            drew, redistributed uniformly across windows.
        candidates-matched (vs B): fixed top-k = round(mean(m_w)) over
            that trial's windows (an approximation to B's non-uniform
            per-window allocation, stated as such).
        beam-matched (vs C): fixed max_paths = C's realized P' applied
            unconditionally (exact, since P' is already trial-level).
Inputs
    None -- instances and grid hardcoded below, reusing phase6_ablation's
    GRID/ARMS/WINDOW_SIZE/seeding exactly so held-out grids are structurally
    identical to the canonical one.
Outputs
    Data/failure_study/phase6_ablation_diagnostics.csv (canonical grid
      re-run with diagnostics; success_* must equal phase6_ablation.csv
      exactly -- checked by verify_diagnostics_match_canonical()).
    Data/failure_study/phase7_heldout_N{N}_a{A}.csv per held-out instance.
    Data/failure_study/phase7_matched_resource.csv.
    Data/failure_study/summary_tables/phase7_heldout_summary.csv and
      phase7_matched_resource_summary.csv.
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

import pandas as pd

from Experiments import harness
from Experiments.phase6_ablation import ARMS, REPEATS, WINDOW_SIZE, all_combos, mcnemar_test
from config import DEFAULT_ADAPTIVE, ShorConfig, WindowConfig

# Pre-registered before any trial in this script ran (see module docstring).
HELD_OUT_INSTANCES = [
    (15, 2, 4),
    (21, 8, 2),
    (21, 2, 6),
]


# --------------------------------------------------------------------------
# Instrumented grid runner (shared by the canonical re-run and the
# held-out instances)
# --------------------------------------------------------------------------

def run_one_trial_with_diagnostics(N: int, A: int, combo: tuple, repeat: int) -> dict:
    """Sample the base window counts once per trial and share them across
    all five arms: given the shared seed, every arm would recompute
    bit-identical base counts anyway (see reconstruct_adaptive's
    `precomputed` docstring), so this only removes 5x redundant Aer
    transpile+execute -- it cannot change any success_<arm> outcome."""
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
    base_sample = harness.sample_window_counts(shor, window)
    for arm_name, flags in ARMS.items():
        _stitched, recovered, diag = harness.reconstruct_adaptive(
            shor, window,
            delta_cov=DEFAULT_ADAPTIVE.delta_cov, beam_max=DEFAULT_ADAPTIVE.beam_max,
            epsilon=DEFAULT_ADAPTIVE.epsilon, s_max=DEFAULT_ADAPTIVE.s_max,
            return_diagnostics=True,
            precomputed=base_sample,
            **flags,
        )
        row[f"success_{arm_name}"] = recovered is not None and recovered.factors is not None
        row[f"diagnostics_{arm_name}"] = json.dumps(diag)
    return row


def run_instance_grid(N: int, A: int, csv_path: Path, resume: bool = True) -> pd.DataFrame:
    """Run the full 72-combo x 2-repeat x 5-arm grid for one (N,A) instance,
    with per-arm diagnostics attached to every row. Structurally identical
    to phase6_ablation.run_ablation_grid except for the diagnostics column
    and the output path, so held-out grids and the canonical re-run share
    exactly one code path."""
    done_seeds: set[int] = set()
    rows: list[dict] = []
    if resume and csv_path.exists():
        existing = pd.read_csv(csv_path)
        rows = existing.to_dict("records")
        done_seeds = set(existing["seed"])

    work_items = [(combo, repeat) for combo in all_combos() for repeat in range(REPEATS)]
    total = len(work_items)
    print(f"  N={N} a={A}: {len(done_seeds)}/{total} already done, resuming", flush=True)
    for i, (combo, repeat) in enumerate(work_items):
        seed = hash((*combo, repeat)) % (2**31)
        if seed in done_seeds:
            continue
        rows.append(run_one_trial_with_diagnostics(N, A, combo, repeat))
        pd.DataFrame(rows).to_csv(csv_path, index=False)
        print(f"  N={N} a={A} progress: {i + 1}/{total}", flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(csv_path, index=False)
    return frame


def verify_diagnostics_match_canonical(instrumented: pd.DataFrame) -> None:
    """Gate: the instrumented canonical (21,19) re-run must reproduce the
    already-published phase6_ablation.csv success_* columns bit-for-bit.
    If this fails, stop -- it means adding diagnostics perturbed the
    result the published Table 2 depends on, or the seed-paired
    determinism the whole paper relies on does not actually hold."""
    published_path = harness.FAILURE_STUDY_DIR / "phase6_ablation.csv"
    published = pd.read_csv(published_path)
    success_cols = [c for c in published.columns if c.startswith("success_")]

    pub = published.set_index("seed")[success_cols].sort_index()
    new = instrumented.set_index("seed")[success_cols].sort_index()
    if not pub.index.equals(new.index):
        raise RuntimeError(
            f"seed sets differ between published ({len(pub)}) and instrumented re-run ({len(new)}); "
            "cannot verify bit-identical outcomes."
        )
    mismatches = (pub != new)
    if mismatches.any().any():
        bad = mismatches[mismatches.any(axis=1)]
        raise RuntimeError(
            f"instrumented re-run diverges from published phase6_ablation.csv on {len(bad)} trial(s): "
            f"seeds {list(bad.index)[:10]}. Stopping -- see plan Part 6 step 2 gate."
        )
    print(f"  verified: instrumented re-run matches published phase6_ablation.csv on all {len(pub)} trials.")


# --------------------------------------------------------------------------
# Held-out instance summaries
# --------------------------------------------------------------------------

def summarize_instance(frame: pd.DataFrame, label: str) -> pd.DataFrame:
    rows = []
    for arm_name in ARMS:
        if arm_name == "baseline":
            continue
        stats = mcnemar_test(frame["success_baseline"], frame[f"success_{arm_name}"])
        rows.append({"instance": label, "arm": arm_name, **stats})
    b_vs_c = mcnemar_test(frame["success_B"], frame["success_C"])
    rows.append({"instance": label, "arm": "B_vs_C", **b_vs_c})
    return pd.DataFrame(rows)


def run_held_out_instances() -> pd.DataFrame:
    summaries = []
    for N, A, order in HELD_OUT_INSTANCES:
        label = f"N={N},a={A},r={order}"
        csv_path = harness.FAILURE_STUDY_DIR / f"phase7_heldout_N{N}_a{A}.csv"
        print(f"held-out instance {label} -> {csv_path.name}")
        frame = run_instance_grid(N, A, csv_path)
        summaries.append(summarize_instance(frame, label))
    combined = pd.concat(summaries, ignore_index=True)
    combined.to_csv(harness.FAILURE_STUDY_DIR / "summary_tables" / "phase7_heldout_summary.csv", index=False)
    return combined


# --------------------------------------------------------------------------
# Matched-resource control arms (built from the canonical grid's
# diagnostics; W2)
# --------------------------------------------------------------------------

def _diag(row: dict, arm: str) -> dict:
    return json.loads(row[f"diagnostics_{arm}"])


def run_matched_resource_trial(row: dict) -> dict:
    """For one already-instrumented canonical-grid trial row, build the
    three matched-resource controls and evaluate each with
    reconstruct_adaptive's Baseline path (all flags False) at the
    post-hoc-derived uniform parameters. Reuses the trial's own seed, so
    each matched arm stays paired with the same trial's Baseline/D/B/C/Full
    outcomes already in the row."""
    N, A = 21, 19
    combo = (row["candidate_count"], row["max_paths"], row["shots"], row["total_precision"], row["overlap"])
    seed = row["seed"]
    base_window = WindowConfig(
        total_precision=row["total_precision"], window_size=WINDOW_SIZE, overlap=row["overlap"],
        candidate_count=row["candidate_count"], max_paths=row["max_paths"],
    )
    kwargs = dict(delta_cov=DEFAULT_ADAPTIVE.delta_cov, beam_max=DEFAULT_ADAPTIVE.beam_max,
                  epsilon=DEFAULT_ADAPTIVE.epsilon, s_max=DEFAULT_ADAPTIVE.s_max)

    out = {"candidate_count": combo[0], "max_paths": combo[1], "shots": combo[2],
           "total_precision": combo[3], "overlap": combo[4], "repeat": row["repeat"], "seed": seed,
           "success_baseline": row["success_baseline"], "success_B": row["success_B"],
           "success_C": row["success_C"], "success_D": row["success_D"], "success_full": row["success_full"]}

    num_windows = len(_diag(row, "baseline")["s_w"])
    baseline_shots_total = row["shots"] * num_windows

    # Shots-matched vs D and vs Full.
    for arm, out_key in (("D", "shots_matched_vs_D"), ("full", "shots_matched_vs_full")):
        s_w = _diag(row, arm)["s_w"]
        extra_total = sum(s_w) - baseline_shots_total
        extra_per_window = max(0, extra_total) // num_windows
        boosted_shots = row["shots"] + extra_per_window
        shor = ShorConfig(N=N, a=A, phase_qubits=row["total_precision"], shots=boosted_shots, random_seed=seed)
        _s, recovered = harness.reconstruct_adaptive(shor, base_window, **kwargs)
        out[out_key] = recovered is not None and recovered.factors is not None
        out[f"{out_key}_extra_per_window"] = int(extra_per_window)

    # Candidates-matched vs B.
    m_w = _diag(row, "B")["m_w"]
    k_matched = max(1, round(sum(m_w) / len(m_w)))
    window_k = WindowConfig(total_precision=row["total_precision"], window_size=WINDOW_SIZE, overlap=row["overlap"],
                             candidate_count=k_matched, max_paths=row["max_paths"])
    shor_base = ShorConfig(N=N, a=A, phase_qubits=row["total_precision"], shots=row["shots"], random_seed=seed)
    _s, recovered = harness.reconstruct_adaptive(shor_base, window_k, **kwargs)
    out["candidates_matched_vs_B"] = recovered is not None and recovered.factors is not None
    out["candidates_matched_vs_B_k"] = int(k_matched)

    # Beam-matched vs C.
    p_prime = _diag(row, "C")["max_paths"]
    window_p = WindowConfig(total_precision=row["total_precision"], window_size=WINDOW_SIZE, overlap=row["overlap"],
                             candidate_count=row["candidate_count"], max_paths=p_prime)
    _s, recovered = harness.reconstruct_adaptive(shor_base, window_p, **kwargs)
    out["beam_matched_vs_C"] = recovered is not None and recovered.factors is not None
    out["beam_matched_vs_C_p_prime"] = int(p_prime)

    return out


def run_matched_resource_grid(instrumented: pd.DataFrame, csv_path: Path, resume: bool = True) -> pd.DataFrame:
    done_seeds: set[int] = set()
    rows: list[dict] = []
    if resume and csv_path.exists():
        existing = pd.read_csv(csv_path)
        rows = existing.to_dict("records")
        done_seeds = set(existing["seed"])

    trial_rows = instrumented.to_dict("records")
    total = len(trial_rows)
    print(f"  matched-resource: {len(done_seeds)}/{total} already done, resuming", flush=True)
    for i, row in enumerate(trial_rows):
        if row["seed"] in done_seeds:
            continue
        rows.append(run_matched_resource_trial(row))
        pd.DataFrame(rows).to_csv(csv_path, index=False)
        print(f"  matched-resource progress: {i + 1}/{total}", flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(csv_path, index=False)
    return frame


def summarize_matched_resource(frame: pd.DataFrame) -> pd.DataFrame:
    pairs = [
        ("D", "shots_matched_vs_D"),
        ("full", "shots_matched_vs_full"),
        ("B", "candidates_matched_vs_B"),
        ("C", "beam_matched_vs_C"),
    ]
    rows = []
    for adaptive_arm, matched_col in pairs:
        vs_baseline = mcnemar_test(frame["success_baseline"], frame[matched_col])
        rows.append({"comparison": f"{matched_col}_vs_baseline", **vs_baseline})
        vs_adaptive = mcnemar_test(frame[f"success_{adaptive_arm}"], frame[matched_col])
        rows.append({"comparison": f"{matched_col}_vs_{adaptive_arm}", **vs_adaptive})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------

def main() -> None:
    print("Step 1: instrumented re-run of canonical (21,19) grid")
    canonical_diag_path = harness.FAILURE_STUDY_DIR / "phase6_ablation_diagnostics.csv"
    canonical = run_instance_grid(21, 19, canonical_diag_path)
    verify_diagnostics_match_canonical(canonical)

    print("Step 2: held-out instances")
    heldout_summary = run_held_out_instances()
    print(heldout_summary.to_string(index=False))

    print("Step 3: matched-resource control arms")
    matched_path = harness.FAILURE_STUDY_DIR / "phase7_matched_resource.csv"
    matched = run_matched_resource_grid(canonical, matched_path)
    matched_summary = summarize_matched_resource(matched)
    matched_summary.to_csv(harness.FAILURE_STUDY_DIR / "summary_tables" / "phase7_matched_resource_summary.csv", index=False)
    print(matched_summary.to_string(index=False))


if __name__ == "__main__":
    main()
