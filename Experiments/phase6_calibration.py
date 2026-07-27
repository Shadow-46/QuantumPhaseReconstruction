"""
Purpose
    Calibration check for the confidence quantity C_w against ground-truth
    window correctness (Data/failure_study/FORMULATION.md Section 6).
Theory
    C_w estimates Pr[p_w(b1) > p_w(b2)] from a finite sample. The correct
    ground truth to check it against is therefore whether the *true*
    (noiseless, infinite-shot) probability of the observed top-1 pattern
    really exceeds that of the observed top-2 pattern -- not an externally
    chosen reference phase, since the standard Shor |1> initial state is an
    equal superposition over all order eigenphases, so there is no single
    "true bit pattern" to bit-slice per window. The exact per-window marginal
    distribution is computed once per (N, a, window spec) via statevector
    (Circuits.windowed_qpe.build_windowed_qpe_circuit +
    _phase_register_probabilities), independent of shots/seed, and cached.
Inputs
    Data/failure_study/raw_windows/*.json (from recapture_window_counts.py).
Outputs
    Data/failure_study/phase6_calibration.csv (per-window C_w + correctness)
    and Data/failure_study/plots/calibration_reliability.png.
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

from Circuits.modular_multiplication import ModularMultiplicationOperator
from Circuits.windowed_qpe import WindowSpec, build_windowed_qpe_circuit, _phase_register_probabilities
from Experiments import harness
from Reconstruction.confidence import top1_vs_top2_confidence

RAW_WINDOWS_DIR = harness.FAILURE_STUDY_DIR / "raw_windows"
_TRUE_DIST_CACHE: dict[tuple, dict[str, float]] = {}


def true_window_distribution(a: int, N: int, spec: WindowSpec) -> dict[str, float]:
    """Return the exact (noiseless, infinite-shot) marginal probability
    distribution over one window's bit patterns, cached by (a, N, spec)."""
    key = (a, N, spec.start, spec.width, spec.total_precision)
    if key not in _TRUE_DIST_CACHE:
        op = ModularMultiplicationOperator(a, N)
        circuit = build_windowed_qpe_circuit(op, eigenstate=1, spec=spec)
        _TRUE_DIST_CACHE[key] = _phase_register_probabilities(circuit, spec.width)
    return _TRUE_DIST_CACHE[key]


def per_window_rows(record: dict) -> list[dict]:
    """Return one row per window: observed C_w and whether the observed
    top-1 pattern's true probability really exceeds the top-2's."""
    spec_dict = record["spec"]
    a, N = spec_dict["a"], spec_dict["N"]
    rows = []
    for window_meta, counts in zip(record["windows"], record["window_counts"]):
        spec = WindowSpec(window_meta["start"], window_meta["width"], window_meta["total_precision"])
        ranked = sorted(counts.items(), key=lambda item: item[1], reverse=True)
        if len(ranked) < 2:
            continue
        (b1, n1), (b2, n2) = ranked[0], ranked[1]
        s_w = sum(counts.values())
        c_w = top1_vs_top2_confidence(n1, n2, s_w)
        true_dist = true_window_distribution(a, N, spec)
        correct = true_dist.get(b1, 0.0) > true_dist.get(b2, 0.0)
        rows.append({
            "trial_id": record["trial_id"],
            "phase": record["phase"],
            "window_start": spec.start,
            "window_width": spec.width,
            "n1": n1,
            "n2": n2,
            "s_w": s_w,
            "c_w": c_w,
            "correct": correct,
        })
    return rows


def run_calibration() -> pd.DataFrame:
    records = [json.loads(p.read_text()) for p in RAW_WINDOWS_DIR.glob("*.json")]
    rows = [row for record in records for row in per_window_rows(record)]
    frame = pd.DataFrame(rows)
    frame.to_csv(harness.FAILURE_STUDY_DIR / "phase6_calibration.csv", index=False)
    return frame


def reliability_diagram(frame: pd.DataFrame, n_bins: int = 10) -> tuple[pd.DataFrame, float]:
    """Return (per-bin mean predicted vs. observed correctness, Brier score)."""
    import matplotlib.pyplot as plt
    import numpy as np

    brier = float(((frame["c_w"] - frame["correct"].astype(float)) ** 2).mean())

    bins = pd.cut(frame["c_w"], bins=np.linspace(0.0, 1.0, n_bins + 1), include_lowest=True)
    agg = frame.groupby(bins, observed=True).agg(
        mean_predicted=("c_w", "mean"),
        observed_rate=("correct", "mean"),
        count=("correct", "size"),
    ).dropna()

    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax.plot([0, 1], [0, 1], "--", color="gray", label="perfect calibration")
    ax.plot(agg["mean_predicted"], agg["observed_rate"], "o-", color="#2b6cb0", label="observed")
    ax.set_xlabel("predicted $C_w$")
    ax.set_ylabel("observed correctness rate")
    ax.set_title(f"C_w reliability diagram (Brier score = {brier:.4f})")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(harness.PLOTS_DIR / "calibration_reliability.pdf")
    plt.close(fig)

    return agg, brier


def fit_isotonic_recalibration(frame: pd.DataFrame, n_bins: int = 10):
    """Fit an isotonic regression mapping raw C_w to a recalibrated
    confidence, built only because the reliability diagram showed systematic
    overconfidence (FORMULATION.md Section 6's stated contingency, not
    built preemptively). Returns (fitted model, recalibrated Brier score,
    recalibrated per-bin reliability table); also plots the recalibrated
    diagram alongside the raw one."""
    import matplotlib.pyplot as plt
    import numpy as np
    from sklearn.isotonic import IsotonicRegression

    model = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
    calibrated = model.fit_transform(frame["c_w"], frame["correct"].astype(float))
    brier = float(((calibrated - frame["correct"].astype(float)) ** 2).mean())

    bins = pd.cut(frame["c_w"], bins=np.linspace(0.0, 1.0, n_bins + 1), include_lowest=True)
    agg = frame.assign(calibrated=calibrated).groupby(bins, observed=True).agg(
        mean_calibrated=("calibrated", "mean"),
        observed_rate=("correct", "mean"),
        count=("correct", "size"),
    ).dropna()

    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax.plot([0, 1], [0, 1], "--", color="gray", label="perfect calibration")
    ax.plot(agg["mean_calibrated"], agg["observed_rate"], "o-", color="#2f855a", label="recalibrated")
    ax.set_xlabel("recalibrated $C_w$")
    ax.set_ylabel("observed correctness rate")
    ax.set_title(f"C_w reliability diagram, isotonic-recalibrated (Brier = {brier:.4f})")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(harness.PLOTS_DIR / "calibration_reliability_recalibrated.pdf")
    plt.close(fig)

    return model, brier, agg


def main() -> None:
    frame = run_calibration()
    print(f"{len(frame)} per-window rows from {frame['trial_id'].nunique()} trials")
    agg, brier = reliability_diagram(frame)
    print(f"Raw Brier score: {brier:.4f}")
    print(agg)

    # Systematic miscalibration triggers the stated contingency
    # (FORMULATION.md Section 6): build recalibration only now that it's
    # actually shown to be needed.
    is_miscalibrated = ((agg["mean_predicted"] - agg["observed_rate"]).abs() > 0.1).any()
    if is_miscalibrated:
        print("Reliability diagram shows systematic miscalibration; fitting isotonic recalibration.")
        _model, recal_brier, recal_agg = fit_isotonic_recalibration(frame)
        print(f"Recalibrated Brier score: {recal_brier:.4f}")
        print(recal_agg)


if __name__ == "__main__":
    main()
