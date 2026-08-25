"""
Purpose
    Window-level validation of the confidence-guided reconstruction method,
    independent of end-to-end Shor factorisation outcomes. Produces the
    Method validation evidence for the method article, whose unit of
    analysis is the window rather than the trial.
Theory
    Module B replaces a fixed top-k candidate retention with a rule that
    retains until observed cumulative mass reaches 1 - delta_cov. The
    validation question is whether that rule, stated on *observed* counts,
    actually captures the window's *true* probability mass.

    Ground truth is the exact noiseless per-window marginal distribution
    obtained by statevector simulation, not a single "true bit pattern":
    the standard Shor |1> initial state is an equal superposition over all
    order eigenphases, so no single correct per-window value exists (see
    phase6_calibration.py, which establishes this ground-truth convention
    and which this module reuses via true_window_distribution).

    Two metrics per window, both well posed under that convention:
      true mass retained  -- sum of true marginal probability over the
                             retained candidate set
      true argmax kept    -- whether the highest-true-probability pattern
                             survives retention
    Both are reported for the fixed top-k baseline and for module B, over
    the same observed counts, so the comparison is paired per window.
Inputs
    Data/failure_study/raw_windows/*.json   (from recapture_window_counts.py)
    Data/failure_study/phase6_calibration.csv (for C_w discrimination)
Outputs
    Data/failure_study/paperA_window_validation.csv
    Data/failure_study/summary_tables/paperA_window_validation_summary.csv
    Data/failure_study/plots/paper_figA_window_retention.png
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

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from Circuits.windowed_qpe import WindowSpec
from Experiments import harness
from Experiments.phase6_calibration import true_window_distribution
from Reconstruction.candidate_generation import (
    generate_window_candidates,
    generate_window_candidates_by_coverage,
)

RAW_WINDOWS_DIR = harness.FAILURE_STUDY_DIR / "raw_windows"
SUMMARY_DIR = harness.FAILURE_STUDY_DIR / "summary_tables"
DELTA_COV = 0.05

# Established elsewhere in the campaign; asserted here so this script cannot
# silently drift from the numbers the companion documents report.
EXPECTED_CALIBRATION_ROWS = 2687
EXPECTED_RAW_BRIER = 0.2239


def per_window_rows(record: dict) -> list[dict]:
    """One row per window: true mass and true argmax retained, under the
    fixed top-k baseline and under module B, from the same observed counts."""
    spec_dict = record["spec"]
    a, N = spec_dict["a"], spec_dict["N"]
    candidate_count = spec_dict["candidate_count"]
    rows = []
    for window_meta, counts in zip(record["windows"], record["window_counts"]):
        spec = WindowSpec(window_meta["start"], window_meta["width"], window_meta["total_precision"])
        true_dist = true_window_distribution(a, N, spec)
        if not true_dist:
            continue
        true_argmax = max(true_dist, key=true_dist.get)

        baseline = generate_window_candidates(
            counts, spec.start, spec.width, spec.total_precision, candidate_count
        )
        coverage = generate_window_candidates_by_coverage(
            counts, spec.start, spec.width, spec.total_precision, DELTA_COV
        )
        bits_baseline = {c.local_bits for c in baseline}
        bits_coverage = {c.local_bits for c in coverage}

        rows.append({
            "trial_id": record["trial_id"],
            "phase": record["phase"],
            "N": N,
            "a": a,
            "window_start": spec.start,
            "window_width": spec.width,
            "total_precision": spec.total_precision,
            "candidate_count": candidate_count,
            "shots": sum(counts.values()),
            "patterns_observed": len(counts),
            "m_baseline": len(bits_baseline),
            "m_coverage": len(bits_coverage),
            "true_mass_baseline": sum(true_dist.get(b, 0.0) for b in bits_baseline),
            "true_mass_coverage": sum(true_dist.get(b, 0.0) for b in bits_coverage),
            "argmax_kept_baseline": true_argmax in bits_baseline,
            "argmax_kept_coverage": true_argmax in bits_coverage,
        })
    return rows


def build_frame() -> pd.DataFrame:
    records = [json.loads(p.read_text()) for p in sorted(RAW_WINDOWS_DIR.glob("*.json"))]
    rows = [row for record in records for row in per_window_rows(record)]
    frame = pd.DataFrame(rows)
    frame.to_csv(harness.FAILURE_STUDY_DIR / "paperA_window_validation.csv", index=False)
    return frame


def cw_discrimination() -> dict[str, float]:
    """Ranking quality of C_w on the already-frozen calibration corpus.

    Discrimination is what module D actually relies on -- only the ordering
    of windows is load-bearing -- whereas the Brier score reported alongside
    it measures the calibration that module D does not rely on."""
    frame = pd.read_csv(harness.FAILURE_STUDY_DIR / "phase6_calibration.csv")
    assert len(frame) == EXPECTED_CALIBRATION_ROWS, (
        f"calibration corpus changed: {len(frame)} rows, expected {EXPECTED_CALIBRATION_ROWS}"
    )
    correct = frame["correct"].astype(bool)
    brier = float(((frame["c_w"] - correct.astype(float)) ** 2).mean())
    assert round(brier, 4) == EXPECTED_RAW_BRIER, (
        f"raw Brier changed: {brier:.4f}, expected {EXPECTED_RAW_BRIER}"
    )

    # AUC via the Mann-Whitney U identity, so this adds no new dependency.
    ranks = frame["c_w"].rank()
    n_pos = int(correct.sum())
    n_neg = int((~correct).sum())
    auc = float((ranks[correct].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))
    return {"n": len(frame), "n_correct": n_pos, "auc": auc, "brier": brier}


def summarise(frame: pd.DataFrame) -> pd.DataFrame:
    def block(label: str, sub: pd.DataFrame) -> dict:
        return {
            "group": label,
            "windows": len(sub),
            "true_mass_baseline": sub["true_mass_baseline"].mean(),
            "true_mass_coverage": sub["true_mass_coverage"].mean(),
            "mass_gain": sub["true_mass_coverage"].mean() - sub["true_mass_baseline"].mean(),
            "argmax_kept_baseline": sub["argmax_kept_baseline"].mean(),
            "argmax_kept_coverage": sub["argmax_kept_coverage"].mean(),
            "m_baseline": sub["m_baseline"].mean(),
            "m_coverage": sub["m_coverage"].mean(),
        }

    rows = [block("all windows", frame)]
    for k, sub in frame.groupby("candidate_count"):
        rows.append(block(f"candidate_count = {k}", sub))
    for (N, a), sub in frame.groupby(["N", "a"]):
        rows.append(block(f"instance N={N}, a={a}", sub))
    for phase, sub in frame.groupby("phase"):
        rows.append(block(f"corpus {phase}", sub))

    summary = pd.DataFrame(rows)
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    summary.to_csv(SUMMARY_DIR / "paperA_window_validation_summary.csv", index=False)
    return summary


def figure(frame: pd.DataFrame) -> Path:
    by_k = frame.groupby("candidate_count").agg(
        mass_baseline=("true_mass_baseline", "mean"),
        mass_coverage=("true_mass_coverage", "mean"),
        m_baseline=("m_baseline", "mean"),
        m_coverage=("m_coverage", "mean"),
    )
    ks = by_k.index.to_list()

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    ax = axes[0]
    ax.plot(ks, by_k["mass_baseline"], "o-", color="#c05621", label="fixed top-$k$ (baseline)")
    ax.plot(ks, by_k["mass_coverage"], "s-", color="#2b6cb0", label="coverage rule (module B)")
    ax.axhline(1 - DELTA_COV, ls=":", color="gray", lw=1)
    ax.annotate(f"target mass $1-\\delta_{{cov}} = {1 - DELTA_COV}$",
                (ks[0], 1 - DELTA_COV), textcoords="offset points",
                xytext=(4, 5), fontsize=8, color="gray")
    ax.set_xscale("log", base=2)
    ax.set_xticks(ks)
    ax.set_xticklabels(ks)
    ax.set_xlabel("configured candidate count $k$")
    ax.set_ylabel("mean true probability mass retained")
    ax.set_title("(a) What the retention rule actually captures")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.plot(ks, by_k["m_baseline"], "o-", color="#c05621", label="fixed top-$k$ (baseline)")
    ax.plot(ks, by_k["m_coverage"], "s-", color="#2b6cb0", label="coverage rule (module B)")
    ax.set_xscale("log", base=2)
    ax.set_xticks(ks)
    ax.set_xticklabels(ks)
    ax.set_xlabel("configured candidate count $k$")
    ax.set_ylabel("mean candidates retained $m_w$")
    ax.set_title("(b) What it costs")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    path = harness.PLOTS_DIR / "paper_figA_window_retention.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path


def main() -> None:
    frame = build_frame()
    print(f"{len(frame)} windows from {frame['trial_id'].nunique()} trials")

    summary = summarise(frame)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 20)
    print("\n" + summary.round(4).to_string(index=False))

    disc = cw_discrimination()
    print(f"\nC_w discrimination on {disc['n']} calibration windows "
          f"({disc['n_correct']} correct): AUC = {disc['auc']:.4f}, "
          f"raw Brier = {disc['brier']:.4f}")

    print(f"\nfigure: {figure(frame)}")


if __name__ == "__main__":
    main()
