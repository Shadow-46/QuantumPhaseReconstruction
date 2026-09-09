"""
Purpose
    Quantify the cost of assumption (A2) -- the working independence of the
    top-2 marginals -- by recomputing the top-2 comparison under the exact
    joint posterior on the frozen calibration corpus, and reporting what
    changes: the induced ranking, discrimination, calibration, and how many
    windows change side of module D's stopping threshold.
Theory
    The shipped statistic C_w models the two marginals as independent
    Binomials. They are not: they are negatively correlated coordinates of one
    K-way multinomial, so the exact posterior is a Dirichlet marginal
    comparison. Reconstruction.confidence.top1_vs_top2_confidence_dirichlet
    evaluates that exact comparison in closed form (derivation in its
    docstring). This script measures the difference on the same 2,687
    observations the calibration study already reports.

    The distinction that matters for the method is ordinal vs. numeric. Module
    D consumes only the ORDERING of windows plus one threshold, so a change
    that leaves the ranking intact does not change what D knows -- it changes
    where the threshold falls. Both are measured separately here.

    This is diagnostic. Nothing in this script is wired into a reconstruction
    path, and phase6_calibration.csv is not modified: the shipped statistic
    remains the one every frozen result in both articles was produced with.
Inputs
    Data/failure_study/phase6_calibration.csv (frozen; read-only here)
Outputs
    Data/failure_study/phase6_calibration_dirichlet.csv
    Data/failure_study/summary_tables/paperA_dirichlet_comparison.csv
    Data/failure_study/summary_tables/paperA_dirichlet_reliability.csv
    Data/failure_study/plots/paperA_figS_dirichlet.png
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr
from sklearn.isotonic import IsotonicRegression

from Experiments import harness
from Reconstruction.confidence import top1_vs_top2_confidence_dirichlet

SUMMARY_DIR = harness.FAILURE_STUDY_DIR / "summary_tables"

# Module D's operating point, and the low-confidence band the calibration
# section reports alongside it. Both come from the shipped configuration.
STOP_THRESHOLD = 0.95
LOW_BAND = 0.60

# Frozen elsewhere in the campaign; asserted so this script cannot silently
# drift from the numbers both articles report.
EXPECTED_CALIBRATION_ROWS = 2687
EXPECTED_RAW_BRIER = 0.2239


def auc_mann_whitney(scores: pd.Series, correct: pd.Series) -> float:
    """AUC via the Mann-Whitney U identity, matching the convention already
    used in paperA_window_validation.cw_discrimination. Average ranks are
    what makes this correct under ties, which matter here: the Dirichlet
    statistic maps every n1 == n2 window to exactly 0.5."""
    ranks = scores.rank()
    n_pos = int(correct.sum())
    n_neg = int((~correct).sum())
    return float((ranks[correct].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def brier(scores: pd.Series, correct: pd.Series) -> float:
    return float(((scores - correct.astype(float)) ** 2).mean())


def isotonic_brier(scores: pd.Series, correct: pd.Series) -> float:
    """In-sample isotonic fit. This bounds what recalibration could achieve on
    this corpus; it does not estimate out-of-sample performance, and it is not
    wired into any decision rule. Computed inline rather than via
    phase6_calibration.fit_isotonic_recalibration because that function writes
    a plot file belonging to the shipped statistic."""
    model = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
    fitted = model.fit_transform(scores, correct.astype(float))
    return float(((fitted - correct.astype(float)) ** 2).mean())


def reliability_bins(frame: pd.DataFrame, column: str, label: str, n_bins: int = 10) -> pd.DataFrame:
    bins = pd.cut(frame[column], bins=np.linspace(0.0, 1.0, n_bins + 1), include_lowest=True)
    agg = (
        frame.groupby(bins, observed=True)
        .agg(
            mean_predicted=(column, "mean"),
            observed_rate=("correct", "mean"),
            count=("correct", "size"),
        )
        .dropna()
    )
    # Keep the interval edges: empty bins are dropped, so a positional index
    # would not say which decile a row describes.
    edges = agg.index.categories[agg.index.codes] if hasattr(agg.index, "codes") else agg.index
    agg = agg.reset_index(drop=True)
    agg.insert(0, "bin_right", [float(iv.right) for iv in edges])
    agg.insert(0, "bin_left", [float(iv.left) for iv in edges])
    agg.insert(0, "statistic", label)
    return agg


def build_frame() -> pd.DataFrame:
    frame = pd.read_csv(harness.FAILURE_STUDY_DIR / "phase6_calibration.csv")
    assert len(frame) == EXPECTED_CALIBRATION_ROWS, (
        f"calibration corpus changed: {len(frame)} rows, expected {EXPECTED_CALIBRATION_ROWS}"
    )
    frame["correct"] = frame["correct"].astype(bool)
    raw_brier = brier(frame["c_w"], frame["correct"])
    assert round(raw_brier, 4) == EXPECTED_RAW_BRIER, (
        f"raw Brier changed: {raw_brier:.4f}, expected {EXPECTED_RAW_BRIER}"
    )

    frame["k_patterns"] = 2 ** frame["window_width"]
    frame["residual"] = frame["s_w"] - frame["n1"] - frame["n2"]
    frame["c_w_indep"] = frame["c_w"]
    frame["c_w_dirichlet"] = [
        top1_vs_top2_confidence_dirichlet(int(n1), int(n2))
        for n1, n2 in zip(frame["n1"], frame["n2"])
    ]
    frame["delta"] = frame["c_w_indep"] - frame["c_w_dirichlet"]
    above_indep = frame["c_w_indep"] >= STOP_THRESHOLD
    above_dir = frame["c_w_dirichlet"] >= STOP_THRESHOLD
    frame["crosses_below_threshold"] = above_indep & ~above_dir
    frame["crosses_above_threshold"] = ~above_indep & above_dir

    columns = [
        "trial_id", "phase", "window_start", "window_width", "k_patterns",
        "n1", "n2", "s_w", "residual",
        "c_w_indep", "c_w_dirichlet", "delta",
        "crosses_below_threshold", "crosses_above_threshold", "correct",
    ]
    out = frame[columns]
    out.to_csv(harness.FAILURE_STUDY_DIR / "phase6_calibration_dirichlet.csv", index=False)
    return frame


def comparison_table(frame: pd.DataFrame) -> pd.DataFrame:
    correct = frame["correct"]
    rows = []
    for label, column in (("independent Beta (shipped)", "c_w_indep"), ("Dirichlet joint", "c_w_dirichlet")):
        scores = frame[column]
        above = scores >= STOP_THRESHOLD
        below = scores < LOW_BAND
        rows.append({
            "statistic": label,
            "n": len(frame),
            "auc": auc_mann_whitney(scores, correct),
            "brier": brier(scores, correct),
            "brier_isotonic": isotonic_brier(scores, correct),
            "n_above_threshold": int(above.sum()),
            "acc_above_threshold": float(correct[above].mean()) if above.any() else float("nan"),
            # The quantity epsilon actually names: the nominal stopping rate a
            # practitioner sets, minus what halting windows deliver. This is
            # what makes epsilon a tuning knob rather than an error target.
            "gap_vs_nominal": (
                float(STOP_THRESHOLD - correct[above].mean()) if above.any() else float("nan")
            ),
            # The same shortfall measured against the mean confidence actually
            # asserted in that set, which exceeds the threshold.
            "gap_vs_asserted": (
                float(scores[above].mean() - correct[above].mean()) if above.any() else float("nan")
            ),
            "n_below_band": int(below.sum()),
            "acc_below_band": float(correct[below].mean()) if below.any() else float("nan"),
        })

    summary = pd.DataFrame(rows)
    # Rank agreement and threshold movement are properties of the pair, so
    # they are attached to the Dirichlet row rather than duplicated.
    summary["spearman_vs_indep"] = [
        float("nan"), float(spearmanr(frame["c_w_indep"], frame["c_w_dirichlet"]).statistic)
    ]
    summary["kendall_vs_indep"] = [
        float("nan"), float(kendalltau(frame["c_w_indep"], frame["c_w_dirichlet"]).statistic)
    ]
    summary["mean_delta"] = [float("nan"), float(frame["delta"].mean())]
    summary["max_delta"] = [float("nan"), float(frame["delta"].max())]
    summary["min_delta"] = [float("nan"), float(frame["delta"].min())]
    summary["n_cross_down"] = [float("nan"), int(frame["crosses_below_threshold"].sum())]
    summary["n_cross_up"] = [float("nan"), int(frame["crosses_above_threshold"].sum())]

    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    summary.to_csv(SUMMARY_DIR / "paperA_dirichlet_comparison.csv", index=False)
    return summary


def reliability_table(frame: pd.DataFrame) -> pd.DataFrame:
    table = pd.concat([
        reliability_bins(frame, "c_w_indep", "independent Beta (shipped)"),
        reliability_bins(frame, "c_w_dirichlet", "Dirichlet joint"),
    ], ignore_index=True)
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    table.to_csv(SUMMARY_DIR / "paperA_dirichlet_reliability.csv", index=False)
    return table


def figure(frame: pd.DataFrame, reliability: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))

    ax = axes[0]
    scatter = ax.scatter(
        frame["c_w_indep"], frame["c_w_dirichlet"],
        c=np.log2(frame["s_w"]), cmap="viridis", s=6, alpha=0.45, linewidths=0,
    )
    ax.plot([0, 1], [0, 1], "--", color="gray", lw=1, label="agreement")
    ax.axvline(STOP_THRESHOLD, color="#805ad5", ls="-.", lw=1)
    ax.axhline(STOP_THRESHOLD, color="#805ad5", ls="-.", lw=1)
    ax.set_xlabel("shipped $C_w$ (independent Beta)")
    ax.set_ylabel("exact $C_w$ (Dirichlet joint)")
    ax.set_title("(a) Where the independence assumption costs")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.3)
    bar = fig.colorbar(scatter, ax=ax)
    bar.set_label("$\\log_2 S_w$", fontsize=8)

    ax = axes[1]
    ax.plot([0, 1], [0, 1], "--", color="gray", lw=1, label="perfect calibration")
    for label, colour, marker in (
        ("independent Beta (shipped)", "#c05621", "o"),
        ("Dirichlet joint", "#2b6cb0", "s"),
    ):
        sub = reliability[reliability["statistic"] == label]
        ax.plot(sub["mean_predicted"], sub["observed_rate"], marker=marker, ls="-",
                color=colour, label=label, ms=5)
    ax.axvline(STOP_THRESHOLD, color="#805ad5", ls="-.", lw=1)
    ax.annotate("module D stops here", (STOP_THRESHOLD, 0.52), rotation=90,
                fontsize=7, color="#805ad5", ha="right")
    ax.set_xlabel("predicted confidence")
    ax.set_ylabel("observed fraction correct")
    ax.set_title("(b) Reliability, before and after the correction")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    path = harness.PLOTS_DIR / "paperA_figS_dirichlet.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path


def self_test() -> None:
    """Verify the metric helpers on constructed inputs, independent of the
    corpus, so a failure here is not confused with a data change."""
    perfect = pd.Series([0.9, 0.8, 0.2, 0.1])
    labels = pd.Series([True, True, False, False])
    assert abs(auc_mann_whitney(perfect, labels) - 1.0) < 1e-12
    inverted = pd.Series([0.1, 0.2, 0.8, 0.9])
    assert abs(auc_mann_whitney(inverted, labels) - 0.0) < 1e-12
    # All-ties must give a coin-flip AUC, the case average ranks exist for.
    assert abs(auc_mann_whitney(pd.Series([0.5] * 4), labels) - 0.5) < 1e-12
    # A constant 0.5 forecast is the Brier reference the article quotes.
    assert abs(brier(pd.Series([0.5] * 4), labels) - 0.25) < 1e-12
    # Isotonic can only improve an in-sample fit.
    assert isotonic_brier(perfect, labels) <= brier(perfect, labels) + 1e-12


def main() -> None:
    frame = build_frame()
    print(f"{len(frame)} windows from {frame['trial_id'].nunique()} trials")

    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 30)

    summary = comparison_table(frame)
    print("\n" + summary.round(4).to_string(index=False))

    reliability = reliability_table(frame)
    print("\nReliability bins:")
    print(reliability.round(4).to_string(index=False))

    down = int(frame["crosses_below_threshold"].sum())
    up = int(frame["crosses_above_threshold"].sum())
    print(f"\nThreshold movement at {STOP_THRESHOLD}: {down} windows fall below, {up} rise above.")
    print("Module D would resample strictly more often under the exact statistic."
          if up == 0 else "Threshold movement is two-sided.")

    print(f"\nfigure: {figure(frame, reliability)}")


if __name__ == "__main__":
    self_test()
    main()
