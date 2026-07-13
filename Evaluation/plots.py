"""
Purpose
    Generate publication-oriented plots from experiment logs.
Theory
    Windowed reconstruction quality is visualized as success rate and runtime
    across geometry, precision, and candidate-budget sweeps.
Inputs
    Pandas DataFrames produced by Evaluation.experiments.
Outputs
    Matplotlib PNG figures under Data/.
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib.pyplot as plt
import pandas as pd


def plot_success_by_window(frame: pd.DataFrame, output_dir: Path) -> Path:
    """Plot mean success rate versus window size."""
    output_dir.mkdir(parents=True, exist_ok=True)
    grouped = frame.groupby("window_size", as_index=False)["success"].mean()
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(grouped["window_size"], grouped["success"], marker="o")
    ax.set_xlabel("Window size")
    ax.set_ylabel("Success rate")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.3)
    path = output_dir / "success_by_window.png"
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def plot_runtime_by_precision(frame: pd.DataFrame, output_dir: Path) -> Path:
    """Plot mean runtime versus phase precision."""
    output_dir.mkdir(parents=True, exist_ok=True)
    grouped = frame.groupby("precision", as_index=False)["runtime_seconds"].mean()
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(grouped["precision"].astype(str), grouped["runtime_seconds"])
    ax.set_xlabel("Phase precision")
    ax.set_ylabel("Runtime (s)")
    ax.grid(True, axis="y", alpha=0.3)
    path = output_dir / "runtime_by_precision.png"
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def generate_standard_plots(frame: pd.DataFrame, output_dir: Path) -> list[Path]:
    """Generate the standard experiment plot set."""
    if frame.empty:
        raise ValueError("cannot plot an empty DataFrame.")
    return [plot_success_by_window(frame, output_dir), plot_runtime_by_precision(frame, output_dir)]


def self_test() -> None:
    """Verify plotting on a minimal DataFrame."""
    frame = pd.DataFrame({"window_size": [3, 4], "success": [True, False], "precision": [6, 8], "runtime_seconds": [0.1, 0.2]})
    paths = generate_standard_plots(frame, Path("Data"))
    assert all(path.exists() for path in paths)


if __name__ == "__main__":
    self_test()
    print("plots self-test passed")
