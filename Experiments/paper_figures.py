"""
Purpose
    Render the four data figures used by Paper/main.tex at publication
    quality. Plotting only: this script executes no trial, samples no
    circuit, and recomputes no experimental outcome. Every value it draws
    is read from a frozen result CSV.
Theory
    The originally generated exploratory plots were written for the
    failure-study campaign, not for a manuscript, and carry the problems
    that implies -- raw identifier strings as tick labels, rotated text,
    no cumulative Pareto line, no annotation of the thresholds the paper
    argues about. These renderings answer one stated scientific question
    each:
      fig2  Which failure mechanism should a fix target?
      fig3  Does the canonical ablation result survive on instances the
            modules were never tuned against?
      fig4  Is the gain attributable to confidence guidance, or to the
            size of the budget the guidance happens to request?
      fig5  Does C_w's numerical value behave like a probability?
    The single isotonic refit in fig5 reproduces the recalibrated Brier
    score of 0.1988 already reported by phase6_calibration.py; it is
    recomputed here only so the recalibrated curve can be drawn on the
    same axes as the raw one.
Inputs
    Data/failure_study/phase5_attribution.csv
    Data/failure_study/phase6_ablation_diagnostics.csv
    Data/failure_study/phase7_heldout_N{15_a2,21_a8,21_a2}.csv
    Data/failure_study/phase7_matched_resource.csv
    Data/failure_study/phase6_calibration.csv
Outputs
    Data/failure_study/plots/paper_fig2_failure_attribution.png
    Data/failure_study/plots/paper_fig3_heldout.png
    Data/failure_study/plots/paper_fig4_matched_resource.png
    Data/failure_study/plots/paper_fig5_calibration.png
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

DATA_DIR = PROJECT_ROOT / "Data" / "failure_study"
PLOT_DIR = DATA_DIR / "plots"

plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 9.5,
    "axes.labelsize": 9,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 300,
})

BUDGET = "#3d6fb4"      # classical-budget starvation
NOISE = "#c98b3a"       # noise
OTHER = "#9aa0a6"       # unattributed
ADAPTIVE = "#3d6fb4"
CONTROL = "#ffffff"


# ---------------------------------------------------------------------
# Figure 2 -- failure attribution
# ---------------------------------------------------------------------

ATTRIBUTION_LABEL = {
    "quantum_estimation_shot_limited":
        "Shot limitation\n(resolved by more shots alone)",
    "unknown":
        "Compounded budget interaction\n(needs candidates $+$ beam together)",
    "noise_sensitivity":
        "Noise sensitivity\n(resolved by disabling noise)",
    "unknown_attribution_timeout":
        "Attribution timed out",
}
ATTRIBUTION_COLOR = {
    "quantum_estimation_shot_limited": BUDGET,
    "unknown": BUDGET,
    "noise_sensitivity": NOISE,
    "unknown_attribution_timeout": OTHER,
}
# Ordered so the two classical-budget mechanisms are adjacent, which is
# what makes the 72% bracket legible.
ATTRIBUTION_ORDER = ["quantum_estimation_shot_limited", "unknown",
                     "noise_sensitivity", "unknown_attribution_timeout"]


def fig2_failure_attribution() -> Path:
    frame = pd.read_csv(DATA_DIR / "phase5_attribution.csv")
    counts = frame["failure_category"].value_counts()
    values = [int(counts.get(k, 0)) for k in ATTRIBUTION_ORDER]
    total = sum(values)

    fig, ax = plt.subplots(figsize=(6.6, 2.9))
    ypos = np.arange(len(values))[::-1]
    bars = ax.barh(ypos, values, height=0.62,
                   color=[ATTRIBUTION_COLOR[k] for k in ATTRIBUTION_ORDER],
                   edgecolor="black", linewidth=0.5)
    for y, v in zip(ypos, values):
        ax.text(v + 0.9, y, f"{v}  ({v/total*100:.0f}%)",
                va="center", ha="left", fontsize=8.5)

    ax.set_yticks(ypos)
    ax.set_yticklabels([ATTRIBUTION_LABEL[k] for k in ATTRIBUTION_ORDER])
    ax.set_xlabel(f"Failures attributed by the counterfactual ladder (n = {total})")
    ax.set_xlim(0, max(values) * 1.62)

    # Bracket spanning the two classical-budget mechanisms.
    budget_share = (values[0] + values[1]) / total * 100.0
    x = max(values) * 1.36
    ax.plot([x, x], [ypos[1] - 0.34, ypos[0] + 0.34], color="black", lw=1.0)
    ax.plot([x - 0.9, x], [ypos[0] + 0.34, ypos[0] + 0.34], color="black", lw=1.0)
    ax.plot([x - 0.9, x], [ypos[1] - 0.34, ypos[1] - 0.34], color="black", lw=1.0)
    ax.text(x + 1.4, (ypos[0] + ypos[1]) / 2,
            f"classical budget\nstarvation: {budget_share:.0f}%",
            va="center", ha="left", fontsize=8.5, fontweight="bold")

    fig.tight_layout()
    out = PLOT_DIR / "paper_fig2_failure_attribution.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------
# Figure 3 -- held-out generalisation
# ---------------------------------------------------------------------

ARM_ORDER = ["baseline", "B", "C", "D", "full"]
ARM_LABEL = {"baseline": "Baseline (published)", "B": "B  candidate coverage",
             "C": "C  beam width", "D": "D  shot stopping",
             "full": "Full  (B+C+D)"}
ARM_COLOR = {"baseline": "#9aa0a6", "B": "#3d6fb4", "C": "#c98b3a",
             "D": "#4f9d69", "full": "#b5474d"}

INSTANCES = [
    ("phase6_ablation_diagnostics.csv", "$N{=}21,\\ a{=}19,\\ r{=}6$", "canonical\n(design grid)"),
    ("phase7_heldout_N21_a2.csv", "$N{=}21,\\ a{=}2,\\ r{=}6$", "held-out\n(new base)"),
    ("phase7_heldout_N15_a2.csv", "$N{=}15,\\ a{=}2,\\ r{=}4$", "held-out\n(new modulus)"),
    ("phase7_heldout_N21_a8.csv", "$N{=}21,\\ a{=}8,\\ r{=}2$", "held-out\n(no headroom)"),
]


def fig3_heldout() -> Path:
    fig, axes = plt.subplots(1, 4, figsize=(7.4, 3.3), sharey=True)
    for ax, (filename, title, subtitle) in zip(axes, INSTANCES):
        frame = pd.read_csv(DATA_DIR / filename)
        rates = [frame[f"success_{a}"].mean() * 100.0 for a in ARM_ORDER]
        base = rates[0]
        bars = ax.bar(np.arange(len(ARM_ORDER)), rates, width=0.74,
                      color=[ARM_COLOR[a] for a in ARM_ORDER],
                      edgecolor="black", linewidth=0.45)
        if max(rates) - min(rates) < 1e-9:
            # Ceiling panel: five identical labels would collide and say
            # nothing five times.
            ax.text(2.0, max(rates) + 2.0, f"all arms {rates[0]:.1f}",
                    ha="center", va="bottom", fontsize=7.6)
        else:
            for x, v in zip(np.arange(len(ARM_ORDER)), rates):
                ax.text(x, v + 2.0, f"{v:.1f}", ha="center", va="bottom", fontsize=7.4)
        ax.axhline(base, color="#9aa0a6", lw=0.9, ls="--", zorder=0)
        ax.axhline(100, color="black", lw=0.6, ls=":", zorder=0)
        ax.set_title(title + "\n" + subtitle, fontsize=8.2, linespacing=1.25)
        ax.set_xticks(np.arange(len(ARM_ORDER)))
        ax.set_xticklabels(["Base", "B", "C", "D", "Full"], fontsize=8)
        ax.set_ylim(0, 118)
        ax.set_yticks([0, 20, 40, 60, 80, 100])
        # Flag the ceiling case, which is why no arm can move it.
        if base >= 99.9:
            ax.text(2.0, 50, "baseline already at\nceiling: no headroom\nfor any arm to recover",
                    ha="center", va="center", fontsize=7.4, style="italic",
                    color="#333333",
                    bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                              edgecolor="#999999", linewidth=0.5, alpha=0.95))
    axes[0].set_ylabel("End-to-end success rate (%)")

    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=ARM_COLOR[a],
                             edgecolor="black", linewidth=0.45) for a in ARM_ORDER]
    fig.legend(handles, [ARM_LABEL[a] for a in ARM_ORDER], ncol=5,
               loc="lower center", frameon=False, bbox_to_anchor=(0.5, -0.02),
               fontsize=7.8, columnspacing=1.2, handlelength=1.1)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    out = PLOT_DIR / "paper_fig3_heldout.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------
# Figure 4 -- matched-resource controls
# ---------------------------------------------------------------------

MATCHED = [
    ("B", "candidates_matched_vs_B", "B\ncandidate coverage", "$p = 1.00$"),
    ("C", "beam_matched_vs_C", "C\nbeam width", "$p = 1.00$"),
    ("D", "shots_matched_vs_D", "D\nshot stopping", "$p = 0.67$"),
    ("full", "shots_matched_vs_full", "Full\n(B+C+D)", "$p = 7{\\times}10^{-12}$"),
]


def fig4_matched_resource() -> Path:
    frame = pd.read_csv(DATA_DIR / "phase7_matched_resource.csv")
    base = frame["success_baseline"].mean() * 100.0
    adaptive = [frame[f"success_{a}"].mean() * 100.0 for a, _, _, _ in MATCHED]
    control = [frame[c].mean() * 100.0 for _, c, _, _ in MATCHED]

    fig, ax = plt.subplots(figsize=(6.8, 3.6))
    width = 0.33
    pos = np.arange(len(MATCHED))
    b1 = ax.bar(pos - width/2, adaptive, width, label="Adaptive arm (confidence-guided)",
                color=ADAPTIVE, edgecolor="black", linewidth=0.5)
    b2 = ax.bar(pos + width/2, control, width,
                label="Matched-resource control (same budget, uniform, no guidance)",
                color=CONTROL, edgecolor="black", linewidth=0.9, hatch="////")
    for bars in (b1, b2):
        for bar in bars:
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1.1,
                    f"{bar.get_height():.1f}", ha="center", va="bottom", fontsize=8)

    ax.axhline(base, color="#9aa0a6", lw=1.1, ls="--", zorder=0,
               label=f"unmodified baseline ({base:.1f}%)")

    for x, (_, _, _, ptxt) in zip(pos, MATCHED):
        top = max(adaptive[int(x)], control[int(x)])
        significant = "times" in ptxt
        ax.annotate("", xy=(x - width/2, top + 9.0), xytext=(x + width/2, top + 9.0),
                    arrowprops=dict(arrowstyle="-", lw=0.8, color="black"))
        ax.text(x, top + 10.2, ("significant\n" if significant else "not significant\n") + ptxt,
                ha="center", va="bottom", fontsize=7.6,
                fontweight="bold" if significant else "normal")

    ax.set_xticks(pos)
    ax.set_xticklabels([lbl for _, _, lbl, _ in MATCHED], fontsize=8.4)
    ax.set_ylabel("End-to-end success rate (%)")
    ax.set_ylim(0, 140)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.legend(loc="lower center", frameon=False, ncol=1,
              bbox_to_anchor=(0.5, 1.10), fontsize=7.8)
    fig.tight_layout()
    out = PLOT_DIR / "paper_fig4_matched_resource.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------
# Figure 5 -- calibration
# ---------------------------------------------------------------------

def fig5_calibration() -> Path:
    from sklearn.isotonic import IsotonicRegression

    frame = pd.read_csv(DATA_DIR / "phase6_calibration.csv")
    truth = frame["correct"].astype(float)
    raw_brier = float(((frame["c_w"] - truth) ** 2).mean())
    model = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
    calibrated = model.fit_transform(frame["c_w"], truth)
    cal_brier = float(((calibrated - truth) ** 2).mean())

    edges = np.linspace(0.0, 1.0, 11)

    def binned(pred):
        d = pd.DataFrame({"p": pred, "y": truth})
        d["bin"] = pd.cut(d["p"], bins=edges, include_lowest=True)
        g = d.groupby("bin", observed=True).agg(n=("p", "size"), pred=("p", "mean"),
                                                obs=("y", "mean"))
        return g

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.5), sharey=True)
    for ax, (pred, brier, title) in zip(axes, [
            (frame["c_w"], raw_brier, "Raw $C_w$"),
            (pd.Series(calibrated), cal_brier, "After isotonic recalibration")]):
        g = binned(pred)
        ax.plot([0.4, 1.0], [0.4, 1.0], color="black", lw=0.9, ls=":",
                label="perfect calibration", zorder=1)
        ax.fill_between(g["pred"], g["obs"], g["pred"], color="#c0504d", alpha=0.18,
                        zorder=2, label="overconfidence gap")
        sizes = 12 + 90 * (g["n"] / g["n"].max())
        ax.plot(g["pred"], g["obs"], color="#3d6fb4", lw=1.3, zorder=3)
        ax.scatter(g["pred"], g["obs"], s=sizes, color="#3d6fb4",
                   edgecolor="black", linewidth=0.5, zorder=4,
                   label="observed (marker area $\\propto$ windows in bin)")
        ax.set_xlim(0.42, 1.02); ax.set_ylim(0.42, 1.02)
        ax.set_xlabel("predicted confidence")
        ax.set_title(f"{title}  —  Brier $=$ {brier:.4f}")
        ax.set_aspect("equal", adjustable="box")

    axes[0].set_ylabel("observed fraction correct")
    # Mark module D's stopping threshold on the raw panel only: it is the
    # raw statistic that the shipped rule thresholds on.
    axes[0].axvline(0.95, color="#7a3b8f", lw=1.0, ls="-.", zorder=5)
    axes[0].text(0.945, 0.47, "module D stops here\n$1-\\varepsilon = 0.95$",
                 rotation=90, ha="right", va="bottom", fontsize=7.2, color="#7a3b8f")
    axes[0].annotate("predicted 0.96,\nobserved 0.79", xy=(0.964, 0.792),
                     xytext=(0.60, 0.90), fontsize=7.6,
                     arrowprops=dict(arrowstyle="->", lw=0.8, color="#c0504d"))
    axes[0].legend(loc="lower right", frameon=False, fontsize=6.9)
    fig.tight_layout()
    out = PLOT_DIR / "paper_fig5_calibration.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"   raw Brier {raw_brier:.4f} / recalibrated {cal_brier:.4f}"
          "  (must match 0.2239 / 0.1988)")
    assert f"{raw_brier:.4f}" == "0.2239" and f"{cal_brier:.4f}" == "0.1988"
    return out


def main() -> None:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    for fn in (fig2_failure_attribution, fig3_heldout,
               fig4_matched_resource, fig5_calibration):
        print(f"wrote {fn()}")


if __name__ == "__main__":
    main()
