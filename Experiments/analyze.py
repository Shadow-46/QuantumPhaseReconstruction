"""
Purpose
    Statistical analysis and publication-quality plotting for the
    failure-mode study: Wilson CIs on every sensitivity curve, a
    multivariable logistic regression for parameter effect-size ranking,
    noise-interaction tests, log-log scaling fits, and a failure-mode
    Pareto chart.
Inputs
    The CSVs written by Phases 1-5 under Data/failure_study/.
Outputs
    Data/failure_study/plots/*.png and summary_tables/*.csv.
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Experiments import harness

FS = harness.FAILURE_STUDY_DIR
PLOTS = harness.PLOTS_DIR
TABLES = harness.TABLES_DIR


def _load(name: str) -> pd.DataFrame | None:
    path = FS / f"{name}.csv"
    return pd.read_csv(path) if path.exists() else None


def ofat_sensitivity_plots() -> None:
    frame = _load("phase1_ofat")
    if frame is None:
        print("phase1_ofat.csv missing, skipping")
        return
    params = frame["param_swept"].dropna().unique()
    summary_rows = []
    for param in params:
        sub = frame[frame["param_swept"] == param]
        grouped = sub.groupby("spec_sweep_value")["result_success"]
        xs, ps, los, his, ns = [], [], [], [], []
        for value, group in grouped:
            n = len(group)
            successes = int(group.sum())
            p, lo, hi = harness.wilson_ci(successes, n)
            xs.append(value)
            ps.append(p)
            los.append(lo)
            his.append(hi)
            ns.append(n)
            summary_rows.append({"param": param, "value": value, "n": n, "success_rate": p, "ci_low": lo, "ci_high": hi})
        order = np.argsort(xs)
        xs, ps, los, his = (np.array(a)[order] for a in (xs, ps, los, his))
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(xs, ps, marker="o", color="#2b6cb0")
        ax.fill_between(xs, los, his, alpha=0.25, color="#2b6cb0")
        ax.set_xlabel(param)
        ax.set_ylabel("success rate (Wilson 95% CI)")
        ax.set_ylim(-0.02, 1.02)
        if param == "shots":
            ax.set_xscale("log")
        ax.set_title(f"Success rate vs {param}")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(PLOTS / f"ofat_{param}.png", dpi=160)
        plt.close(fig)
    pd.DataFrame(summary_rows).to_csv(TABLES / "phase1_ofat_summary.csv", index=False)
    print("wrote OFAT sensitivity plots")


def noise_plots() -> None:
    single = _load("phase2_noise_single")
    if single is not None:
        fig, ax = plt.subplots(figsize=(6, 4))
        for channel, group in single.groupby("param_swept"):
            agg = group.groupby("spec_sweep_value")["result_success"].mean().sort_index()
            ax.plot(agg.index, agg.values, marker="o", label=channel)
        ax.set_xscale("symlog", linthresh=1e-3)
        ax.set_xlabel("noise magnitude")
        ax.set_ylabel("success rate")
        ax.set_title("Success rate vs single-source noise magnitude")
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(PLOTS / "noise_single_source.png", dpi=160)
        plt.close(fig)

    combined = _load("phase2_noise_combined")
    if combined is not None:
        order = [lvl for lvl in ["none", "low", "med", "high"] if lvl in set(combined["spec_level_2q"]) | set(combined["spec_level_ro"])]
        pivot = combined.pivot_table(index="spec_level_2q", columns="spec_level_ro", values="result_success", aggfunc="mean")
        pivot = pivot.reindex(index=order, columns=order)
        fig, ax = plt.subplots(figsize=(5.5, 5))
        im = ax.imshow(pivot.values, cmap="RdYlGn", vmin=0, vmax=1)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order)
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels(order)
        ax.set_xlabel("readout error level")
        ax.set_ylabel("2-qubit depolarizing level")
        ax.set_title("Success rate: 2q-depolarizing x readout (averaged over 1q levels)")
        fig.colorbar(im, ax=ax, label="success rate")
        fig.tight_layout()
        fig.savefig(PLOTS / "noise_combined_heatmap.png", dpi=160)
        plt.close(fig)

        # interaction check: compare observed combined degradation to the
        # product of the two marginal single-source success rates (a simple,
        # assumption-light additivity baseline in odds-ratio terms).
        marginal_2q = combined.groupby("spec_level_2q")["result_success"].mean()
        marginal_ro = combined.groupby("spec_level_ro")["result_success"].mean()
        rows = []
        for l2q in order:
            for lro in order:
                observed = combined[(combined.spec_level_2q == l2q) & (combined.spec_level_ro == lro)]["result_success"].mean()
                independence_baseline = marginal_2q.get(l2q, np.nan) * marginal_ro.get(lro, np.nan)
                rows.append({"spec_level_2q": l2q, "spec_level_ro": lro, "observed": observed, "independence_baseline": independence_baseline, "excess": observed - independence_baseline})
        pd.DataFrame(rows).to_csv(TABLES / "phase2_noise_interaction.csv", index=False)
    print("wrote noise plots")


def adversarial_plots() -> None:
    adv = _load("phase3_adversarial")
    if adv is not None:
        agg = adv.groupby("condition")["success"].agg(["mean", "count"])
        cis = adv.groupby("condition")["success"].apply(lambda s: harness.wilson_ci(int(s.sum()), len(s)))
        agg["ci_low"] = [c[1] for c in cis]
        agg["ci_high"] = [c[2] for c in cis]
        agg = agg.sort_values("mean")
        fig, ax = plt.subplots(figsize=(7, 4.5))
        yerr = [agg["mean"] - agg["ci_low"], agg["ci_high"] - agg["mean"]]
        ax.bar(agg.index, agg["mean"], yerr=yerr, capsize=4, color="#c05621")
        ax.set_ylabel("success rate")
        ax.set_title("Success rate under adversarial corruption conditions")
        ax.set_ylim(0, 1.05)
        plt.setp(ax.get_xticklabels(), rotation=25, ha="right")
        ax.grid(True, axis="y", alpha=0.3)
        fig.tight_layout()
        fig.savefig(PLOTS / "adversarial_conditions.png", dpi=160)
        plt.close(fig)

    carry = _load("phase3_carry_truth_table")
    if carry is not None:
        counts = carry["outcome"].value_counts()
        fig, ax = plt.subplots(figsize=(5.5, 4))
        colors = {"true_accept": "#2f855a", "true_reject": "#2b6cb0", "false_accept": "#c53030", "false_reject": "#dd6b20"}
        ax.bar(counts.index, counts.values, color=[colors.get(k, "gray") for k in counts.index])
        ax.set_ylabel("count")
        ax.set_title("Carry-check truth table (exhaustive enumeration)")
        fig.tight_layout()
        fig.savefig(PLOTS / "carry_truth_table.png", dpi=160)
        plt.close(fig)

    ambig = _load("phase3_ambiguity_boundary")
    if ambig is not None:
        agg = ambig.groupby("delta_from_boundary").agg(
            ambiguous_rate=("ambiguous_by_0.9_threshold", "mean"),
            flip_rate=("winner_disagreement_across_repeats", "mean"),
        ).sort_index(ascending=False)
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(agg.index, agg["ambiguous_rate"], marker="o", label="flagged ambiguous (top2/top1 > 0.9)")
        ax.plot(agg.index, agg["flip_rate"], marker="s", label="winner disagrees across repeats")
        ax.set_xlabel("distance from 0.5 bucket boundary (delta)")
        ax.invert_xaxis()
        ax.set_ylabel("rate")
        ax.set_title("Measurement-ambiguity boundary sweep")
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(PLOTS / "ambiguity_boundary.png", dpi=160)
        plt.close(fig)
    print("wrote adversarial plots")


def scaling_plots() -> None:
    by_n = _load("phase4_scaling_by_N")
    stats = _load("phase4_circuit_stats_by_N")
    if by_n is not None and stats is not None:
        runtime = by_n.groupby("spec_N")["runtime_seconds"].mean().reset_index()
        merged = runtime.merge(stats, left_on="spec_N", right_on="N")
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
        ax1.scatter(merged["work_qubits"], merged["runtime_seconds"], color="#2b6cb0")
        if (merged["runtime_seconds"] > 0).all():
            coeffs = np.polyfit(merged["work_qubits"], np.log(merged["runtime_seconds"]), 1)
            xs = np.linspace(merged["work_qubits"].min(), merged["work_qubits"].max(), 50)
            ax1.plot(xs, np.exp(coeffs[1] + coeffs[0] * xs), "--", color="gray", label=f"fit: log(t) ~ {coeffs[0]:.3f}*qubits")
            ax1.legend()
        ax1.set_yscale("log")
        ax1.set_xlabel("work-register qubits (ceil(log2 N))")
        ax1.set_ylabel("mean runtime (s, log scale)")
        ax1.set_title("Runtime vs problem size")
        ax1.grid(True, alpha=0.3)

        peak = by_n.groupby("spec_N")["peak_memory_bytes"].mean().reset_index()
        merged2 = peak.merge(stats, left_on="spec_N", right_on="N")
        ax2.scatter(merged2["max_circuit_qubits"], merged2["peak_memory_bytes"], color="#c05621")
        ax2.set_yscale("log")
        ax2.set_xlabel("max per-block circuit qubits (window+work)")
        ax2.set_ylabel("mean peak traced memory (bytes, log scale)")
        ax2.set_title("Memory vs per-block qubit count")
        ax2.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(PLOTS / "scaling_by_N.png", dpi=160)
        plt.close(fig)

    by_prec = _load("phase4_scaling_by_precision")
    by_win = _load("phase4_scaling_by_window_size")
    if by_prec is not None and by_win is not None:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
        agg1 = by_prec.groupby("spec_total_precision")["runtime_seconds"].mean()
        ax1.plot(agg1.index, agg1.values, marker="o", color="#2b6cb0")
        ax1.set_xlabel("total_precision (window_size fixed at 4)")
        ax1.set_ylabel("mean runtime (s)")
        ax1.set_title("Cost of adding more windows")
        ax1.grid(True, alpha=0.3)

        agg2 = by_win.groupby("spec_window_size")["runtime_seconds"].mean()
        ax2.plot(agg2.index, agg2.values, marker="s", color="#c05621")
        ax2.set_xlabel("window_size (total_precision fixed at 16)")
        ax2.set_ylabel("mean runtime (s)")
        ax2.set_title("Cost of widening each block")
        ax2.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(PLOTS / "scaling_decoupling.png", dpi=160)
        plt.close(fig)
    print("wrote scaling plots")


def logistic_regression_effects(n_bootstrap: int = 300) -> None:
    """Fit a standardized multivariable logistic regression (sklearn) of
    success on every OFAT parameter, with bootstrap 95% CIs on each
    coefficient, to rank parameters by effect size on the log-odds of
    success. Avoids a statsmodels dependency by using sklearn (already a
    transitive dependency here) plus a straightforward nonparametric
    bootstrap for interval estimates."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import log_loss

    # phase1_ofat showed zero-variance success (100%) under the qubit-capped
    # small-N regime this campaign was forced into (see REPORT.md scope
    # note), so a logistic fit on it is undefined (single class). Fall back
    # to whichever available dataset actually has outcome variance, noting
    # the swap plainly rather than silently reporting nothing.
    candidates = [
        ("phase1_ofat", ["spec_total_precision", "spec_window_size", "spec_overlap", "spec_candidate_count", "spec_max_paths", "log_shots"]),
        ("phase1b_combined_stress", ["spec_total_precision", "spec_overlap", "spec_candidate_count", "spec_max_paths", "log_shots"]),
        ("phase2_noise_combined", ["spec_noise_1q", "spec_noise_2q", "spec_noise_ro"]),
    ]
    frame = None
    predictors: list[str] = []
    source_name = None
    for name, preds in candidates:
        candidate_frame = _load(name)
        if candidate_frame is None:
            continue
        candidate_frame = candidate_frame.copy()
        if "spec_shots" in candidate_frame.columns and "log_shots" in preds:
            candidate_frame["log_shots"] = np.log(candidate_frame["spec_shots"])
        if candidate_frame["result_success"].nunique() > 1:
            frame, predictors, source_name = candidate_frame, preds, name
            break
    if frame is None:
        print("no dataset with outcome variance available; skipping logistic regression")
        return
    print(f"logistic regression effect sizes computed from: {source_name}")

    df = frame
    df = df.dropna(subset=predictors + ["result_success"])
    X_raw = df[predictors].astype(float).values
    mean, std = X_raw.mean(axis=0), X_raw.std(axis=0)
    std[std == 0] = 1.0
    X = (X_raw - mean) / std
    y = df["result_success"].astype(int).values

    def fit_coefs(Xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
        model = LogisticRegression(penalty=None, max_iter=2000)
        model.fit(Xs, ys)
        return model.coef_.ravel()

    point = fit_coefs(X, y)
    rng = np.random.default_rng(0)
    boot = np.empty((n_bootstrap, len(predictors)))
    for i in range(n_bootstrap):
        idx = rng.integers(0, len(y), len(y))
        try:
            boot[i] = fit_coefs(X[idx], y[idx])
        except Exception:  # noqa: BLE001 - a degenerate resample (all one class) is skipped
            boot[i] = np.nan
    lo = np.nanpercentile(boot, 2.5, axis=0)
    hi = np.nanpercentile(boot, 97.5, axis=0)

    base_model = LogisticRegression(penalty=None, max_iter=2000).fit(X, y)
    ll_model = -log_loss(y, base_model.predict_proba(X), normalize=False)
    p_null = y.mean()
    ll_null = -log_loss(y, np.full_like(y, p_null, dtype=float), normalize=False)
    pseudo_r2 = 1 - ll_model / ll_null if ll_null != 0 else float("nan")

    order = np.argsort(-np.abs(point))
    names = [predictors[i].replace("spec_", "") for i in order]
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.errorbar(point[order], range(len(order)), xerr=[point[order] - lo[order], hi[order] - point[order]], fmt="o", color="#2b6cb0")
    ax.axvline(0, color="gray", linestyle="--")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(names)
    ax.set_xlabel("standardized logistic-regression coefficient (log-odds of success)")
    ax.set_title(f"Parameter effect-size ranking (McFadden pseudo-R2={pseudo_r2:.2f})")
    ax.grid(True, axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "logistic_effect_sizes.png", dpi=160)
    plt.close(fig)
    pd.DataFrame({"parameter": names, "coef": point[order], "ci_low": lo[order], "ci_high": hi[order]}).to_csv(
        TABLES / "phase1_logistic_regression.csv", index=False
    )
    print("wrote logistic regression effect-size plot")


def failure_mode_pareto() -> None:
    attribution = _load("phase5_attribution")
    if attribution is None:
        return
    counts = attribution["failure_category"].value_counts()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(counts.index, counts.values, color="#805ad5")
    ax.set_ylabel("count of attributed failures")
    ax.set_title("Failure-mode Pareto chart")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "failure_mode_pareto.pdf")
    plt.close(fig)
    counts.to_csv(TABLES / "phase5_failure_mode_counts.csv")
    print("wrote failure-mode Pareto chart")


def main() -> None:
    ofat_sensitivity_plots()
    noise_plots()
    adversarial_plots()
    scaling_plots()
    logistic_regression_effects()
    failure_mode_pareto()


if __name__ == "__main__":
    main()
