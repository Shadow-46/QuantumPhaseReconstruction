"""
Purpose
    Window-level validation of module D (confidence-gated shot stopping),
    against a shots-matched uniform control, swept over three values of the
    shot cap S_max. Supplies Paper A's own evidence for module D, at the
    window, complementing the companion study's trial-level evaluation.
Theory
    Module D spends extra shots on windows whose evidence is thin. Two
    questions are separable and are separated here:

      (1) does spending help?          D vs. an unspent baseline
      (2) does spending THERE help?    D vs. the same total budget spread
                                       uniformly, with no confidence
                                       computation anywhere in the control

    Only (2) is reviewer-facing. (1) is near-uninformative on its own because
    conditioning on C_w < 1 - epsilon selects the windows with the thinnest
    evidence and hence the most headroom, so every spending arm beats
    baseline. The control is an oracle, not an algorithm: its budget can only
    be known after watching module D choose it.

    Primary outcome is TOP-1 IDENTIFICATION -- whether the observed argmax
    lies in the true marginal's argmax set. It is a fixed target, identical
    across arms, so the pairing is clean. The phase6 `correct` label is not
    usable as a primary outcome here because its target (the observed top-2
    pair) moves between arms; retained mass is module B's outcome and would
    confound the two modules. Both are reported as secondaries, the latter
    under a FIXED top-k rule so only the counts vary.

    Why the cap sweep is core. The pilot (Experiments/replay.py) found that at
    the shipped S_max = 64, 79.4% of firing windows run all the way to the
    cap. Where that happens the cap, not the confidence rule, sets the budget,
    and D's allocation approaches uniform by construction. Sweeping
    S_max in {32, 64, 256} is what separates "the statistic carries no
    placement information" from "the cap binds before it can act". Without the
    sweep a null result at the shipped cap is uninterpretable.

    Scope. The replay corpus is heavily concentrated: the overwhelming
    majority of D-eligible windows come from one instance (N=21, a=19) at one
    window width. Concentration is measured and reported by this script rather
    than assumed; see the `concentration` block in the summary output. No
    claim here generalises beyond the sampled corpus.
Inputs
    Data/failure_study/raw_windows/*.json (from recapture_window_counts.py)
Outputs
    Data/failure_study/paperA_moduleD_window_validation.csv.gz
    Data/failure_study/summary_tables/paperA_moduleD_summary.csv
    Data/failure_study/summary_tables/paperA_moduleD_tests.csv
    Data/failure_study/plots/paperA_figA_moduleD_window.png
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
from functools import lru_cache

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.special import betainc, betaln
from scipy.stats import binomtest

from Circuits.windowed_qpe import WindowSpec
from Experiments import harness
from Experiments.phase6_calibration import true_window_distribution
from Experiments.replay import (
    RAW_WINDOWS_DIR,
    SALT,
    apply_prefix,
    crn_stream,
    marginal_vector,
    replay_module_d,
)
from Reconstruction.candidate_generation import generate_window_candidates
from Reconstruction.confidence import top1_vs_top2_confidence

SUMMARY_DIR = harness.FAILURE_STUDY_DIR / "summary_tables"

EPSILON = 0.05
CAPS = (32, 64, 256)
SHIPPED_CAP = 64
REPLICATES = 20
BOOTSTRAP_RESAMPLES = 2000
BOOTSTRAP_SEED = 0

# The corpus is frozen upstream; asserted so this script cannot silently drift.
EXPECTED_WINDOWS = 3354
EXPECTED_TRIALS = 579

ARMS = ("baseline", "D", "uniform")


# Quantities that are invariant across replicates and caps. The exact
# marginal depends only on the window key, and the baseline arm never spends,
# so both would otherwise be recomputed once per replicate per cap.
_MARGINAL_CACHE: dict[tuple, tuple[tuple[str, ...], np.ndarray]] = {}
_BASELINE_CACHE: dict[tuple, dict] = {}


def cached_marginal(a: int, N: int, spec: WindowSpec) -> tuple[tuple[str, ...], np.ndarray]:
    key = (a, N, spec.start, spec.width, spec.total_precision)
    if key not in _MARGINAL_CACHE:
        _MARGINAL_CACHE[key] = marginal_vector(a, N, spec)
    return _MARGINAL_CACHE[key]


def _fast_confidence(n1: int, n2: int, s_w: int, alpha0: float = 0.5) -> float:
    """Numerically identical restatement of the shipped C_w, for replay only.

    The shipped implementation integrates `beta.pdf(x,a1,b1)*beta.cdf(x,a2,b2)`
    via scipy's rv_continuous objects. Profiling the replay showed 99.2% of
    total runtime inside that call: each quadrature evaluates the integrand
    ~170 times, and every evaluation pays the generic distribution
    dispatcher's argument validation and broadcasting. The arithmetic itself
    is negligible.

    This evaluates THE SAME integrand through scipy.special directly -- the
    Beta density in log form, and betainc for the CDF -- with the same
    quadrature call, the same limit, and the same clamp. Reconstruction/
    confidence.py is deliberately left untouched: the shipped statistic is
    what module D thresholds and what every frozen result was produced with.
    Agreement is asserted, not assumed, by verify_confidence_equivalence()."""
    a1, b1 = alpha0 + n1, alpha0 + s_w - n1
    a2, b2 = alpha0 + n2, alpha0 + s_w - n2

    def integrand(x: float) -> float:
        density = np.exp((a1 - 1.0) * np.log(x) + (b1 - 1.0) * np.log1p(-x) - betaln(a1, b1))
        return density * betainc(a2, b2, x)

    value, _ = quad(integrand, 0.0, 1.0, limit=200)
    return float(min(1.0, max(0.0, value)))


# Every distinct argument triple the replay evaluates, so the equivalence
# guard in main() can re-check exactly what was used rather than a sample.
_CONSUMED_TRIPLES: set[tuple[int, int, int]] = set()


@lru_cache(maxsize=None)
def _confidence_cached(n1: int, n2: int, s_w: int) -> float:
    return _fast_confidence(n1, n2, s_w)


def confidence(n1: int, n2: int, s_w: int) -> float:
    """Memoised C_w. Pure function of three small integers, so caching is
    exact. Routes to the fast restatement above; see its docstring."""
    _CONSUMED_TRIPLES.add((n1, n2, s_w))
    return _confidence_cached(n1, n2, s_w)


def verify_confidence_equivalence(triples: list[tuple[int, int, int]],
                                  tolerance: float = 1e-9) -> tuple[float, int]:
    """Assert the fast restatement reproduces the shipped statistic, on the
    triples actually consumed. Two things are checked: agreement to within
    `tolerance`, and -- what actually matters for module D -- that no window
    lands on a different side of the stopping threshold under the two
    implementations. Returns (max absolute difference, threshold
    disagreements)."""
    worst = 0.0
    disagreements = 0
    for n1, n2, s_w in triples:
        shipped = top1_vs_top2_confidence(n1, n2, s_w)
        fast = _fast_confidence(n1, n2, s_w)
        worst = max(worst, abs(shipped - fast))
        if (shipped >= 1.0 - EPSILON) != (fast >= 1.0 - EPSILON):
            disagreements += 1
    assert worst < tolerance, f"fast C_w differs from shipped by {worst:.3e}"
    assert disagreements == 0, f"{disagreements} windows change stopping decision"
    return worst, disagreements


def top_two(counts: dict[str, int]) -> tuple[int, int]:
    ranked = sorted(counts.values(), reverse=True)
    return (ranked[0] if ranked else 0, ranked[1] if len(ranked) > 1 else 0)


def observed_argmax(counts: dict[str, int]) -> str:
    """Deterministic observed top-1: ties broken lexicographically, so the
    result cannot depend on dict insertion order across arms."""
    top = max(counts.values())
    return min(bits for bits, count in counts.items() if count == top)


def true_argmax_set(true_dist: dict[str, float], tol: float = 1e-12) -> set[str]:
    top = max(true_dist.values())
    return {bits for bits, p in true_dist.items() if p >= top - tol}


def true_margin(true_dist: dict[str, float]) -> float:
    ordered = sorted(true_dist.values(), reverse=True)
    return float(ordered[0] - ordered[1]) if len(ordered) > 1 else float(ordered[0])


def total_variation(counts: dict[str, int], true_dist: dict[str, float]) -> float:
    s_w = sum(counts.values())
    if s_w == 0:
        return float("nan")
    keys = set(counts) | set(true_dist)
    return 0.5 * sum(abs(counts.get(b, 0) / s_w - true_dist.get(b, 0.0)) for b in keys)


def largest_remainder(total: int, parts: int) -> list[int]:
    """Split `total` into `parts` integers summing EXACTLY to total. The
    companion study's trial-level control used floor division, which
    under-matched the budget it was supposed to equal; at the window level
    with few windows per trial that residual is proportionally larger, so the
    match is made exact here."""
    if parts <= 0:
        return []
    base, remainder = divmod(total, parts)
    return [base + (1 if i < remainder else 0) for i in range(parts)]


def outcomes(counts: dict[str, int], true_dist: dict[str, float], argmax_set: set[str],
             spec: WindowSpec, candidate_count: int, suffix: str) -> dict:
    """The three outcome families for one arm's final counts."""
    kept = {c.local_bits for c in generate_window_candidates(
        counts, spec.start, spec.width, spec.total_precision, candidate_count)}
    return {
        f"argmax_correct_{suffix}": observed_argmax(counts) in argmax_set,
        f"tv_{suffix}": total_variation(counts, true_dist),
        f"true_mass_topk_{suffix}": sum(true_dist.get(b, 0.0) for b in kept),
        f"argmax_kept_topk_{suffix}": bool(argmax_set & kept),
    }


def replay_trial(record: dict, ordinal: int, replicate: int, s_max: int) -> list[dict]:
    """Replay one trial under all three arms at one cap, one replicate.

    Arms share a common random number stream per window: each consumes a
    PREFIX of the same draws sized to its own budget, so the arms differ only
    in how many shots a window receives, never in which shots."""
    spec_dict = record["spec"]
    a, N = spec_dict["a"], spec_dict["N"]
    candidate_count = spec_dict["candidate_count"]

    per_window = []
    for window_meta, counts in zip(record["windows"], record["window_counts"]):
        spec = WindowSpec(window_meta["start"], window_meta["width"], window_meta["total_precision"])
        true_dist = true_window_distribution(a, N, spec)
        if not true_dist:
            continue
        patterns, probs = cached_marginal(a, N, spec)
        s_w = sum(counts.values())
        stream = crn_stream(
            probs, max(0, s_max - s_w),
            [SALT, ordinal, spec.start, spec.width, s_w, replicate],
        )
        final_d, extra_d, rounds_d = replay_module_d(
            counts, patterns, stream, epsilon=EPSILON, s_max=s_max,
            confidence_fn=confidence)
        per_window.append({
            "spec": spec, "counts": counts, "true_dist": true_dist,
            "patterns": patterns, "stream": stream, "s_w": s_w,
            "final_d": final_d, "extra_d": extra_d, "rounds_d": rounds_d,
        })

    if not per_window:
        return []

    # The control receives module D's OWN realised total, spread uniformly
    # over every window of the trial -- including windows D declined to fund.
    total_extra = sum(w["extra_d"] for w in per_window)
    allocation = largest_remainder(total_extra, len(per_window))
    assert sum(allocation) == total_extra, "control budget must match D exactly"
    trial_fired = total_extra > 0
    allocation_uniform = len({w["extra_d"] for w in per_window}) == 1

    rows = []
    for w, extra_uniform in zip(per_window, allocation):
        spec = w["spec"]
        true_dist = w["true_dist"]
        argmax_set = true_argmax_set(true_dist)
        final_uniform = apply_prefix(w["counts"], w["patterns"], w["stream"], extra_uniform)
        n1_final, n2_final = top_two(w["final_d"])
        n1_init, n2_init = top_two(w["counts"])
        row = {
            "trial_id": record["trial_id"],
            "phase": record["phase"],
            "N": N,
            "a": a,
            "window_start": spec.start,
            "window_width": spec.width,
            "total_precision": spec.total_precision,
            "candidate_count": candidate_count,
            "s_w_initial": w["s_w"],
            "s_max": s_max,
            "replicate": replicate,
            "k_patterns": 2 ** spec.width,
            "d_eligible": w["extra_d"] > 0,
            "trial_fired": trial_fired,
            "d_allocation_uniform": allocation_uniform,
            "true_argmax_tied": len(argmax_set) > 1,
            "true_margin": true_margin(true_dist),
            "c_w_initial": confidence(n1_init, n2_init, w["s_w"]),
            "s_final_D": w["s_w"] + w["extra_d"],
            "extra_D": w["extra_d"],
            "rounds_D": w["rounds_d"],
            "hit_cap_D": bool(w["extra_d"] > 0 and w["s_w"] + w["extra_d"] >= s_max),
            "c_w_final_D": confidence(n1_final, n2_final, w["s_w"] + w["extra_d"]),
            "extra_uniform": extra_uniform,
            # Did the CONTROL's own allocation push a window past the cap?
            # Windows already configured above the cap receive nothing and are
            # not counted here -- they were never below it to begin with.
            "control_exceeds_smax": bool(extra_uniform > 0 and w["s_w"] + extra_uniform > s_max),
        }
        baseline_key = (record["trial_id"], spec.start)
        if baseline_key not in _BASELINE_CACHE:
            _BASELINE_CACHE[baseline_key] = outcomes(
                w["counts"], true_dist, argmax_set, spec, candidate_count, "baseline")
        row.update(_BASELINE_CACHE[baseline_key])
        row.update(outcomes(w["final_d"], true_dist, argmax_set, spec, candidate_count, "D"))
        row.update(outcomes(final_uniform, true_dist, argmax_set, spec, candidate_count, "uniform"))
        rows.append(row)
    return rows


def build_frame(caps: tuple[int, ...] = CAPS, replicates: int = REPLICATES) -> pd.DataFrame:
    paths = sorted(RAW_WINDOWS_DIR.glob("*.json"))
    records = [(ordinal, json.loads(p.read_text())) for ordinal, p in enumerate(paths)]
    assert len(records) == EXPECTED_TRIALS, (
        f"corpus changed: {len(records)} trials, expected {EXPECTED_TRIALS}")

    rows = []
    for s_max in caps:
        for replicate in range(replicates):
            for ordinal, record in records:
                rows.extend(replay_trial(record, ordinal, replicate, s_max))
        print(f"  cap {s_max}: done")

    frame = pd.DataFrame(rows)
    windows_per_cap = len(frame) // (len(caps) * replicates)
    assert windows_per_cap == EXPECTED_WINDOWS, (
        f"window count changed: {windows_per_cap}, expected {EXPECTED_WINDOWS}")
    # Gzipped: 201,240 rows x 37 columns is 62 MB plain, 6.4 MB compressed.
    # pandas reads it back transparently from the .gz suffix, so the full
    # 20-replicate design is released intact rather than thinned.
    frame.to_csv(harness.FAILURE_STUDY_DIR / "paperA_moduleD_window_validation.csv.gz",
                 index=False, compression="gzip")
    return frame


# --------------------------------------------------------------------------
# Statistics
# --------------------------------------------------------------------------

def clustered_bootstrap(per_window: pd.DataFrame, value: str, cluster: str = "cluster",
                        resamples: int = BOOTSTRAP_RESAMPLES) -> tuple[float, float, float]:
    """Percentile bootstrap CI on the mean of `value`, resampling CLUSTERS.

    Windows within a trial share a seed, a configuration and an instance, so a
    window-level bootstrap would be badly anti-conservative: the effective
    sample size is the number of trials, not the number of windows. The
    companion study's helper resamples rows and is deliberately not reused."""
    groups = [g[value].to_numpy() for _, g in per_window.groupby(cluster, sort=True)]
    if not groups:
        return float("nan"), float("nan"), float("nan")
    point = float(np.concatenate(groups).mean())
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    n = len(groups)
    means = np.empty(resamples)
    for i in range(resamples):
        pick = rng.integers(0, n, size=n)
        means[i] = np.concatenate([groups[j] for j in pick]).mean()
    lo, hi = np.percentile(means, [2.5, 97.5])
    return point, float(lo), float(hi)


def mcnemar_exact(reference: pd.Series, arm: pd.Series) -> tuple[int, int, float]:
    """Exact two-sided binomial McNemar, the convention the companion study
    uses throughout. b = reference succeeded and arm failed; c = the reverse."""
    b = int((reference & ~arm).sum())
    c = int((~reference & arm).sum())
    discordant = b + c
    p = 1.0 if discordant == 0 else binomtest(min(b, c), discordant, 0.5).pvalue
    return b, c, float(p)


def holm(pvalues: list[float]) -> list[float]:
    """Holm sequentially rejective correction, returned as adjusted p-values."""
    order = np.argsort(pvalues)
    m = len(pvalues)
    adjusted = [0.0] * m
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, min(1.0, (m - rank) * pvalues[idx]))
        adjusted[idx] = running
    return adjusted


def comparison(frame: pd.DataFrame, cap: int, reference: str, arm: str,
               outcome: str, subset: str, label: str) -> dict:
    """One pre-declared comparison. Effect size is the clustered bootstrap on
    replicate-averaged per-window differences; the formal test is exact
    McNemar on replicate 0, so the paired binary structure is preserved."""
    sub = frame[frame["s_max"] == cap]
    if subset == "firing trials":
        sub = sub[sub["trial_fired"]]
    elif subset == "non-uniform allocation":
        sub = sub[sub["trial_fired"] & ~sub["d_allocation_uniform"]]
    elif subset == "D-eligible windows":
        sub = sub[sub["d_eligible"]]

    ref_col, arm_col = f"{outcome}_{reference}", f"{outcome}_{arm}"

    averaged = (
        sub.groupby(["trial_id", "window_start"], sort=True)
        .agg(ref=(ref_col, "mean"), arm=(arm_col, "mean"))
        .reset_index()
    )
    averaged["diff"] = averaged["arm"] - averaged["ref"]
    averaged["cluster"] = averaged["trial_id"]
    point, lo, hi = clustered_bootstrap(averaged, "diff")

    first = sub[sub["replicate"] == 0]
    b, c, p = mcnemar_exact(first[ref_col].astype(bool), first[arm_col].astype(bool))

    return {
        "cap": cap,
        "label": label,
        "outcome": outcome,
        "reference": reference,
        "arm": arm,
        "subset": subset,
        "n_windows": int(averaged.shape[0]),
        "n_trials": int(averaged["cluster"].nunique()),
        "rate_reference": float(averaged["ref"].mean()),
        "rate_arm": float(averaged["arm"].mean()),
        "diff_pp": 100.0 * point,
        "ci_lo_pp": 100.0 * lo,
        "ci_hi_pp": 100.0 * hi,
        "b": b,
        "c": c,
        "mcnemar_p": p,
    }


def run_tests(frame: pd.DataFrame) -> pd.DataFrame:
    """The confirmatory family is three comparisons on the primary outcome at
    the SHIPPED cap. The cap sweep is exploratory and reported uncorrected.
    This family is disjoint from the companion study's family of nine."""
    confirmatory = [
        comparison(frame, SHIPPED_CAP, "baseline", "D", "argmax_correct",
                   "firing trials", "C1  D vs unspent baseline"),
        comparison(frame, SHIPPED_CAP, "uniform", "D", "argmax_correct",
                   "firing trials", "C2  D vs shots-matched uniform control"),
        comparison(frame, SHIPPED_CAP, "baseline", "D", "argmax_correct",
                   "non-uniform allocation", "C3  D vs baseline, non-degenerate subset"),
    ]
    adjusted = holm([row["mcnemar_p"] for row in confirmatory])
    for row, value in zip(confirmatory, adjusted):
        row["holm_p"] = value
        row["family"] = "confirmatory"

    exploratory = []
    for cap in CAPS:
        for reference, arm, label in (
            ("baseline", "D", "sweep  D vs baseline"),
            ("uniform", "D", "sweep  D vs matched control"),
        ):
            row = comparison(frame, cap, reference, arm, "argmax_correct",
                             "firing trials", label)
            row["holm_p"] = float("nan")
            row["family"] = "exploratory (cap sweep)"
            exploratory.append(row)
    for cap in CAPS:
        row = comparison(frame, cap, "uniform", "D", "argmax_kept_topk",
                         "firing trials", "sweep  D vs control, argmax kept under fixed top-k")
        row["holm_p"] = float("nan")
        row["family"] = "exploratory (secondary outcome)"
        exploratory.append(row)

    tests = pd.DataFrame(confirmatory + exploratory)
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    tests.to_csv(SUMMARY_DIR / "paperA_moduleD_tests.csv", index=False)
    return tests


def summarise(frame: pd.DataFrame) -> pd.DataFrame:
    replicates = frame["replicate"].nunique()
    rows = []
    for cap in sorted(frame["s_max"].unique()):
        sub = frame[frame["s_max"] == cap]
        firing = sub[sub["trial_fired"]]
        eligible = sub[sub["d_eligible"]]
        block = {
            "s_max": cap,
            "windows_total": int(len(sub) / replicates),
            "windows_eligible_mean": float(len(eligible) / replicates),
            "eligible_fraction": float(sub["d_eligible"].mean()),
            "cap_hit_fraction_of_eligible": float(eligible["hit_cap_D"].mean()) if len(eligible) else float("nan"),
            "allocation_uniform_fraction": float(firing["d_allocation_uniform"].mean()) if len(firing) else float("nan"),
            "control_exceeds_cap": int(sub["control_exceeds_smax"].sum()),
            "mean_extra_D": float(firing["extra_D"].mean()) if len(firing) else float("nan"),
            "mean_extra_uniform": float(firing["extra_uniform"].mean()) if len(firing) else float("nan"),
            "mean_rounds_D": float(eligible["rounds_D"].mean()) if len(eligible) else float("nan"),
        }
        for arm in ARMS:
            block[f"argmax_correct_{arm}"] = float(firing[f"argmax_correct_{arm}"].mean()) if len(firing) else float("nan")
            block[f"tv_{arm}"] = float(firing[f"tv_{arm}"].mean()) if len(firing) else float("nan")
            block[f"mass_topk_{arm}"] = float(firing[f"true_mass_topk_{arm}"].mean()) if len(firing) else float("nan")
        rows.append(block)
    summary = pd.DataFrame(rows)
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    summary.to_csv(SUMMARY_DIR / "paperA_moduleD_summary.csv", index=False)
    return summary


def concentration(frame: pd.DataFrame) -> pd.DataFrame:
    """Quantify how narrow the D-eligible population is. This is a limit on
    what the replay can support, and is measured rather than assumed."""
    eligible = frame[(frame["s_max"] == SHIPPED_CAP) & frame["d_eligible"] & (frame["replicate"] == 0)]
    if eligible.empty:
        return pd.DataFrame()
    rows = []
    for keys, label in ((["N", "a"], "instance"), (["window_width"], "window width"),
                        (["s_w_initial"], "initial shots")):
        for key, count in eligible.groupby(keys).size().sort_values(ascending=False).items():
            rows.append({
                "dimension": label,
                "value": str(key),
                "eligible_windows": int(count),
                "share": float(count / len(eligible)),
            })
    return pd.DataFrame(rows)


def figure(frame: pd.DataFrame, summary: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    caps = summary["s_max"].to_list()
    x = np.arange(len(caps))

    ax = axes[0]
    ax.plot(x, summary["cap_hit_fraction_of_eligible"], "o-", color="#c05621",
            label="firing windows that run to the cap")
    ax.plot(x, summary["allocation_uniform_fraction"], "s-", color="#805ad5",
            label="trials whose allocation is uniform")
    ax.set_xticks(x)
    ax.set_xticklabels(caps)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("shot cap $S_{\\max}$")
    ax.set_ylabel("fraction")
    ax.set_title("(a) Does the cap or the statistic set the budget?")
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    styles = {"baseline": ("#718096", "o", "baseline (no extra shots)"),
              "uniform": ("#2f855a", "^", "shots-matched uniform control"),
              "D": ("#2b6cb0", "s", "module D (confidence-gated)")}
    for arm, (colour, marker, label) in styles.items():
        ax.plot(x, summary[f"argmax_correct_{arm}"], marker=marker, ls="-",
                color=colour, label=label, ms=6)
    ax.set_xticks(x)
    ax.set_xticklabels(caps)
    ax.set_xlabel("shot cap $S_{\\max}$")
    ax.set_ylabel("top-1 identification rate")
    ax.set_title("(b) Primary outcome, by arm and cap")
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    path = harness.PLOTS_DIR / "paperA_figA_moduleD_window.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path


def self_test() -> None:
    """Verify the statistical and allocation helpers on constructed inputs."""
    # Exact budget matching, including the indivisible case.
    assert largest_remainder(10, 4) == [3, 3, 2, 2]
    assert sum(largest_remainder(37, 5)) == 37
    assert largest_remainder(0, 3) == [0, 0, 0]

    # Deterministic argmax, independent of insertion order.
    assert observed_argmax({"01": 5, "00": 5}) == "00"
    assert observed_argmax({"00": 5, "01": 5}) == "00"
    assert observed_argmax({"10": 7, "00": 5}) == "10"

    # Top-two extraction, including the single-pattern case module D must
    # handle as n2 = 0.
    assert top_two({"00": 7, "01": 3, "10": 1}) == (7, 3)
    assert top_two({"00": 4}) == (4, 0)

    # Argmax sets and margins, including an exact tie.
    assert true_argmax_set({"00": 0.5, "01": 0.5}) == {"00", "01"}
    assert true_argmax_set({"00": 0.6, "01": 0.4}) == {"00"}
    assert abs(true_margin({"00": 0.6, "01": 0.4}) - 0.2) < 1e-12

    # Total variation: zero against its own law, one against disjoint support.
    assert abs(total_variation({"00": 5, "01": 5}, {"00": 0.5, "01": 0.5})) < 1e-12
    assert abs(total_variation({"00": 4}, {"01": 1.0}) - 1.0) < 1e-12

    # McNemar: discordance drives the test, concordance gives p = 1.
    ref = pd.Series([True, True, False, False])
    assert mcnemar_exact(ref, pd.Series([True, True, False, False])) == (0, 0, 1.0)
    b, c, p = mcnemar_exact(pd.Series([True] * 10), pd.Series([False] * 10))
    assert (b, c) == (10, 0) and p < 0.01

    # Holm is monotone, never reduces a p-value, and caps at 1.
    adjusted = holm([0.01, 0.02, 0.5])
    assert adjusted[0] <= adjusted[1] <= adjusted[2]
    assert all(a >= p for a, p in zip(adjusted, [0.01, 0.02, 0.5]))
    assert all(a <= 1.0 for a in adjusted)

    # Clustered bootstrap: a constant difference has a degenerate interval.
    constant = pd.DataFrame({"cluster": ["t1", "t1", "t2", "t2"], "diff": [0.25] * 4})
    point, lo, hi = clustered_bootstrap(constant, "diff", resamples=200)
    assert abs(point - 0.25) < 1e-12 and abs(hi - lo) < 1e-12

    # The fast restatement must reproduce the shipped statistic across the
    # shot budgets and count splits the replay can reach, ties and extremes
    # included. The full consumed set is re-verified in main().
    spread = [(n1, n2, s_w)
              for s_w in (4, 8, 16, 32, 64, 128, 256)
              for n1 in (s_w, s_w // 2, s_w // 2 + 1, s_w // 4, 1)
              for n2 in (0, 1, s_w // 4, s_w // 2)
              if 0 <= n2 <= n1 <= s_w]
    verify_confidence_equivalence(spread)


def main() -> None:
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)

    print(f"replaying {EXPECTED_TRIALS} trials x {REPLICATES} replicates x {len(CAPS)} caps")
    frame = build_frame()
    print(f"{len(frame)} rows written")

    # Every (n1, n2, s_w) the replay actually consumed, re-evaluated through
    # the shipped statistic. This is the guard that the speed-up above did not
    # change a single stopping decision.
    consumed = sorted(_CONSUMED_TRIPLES)
    worst, disagreements = verify_confidence_equivalence(consumed)
    print(f"C_w equivalence: {len(consumed)} distinct triples consumed, "
          f"max |shipped - fast| = {worst:.3e}, stopping disagreements = {disagreements}\n")

    summary = summarise(frame)
    print("Per-cap summary (rates over all windows of firing trials):")
    print(summary.round(4).to_string(index=False))

    conc = concentration(frame)
    print("\nConcentration of the D-eligible population at the shipped cap:")
    print(conc.round(4).to_string(index=False))

    tests = run_tests(frame)
    print("\nStatistical tests:")
    columns = ["family", "label", "cap", "subset", "n_windows", "n_trials",
               "rate_reference", "rate_arm", "diff_pp", "ci_lo_pp", "ci_hi_pp",
               "b", "c", "mcnemar_p", "holm_p"]
    print(tests[columns].round(4).to_string(index=False))

    print(f"\nfigure: {figure(frame, summary)}")


if __name__ == "__main__":
    self_test()
    main()
