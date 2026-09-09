"""
Purpose
    Two synthetic-grid analyses that bound what the released-corpus replay
    (paperA_moduleD_validation.py) can and cannot say about module D.

      Family A -- REALISABLE. All 39 observed trial configurations, across
      all three instances and window widths 2/3/4/6, run from freshly drawn
      counts. Asks whether the replay's conclusion generalises beyond the
      dominant instance and width.

      Family B -- CONSTRUCTED MECHANISM PROBE. Ensembles assembled from the
      83 resolvable-ambiguous marginals of (N=21, a=19) and the 46 decisive
      marginals of (N=15, a=4) and (N=21, a=8), at controlled ambiguity
      concentration. Asks whether placement information can matter at all
      when ambiguity is concentrated in a minority of windows.
Theory
    Inventorying all 156 exact marginals shows ambiguity and instance are
    almost perfectly confounded. Margins take three values: 0.0 (27 windows,
    EXACT ties -- both patterns sit in the true argmax set, so a top-1
    endpoint is unconditionally satisfied and extra shots cannot help),
    ~0.106-0.125 (83 windows, resolvable, ALL from N=21 a=19), and 1.0 (46
    windows, decisive, ALL from the other two instances). There is therefore
    no within-instance way to vary ambiguity concentration, and the exact
    ties cannot serve as an ambiguous pool.

    Family B consequently mixes instances. It is NOT a physical trial, NOT an
    end-to-end experiment, and NOT a general validation of module D: no
    quantum circuit produces such a window set. It is a probe of the
    allocation rule's behaviour under a controlled ambiguity profile.

    A tautology warning belongs with Family B's headline number. At low
    ambiguity fractions the ensemble is "windows that need shots" plus
    "windows that need none", so module D beats a uniform split largely
    because the uniform split wastes budget on windows already at probability
    one. The informative content is the SHAPE of the effect against
    concentration -- where the crossover sits and how steep it is -- not the
    value at the extreme.

    Nothing here is confirmatory. No p-values are computed and no result may
    revise the replay study's pre-declared conclusion, which stands as
    measured: module D does not significantly outperform its shots-matched
    uniform control at the shipped cap, and the cap sweep does not reveal a
    hidden placement advantage.
Inputs
    Data/failure_study/raw_windows/*.json (window geometry and configurations)
Outputs
    Data/failure_study/paperA_synthetic_grid_A.csv.gz
    Data/failure_study/paperA_synthetic_grid_B.csv.gz
    Data/failure_study/summary_tables/paperA_synthetic_A_summary.csv
    Data/failure_study/summary_tables/paperA_synthetic_B_summary.csv
    Data/failure_study/plots/paperA_figS_synthetic.png
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
import numpy as np
import pandas as pd

from Circuits.windowed_qpe import WindowSpec
from Experiments import harness
from Experiments.phase6_calibration import true_window_distribution
from Experiments.paperA_moduleD_validation import (
    EPSILON,
    cached_marginal,
    clustered_bootstrap,
    confidence,
    largest_remainder,
    outcomes,
    true_argmax_set,
)
from Experiments.replay import RAW_WINDOWS_DIR, apply_prefix, crn_stream, replay_module_d

SUMMARY_DIR = harness.FAILURE_STUDY_DIR / "summary_tables"

# Distinct salts so the two families never share a stream.
SALT_A = 0x0A17E5
SALT_B = 0x0B17E5

CAP = 64                          # shipped S_max
CAP_SENSITIVITY = 256
S0_FAMILY_A = (4, 8, 16, 32)
S0_FAMILY_B = (4, 8, 16)
REPLICATES = 20

ENSEMBLE_SIZE = 8                 # windows per constructed ensemble
AMBIGUOUS_COUNTS = (1, 2, 4, 8)   # concentration levels out of ENSEMBLE_SIZE
ENSEMBLES_PER_CELL = 30

# Margin class boundaries, read off the measured trimodal distribution
# (exactly 0.0, ~0.106-0.125, exactly 1.0) rather than tuned.
TIE_TOL = 0.01
DECISIVE_MIN = 0.5

# The synthetic grid has no configured k, so the fixed-top-k secondary uses
# one value throughout, making it comparable across every cell.
K_FIXED = 2

EXPECTED_MARGINALS = 156
EXPECTED_CONFIGS = 39
EXPECTED_RESOLVABLE = 83
EXPECTED_DECISIVE = 46
EXPECTED_TIES = 27


# --------------------------------------------------------------------------
# Marginal inventory and pools
# --------------------------------------------------------------------------

def inventory() -> tuple[dict, list]:
    """Return (configuration -> window list, sorted marginal keys)."""
    configs: dict[tuple, list] = {}
    keys: dict[tuple, None] = {}
    for path in sorted(RAW_WINDOWS_DIR.glob("*.json")):
        record = json.loads(path.read_text())
        s = record["spec"]
        cfg = (s["N"], s["a"], s["total_precision"], s["window_size"], s["overlap"])
        if cfg not in configs:
            configs[cfg] = [(w["start"], w["width"], w["total_precision"])
                            for w in record["windows"]]
        for w in record["windows"]:
            keys[(s["a"], s["N"], w["start"], w["width"], w["total_precision"])] = None
    return configs, sorted(keys)


def classify(key: tuple) -> dict:
    a, N, start, width, precision = key
    dist = true_window_distribution(a, N, WindowSpec(start, width, precision))
    ordered = sorted(dist.values(), reverse=True)
    top1 = ordered[0]
    top2 = ordered[1] if len(ordered) > 1 else 0.0
    margin = top1 - top2
    if margin < TIE_TOL:
        cls = "tie"
    elif margin >= DECISIVE_MIN:
        cls = "decisive"
    else:
        cls = "resolvable"
    return {"key": key, "a": a, "N": N, "start": start, "width": width,
            "precision": precision, "margin": margin, "top1": top1,
            "class": cls, "support": len(dist)}


def build_pools() -> tuple[dict, dict, list, list, list]:
    configs, keys = inventory()
    assert len(keys) == EXPECTED_MARGINALS, f"{len(keys)} marginals, expected {EXPECTED_MARGINALS}"
    assert len(configs) == EXPECTED_CONFIGS, f"{len(configs)} configs, expected {EXPECTED_CONFIGS}"
    scored = {k: classify(k) for k in keys}
    resolvable = [r for r in scored.values() if r["class"] == "resolvable"]
    decisive = [r for r in scored.values() if r["class"] == "decisive"]
    ties = [r for r in scored.values() if r["class"] == "tie"]
    assert len(resolvable) == EXPECTED_RESOLVABLE, f"{len(resolvable)} resolvable"
    assert len(decisive) == EXPECTED_DECISIVE, f"{len(decisive)} decisive"
    assert len(ties) == EXPECTED_TIES, f"{len(ties)} ties"
    # The confound this design is built around: asserted, not assumed, so a
    # corpus change cannot silently invalidate the framing.
    assert {(r["N"], r["a"]) for r in resolvable} == {(21, 19)}
    assert {(r["N"], r["a"]) for r in decisive} == {(15, 4), (21, 8)}
    return configs, scored, resolvable, decisive, ties


# --------------------------------------------------------------------------
# One window under all arms
# --------------------------------------------------------------------------

def prepare_window(entry: dict, s0: int, cap: int, entropy: list[int]) -> dict:
    """Draw a window's initial counts and module D's response from one common
    random number stream. The stream's first s0 draws are the baseline counts;
    every arm then consumes a prefix of the remainder, so the arms differ only
    in how many extra shots a window receives, never in which."""
    a, N, start, width, precision = entry["key"]
    spec = WindowSpec(start, width, precision)
    true_dist = true_window_distribution(a, N, spec)
    patterns, probs = cached_marginal(a, N, spec)
    stream = crn_stream(probs, cap, entropy)
    base_counts = apply_prefix({}, patterns, stream, s0)
    extras_stream = stream[s0:]
    final_d, extra_d, rounds_d = replay_module_d(
        base_counts, patterns, extras_stream, epsilon=EPSILON, s_max=cap,
        confidence_fn=confidence)
    return {"entry": entry, "spec": spec, "true_dist": true_dist,
            "patterns": patterns, "base_counts": base_counts,
            "extras_stream": extras_stream, "final_d": final_d,
            "extra_d": extra_d, "rounds_d": rounds_d, "s0": s0}


def score_arms(prepared: list[dict], allocations: dict[str, list[int]],
               cap: int) -> list[dict]:
    """Score every arm on every window of one ensemble or configuration."""
    rows = []
    for index, w in enumerate(prepared):
        true_dist = w["true_dist"]
        argmax_set = true_argmax_set(true_dist)
        row = {
            "window_index": index,
            "N": w["entry"]["N"], "a": w["entry"]["a"],
            "window_start": w["entry"]["start"],
            "window_width": w["entry"]["width"],
            "total_precision": w["entry"]["precision"],
            "margin_class": w["entry"]["class"],
            "true_margin": w["entry"]["margin"],
            "true_argmax_tied": len(argmax_set) > 1,
            "s0": w["s0"], "s_max": cap,
            "extra_D": w["extra_d"], "rounds_D": w["rounds_d"],
            "d_eligible": w["extra_d"] > 0,
            "hit_cap_D": bool(w["extra_d"] > 0 and w["s0"] + w["extra_d"] >= cap),
        }
        row.update(outcomes(w["base_counts"], true_dist, argmax_set,
                            w["spec"], K_FIXED, "baseline"))
        row.update(outcomes(w["final_d"], true_dist, argmax_set,
                            w["spec"], K_FIXED, "D"))
        for name, alloc in allocations.items():
            counts = apply_prefix(w["base_counts"], w["patterns"],
                                  w["extras_stream"], alloc[index])
            row[f"extra_{name}"] = alloc[index]
            row.update(outcomes(counts, true_dist, argmax_set,
                                w["spec"], K_FIXED, name))
        rows.append(row)
    return rows


def uniform_allocation(prepared: list[dict]) -> list[int]:
    total = sum(w["extra_d"] for w in prepared)
    return largest_remainder(total, len(prepared))


def anti_allocation(prepared: list[dict], cap: int) -> tuple[list[int], bool]:
    """Spend module D's own budget on the windows D funded LEAST. This is the
    worst-placement bracket: with D above and this below, the pair bounds how
    much the placement decision can be worth at all. A bracket, not a result."""
    total = sum(w["extra_d"] for w in prepared)
    size = len(prepared)
    order = sorted(range(size), key=lambda i: (prepared[i]["extra_d"], i))
    bottom = order[: max(1, size // 2)]
    headroom = [cap - prepared[i]["s0"] for i in bottom]
    capacity = sum(headroom)
    saturated = total > capacity
    spend = min(total, capacity)
    alloc = [0] * size
    give = largest_remainder(spend, len(bottom))
    for slot, (i, amount) in enumerate(zip(bottom, give)):
        alloc[i] = min(amount, headroom[slot])
    return alloc, saturated


# --------------------------------------------------------------------------
# Family A -- realisable configurations
# --------------------------------------------------------------------------

def run_family_a(configs: dict, scored: dict,
                 caps: tuple[int, ...] = (CAP, CAP_SENSITIVITY)) -> pd.DataFrame:
    rows = []
    for cap in caps:
        for s0 in S0_FAMILY_A:
            if s0 >= cap:
                continue
            for cfg_index, (cfg, windows) in enumerate(sorted(configs.items())):
                N, a, precision, window_size, overlap = cfg
                for replicate in range(REPLICATES):
                    prepared = [
                        prepare_window(scored[(a, N, start, width, prec)], s0, cap,
                                       [SALT_A, cfg_index, w_index, s0, cap, replicate])
                        for w_index, (start, width, prec) in enumerate(windows)
                    ]
                    allocations = {"uniform": uniform_allocation(prepared)}
                    total_extra = sum(w["extra_d"] for w in prepared)
                    for row in score_arms(prepared, allocations, cap):
                        row.update({
                            "family": "A",
                            "config_id": f"N{N}_a{a}_p{precision}_w{window_size}_o{overlap}",
                            "config_windows": len(windows),
                            "overlap": overlap,
                            "replicate": replicate,
                            "trial_fired": total_extra > 0,
                            "allocation_uniform": len({w["extra_d"] for w in prepared}) == 1,
                        })
                        rows.append(row)
        print(f"  family A: cap {cap} done")
    frame = pd.DataFrame(rows)
    frame.to_csv(harness.FAILURE_STUDY_DIR / "paperA_synthetic_grid_A.csv.gz",
                 index=False, compression="gzip")
    return frame


# --------------------------------------------------------------------------
# Family B -- constructed concentration ensembles
# --------------------------------------------------------------------------

def run_family_b(resolvable: list, decisive: list) -> pd.DataFrame:
    rows = []
    for n_ambiguous in AMBIGUOUS_COUNTS:
        for s0 in S0_FAMILY_B:
            for ensemble in range(ENSEMBLES_PER_CELL):
                picker = np.random.default_rng(
                    np.random.SeedSequence([SALT_B, n_ambiguous, s0, ensemble]))
                amb_idx = picker.choice(len(resolvable), size=n_ambiguous, replace=False)
                dec_idx = picker.choice(len(decisive),
                                        size=ENSEMBLE_SIZE - n_ambiguous, replace=False)
                members = ([resolvable[i] for i in amb_idx]
                           + [decisive[i] for i in dec_idx])
                for replicate in range(REPLICATES):
                    prepared = [
                        prepare_window(entry, s0, CAP,
                                       [SALT_B, n_ambiguous, s0, ensemble, w_index, replicate])
                        for w_index, entry in enumerate(members)
                    ]
                    anti, saturated = anti_allocation(prepared, CAP)
                    allocations = {"uniform": uniform_allocation(prepared), "anti": anti}
                    total_extra = sum(w["extra_d"] for w in prepared)
                    for row in score_arms(prepared, allocations, CAP):
                        row.update({
                            "family": "B",
                            "ensemble_id": f"f{n_ambiguous}_s{s0}_e{ensemble}",
                            "n_ambiguous": n_ambiguous,
                            "ambiguous_fraction": n_ambiguous / ENSEMBLE_SIZE,
                            "replicate": replicate,
                            "trial_fired": total_extra > 0,
                            "anti_saturated": saturated,
                            "allocation_uniform": len({w["extra_d"] for w in prepared}) == 1,
                        })
                        rows.append(row)
        print(f"  family B: concentration {n_ambiguous}/{ENSEMBLE_SIZE} done")
    frame = pd.DataFrame(rows)
    frame.to_csv(harness.FAILURE_STUDY_DIR / "paperA_synthetic_grid_B.csv.gz",
                 index=False, compression="gzip")
    return frame


# --------------------------------------------------------------------------
# Summaries -- effect sizes with clustered intervals, no hypothesis tests
# --------------------------------------------------------------------------

def paired_effect(frame: pd.DataFrame, reference: str, arm: str, cluster: str,
                  outcome: str = "argmax_correct") -> tuple[float, float, float]:
    """Replicate-averaged paired difference with a bootstrap interval over
    clusters. Reported as an effect size only: this family is characterisation,
    so no p-value is computed and no hypothesis is tested."""
    averaged = (frame.groupby([cluster, "window_index"], sort=True)
                .agg(ref=(f"{outcome}_{reference}", "mean"),
                     arm=(f"{outcome}_{arm}", "mean"))
                .reset_index())
    averaged["diff"] = averaged["arm"] - averaged["ref"]
    averaged["cluster"] = averaged[cluster]
    return clustered_bootstrap(averaged, "diff")


def summarise_family_a(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (cap, s0), sub in frame.groupby(["s_max", "s0"]):
        blocks = [("all instances", sub)]
        blocks += [(f"N={n} a={a}", g) for (n, a), g in sub.groupby(["N", "a"])]
        blocks += [(f"width {w}", g) for w, g in sub.groupby("window_width")]
        for label, group in blocks:
            block = {
                "s_max": cap, "s0": s0, "group": label,
                "windows": len(group) // REPLICATES,
                "eligible_fraction": float(group["d_eligible"].mean()),
                "mean_extra_D": float(group["extra_D"].mean()),
                "baseline": float(group["argmax_correct_baseline"].mean()),
                "D": float(group["argmax_correct_D"].mean()),
                "uniform": float(group["argmax_correct_uniform"].mean()),
                "tv_baseline": float(group["tv_baseline"].mean()),
                "tv_D": float(group["tv_D"].mean()),
                "tv_uniform": float(group["tv_uniform"].mean()),
            }
            fired = group[group["trial_fired"]]
            if len(fired):
                for ref, arm, name in (("baseline", "D", "D_vs_baseline"),
                                       ("uniform", "D", "D_vs_uniform"),
                                       ("baseline", "uniform", "uniform_vs_baseline")):
                    point, lo, hi = paired_effect(fired, ref, arm, "config_id")
                    block[f"{name}_pp"] = 100 * point
                    block[f"{name}_lo"] = 100 * lo
                    block[f"{name}_hi"] = 100 * hi
            rows.append(block)
    summary = pd.DataFrame(rows)
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    summary.to_csv(SUMMARY_DIR / "paperA_synthetic_A_summary.csv", index=False)
    return summary


def summarise_family_b(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (n_amb, s0), sub in frame.groupby(["n_ambiguous", "s0"]):
        ambiguous_only = sub[sub["margin_class"] == "resolvable"]
        block = {
            "n_ambiguous": n_amb,
            "ambiguous_fraction": n_amb / ENSEMBLE_SIZE,
            "s0": s0,
            "ensembles": sub["ensemble_id"].nunique(),
            "eligible_fraction": float(sub["d_eligible"].mean()),
            "mean_extra_D": float(sub["extra_D"].mean()),
            "mean_extra_uniform": float(sub["extra_uniform"].mean()),
            "anti_saturated": float(sub["anti_saturated"].mean()),
        }
        for arm in ("baseline", "D", "uniform", "anti"):
            block[f"ensemble_{arm}"] = float(sub[f"argmax_correct_{arm}"].mean())
            block[f"ambiguous_{arm}"] = (float(ambiguous_only[f"argmax_correct_{arm}"].mean())
                                         if len(ambiguous_only) else float("nan"))
        for ref, arm, name in (("baseline", "D", "D_vs_baseline"),
                               ("uniform", "D", "D_vs_uniform"),
                               ("baseline", "uniform", "uniform_vs_baseline"),
                               ("anti", "D", "D_vs_anti")):
            point, lo, hi = paired_effect(sub, ref, arm, "ensemble_id")
            block[f"{name}_pp"] = 100 * point
            block[f"{name}_lo"] = 100 * lo
            block[f"{name}_hi"] = 100 * hi
        rows.append(block)
    summary = pd.DataFrame(rows)
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    summary.to_csv(SUMMARY_DIR / "paperA_synthetic_B_summary.csv", index=False)
    return summary


def figure(summary_a: pd.DataFrame, summary_b: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

    ax = axes[0]
    sub = summary_a[(summary_a["s_max"] == CAP) & summary_a["group"].str.startswith("N=")]
    for group in sorted(sub["group"].unique()):
        g = sub[sub["group"] == group].sort_values("s0")
        ax.plot(np.arange(len(g)), 100 * (g["D"] - g["uniform"]), "o-", label=group, ms=5)
    ax.axhline(0, color="gray", ls="--", lw=1)
    ax.set_xticks(np.arange(len(S0_FAMILY_A)))
    ax.set_xticklabels(S0_FAMILY_A)
    ax.set_xlabel("initial shots per window $S_0$")
    ax.set_ylabel("D minus uniform control (pp)")
    ax.set_title("(a) Family A: realisable configurations")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    for s0 in S0_FAMILY_B:
        g = summary_b[summary_b["s0"] == s0].sort_values("ambiguous_fraction")
        ax.errorbar(g["ambiguous_fraction"], g["D_vs_uniform_pp"],
                    yerr=[g["D_vs_uniform_pp"] - g["D_vs_uniform_lo"],
                          g["D_vs_uniform_hi"] - g["D_vs_uniform_pp"]],
                    marker="o", capsize=3, label=f"$S_0$ = {s0}", ms=5)
    ax.axhline(0, color="gray", ls="--", lw=1)
    ax.set_xlabel("fraction of windows that are ambiguous")
    ax.set_ylabel("D minus uniform control (pp)")
    ax.set_title("(b) Family B: constructed concentration probe")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    path = harness.PLOTS_DIR / "paperA_figS_synthetic.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path


def self_test() -> None:
    """Verify the allocation helpers on constructed inputs."""
    fake = [{"extra_d": e, "s0": 4} for e in (60, 0, 0, 12)]
    assert sum(uniform_allocation(fake)) == 72
    assert uniform_allocation(fake) == [18, 18, 18, 18]

    # The anti-control spends D's total on the windows D funded least, and
    # never on the window D funded most.
    alloc, saturated = anti_allocation(fake, CAP)
    assert sum(alloc) == 72 and not saturated
    assert alloc[0] == 0, "the most-funded window must receive nothing"
    assert alloc[1] > 0 and alloc[2] > 0

    # Saturation is detected rather than silently over-spending the cap.
    tight = [{"extra_d": 60, "s0": 4} for _ in range(4)]
    alloc, saturated = anti_allocation(tight, CAP)
    assert saturated
    assert all(a <= CAP - 4 for a in alloc)

    # Class boundaries bracket the measured trimodal margin values.
    assert TIE_TOL < 0.106 < DECISIVE_MIN < 1.0


def main() -> None:
    pd.set_option("display.width", 260)
    pd.set_option("display.max_columns", 40)

    configs, scored, resolvable, decisive, ties = build_pools()
    print(f"pools: {len(resolvable)} resolvable (ambiguous), {len(decisive)} decisive, "
          f"{len(ties)} exact ties excluded from both")
    print(f"       {len(configs)} realisable configurations\n")

    frame_a = run_family_a(configs, scored)
    print(f"family A: {len(frame_a):,} rows\n")
    summary_a = summarise_family_a(frame_a)
    print("Family A -- realisable configurations (primary endpoint, cap 64):")
    cols_a = ["s0", "group", "windows", "eligible_fraction", "mean_extra_D",
              "baseline", "uniform", "D", "D_vs_uniform_pp", "D_vs_uniform_lo",
              "D_vs_uniform_hi", "uniform_vs_baseline_pp"]
    print(summary_a[summary_a["s_max"] == CAP][cols_a].round(4).to_string(index=False))

    frame_b = run_family_b(resolvable, decisive)
    print(f"\nfamily B: {len(frame_b):,} rows\n")
    summary_b = summarise_family_b(frame_b)
    print("Family B -- constructed concentration probe (NOT a physical trial):")
    cols_b = ["ambiguous_fraction", "s0", "ensembles", "eligible_fraction",
              "mean_extra_D", "ensemble_baseline", "ensemble_uniform", "ensemble_D",
              "ensemble_anti", "ambiguous_D", "ambiguous_uniform",
              "D_vs_uniform_pp", "D_vs_uniform_lo", "D_vs_uniform_hi",
              "D_vs_anti_pp", "uniform_vs_baseline_pp"]
    print(summary_b[cols_b].round(4).to_string(index=False))

    print(f"\nfigure: {figure(summary_a, summary_b)}")


if __name__ == "__main__":
    self_test()
    main()
