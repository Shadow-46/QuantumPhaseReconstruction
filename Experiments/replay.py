"""
Purpose
    Circuit-free replay of module D's shot-stopping loop against the exact
    per-window marginal, plus the pilot that decides whether a window-level
    module D study can say anything on the released corpus.
Theory
    Every trial in raw_windows/ is noiseless, so the exact statevector
    marginal returned by phase6_calibration.true_window_distribution IS the
    sampling law of the stored counts. Drawing extra shots is therefore a
    multinomial draw from that marginal, distributionally identical to what
    Aer would return, at zero circuit cost. The replay reproduces the LAW, not
    the counts: module D's real loop calls Aer with a fresh seed per round
    (Algorithms/adaptive_reconstruction.py), so no replayed trial is
    bit-identical to a seeded re-run. Every statement built on this module is
    distributional.

    Common random numbers. Each (window, replicate) gets ONE categorical
    stream drawn once; every arm consumes a PREFIX of it sized to its own
    budget. Two consequences: the doubling schedule is internally consistent
    (the counts at S=16 extend those at S=8 rather than being an independent
    draw), and module D differs from a shots-matched control only in how many
    shots each window gets, never in which shots. That removes the dominant
    Monte Carlo term from the comparison the study exists to make.

    What the pilot answers. Module D can only allocate non-uniformly by
    stopping different windows at different rounds. Every window in a trial
    starts at the same S_w, and the shipped cap allows at most a few doubling
    rounds, so if most windows run to the cap then D's realised allocation is
    uniform and a shots-matched uniform control is IDENTICAL to D by
    construction -- the degeneracy that produced (b=0, c=0) for module C in
    the companion study. The pilot measures that fraction before any study is
    built on top of it.
Inputs
    Data/failure_study/raw_windows/*.json (from recapture_window_counts.py)
Outputs
    None. The pilot prints; it writes no file.
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

import numpy as np

from Circuits.windowed_qpe import WindowSpec
from Experiments import harness
from Experiments.phase6_calibration import true_window_distribution
from Reconstruction.confidence import top1_vs_top2_confidence

RAW_WINDOWS_DIR = harness.FAILURE_STUDY_DIR / "raw_windows"

# Shipped configuration (config.AdaptiveReconstructionConfig).
EPSILON = 0.05
S_MAX = 64

# Distinguishes this module's seed derivation from every other consumer of
# grid coordinates in the campaign.
SALT = 0x0DEEDBEE


def marginal_vector(a: int, N: int, spec: WindowSpec) -> tuple[tuple[str, ...], np.ndarray]:
    """Return (patterns, probabilities) for one window's exact marginal, with
    patterns in a fixed sorted order so a stream drawn against it is
    reproducible independent of dict iteration order."""
    dist = true_window_distribution(a, N, spec)
    patterns = tuple(sorted(dist))
    probs = np.array([dist[b] for b in patterns], dtype=float)
    total = probs.sum()
    if total <= 0:
        return patterns, probs
    return patterns, probs / total


def crn_stream(probs: np.ndarray, length: int, entropy: list[int]) -> np.ndarray:
    """One categorical stream of pattern indices, drawn once per (window,
    replicate). Arms consume prefixes of it -- see the module docstring."""
    if length <= 0 or probs.size == 0:
        return np.empty(0, dtype=np.int64)
    rng = np.random.default_rng(np.random.SeedSequence(entropy))
    return rng.choice(probs.size, size=length, p=probs)


def apply_prefix(counts: dict[str, int], patterns: tuple[str, ...],
                 stream: np.ndarray, take: int) -> dict[str, int]:
    """Return `counts` plus the first `take` draws of `stream`."""
    merged = dict(counts)
    if take <= 0:
        return merged
    idx, hits = np.unique(stream[:take], return_counts=True)
    for i, hit in zip(idx, hits):
        bits = patterns[int(i)]
        merged[bits] = merged.get(bits, 0) + int(hit)
    return merged


def replay_module_d(counts: dict[str, int], patterns: tuple[str, ...], stream: np.ndarray,
                    epsilon: float = EPSILON, s_max: int = S_MAX,
                    confidence_fn=top1_vs_top2_confidence) -> tuple[dict[str, int], int, int]:
    """Replay module D on one window. Mirrors _resample_shot_stopping exactly:
    the same doubling schedule `min(s_w, s_max - s_w)`, the same raw-C_w
    threshold, the same n2 = 0 handling when only one pattern is observed, and
    re-ranking after every merge. Returns (final counts, extra shots drawn,
    rounds).

    `confidence_fn` defaults to the shipped statistic and exists only so a
    caller running many replicates can inject a memoised or numerically
    identical restatement of it. Any substitute must be verified against the
    shipped function; see
    paperA_moduleD_validation.verify_confidence_equivalence."""
    current = dict(counts)
    s_w = sum(current.values())
    used = 0
    rounds = 0
    while True:
        ranked = sorted(current.values(), reverse=True)
        n1 = ranked[0] if ranked else 0
        n2 = ranked[1] if len(ranked) > 1 else 0
        if confidence_fn(n1, n2, s_w) >= 1.0 - epsilon:
            break
        extra = min(s_w, s_max - s_w)
        if extra <= 0:
            break
        current = apply_prefix(counts, patterns, stream, used + extra)
        used += extra
        s_w += extra
        rounds += 1
    return current, used, rounds


def trial_allocation(record: dict, trial_ordinal: int, replicate: int,
                     s_max: int = S_MAX) -> list[int]:
    """Extra shots module D draws for each window of one trial, one replicate."""
    spec_dict = record["spec"]
    a, N = spec_dict["a"], spec_dict["N"]
    extras = []
    for window_meta, counts in zip(record["windows"], record["window_counts"]):
        spec = WindowSpec(window_meta["start"], window_meta["width"], window_meta["total_precision"])
        patterns, probs = marginal_vector(a, N, spec)
        s_w = sum(counts.values())
        stream = crn_stream(
            probs, max(0, s_max - s_w),
            [SALT, trial_ordinal, spec.start, spec.width, s_w, replicate],
        )
        _final, used, _rounds = replay_module_d(counts, patterns, stream, s_max=s_max)
        extras.append(used)
    return extras


def pilot(n_trials: int = 20, replicates: int = 3, s_max: int = S_MAX) -> dict:
    """Report the one number that decides the module D study: among trials
    where D actually fires, the fraction whose realised allocation is uniform
    across windows. Where it is uniform, a shots-matched uniform control is
    identical to D by construction and the comparison is empty."""
    paths = sorted(RAW_WINDOWS_DIR.glob("*.json"))
    # Only trials configured below the cap can draw anything at all.
    candidates = []
    for ordinal, path in enumerate(paths):
        record = json.loads(path.read_text())
        if record["spec"]["shots"] < s_max:
            candidates.append((ordinal, record))
    if not candidates:
        return {"eligible_trials": 0}

    step = max(1, len(candidates) // n_trials)
    sample = candidates[::step][:n_trials]

    fired = 0
    uniform = 0
    rows = []
    windows_firing_total = 0
    windows_at_cap = 0
    for ordinal, record in sample:
        shots = record["spec"]["shots"]
        # A window that spends every doubling round without ever clearing the
        # threshold lands exactly on the cap; that is the cap binding, not the
        # confidence rule deciding.
        cap_extra = s_max - shots
        for replicate in range(replicates):
            extras = trial_allocation(record, ordinal, replicate, s_max=s_max)
            active = [e for e in extras if e > 0]
            if not active:
                continue
            fired += 1
            is_uniform = len(set(extras)) == 1
            uniform += int(is_uniform)
            windows_firing_total += len(active)
            windows_at_cap += sum(1 for e in active if e == cap_extra)
            rows.append({
                "trial_id": record["trial_id"],
                "shots": shots,
                "windows": len(extras),
                "windows_firing": len(active),
                "extras": extras,
                "uniform": is_uniform,
            })

    return {
        "trials_sampled": len(sample),
        "replicates": replicates,
        "s_max": s_max,
        "trial_replicates_firing": fired,
        "uniform_allocations": uniform,
        "uniform_fraction": (uniform / fired) if fired else float("nan"),
        "windows_firing": windows_firing_total,
        "windows_at_cap": windows_at_cap,
        "cap_hit_fraction": (windows_at_cap / windows_firing_total) if windows_firing_total else float("nan"),
        "rows": rows,
    }


def self_test() -> None:
    """Verify the replay primitives on constructed inputs."""
    patterns = ("00", "01", "10", "11")
    probs = np.array([0.7, 0.2, 0.06, 0.04])

    # A stream is reproducible from its entropy and only from its entropy.
    a = crn_stream(probs, 50, [SALT, 1, 0, 2, 4, 0])
    b = crn_stream(probs, 50, [SALT, 1, 0, 2, 4, 0])
    c = crn_stream(probs, 50, [SALT, 1, 0, 2, 4, 1])
    assert np.array_equal(a, b)
    assert not np.array_equal(a, c)

    # Prefix semantics: a longer take extends a shorter one, never replaces it.
    base = {"00": 3, "01": 1}
    short = apply_prefix(base, patterns, a, 8)
    long = apply_prefix(base, patterns, a, 16)
    assert sum(short.values()) == sum(base.values()) + 8
    assert sum(long.values()) == sum(base.values()) + 16
    for bits, count in short.items():
        assert long.get(bits, 0) >= count, "prefix growth must be monotone per pattern"

    # A decisive window stops immediately and spends nothing.
    decisive = {"00": 60, "01": 1}
    _f, used, rounds = replay_module_d(decisive, patterns, crn_stream(probs, 3, [SALT, 0, 0, 2, 61, 0]))
    assert used == 0 and rounds == 0

    # A starved near-tie spends, respects the cap, and follows the doubling
    # schedule 4 -> 8 -> 16 -> 32 -> 64, i.e. at most 60 extra shots.
    tie = {"00": 2, "01": 2}
    stream = crn_stream(np.array([0.5, 0.5, 0.0, 0.0]), 60, [SALT, 0, 0, 2, 4, 0])
    final, used, rounds = replay_module_d(tie, patterns, stream)
    assert 0 < used <= 60
    assert sum(final.values()) == 4 + used <= 64
    assert rounds >= 1

    # Zero headroom is a no-op even when confidence is low.
    at_cap = {"00": 32, "01": 32}
    _f, used, _r = replay_module_d(at_cap, patterns, crn_stream(probs, 1, [SALT, 0, 0, 2, 64, 0]))
    assert used == 0


def main() -> None:
    result = pilot()
    print(f"pilot: {result['trials_sampled']} trials x {result['replicates']} replicates "
          f"at S_max = {result['s_max']}")
    print(f"trial-replicates where module D fired: {result['trial_replicates_firing']}")
    print(f"of those, allocation uniform across windows: {result['uniform_allocations']} "
          f"({result['uniform_fraction']:.1%})")
    print(f"firing windows that ran to the cap: {result['windows_at_cap']}/"
          f"{result['windows_firing']} ({result['cap_hit_fraction']:.1%}) "
          f"-- where this is high, the cap decides the budget, not the confidence rule")
    print("\nper trial-replicate:")
    for row in result["rows"]:
        print(f"  shots={row['shots']:>2}  windows={row['windows']:>2}  "
              f"firing={row['windows_firing']:>2}  uniform={str(row['uniform']):>5}  "
              f"extras={row['extras']}")


if __name__ == "__main__":
    self_test()
    main()
