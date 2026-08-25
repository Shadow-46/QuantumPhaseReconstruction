"""
Purpose
    The confidence-guided adaptive reconstruction method: the posterior
    comparison statistic C_w driving three independently switchable modules
    -- B (coverage-based candidate retention), C (ambiguity-scaled beam
    width), and D (confidence-gated shot stopping).
Theory
    This is the method itself, not campaign tooling, and it lives beside
    paper_algorithm.py because it is the adaptive counterpart of that
    module's fixed-budget reconstruction. It deliberately does NOT live
    under Reconstruction/: module D calls back into the sampler, and the
    published adoption protocol requires the shot-stopping loop to sit in a
    layer above the reconstruction functions so those stay pure and
    side-effect-free. Reconstruction/ therefore remains free of any
    Circuits/ or Algorithms/ dependency.

    With use_coverage, use_beam_rule and use_shot_stopping all False,
    reconstruct_adaptive executes a call sequence structurally identical to
    the published fixed-budget pipeline.
Inputs
    ShorConfig and WindowConfig, optional noise model, the three module
    flags, and the parameters of config.AdaptiveReconstructionConfig.
Outputs
    (stitched phases, OrderRecoveryResult | None), plus a realised-resource
    diagnostics dict when return_diagnostics is True.
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Algorithms.paper_algorithm import sample_window_counts
from Circuits.modular_multiplication import ModularMultiplicationOperator
from Circuits.windowed_qpe import run_windowed_qpe_block
from Reconstruction.candidate_generation import (
    generate_window_candidates,
    generate_window_candidates_by_coverage,
    normalize_counts,
)
from Reconstruction.confidence import top1_vs_top2_confidence
from Reconstruction.continued_fraction import recover_order_and_factors
from Reconstruction.stitching import stitch_candidates
from config import DEFAULT_ADAPTIVE, ShorConfig, WindowConfig


def _windows_needing_expansion(window_counts: list[dict[str, int]], candidate_count: int, delta_cov: float) -> int:
    """Return how many windows' fixed top-`candidate_count` truncation covers
    less than `1 - delta_cov` of the observed probability mass -- the signal
    module C's beam-width rule scales on, computed independently of whether
    module B's coverage-based candidate generation is actually enabled (so B
    and C remain separately ablatable, FORMULATION.md Section 5)."""
    target = 1.0 - delta_cov
    expanded = 0
    for counts in window_counts:
        ranked = sorted(normalize_counts(counts).values(), reverse=True)
        if sum(ranked[:candidate_count]) < target:
            expanded += 1
    return expanded


def _resample_shot_stopping(
    shor: ShorConfig,
    window: WindowConfig,
    window_counts: list[dict[str, int]],
    specs,
    epsilon: float,
    s_max: int,
    noise_model=None,
) -> list[dict[str, int]]:
    """Resample any window below the confidence target `1 - epsilon` in
    doubling increments up to a hard cap `s_max` (module D, FORMULATION.md
    Section 4.3). Uses the same operator/eigenstate as the original sample;
    each resample round uses a distinct seed so repeated calls are
    deterministic but not degenerate repeats of the same shots.

    Thresholds on the RAW top1_vs_top2_confidence (C_w), not the isotonic
    recalibration fit in Experiments/phase6_calibration.py. That
    recalibration is known to correct real overconfidence in C_w (Brier
    0.224 -> 0.199, FORMULATION.md Section 6), but is not applied here: this
    is a deliberate scope decision (documented in FORMULATION.md Section 6),
    not an oversight -- swapping in the recalibrated value would change
    every ablation arm that uses D (D, Full), invalidating the already-run
    phase6_ablation.py comparison, which was executed entirely against raw
    C_w. Wiring recalibration into this function is a valid follow-up but
    requires re-running the ablation to attribute any resulting change
    correctly."""
    op = ModularMultiplicationOperator(shor.a, shor.N)
    updated = [dict(c) for c in window_counts]
    for i, spec in enumerate(specs):
        round_index = 0
        while True:
            ranked = sorted(updated[i].values(), reverse=True)
            n1 = ranked[0] if ranked else 0
            n2 = ranked[1] if len(ranked) > 1 else 0
            s_w = sum(updated[i].values())
            if top1_vs_top2_confidence(n1, n2, s_w) >= 1.0 - epsilon:
                break
            extra_shots = min(s_w, s_max - s_w)
            if extra_shots <= 0:
                break
            round_index += 1
            extra_counts = run_windowed_qpe_block(
                op, eigenstate=1, spec=spec, shots=extra_shots,
                method=shor.simulator_method, seed=shor.random_seed + round_index,
                noise_model=noise_model,
            )
            for bits, count in extra_counts.items():
                updated[i][bits] = updated[i].get(bits, 0) + count
    return updated


def reconstruct_adaptive(
    shor: ShorConfig,
    window: WindowConfig,
    noise_model=None,
    use_coverage: bool = False,
    use_beam_rule: bool = False,
    use_shot_stopping: bool = False,
    delta_cov: float = DEFAULT_ADAPTIVE.delta_cov,
    beam_max: int = DEFAULT_ADAPTIVE.beam_max,
    epsilon: float = DEFAULT_ADAPTIVE.epsilon,
    s_max: int = DEFAULT_ADAPTIVE.s_max,
    return_diagnostics: bool = False,
    precomputed: tuple[list[dict[str, int]], list] | None = None,
):
    """Run the classical reconstruction pipeline with modules B (candidate
    coverage), C (beam width), and D (shot stopping) independently
    switchable, for the ablation study (FORMULATION.md Section 5). With all
    three flags False this reproduces Baseline exactly.

    When return_diagnostics is True, also returns a dict of realized
    per-trial resource usage -- m_w (candidates kept per window), max_paths
    (realized beam width P'), and s_w (final shot count per window after
    any D resampling) -- read off values already computed on the existing
    control-flow path. Purely additive: no RNG draw, branch, or existing
    return value changes when this flag is False.

    precomputed, if given, is a (window_counts, specs) pair to reuse instead
    of calling sample_window_counts again. Safe to share the same base
    sample across multiple arms of the same (shor, window) trial: given a
    fixed seed, sample_window_counts is deterministic, so every arm would
    recompute bit-identical base counts anyway (Aer's own seeded sampling,
    Algorithms/paper_algorithm.py:76-92) -- this only skips redundant
    circuit transpile+execute, it does not change any arm's result.
    _resample_shot_stopping already defensively copies
    (`[dict(c) for c in window_counts]`) before mutating, so reuse across
    arms cannot leak D's resampled counts back into a shared base."""
    if precomputed is not None:
        window_counts, specs = precomputed
    else:
        window_counts, specs = sample_window_counts(shor, window, noise_model=noise_model)

    if use_shot_stopping:
        window_counts = _resample_shot_stopping(shor, window, window_counts, specs, epsilon, s_max, noise_model)

    if use_coverage:
        groups = [
            generate_window_candidates_by_coverage(counts, spec.start, spec.width, window.total_precision, delta_cov)
            for counts, spec in zip(window_counts, specs)
        ]
    else:
        groups = [
            generate_window_candidates(counts, spec.start, spec.width, window.total_precision, window.candidate_count)
            for counts, spec in zip(window_counts, specs)
        ]

    if use_beam_rule:
        expanded = _windows_needing_expansion(window_counts, window.candidate_count, delta_cov)
        max_paths = min(beam_max, window.max_paths * (1 + expanded))
    else:
        max_paths = window.max_paths

    stitched = stitch_candidates(groups, max_paths=max_paths)
    result = None
    for phase in stitched:
        recovered = recover_order_and_factors(phase.phase, shor.a, shor.N)
        if recovered.factors is not None:
            result = (stitched, recovered)
            break
    if result is None:
        result = (stitched, None)

    if not return_diagnostics:
        return result

    diagnostics = {
        "m_w": [len(g) for g in groups],
        "max_paths": max_paths,
        "s_w": [sum(c.values()) for c in window_counts],
    }
    return result[0], result[1], diagnostics
