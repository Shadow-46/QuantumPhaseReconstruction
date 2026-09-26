"""Verification of the arXiv-v3 infinite-shot failure floor (PRELIMINARY).

Tasks (one experiment, shards tagged by `task`):
  crosscheck  primary vectorised decoder vs independent string decoder on the
              full infinite-shot grid (both deterministic lowest-index ties).
  sweep       dense epsilon sweep with the primary decoder; observed failures
              vs the per-phase theoretical predicate (TP/FP/FN), for the
              faithful decoder and the special-chunk ablation.
  failures    one row per infinite-shot failure with full block diagnostics.
  variants    every combination of the ambiguity readings in
              verification/awqpe_strings.Variant (216), several epsilons.
  golden      whether each variant still reproduces the paper's Sec. 6 examples.

    python -m research.awqpe.run_verify_epsilon [--pilot] [--max-workers N]
"""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd

from research.awqpe.baseline.awqpe import awqpe_vectorised
from research.awqpe.blocks.geometry import partition_blocks
from research.awqpe.evaluation.metrics import best_nbit_integer, circular_error, exact_success
from research.awqpe.model.kernel import block_probabilities
from research.awqpe.runner.core import common_arguments, resolve_config, run_sharded
from research.awqpe.verification import awqpe_strings as S
from research.awqpe.verification.theory import partition_epsilon_star, predicted_failure

GOLDEN = [
    (0.8203125, [3, 2, 3], "11110010", "11010010"),
    (0.3, [2, 2], "0101", "0101"),
    (float(np.pi / 6), [3, 2, 2, 3], "1000111000", "1000011000"),
    (0.671875, [4, 4], "10111100", "10101100"),
    (float(1 / np.sqrt(2)), [3] * 10, "110101010000010100110011010101", "101101010000010011110011001101"),
    (float(np.sin(np.pi / 12)), [5, 6, 7, 4], "0100001001000010001110", "0100001001000001111110"),
]


def wkey(w) -> str:
    return "-".join(map(str, w))


def grid(widths, refine: int, offset: float) -> np.ndarray:
    G = 1 << (int(sum(widths)) + refine)
    return (np.arange(G) + offset) / G


def probs(phis, widths):
    return [block_probabilities(phis, s.offset, s.width) for s in partition_blocks(widths)]


def all_variants():
    keys = list(S.OPTIONS)
    return [S.Variant(**dict(zip(keys, combo))) for combo in itertools.product(*(S.OPTIONS[k] for k in keys))]


def task_crosscheck(cfg, spec):
    w, off = spec["widths"], spec["offset"]
    n = sum(w)
    phis = grid(w, cfg["grid_refine"], off)
    P = probs(phis, w)
    rows = []
    for eps in cfg["crosscheck_eps"]:
        prim = awqpe_vectorised(P, w, eps)["estimate"]  # zero jitter: lowest index wins exact ties
        for t, phi in enumerate(phis):
            r = S.decode([p[t] for p in P], w, eps, S.Variant(ties="low"))
            rows.append((phi, eps, int(prim[t]), int(r["est"], 2)))
    df = pd.DataFrame(rows, columns=["phi", "epsilon", "est_primary", "est_independent"])
    df["agree"] = df.est_primary == df.est_independent
    df["widths"], df["grid_offset"], df["n"], df["task"] = wkey(w), off, n, "crosscheck"
    return df


def task_sweep(cfg, spec):
    w, off = spec["widths"], spec["offset"]
    n = sum(w)
    phis = grid(w, cfg["grid_refine_sweep"].get(str(n), cfg["grid_refine"]), off)
    P = probs(phis, w)
    rows = []
    for eps in np.round(np.arange(cfg["eps_lo"], cfg["eps_hi"] + 1e-9, cfg["eps_step"]), 4):
        for name, rule in (("awqpe", True), ("awqpe_ablate_special", False)):
            est = awqpe_vectorised(P, w, float(eps), special_chunk_rule=rule)["estimate"]
            obs = ~exact_success(est, phis, n)
            pred, _ = predicted_failure(phis, w, float(eps), special_rule=rule)
            rows.append({"decoder": name, "epsilon": float(eps), "n_phases": len(phis), "observed_fail": int(obs.sum()),
                         "predicted_fail": int(pred.sum()), "tp": int((obs & pred).sum()), "fp": int((~obs & pred).sum()),
                         "fn": int((obs & ~pred).sum()), "max_error_lsb": float((circular_error(est / 2**n, phis) * 2**n).max())})
    df = pd.DataFrame(rows)
    df["widths"], df["grid_offset"], df["n"], df["task"] = wkey(w), off, n, "sweep"
    df["partition_eps_star"] = partition_epsilon_star(w)
    return df


def task_failures(cfg, spec):
    w, off = spec["widths"], spec["offset"]
    n = sum(w)
    phis = grid(w, cfg["grid_refine"], off)
    P = probs(phis, w)
    frames = []
    for eps in cfg["failure_eps"]:
        r = awqpe_vectorised(P, w, eps)
        est = r["estimate"]
        fail = ~exact_success(est, phis, n)
        _, pblock = predicted_failure(phis, w, eps)
        truth = best_nbit_integer(phis, n)
        for t in np.flatnonzero(fail):
            true_bits, est_bits = format(int(truth[t]), f"0{n}b"), format(int(est[t]), f"0{n}b")
            pos, first_bad = 0, 0
            for j, m in enumerate(w, start=1):  # first chunk whose decoded bits differ from the best n-bit value
                if est_bits[pos:pos + m] != true_bits[pos:pos + m]:
                    first_bad = j
                    break
                pos += m
            b = int(pblock[t]) or first_bad
            pb = P[b - 1][t]
            t1, t2 = int(r["top1"][t, b - 1]), int(r["top2"][t, b - 1])
            frames.append({
                "phi": float(phis[t]), "widths": wkey(w), "n": n, "epsilon": eps, "grid_offset": off,
                "affected_block": b, "predicted_block": int(pblock[t]), "first_wrong_chunk": first_bad,
                "block_probabilities": ";".join(f"{v:.6f}" for v in pb),
                "top1": t1, "top2": t2, "p_top1": float(pb[t1]), "p_top2": float(pb[t2]),
                "ratio_c2_c1": float(pb[t2] / pb[t1]), "flag": bool(r["flags"][t, b - 1]),
                "special_index": int(r["special_index"][t]),
                "raw_bits": "".join(format(int(x), f"0{m}b") for x, m in zip(r["raw"][t], w)),
                "decoded_bits": est_bits, "true_nearest_bits": true_bits,
                "true_nearest_phase": int(truth[t]) / 2**n, "reconstructed_phase": int(est[t]) / 2**n,
                "circular_error": float(circular_error(est[t] / 2**n, phis[t])),
                "absolute_error": float(abs(est[t] / 2**n - phis[t])),
                "exact_tie_case": bool(np.any(np.isclose(np.mod(phis[t] * 2.0 ** np.arange(1, n + 2), 1.0), 0.5, atol=1e-12))),
            })
    df = pd.DataFrame(frames)
    if df.empty:
        df = pd.DataFrame({"widths": [wkey(w)], "n": [n]})
    df["task"] = "failures"
    return df


def task_variants(cfg, spec):
    w, off = spec["widths"], spec["offset"]
    n = sum(w)
    phis = grid(w, cfg["grid_refine"], off)
    P = [p.tolist() for p in probs(phis, w)]
    rows = []
    for v in all_variants()[spec["v_lo"]:spec["v_hi"]]:
        for eps in cfg["variant_eps"]:
            fails, big = 0, 0
            for t in range(len(phis)):
                est = int(S.decode([p[t] for p in P], w, eps, v, seed=t)["est"], 2)
                if not exact_success(est, phis[t], n):
                    fails += 1
                    big += int(circular_error(est / 2**n, phis[t]) > 2.0 ** -n)
            rows.append({"variant": v.label, **{f: getattr(v, f) for f in S.OPTIONS}, "epsilon": eps,
                         "n_phases": len(phis), "failures": fails, "failures_beyond_1lsb": big})
    df = pd.DataFrame(rows)
    df["widths"], df["grid_offset"], df["n"], df["task"] = wkey(w), off, n, "variants"
    return df


def task_golden(cfg, spec):
    rows = []
    for v in all_variants():
        ok_all = True
        for phi, w, raw, fin in GOLDEN:
            P = [p.tolist() for p in probs(np.array([phi]), w)]
            r = S.decode([p[0] for p in P], w, 0.9, v)
            ok_all &= (r["raw"] == raw) and (r["est"] == fin)
        rows.append({"variant": v.label, **{f: getattr(v, f) for f in S.OPTIONS}, "golden_all_reproduce": ok_all})
    df = pd.DataFrame(rows)
    df["task"] = "golden"
    return df


TASKS = {"crosscheck": task_crosscheck, "sweep": task_sweep, "failures": task_failures, "variants": task_variants, "golden": task_golden}


def verify_shard(cfg: dict, spec: dict) -> pd.DataFrame:
    return TASKS[spec["task"]](cfg, spec)


def shards(cfg):
    out = []
    for w in cfg["crosscheck_partitions"]:
        for off in cfg["offsets"]:
            out.append({"shard_id": f"cross_w{wkey(w)}_o{off}", "task": "crosscheck", "widths": w, "offset": off})
    for w in cfg["sweep_partitions"]:
        for off in cfg["offsets"]:
            out.append({"shard_id": f"sweep_w{wkey(w)}_o{off}", "task": "sweep", "widths": w, "offset": off})
            out.append({"shard_id": f"fail_w{wkey(w)}_o{off}", "task": "failures", "widths": w, "offset": off})
    nv, step = len(all_variants()), int(cfg["variants_per_shard"])
    for w in cfg["variant_partitions"]:
        for off in cfg["offsets"]:
            for lo in range(0, nv, step):
                out.append({"shard_id": f"var_w{wkey(w)}_o{off}_v{lo}", "task": "variants", "widths": w, "offset": off, "v_lo": lo, "v_hi": lo + step})
    out.append({"shard_id": "golden", "task": "golden"})
    return out


def main() -> None:
    args = common_arguments(__doc__, "verify_epsilon.yaml").parse_args()
    cfg = resolve_config(args)
    run_sharded(args, cfg, shards(cfg), verify_shard, extra_manifest={"status": "PRELIMINARY (arXiv-v3 reading)"})


if __name__ == "__main__":
    main()
