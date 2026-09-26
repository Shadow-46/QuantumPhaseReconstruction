"""Phase 2: faithful AWQPE baseline characterisation (no adaptation).

Two experiments share this entry point:

  --mode limit   p2a_infinite_shot_limit: decode the exact outcome
                 distributions (counts = probabilities) on a dense phase grid
                 to test Theorem 3.7 empirically in the infinite-shot limit.
  --mode mc      p2b_awqpe_baseline_mc: finite-shot Monte Carlo over
                 partitions x shots/block x epsilon x decoder on stratified
                 phases. Counts and tie-breaking jitter are shared by every
                 epsilon and both decoders within a trial (paired design).

Usage:
    python -m research.awqpe.run_baseline --mode limit
    python -m research.awqpe.run_baseline --mode mc --pilot
    python -m research.awqpe.run_baseline --mode mc --max-workers 20
    python -m research.awqpe.run_baseline --mode mc --resume <run_dir>
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from research.awqpe.baseline.awqpe import awqpe_vectorised
from research.awqpe.blocks.geometry import partition_blocks
from research.awqpe.decode.likelihood import likelihood_decode
from research.awqpe.evaluation.metrics import boundary_hardness, circular_error, exact_success, final_residual
from research.awqpe.model.kernel import block_delta, block_probabilities
from research.awqpe.model.noise import NoiseSpec
from research.awqpe.phases import make_phase_table
from research.awqpe.runner.core import common_arguments, resolve_config, run_sharded
from research.awqpe.sim.oracle_sim import batch_block_counts
from research.awqpe.sim.seeds import generator, stable_int


# "awqpe" is the published algorithm. "awqpe_ablate_special" disables only the
# special-chunk borrow suppression (Alg. 2 line 17); it is a diagnostic
# ablation (docs/DECISIONS.md D-007), never reported as AWQPE.
DECODER_VARIANTS = (("awqpe", True), ("awqpe_ablate_special", False))


def widths_key(widths) -> str:
    return "-".join(str(int(m)) for m in widths)


def true_block_outcomes(phis: np.ndarray, widths) -> list[np.ndarray]:
    """T1 per chunk: the most probable outcome floor(M delta + 0.5) mod M (evaluator only)."""
    out = []
    for s in partition_blocks(widths):
        out.append(np.mod(np.floor(block_delta(phis, s.offset) * s.M + 0.5), s.M).astype(np.int64))
    return out


def _bitmask(bool_matrix: np.ndarray) -> np.ndarray:
    weights = 1 << np.arange(bool_matrix.shape[1], dtype=np.int64)
    return (bool_matrix.astype(np.int64) * weights).sum(axis=1)


def decode_all(counts, widths, epsilons, jitter, refine, noise):
    """Run D1 (AWQPE) for every epsilon and D2 (likelihood) once on shared counts."""
    n = int(sum(widths))
    specs = partition_blocks(widths)
    results = []
    for eps in epsilons:
        for name, rule in DECODER_VARIANTS:
            r = awqpe_vectorised(counts, widths, float(eps), jitter=jitter, special_chunk_rule=rule)
            results.append((name, float(eps), r["estimate"], r))
    T = counts[0].shape[0]
    phi_hat = np.concatenate([likelihood_decode(specs, [c[a:a + 2048] for c in counts], n, refine, noise) for a in range(0, T, 2048)])
    est = np.mod(np.floor(phi_hat * (1 << n) + 0.5), 1 << n).astype(np.int64)
    results.append(("likelihood", np.nan, est, {"phi_hat": phi_hat}))
    return results


def _rows(meta: pd.DataFrame, widths, results, t1_true) -> pd.DataFrame:
    n = int(sum(widths))
    B = len(widths)
    frames = []
    for decoder, eps, est, r in results:
        phi_hat = r["phi_hat"] if decoder == "likelihood" else est / float(1 << n)
        df = meta.copy()
        df["decoder"] = decoder
        df["epsilon"] = eps
        df["estimate_int"] = est
        df["phi_hat"] = phi_hat
        df["error"] = circular_error(phi_hat, meta["phi"].to_numpy())
        df["exact"] = exact_success(est, meta["phi"].to_numpy(), n)
        if decoder.startswith("awqpe"):
            df["flags_mask"] = _bitmask(r["flags"])
            df["t1_correct_mask"] = _bitmask(np.stack([r["top1"][:, i] == t1_true[i] for i in range(B)], axis=1))
            df["special_index"] = r["special_index"]
            df["min_ratio"] = r["ratio"].min(axis=1)
            df["max_ratio"] = r["ratio"].max(axis=1)
        else:
            for col in ("flags_mask", "t1_correct_mask", "special_index"):
                df[col] = -1
            df["min_ratio"] = np.nan
            df["max_ratio"] = np.nan
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def _noise(cfg) -> NoiseSpec | None:
    return NoiseSpec(**cfg["noise"]) if cfg.get("noise") else None


# ---------------------------------------------------------------- MC shards
def mc_shard(cfg: dict, spec: dict) -> pd.DataFrame:
    widths, shots, R = spec["widths"], int(spec["shots"]), int(cfg["replicates"])
    n = int(sum(widths))
    phases = pd.DataFrame(spec["phases"])
    meta = phases.loc[phases.index.repeat(R)].reset_index(drop=True)
    meta["replicate_id"] = np.tile(np.arange(R), len(phases))
    phis = meta["phi"].to_numpy()
    rng = generator(cfg["master_seed"], stable_int(cfg["experiment_id"]), stable_int(widths_key(widths)), shots, spec["chunk_index"])
    noise = _noise(cfg)
    specs = partition_blocks(widths)
    counts = [batch_block_counts(phis, s, shots, rng, noise) for s in specs]
    jitter = [rng.random(c.shape) * 0.5 for c in counts]
    results = decode_all(counts, widths, cfg["epsilons"], jitter, int(cfg.get("likelihood_refine", 4)), noise)
    meta["widths"] = widths_key(widths)
    meta["n"] = n
    meta["shots_per_block"] = shots
    meta["shots_total"] = shots * len(specs)
    meta["u_queries_total"] = shots * sum(s.u_queries_per_shot for s in specs)
    meta["noise_model"] = "ideal" if noise is None else repr(noise)
    meta["allocation_policy"] = "uniform_fixed"
    meta["overlap_policy"] = "none"
    meta["stopping_policy"] = "none"
    return _rows(meta, widths, results, true_block_outcomes(phis, widths))


def mc_shards(cfg: dict) -> list[dict]:
    shards = []
    for widths in cfg["partitions"]:
        table = make_phase_table(widths, int(cfg["phases_per_stratum"]), cfg["split"], int(cfg["master_seed"]))
        cols = ["phase_id", "split", "stratum", "phi", "final_residual", "boundary_hardness"]
        records = table[cols].to_dict("records")
        size = int(cfg.get("phase_chunk", 50))
        for shots in cfg["shots_per_block"]:
            for ci, a in enumerate(range(0, len(records), size)):
                shards.append({"shard_id": f"w{widths_key(widths)}_s{shots}_c{ci}", "widths": list(widths), "shots": int(shots), "chunk_index": ci, "phases": records[a:a + size]})
    return shards


# ------------------------------------------------------------- limit shards
def limit_shard(cfg: dict, spec: dict) -> pd.DataFrame:
    """Infinite-shot limit on a dense grid: exact probabilities as counts."""
    widths = spec["widths"]
    n = int(sum(widths))
    G = 1 << (n + int(cfg["grid_refine"]))
    grid = (np.arange(G) + float(spec["grid_offset"])) / G
    rng = generator(cfg["master_seed"], stable_int(widths_key(widths)), int(round(spec["grid_offset"] * 1000)))
    specs = partition_blocks(widths)
    probs = [block_probabilities(grid, s.offset, s.width) for s in specs]
    jitter = [rng.random(p.shape) * 1e-9 for p in probs]  # only breaks exact (floating) ties
    meta = pd.DataFrame({"phase_id": [f"grid{spec['grid_offset']}:{g}" for g in range(G)], "stratum": "grid", "phi": grid})
    meta["final_residual"] = final_residual(grid, n)
    meta["boundary_hardness"] = boundary_hardness(grid, widths)
    meta["replicate_id"] = 0
    meta["widths"] = widths_key(widths)
    meta["n"] = n
    meta["shots_per_block"] = -1  # infinite-shot limit
    meta["grid_offset"] = spec["grid_offset"]
    results = []
    for eps in cfg["epsilons"]:
        for name, rule in DECODER_VARIANTS:
            r = awqpe_vectorised(probs, widths, float(eps), jitter=jitter, special_chunk_rule=rule)
            results.append((name, float(eps), r["estimate"], r))
    return _rows(meta, widths, results, true_block_outcomes(grid, widths))


def limit_shards(cfg: dict) -> list[dict]:
    return [{"shard_id": f"w{widths_key(w)}_o{str(o).replace('.', 'p')}", "widths": list(w), "grid_offset": float(o)}
            for w in cfg["partitions"] for o in cfg["grid_offsets"]]


def main() -> None:
    ap = common_arguments(__doc__, "p2_baseline_mc.yaml")
    ap.add_argument("--mode", choices=["mc", "limit"], default="mc")
    args = ap.parse_args()
    if args.mode == "limit" and args.config.endswith("p2_baseline_mc.yaml"):
        args.config = args.config.replace("p2_baseline_mc.yaml", "p2_infinite_shot_limit.yaml")
    cfg = resolve_config(args)
    if args.mode == "limit":
        run_sharded(args, cfg, limit_shards(cfg), limit_shard)
    else:
        run_sharded(args, cfg, mc_shards(cfg), mc_shard)


if __name__ == "__main__":
    main()
