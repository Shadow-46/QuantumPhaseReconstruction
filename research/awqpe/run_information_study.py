"""Phase 4: do window signals predict the marginal value of extra shots? (PRELIMINARY, dev phases)

Counterfactual design, per trial:
  1. every block b receives S0 shots (counts c0_b);
  2. every signal in information/signals.py is computed from c0 only;
  3. for each block b separately, an extra batch of dS shots (x_b, drawn once
     per block and reused, i.e. common random numbers) is added to b alone,
     and the phase is re-decoded;
  4. improvement = outcome(c0 + x_b on block b) - outcome(c0).
Decoders are never mixed: improvement is measured within D1@0.9, D1@eps_safe
and D2 separately (docs/DECISIONS.md D-005).

Truth is used only after decoding, to score outcomes (evaluator side).

    python -m research.awqpe.run_information_study [--pilot] [--max-workers N]
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from research.awqpe.baseline.awqpe import awqpe_vectorised
from research.awqpe.blocks.geometry import partition_blocks
from research.awqpe.decode.likelihood import GridPosterior
from research.awqpe.evaluation.metrics import circular_error, exact_success
from research.awqpe.information.signals import all_block_signals, global_risk
from research.awqpe.model.kernel import block_delta
from research.awqpe.phases import make_phase_table
from research.awqpe.runner.core import common_arguments, resolve_config, run_sharded
from research.awqpe.sim.oracle_sim import batch_block_counts
from research.awqpe.sim.seeds import generator, stable_int
from research.awqpe.verification.theory import safe_epsilon



def _chunk(n: int, refine: int) -> int:
    """Trials per posterior batch so that a (trials x grid) float64 array stays <= 64 MB."""
    return max(32, (1 << 23) >> (n + refine))


def _d2_estimates(specs, counts, n, refine):
    """MAP estimate from counts (posterior built in trial chunks)."""
    T = counts[0].shape[0]
    est = np.empty(T)
    CHUNK = _chunk(n, refine)
    for a in range(0, T, CHUNK):
        gp = GridPosterior(n, min(CHUNK, T - a), refine)
        for s, c in zip(specs, counts):
            gp.add_counts(s, c[a:a + CHUNK])
        est[a:a + CHUNK] = gp.map_estimate()
    return est


def info_shard(cfg: dict, spec: dict) -> pd.DataFrame:
    widths, S0, dS, R = spec["widths"], int(spec["S0"]), int(spec["dS"]), int(cfg["replicates"])
    n = int(sum(widths))
    B = len(widths)
    specs = partition_blocks(widths)
    refine = int(cfg.get("likelihood_refine", 4))
    phases = pd.DataFrame(spec["phases"])
    meta = phases.loc[phases.index.repeat(R)].reset_index(drop=True)
    meta["replicate_id"] = np.tile(np.arange(R), len(phases))
    phis = meta["phi"].to_numpy()
    T = len(phis)
    rng = generator(cfg["master_seed"], stable_int(cfg["experiment_id"]), stable_int("-".join(map(str, widths))), S0, dS, spec["chunk_index"])
    c0 = [batch_block_counts(phis, s, S0, rng) for s in specs]
    extra = [batch_block_counts(phis, s, dS, rng) for s in specs]
    jitter = [rng.random(c.shape) * 0.5 for c in c0]

    # ---- signals from initial counts only (no truth)
    CHUNK = _chunk(n, refine)
    sig_frames = []
    for a in range(0, T, CHUNK):
        sl = slice(a, a + CHUNK)
        gp = GridPosterior(n, len(phis[sl]), refine)
        for s, c in zip(specs, c0):
            gp.add_counts(s, c[sl])
        sig = all_block_signals(specs, [c[sl] for c in c0], n, gp, with_cw=bool(cfg.get("legacy_cw", True)))
        sig["global_risk"] = np.repeat(global_risk(gp, n)[:, None], B, axis=1)
        sig_frames.append(sig)
    signals = {k: np.concatenate([f[k] for f in sig_frames], axis=0) for k in sig_frames[0]}

    # ---- decoders on c0 and on each single-block counterfactual
    def decode(counts):
        out = {}
        for name, eps in (("awqpe", 0.9), ("awqpe_eps_safe", safe_epsilon(widths))):
            out[name] = awqpe_vectorised(counts, widths, eps, jitter=jitter)["estimate"] / 2.0**n
        out["likelihood"] = _d2_estimates(specs, counts, n, refine)
        return out

    base = decode(c0)
    cf = []
    for b in range(B):
        counts = [c + (extra[b] if i == b else 0) for i, c in enumerate(c0)]
        cf.append(decode(counts))

    # ---- evaluator side: outcomes
    t1_true = [np.mod(np.floor(block_delta(phis, s.offset) * s.M + 0.5), s.M).astype(np.int64) for s in specs]
    t1_obs = [np.argmax(c + j, axis=1) for c, j in zip(c0, jitter)]
    frames = []
    tau_c = 2.0 ** -(n - 2)
    for dec in base:
        phi0 = base[dec]
        err0 = circular_error(phi0, phis)
        exact0 = exact_success(np.mod(np.floor(phi0 * 2**n + 0.5), 2**n), phis, n)
        for b in range(B):
            phi1 = cf[b][dec]
            err1 = circular_error(phi1, phis)
            df = meta.copy()
            df["decoder"], df["block"], df["offset"], df["width"] = dec, b + 1, specs[b].offset, specs[b].width
            df["err0"], df["err1"] = err0, err1
            df["tol0"], df["tol1"] = err0 <= 2.0**-n + 1e-15, err1 <= 2.0**-n + 1e-15
            df["tolc0"], df["tolc1"] = err0 <= tau_c + 1e-15, err1 <= tau_c + 1e-15
            df["exact0"] = exact0
            df["exact1"] = exact_success(np.mod(np.floor(phi1 * 2**n + 0.5), 2**n), phis, n)
            df["t1_correct"] = t1_obs[b] == t1_true[b]
            for k, v in signals.items():
                df[k] = v[:, b]
            frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    out["widths"], out["n"], out["B"], out["S0"], out["dS"] = "-".join(map(str, widths)), n, B, S0, dS
    return out


def shards(cfg):
    out = []
    for widths in cfg["partitions"]:
        table = make_phase_table(widths, int(cfg["phases_per_stratum"]), cfg["split"], int(cfg["master_seed"]))
        records = table[["phase_id", "stratum", "phi", "final_residual", "boundary_hardness"]].to_dict("records")
        size = int(cfg.get("phase_chunk", 40))
        for S0 in cfg["S0"]:
            for mult in cfg["dS_multipliers"]:
                dS = int(S0 * mult)
                for ci, a in enumerate(range(0, len(records), size)):
                    out.append({"shard_id": f"w{'-'.join(map(str, widths))}_S{S0}_d{dS}_c{ci}", "widths": widths, "S0": S0, "dS": dS,
                                "chunk_index": ci, "phases": records[a:a + size]})
    return out


def main() -> None:
    args = common_arguments(__doc__, "p4_information.yaml").parse_args()
    cfg = resolve_config(args)
    run_sharded(args, cfg, shards(cfg), info_shard, extra_manifest={"status": "PRELIMINARY (dev phases)"})


if __name__ == "__main__":
    main()
