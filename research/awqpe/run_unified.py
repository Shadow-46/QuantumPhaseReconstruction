"""Phase 7: unified {shots, overlap, stop} controller vs its ablations (docs/DECISIONS.md D-031).

Arms (D2 decoder, P8 credible stopping at tau, S0 = batch = 4, cap per chunk; CRN streams):
  stop_uniform     round-robin shots + stop                  (P8 arm, run_tolerance.run_stopping)
  stop_eig         information-gain shots + stop             (controller, no overlap candidates)
  unified_bridge   information-gain over chunks + O-bridge (shift 1), <= A overlap batches, + stop
  unified_ext      information-gain over chunks + O-ext (v = 1),      <= A overlap batches, + stop
The overlap variants and A are the P6 dev selections frozen under D-027 (bridge_s1_A2,
ext_v1_eig_A2, for the D2 decoder); nothing is tuned in P7.

    python -m research.awqpe.run_unified --config research/configs/p7_unified_dev.yaml [--pilot]
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from research.awqpe.blocks.geometry import partition_blocks
from research.awqpe.controller import unified_run
from research.awqpe.evaluation.metrics import circular_error
from research.awqpe.overlap.candidates import overlap_candidates
from research.awqpe.phases import make_phase_table
from research.awqpe.run_adaptive_shots import make_streams
from research.awqpe.run_tolerance import run_stopping
from research.awqpe.runner.core import common_arguments, resolve_config, run_sharded
from research.awqpe.sim.seeds import generator, stable_int


def unified_shard(cfg, spec):
    widths, R = list(spec["widths"]), int(cfg["replicates"])
    n = int(sum(widths))
    ph = pd.DataFrame(spec["phases"])
    meta = ph.loc[ph.index.repeat(R)].reset_index(drop=True)
    meta["replicate_id"] = np.tile(np.arange(R), len(ph))
    phis, T = meta.phi.to_numpy(), len(meta)
    wkey = "-".join(map(str, widths))
    root = generator(cfg["master_seed"], stable_int(cfg["experiment_id"]), stable_int(wkey), int(spec["chunk_index"]))
    S0, cap, A, alpha = int(cfg["S0"]), int(cfg["cap_per_block"]), int(cfg["A"]), float(cfg["alpha"])
    chunks = partition_blocks(widths)
    bridge = [c[1] for c in overlap_candidates(widths, "bridge", shift="1")]
    ext = [c[1] for c in overlap_candidates(widths, "ext", v=1)]
    streams = make_streams(phis, chunks + bridge + ext, cap, generator(int(root.integers(0, 2**62)), 1))
    B, nb = len(chunks), len(bridge)
    s_chunk, s_bridge, s_ext = streams[:B], streams[B:B + nb], streams[B + nb:]
    arm_seed = int(root.integers(0, 2**62))
    refine = int(cfg["likelihood_refine"].get(str(n), 3))
    frames = []
    for j in cfg["tau_exponents"]:
        j = n // 2 if j == "half" else n if j == "n" else int(j)
        tau = 2.0**-j
        runs = {}
        est, sh, stp, _ = run_stopping(chunks, s_chunk, n, refine, S0, cap, tau, alpha, "uniform", generator(arm_seed, stable_int("stop_uniform"), j))
        runs["stop_uniform"] = (est, sh, stp, np.zeros(T, dtype=np.int64), chunks)
        for arm, cands, cs in (("stop_eig", [], []), ("unified_bridge", bridge, s_bridge), ("unified_ext", ext, s_ext)):
            est, sh, stp, nov = unified_run(chunks, cands, s_chunk + cs, n, refine, S0, cap, A if cands else 0, tau, alpha,
                                            generator(arm_seed, stable_int(arm), j))
            runs[arm] = (est, sh, stp, nov, chunks + cands)
        for arm, (est, sh, stp, nov, specs) in runs.items():
            u_per = np.array([s.u_queries_per_shot for s in specs])
            err = circular_error(est, phis)
            d = meta.copy()
            d["tau"], d["tau_exp"], d["arm"] = tau, j, arm
            d["shots"], d["overlap_shots"], d["overlap_batches"] = sh.sum(axis=1), sh[:, B:].sum(axis=1), nov
            d["u_queries"], d["stopped"], d["error"] = sh @ u_per, stp, err
            d["covered"] = err <= tau + 1e-15
            frames.append(d)
    out = pd.concat(frames, ignore_index=True)
    out["widths"], out["n"], out["S0"], out["A"], out["cap_per_block"], out["split"] = wkey, n, S0, A, cap, cfg["split"]
    return out


def shards(cfg):
    out, chunk = [], int(cfg.get("phase_chunk", 10))
    for widths in cfg["partitions"]:
        t = make_phase_table(widths, int(cfg["phases_per_stratum"]), cfg["split"], int(cfg["master_seed"]), strata=tuple(cfg["strata"]))
        recs = t[["phase_id", "stratum", "phi"]].to_dict("records")
        for ci, a in enumerate(range(0, len(recs), chunk)):
            out.append({"shard_id": f"w{'-'.join(map(str, widths))}_c{ci}", "widths": widths, "chunk_index": ci, "phases": recs[a:a + chunk]})
    return out


def main() -> None:
    args = common_arguments(__doc__, "p7_unified_dev.yaml").parse_args()
    cfg = resolve_config(args)
    run_sharded(args, cfg, shards(cfg), unified_shard, extra_manifest={"status": cfg.get("status", "PRELIMINARY / DEV")})


if __name__ == "__main__":
    main()
