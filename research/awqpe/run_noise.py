"""Phase 9: noise robustness under the analytic channels of model/noise.py (docs/DECISIONS.md D-030).

These are phenomenological Dirichlet-level channels (readout confusion, per-query global
depolarisation, control-phase jitter), NOT the noise models of the published AWQPE paper,
which could not be checked (docs/PAPER_VERSION_NOTES.md).

For every noise setting (each channel alone at three strengths, one combined, ideal):
  fixed_S<S>  S shots per block, decoded by faithful AWQPE (eps 0.9, eps_safe), by D2 with the
              ideal likelihood (misspecified: "likelihood_ideal") and by D2 with the true noise
              model in the likelihood ("likelihood_aware").
  stop_tau<j> the P8 posterior-credible stopping rule (uniform round-robin batches of S0, cap per
              block) for tau = 2^-j, with the ideal or the noise-aware posterior: does coverage survive?
All arms of a trial share the same noisy outcome streams (common random numbers).

    python -m research.awqpe.run_noise --config research/configs/p9_noise_dev.yaml [--pilot]
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from research.awqpe.baseline.awqpe import awqpe_vectorised
from research.awqpe.blocks.geometry import partition_blocks
from research.awqpe.decode.likelihood import GridPosterior
from research.awqpe.evaluation.metrics import circular_error
from research.awqpe.model.noise import NoiseSpec, noisy_block_probabilities
from research.awqpe.phases import make_phase_table
from research.awqpe.run_adaptive_shots import batch_counts
from research.awqpe.runner.core import common_arguments, resolve_config, run_sharded
from research.awqpe.sim.seeds import generator, stable_int
from research.awqpe.stopping.rules import should_stop
from research.awqpe.verification.theory import safe_epsilon


def noisy_streams(phis, specs, L, noise, rng):
    """(T, L) outcome streams per block by inverse CDF of the noisy block distribution."""
    out = []
    for s in specs:
        cdf = np.cumsum(noisy_block_probabilities(phis, s.offset, s.width, noise), axis=1)
        cdf /= cdf[:, -1:]
        u = rng.random((len(phis), L))
        y = np.empty((len(phis), L), dtype=np.int16)
        for t in range(len(phis)):
            y[t] = np.minimum(np.searchsorted(cdf[t], u[t], side="right"), s.M - 1)
        out.append(y)
    return out


def stop_uniform(specs, streams, n, refine, S0, cap, tau, alpha, noise):
    T, B = streams[0].shape[0], len(specs)
    shots = np.full((T, B), S0, dtype=np.int64)
    gp = GridPosterior(n, T, refine, noise)
    for b, s in enumerate(specs):
        gp.add_counts(s, batch_counts(streams[b], np.zeros(T, dtype=np.int64), S0, s.M))
    stopped, _ = should_stop(gp.loglik, tau, alpha)
    step = 0
    while True:
        b = step % B
        active = ~stopped & (shots[:, b] + S0 <= cap)
        if not (~stopped & (shots + S0 <= cap).any(axis=1)).any():
            break
        if active.any():
            d = batch_counts(streams[b], np.minimum(shots[:, b], streams[b].shape[1] - S0), S0, specs[b].M)
            d[~active] = 0
            shots[active, b] += S0
            gp.add_counts(specs[b], d)
            live = np.flatnonzero(~stopped)
            s_new, _ = should_stop(gp.loglik[live], tau, alpha)
            stopped[live[s_new]] = True
        step += 1
    return gp.map_estimate(), shots, stopped


def noise_shard(cfg, spec):
    widths, R = list(spec["widths"]), int(cfg["replicates"])
    n, specs = int(sum(widths)), partition_blocks(widths)
    noise = NoiseSpec(**spec["noise"])
    ph = pd.DataFrame(spec["phases"])
    meta = ph.loc[ph.index.repeat(R)].reset_index(drop=True)
    meta["replicate_id"] = np.tile(np.arange(R), len(ph))
    phis, T = meta.phi.to_numpy(), len(meta)
    wkey = "-".join(map(str, widths))
    root = generator(cfg["master_seed"], stable_int(cfg["experiment_id"]), stable_int(wkey), stable_int(spec["noise_id"]), int(spec["chunk_index"]))
    cap, S0, alpha = int(cfg["cap_per_block"]), int(cfg["S0"]), float(cfg["alpha"])
    L = max(cap, max(cfg["fixed_S"]))
    streams = noisy_streams(phis, specs, L, noise, generator(int(root.integers(0, 2**62)), 1))
    jit = [generator(int(root.integers(0, 2**62)), 2).random((T, s.M)) * 0.5 for s in specs]
    refine = int(cfg["likelihood_refine"].get(str(n), 3))
    u_per = np.array([s.u_queries_per_shot for s in specs])
    models = {"likelihood_ideal": None, "likelihood_aware": noise}
    rows = []
    for S in cfg["fixed_S"]:
        counts = [batch_counts(streams[b], np.zeros(T, dtype=np.int64), S, s.M) for b, s in enumerate(specs)]
        est = {"awqpe": awqpe_vectorised(counts, widths, 0.9, jitter=jit)["estimate"] / 2.0**n,
               "awqpe_eps_safe": awqpe_vectorised(counts, widths, safe_epsilon(widths), jitter=jit)["estimate"] / 2.0**n}
        for name, model in models.items():
            if name == "likelihood_aware" and noise.is_ideal:
                continue
            gp = GridPosterior(n, T, refine, model)
            for s, c in zip(specs, counts):
                gp.add_counts(s, c)
            est[name] = gp.map_estimate()
        sh = np.full((T, len(specs)), S)
        for dec, e in est.items():
            rows.append((f"fixed_S{S}", dec, e, sh, np.ones(T, bool), np.nan))
    for j in cfg["stop_tau_exponents"]:
        tau = 2.0 ** -(n // 2 if j == "half" else n if j == "n" else int(j))
        for name, model in models.items():
            if name == "likelihood_aware" and noise.is_ideal:
                continue
            e, sh, stp = stop_uniform(specs, streams, n, refine, S0, cap, tau, alpha, model)
            rows.append((f"stop_tau{int(round(-np.log2(tau)))}", name, e, sh, stp, tau))
    frames = []
    for arm, dec, e, sh, stp, tau in rows:
        err = circular_error(e, phis)
        d = meta.copy()
        d["arm"], d["decoder"], d["stopped"], d["tau"] = arm, dec, stp, tau
        d["shots"], d["u_queries"], d["error"] = sh.sum(axis=1), sh @ u_per, err
        d["tol"] = err <= 2.0**-n + 1e-15
        d["tol_half"] = err <= 2.0 ** -(n // 2) + 1e-15
        d["covered"] = err <= tau + 1e-15 if np.isfinite(tau) else d["tol"]
        frames.append(d)
    out = pd.concat(frames, ignore_index=True)
    for k, v in spec["noise"].items():
        out[k] = v
    out["noise_id"], out["channel"], out["widths"], out["n"], out["split"] = spec["noise_id"], spec["channel"], wkey, n, cfg["split"]
    return out


def shards(cfg):
    out, chunk = [], int(cfg.get("phase_chunk", 10))
    for widths in cfg["partitions"]:
        t = make_phase_table(widths, int(cfg["phases_per_stratum"]), cfg["split"], int(cfg["master_seed"]), strata=tuple(cfg["strata"]))
        recs = t[["phase_id", "stratum", "phi"]].to_dict("records")
        for ns in cfg["noise_settings"]:
            for ci, a in enumerate(range(0, len(recs), chunk)):
                out.append({"shard_id": f"w{'-'.join(map(str, widths))}_{ns['id']}_c{ci}", "widths": widths, "chunk_index": ci,
                            "noise_id": ns["id"], "channel": ns["channel"], "noise": ns["spec"], "phases": recs[a:a + chunk]})
    return out


def main() -> None:
    args = common_arguments(__doc__, "p9_noise_dev.yaml").parse_args()
    cfg = resolve_config(args)
    run_sharded(args, cfg, shards(cfg), noise_shard, extra_manifest={"status": cfg.get("status", "PRELIMINARY / DEV")})


if __name__ == "__main__":
    main()
