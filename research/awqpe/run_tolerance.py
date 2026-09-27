"""Phase 8: phase tolerance, Mode B (docs/DECISIONS.md D-028).

For a declared tolerance tau and confidence 1 - alpha: how many shots (and
U-queries) are needed for realised coverage P(|phi_hat - phi|_circle <= tau) >= 1 - alpha,
and do a posterior-credible stopping rule and information-guided allocation reduce them?

Arms (all share per-block common-random-number outcome streams):
  fixed_S<S>        S shots on every block, no adaptivity. Decoded by D2 (grid MAP) and by
                    D1 = faithful AWQPE at eps_safe. The dev split picks the smallest S whose
                    coverage reaches 1 - alpha ("tau-matched fixed" baseline).
  stop_uniform      S0 per block, then round-robin batches of S0; stop as soon as
                    P(|phi - MAP| <= tau | data) >= 1 - alpha (stopping/rules.py). D2 only.
  stop_eig          same rule; each batch goes to the block with the largest expected
                    information gain about the tau-resolution cell (P5 signal, cell bits = tau_bits).
Every arm is run on the full partition and on the tau-matched prefix ("_trunc": the
shortest MSB-first block prefix with at least tau_bits(tau) bits), so savings from stopping
are separated from the trivial saving of measuring fewer bits. A per-block cap bounds
every arm; trials that hit it without stopping are kept and flagged stopped = False.

    python -m research.awqpe.run_tolerance --config research/configs/p8_tolerance_dev.yaml [--pilot]
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from research.awqpe.allocation.eig_cached import eig_cell_support
from research.awqpe.allocation.policies import AllocationState
from research.awqpe.baseline.awqpe import awqpe_vectorised
from research.awqpe.blocks.geometry import partition_blocks
from research.awqpe.decode.likelihood import GridPosterior
from research.awqpe.evaluation.metrics import circular_error
from research.awqpe.phases import make_phase_table
from research.awqpe.run_adaptive_shots import batch_counts, make_streams
from research.awqpe.runner.core import common_arguments, resolve_config, run_sharded
from research.awqpe.sim.seeds import generator, stable_int
from research.awqpe.stopping.rules import credible_mass_at_map, should_stop, tau_bits, truncated_widths
from research.awqpe.verification.theory import safe_epsilon


def eps_safe(widths) -> float:
    return 0.9 if len(widths) < 2 else safe_epsilon(widths)  # one block: no borrow, eps is inert


def run_stopping(specs, streams, n, refine, S0, cap, tau, alpha, policy, rng):
    """Sequential D2 run with credible stopping. Returns MAP (T,), shots (T, B), stopped (T,), mass (T,)."""
    T, B = streams[0].shape[0], len(specs)
    shots = np.full((T, B), S0, dtype=np.int64)
    counts = [batch_counts(streams[b], np.zeros(T, dtype=np.int64), S0, specs[b].M) for b in range(B)]
    gp = GridPosterior(n, T, refine)
    for s, c in zip(specs, counts):
        gp.add_counts(s, c)
    stopped, mass = should_stop(gp.loglik, tau, alpha)
    cell_bits = min(n, tau_bits(tau))
    step = 0
    while True:
        E = shots + S0 <= cap
        if policy == "uniform":
            b_rr = step % B
            chosen = np.full(T, b_rr)
            active = ~stopped & E[:, b_rr]
            if not (~stopped & E.any(axis=1)).any():
                break
        else:
            active = ~stopped & E.any(axis=1)
            if not active.any():
                break
            ia = np.flatnonzero(active)
            st = AllocationState(specs, cell_bits, [c[ia] for c in counts], shots[ia], gp.loglik[ia], gp.G, 0, rng, step)
            sig, _ = eig_cell_support(st)
            sig = np.where(E[ia], sig + rng.random(sig.shape) * 1e-12, -np.inf)
            chosen = np.full(T, -1)
            chosen[ia] = np.argmax(sig, axis=1)
        for b in range(B):
            m = active & (chosen == b)
            if not m.any():
                continue
            d = batch_counts(streams[b], np.minimum(shots[:, b], streams[b].shape[1] - S0), S0, specs[b].M)
            d[~m] = 0
            counts[b] += d
            shots[m, b] += S0
            gp.add_counts(specs[b], d)
        step += 1
        live = ~stopped
        s_new, m_new = should_stop(gp.loglik[live], tau, alpha)
        mass[live] = m_new
        stopped[np.flatnonzero(live)[s_new]] = True
    return gp.map_estimate(), shots, stopped, mass


def tol_shard(cfg, spec):
    widths, S0, R = list(spec["widths"]), int(cfg["S0"]), int(cfg["replicates"])
    n_full = int(sum(widths))
    ph = pd.DataFrame(spec["phases"])
    meta = ph.loc[ph.index.repeat(R)].reset_index(drop=True)
    meta["replicate_id"] = np.tile(np.arange(R), len(ph))
    phis = meta.phi.to_numpy()
    T = len(phis)
    cap, alpha = int(cfg["cap_per_block"]), float(cfg["alpha"])
    wkey = "-".join(map(str, widths))
    root = generator(cfg["master_seed"], stable_int(cfg["experiment_id"]), stable_int(wkey), int(spec["chunk_index"]))
    all_specs = partition_blocks(widths)
    L = max(cap, max(cfg["fixed_S"]))
    streams_all = make_streams(phis, all_specs, L, generator(int(root.integers(0, 2**62)), 1))  # CRN across arms and tau
    jit_all = [generator(int(root.integers(0, 2**62)), 2).random((T, s.M)) * 0.5 for s in all_specs]
    arm_seed = int(root.integers(0, 2**62))
    taus = [2.0 ** -j for j in cfg["tau_exponents"] if j <= n_full]

    fixed_cache = {}

    def fixed(prefix, S):
        key = (tuple(prefix), S)
        if key not in fixed_cache:
            B, n = len(prefix), int(sum(prefix))
            specs = all_specs[:B]
            counts = [batch_counts(streams_all[b], np.zeros(T, dtype=np.int64), S, specs[b].M) for b in range(B)]
            gp = GridPosterior(n, T, int(cfg["likelihood_refine"].get(str(n), 3)))
            for s, c in zip(specs, counts):
                gp.add_counts(s, c)
            d1 = awqpe_vectorised(counts, prefix, eps_safe(prefix), jitter=jit_all[:B])["estimate"] / 2.0**n
            fixed_cache[key] = (gp.map_estimate(), d1, gp.loglik)
        return fixed_cache[key]

    frames = []
    for tau in taus:
        for trunc in (False, True):
            prefix = truncated_widths(widths, tau) if trunc else widths
            if trunc and len(prefix) == len(widths):
                continue  # prefix is the full partition: identical to the untruncated arm
            B, n = len(prefix), int(sum(prefix))
            specs = all_specs[:B]
            u_per = np.array([s.u_queries_per_shot for s in specs])
            refine = int(cfg["likelihood_refine"].get(str(n), 3))
            sfx = "_trunc" if trunc else ""
            results = []
            for S in cfg["fixed_S"]:
                d2, d1, ll = fixed(prefix, S)
                mass, _ = credible_mass_at_map(ll, tau)
                sh = np.full((T, B), S)
                results.append((f"fixed_S{S}{sfx}", "likelihood", d2, sh, np.ones(T, bool), mass, S))
                results.append((f"fixed_S{S}{sfx}", "awqpe_eps_safe", d1, sh, np.ones(T, bool), np.full(T, np.nan), S))
            for policy in cfg["stop_policies"]:
                rng = generator(arm_seed, stable_int(policy), int(round(-np.log2(tau))), int(trunc))
                est, sh, stp, mass = run_stopping(specs, streams_all[:B], n, refine, S0, cap, tau, alpha, policy, rng)
                results.append((f"stop_{policy}{sfx}", "likelihood", est, sh, stp, mass, np.nan))
            for arm, dec, est, sh, stp, mass, S in results:
                err = circular_error(est, phis)
                d = meta.copy()
                d["tau"], d["tau_exp"], d["alpha"] = tau, int(round(-np.log2(tau))), alpha
                d["arm"], d["decoder"], d["fixed_S"], d["truncated"] = arm, dec, float(S), trunc
                d["blocks_used"], d["n_eff"], d["stopped"], d["credible_mass"] = B, n, stp, mass
                d["shots"], d["u_queries"] = sh.sum(axis=1), sh @ u_per
                d["error"], d["covered"] = err, err <= tau + 1e-15
                frames.append(d)
    out = pd.concat(frames, ignore_index=True)
    out["widths"], out["n"], out["S0"], out["cap_per_block"], out["split"] = wkey, n_full, S0, cap, cfg["split"]
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
    args = common_arguments(__doc__, "p8_tolerance_dev.yaml").parse_args()
    cfg = resolve_config(args)
    run_sharded(args, cfg, shards(cfg), tol_shard, extra_manifest={"status": cfg.get("status", "PRELIMINARY / DEV")})


if __name__ == "__main__":
    main()
