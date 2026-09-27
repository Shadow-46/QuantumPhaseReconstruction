"""Phase 5: multi-step global shot allocation under a fixed budget.

Per trial: every block gets S0 shots; then (r - 1) * S0 * B further shots are
handed out in batches of dS, one block per decision, by a policy. Every
block owns a pre-generated outcome stream (common random numbers): shot i of
block b is identical whichever policy requests it, so policies are paired
shot-for-shot. Allocation never depends on the decoder (no practical signal
uses a decoder output), so each trajectory is scored under all three
decoders: AWQPE eps 0.9 (faithful), AWQPE eps_safe (D-011), D2 likelihood.
The greedy oracle is decoder-specific and runs once per decoder.

At every decision the realized benefit of giving the next dS shots to EACH
block is computed under each decoder (evaluator side, truth used only for
scoring), so predicted utility can be compared with realized benefit,
regret, and oracle agreement.

    python -m research.awqpe.run_adaptive_shots [--pilot] [--max-workers N] [--config ...]
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from research.awqpe.allocation.policies import POLICIES, AllocationState
from research.awqpe.baseline.awqpe import awqpe_vectorised
from research.awqpe.blocks.geometry import partition_blocks
from research.awqpe.decode.likelihood import GridPosterior
from research.awqpe.evaluation.metrics import circular_error, exact_success
from research.awqpe.model.kernel import block_probabilities
from research.awqpe.oracle.greedy import greedy_oracle_choice
from research.awqpe.phases import make_phase_table
from research.awqpe.runner.core import common_arguments, resolve_config, run_sharded
from research.awqpe.sim.seeds import generator, stable_int
from research.awqpe.verification.theory import safe_epsilon


def make_streams(phis, specs, L, rng):
    """Pre-generated outcome streams (T, L) per block, by inverse CDF of the kernel."""
    out = []
    for s in specs:
        cdf = np.cumsum(block_probabilities(phis, s.offset, s.width), axis=1)
        cdf /= cdf[:, -1:]
        u = rng.random((len(phis), L))
        y = np.empty((len(phis), L), dtype=np.int16)
        for t in range(len(phis)):
            y[t] = np.minimum(np.searchsorted(cdf[t], u[t], side="right"), s.M - 1)
        out.append(y)
    return out


def batch_counts(stream, start, dS, M):
    """Counts of stream[t, start[t] : start[t] + dS] for every trial -> (T, M)."""
    T = stream.shape[0]
    idx = start[:, None] + np.arange(dS)[None, :]
    ys = np.take_along_axis(stream, idx, axis=1)
    c = np.zeros((T, M), dtype=np.int64)
    np.add.at(c, (np.repeat(np.arange(T), dS), ys.ravel()), 1)
    return c


class Scorer:
    """Evaluator-side decoding and scoring (uses truth only for scoring)."""

    def __init__(self, widths, specs, n, phis, jitter, refine):
        self.widths, self.specs, self.n, self.phis, self.jitter, self.refine = widths, specs, n, phis, jitter, refine
        self.decoders = [("awqpe", 0.9), ("awqpe_eps_safe", safe_epsilon(widths)), ("likelihood", np.nan)]

    def decode(self, counts, loglik):
        out = {}
        for name, eps in self.decoders[:2]:
            out[name] = awqpe_vectorised(counts, self.widths, eps, jitter=self.jitter)["estimate"] / 2.0**self.n
        out["likelihood"] = np.argmax(loglik, axis=1) / loglik.shape[1]
        return out

    def score(self, phi_hat):
        err = circular_error(phi_hat, self.phis)
        return err, err <= 2.0**-self.n + 1e-15


def run_trajectory(policy_name, oracle_decoder, ctx):
    """One allocation trajectory for all trials in the chunk. Returns (final rows, step rows)."""
    specs, n, S0, dS, B, T = ctx["specs"], ctx["n"], ctx["S0"], ctx["dS"], ctx["B"], ctx["T"]
    streams, scorer, G = ctx["streams"], ctx["scorer"], ctx["G"]
    shots = np.full((T, B), S0, dtype=np.int64)
    counts = [batch_counts(streams[b], np.zeros(T, dtype=np.int64), S0, specs[b].M) for b in range(B)]
    gp = GridPosterior(n, T, ctx["refine"])
    for s, c in zip(specs, counts):
        gp.add_counts(s, c)
    policy = None if oracle_decoder else POLICIES[policy_name]()
    prng = generator(ctx["seed_base"], stable_int(policy_name + (oracle_decoder or "")), 7)
    steps = ctx["extra"] // dS
    step_rows = []
    dec_names = [d for d, _ in scorer.decoders]
    rows = np.arange(T)
    for step in range(steps):
        base = scorer.decode(counts, gp.loglik)
        base_scores = {d: scorer.score(base[d]) for d in dec_names}
        # realized counterfactual gains of each candidate block, per decoder (evaluator side only)
        deltas = [batch_counts(streams[b], shots[:, b], dS, specs[b].M) for b in range(B)]
        gain_tol = {d: np.zeros((T, B)) for d in dec_names}
        gain_err = {d: np.zeros((T, B)) for d in dec_names}
        cf_hat = {d: np.zeros((T, B)) for d in dec_names}
        for b in range(B):
            c_b = [c + (deltas[b] if i == b else 0) for i, c in enumerate(counts)]
            tmp = GridPosterior(n, 1, ctx["refine"])
            tmp.loglik = gp.loglik.copy()
            tmp.add_counts(specs[b], deltas[b])
            est = scorer.decode(c_b, tmp.loglik)
            for d in dec_names:
                e1, t1 = scorer.score(est[d])
                e0, t0 = base_scores[d]
                gain_tol[d][:, b] = t1.astype(int) - t0.astype(int)
                gain_err[d][:, b] = (e0 - e1) * 2.0**n
                cf_hat[d][:, b] = est[d]
        state = AllocationState(specs, n, counts, shots, gp.loglik, G, int(ctx["extra"] - step * dS), prng, step)
        if oracle_decoder:
            chosen = greedy_oracle_choice(gain_tol[oracle_decoder], gain_err[oracle_decoder], prng)
            sig = np.zeros((T, B))
        else:
            chosen, sig = policy.choose(state)
        rec = {"trial": np.arange(T), "step": step, "chosen_block": chosen + 1, "remaining_before": state.remaining}
        for b in range(B):
            rec[f"shots_b{b + 1}"] = shots[:, b].copy()
            rec[f"signal_b{b + 1}"] = sig[:, b]
        for d in dec_names:
            gt = gain_tol[d]
            score = gt * 1e6 + gain_err[d]
            best = np.argmax(score, axis=1)
            rec[f"{d}__gain_tol_chosen"] = gt[rows, chosen]
            rec[f"{d}__gain_err_chosen"] = gain_err[d][rows, chosen]
            rec[f"{d}__gain_tol_best"] = gt[rows, best]
            rec[f"{d}__gain_err_best"] = gain_err[d][rows, best]
            rec[f"{d}__chose_best"] = score[rows, chosen] >= score[rows, best] - 1e-12
            rec[f"{d}__tol_before"] = base_scores[d][1]
            rec[f"{d}__changed"] = cf_hat[d][rows, chosen] != base[d]
            for b in range(B):
                rec[f"{d}__gain_err_b{b + 1}"] = gain_err[d][:, b]
                rec[f"{d}__gain_tol_b{b + 1}"] = gt[:, b]
        step_rows.append(pd.DataFrame(rec))
        for b in range(B):  # commit the chosen batch
            m = chosen == b
            if m.any():
                counts[b][m] += deltas[b][m]
                shots[m, b] += dS
                sub = np.zeros_like(deltas[b])
                sub[m] = deltas[b][m]
                gp.add_counts(specs[b], sub)
    final = scorer.decode(counts, gp.loglik)
    frames = []
    for d, eps in scorer.decoders:
        if oracle_decoder and d != oracle_decoder:
            continue
        err, tol = scorer.score(final[d])
        df = pd.DataFrame({"trial": np.arange(T), "decoder": d, "epsilon": eps, "phi_hat": final[d], "error": err, "tol": tol,
                           "exact": exact_success(np.mod(np.floor(final[d] * 2**n + 0.5), 2**n), scorer.phis, n)})
        for b in range(B):
            df[f"shots_b{b + 1}"] = shots[:, b]
        df["total_shots"] = shots.sum(axis=1)
        df["u_queries"] = (shots * np.array([s.u_queries_per_shot for s in specs])[None, :]).sum(axis=1)
        df["n_decisions"] = steps
        frames.append(df)
    steps_df = pd.concat(step_rows, ignore_index=True) if step_rows else pd.DataFrame()
    return pd.concat(frames, ignore_index=True), steps_df


def adaptive_shard(cfg: dict, spec: dict) -> pd.DataFrame:
    widths, S0, r, R = spec["widths"], int(spec["S0"]), int(spec["r"]), int(cfg["replicates"])
    dS = S0 if spec["dS"] == "S0" else int(spec["dS"])
    specs = partition_blocks(widths)
    n, B = int(sum(widths)), len(widths)
    refine = int(cfg["likelihood_refine"].get(str(n), 3))
    extra = (r - 1) * S0 * B
    phases = pd.DataFrame(spec["phases"])
    meta = phases.loc[phases.index.repeat(R)].reset_index(drop=True)
    meta["replicate_id"] = np.tile(np.arange(R), len(phases))
    seed_base = int(generator(cfg["master_seed"], stable_int(cfg["experiment_id"]), stable_int("-".join(map(str, widths))), S0, r,
                              stable_int(str(spec["dS"])), spec["chunk_index"]).integers(0, 2**62))
    rng = generator(seed_base, 1)
    phis_all = meta["phi"].to_numpy()
    streams_all = make_streams(phis_all, specs, S0 + extra, rng)
    jitter_all = [rng.random((len(phis_all), s.M)) * 0.5 for s in specs]
    # decoder-limited labels: does each decoder recover the phase with infinitely many shots? (evaluator side)
    probs = [block_probabilities(phis_all, s.offset, s.width) for s in specs]
    limit_ok = {"likelihood": np.ones(len(phis_all), dtype=bool)}
    for name, eps in (("awqpe", 0.9), ("awqpe_eps_safe", safe_epsilon(widths))):
        est = awqpe_vectorised(probs, widths, eps)["estimate"] / 2.0**n
        limit_ok[name] = circular_error(est, phis_all) <= 2.0**-n + 1e-15
    G = 1 << (n + refine)
    chunk = max(8, (1 << 22) >> (n + refine))
    finals, steps_all = [], []
    for a in range(0, len(phis_all), chunk):
        sl = slice(a, a + chunk)
        phis = phis_all[sl]
        scorer = Scorer(widths, specs, n, phis, [j[sl] for j in jitter_all], refine)
        ctx = {"specs": specs, "n": n, "S0": S0, "dS": dS, "B": B, "T": len(phis), "streams": [s[sl] for s in streams_all],
               "scorer": scorer, "G": G, "refine": refine, "extra": extra, "seed_base": seed_base + a}
        runs = [(p, None) for p in cfg["policies"]] + [("greedy_oracle", d) for d in cfg["oracle_decoders"]]
        for pname, odec in runs:
            fin, stp = run_trajectory(pname, odec, ctx)
            label = pname if odec is None else f"greedy_oracle[{odec}]"
            for df in (fin, stp):
                if len(df):
                    df["trial"] = df["trial"] + a
                    df["policy"] = label
                    df["is_oracle"] = odec is not None
            finals.append(fin)
            if len(stp):
                steps_all.append(stp)
    fin = pd.concat(finals, ignore_index=True).join(meta, on="trial")
    fin["limit_correct"] = [bool(limit_ok[d][t]) for d, t in zip(fin.decoder, fin.trial)]
    fin["row_type"] = "final"
    st = pd.concat(steps_all, ignore_index=True).join(meta[["phase_id", "stratum", "replicate_id"]], on="trial")
    st["row_type"] = "step"
    out = pd.concat([fin, st], ignore_index=True, sort=False)
    out["widths"], out["n"], out["B"], out["S0"], out["r"], out["dS"] = "-".join(map(str, widths)), n, B, S0, r, dS
    out["split"] = spec["split"]
    return out


def shards(cfg):
    out = []
    for widths in cfg["partitions"]:
        table = make_phase_table(widths, int(cfg["phases_per_stratum"]), cfg["split"], int(cfg["master_seed"]))
        records = table[["phase_id", "stratum", "phi", "final_residual", "boundary_hardness"]].to_dict("records")
        size = int(cfg.get("phase_chunk", 20))
        for S0 in cfg["S0"]:
            for r in cfg["budget_multipliers"]:
                for dS in cfg["dS"]:
                    for ci, a in enumerate(range(0, len(records), size)):
                        out.append({"shard_id": f"w{'-'.join(map(str, widths))}_S{S0}_r{r}_d{dS}_c{ci}", "widths": widths, "S0": S0, "r": r,
                                    "dS": dS, "chunk_index": ci, "split": cfg["split"], "phases": records[a:a + size]})
    return out


def main() -> None:
    args = common_arguments(__doc__, "p5_adaptive_dev.yaml").parse_args()
    cfg = resolve_config(args)
    run_sharded(args, cfg, shards(cfg), adaptive_shard, extra_manifest={"status": cfg.get("status", "PRELIMINARY / DEV")})


if __name__ == "__main__":
    main()
