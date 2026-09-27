"""Phase 6: adaptive overlap (O-ext and O-bridge), kept separate from P5.

Every trial has chunk blocks (the AWQPE partition) plus, per mechanism
variant, one overlap candidate per internal boundary (overlap/candidates.py).
All blocks own pre-generated common-random-number outcome streams. Chunk
streams are identical across every arm and variant of a shard (same seed,
same stream length, chunks generated first). An arm is a sequence of
segments (policy, eligible set, length):

  B1 uniform          [(uniform, chunks, steps)]
  B2 p5_eig           [(eig, chunks, steps)]                         (the P5 policy)
  B3 overlap_alone    [(trigger, cands, A), (uniform, chunks, steps - A)]
  B4 eig_plus_overlap [(trigger, chunks + cands (<= A overlap batches), steps)]
  B1u / B2u           uniform / p5_eig on chunks with a per-trial U-query budget
                      equal to what B3 / B4 consumed (equal quantum work; the last
                      batch may overshoot, which favours the baseline)
  rescue diagnostics  continue B2 for A extra batches with greedy-oracle shots on
                      chunks / greedy-oracle overlap batches (capability), and with
                      practical eig shots / practical trigger overlap batches

Overlap triggers (decision side, no truth):
  eig         largest one-shot mutual information with the n-bit cell (P5 signal),
              computed for chunks and overlap candidates alike
  lowerhalf   overlap candidates at boundaries whose decoded lower part (faithful
              eps 0.9 chunk decode of the current counts) is exactly 10..0 are
              taken first; otherwise and for ties, eig

Decoders (identity explicit in every row): awqpe, awqpe_eps_safe (faithful,
ignore overlap blocks), awqpe_ext, awqpe_eps_safe_ext (non-paper O-ext
variants, decode/awqpe_overlap.py; ext variants only), likelihood (D2, fuses
every block).

    python -m research.awqpe.run_adaptive_overlap --config research/configs/p6_overlap_pilot.yaml
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from research.awqpe.allocation.policies import AllocationState, eig_cell_support
from research.awqpe.baseline.awqpe import awqpe_vectorised
from research.awqpe.blocks.geometry import partition_blocks
from research.awqpe.decode.awqpe_overlap import awqpe_ext_decode
from research.awqpe.decode.likelihood import GridPosterior
from research.awqpe.evaluation.metrics import circular_error
from research.awqpe.model.kernel import block_probabilities
from research.awqpe.oracle.greedy import greedy_oracle_choice
from research.awqpe.overlap.candidates import overlap_candidates
from research.awqpe.overlap.triggers import lowerhalf_mask
from research.awqpe.phases import make_phase_table
from research.awqpe.run_adaptive_shots import batch_counts, make_streams
from research.awqpe.runner.core import common_arguments, resolve_config, run_sharded
from research.awqpe.sim.seeds import generator, stable_int
from research.awqpe.verification.theory import safe_epsilon


class Ctx:
    def __init__(self, cfg, widths, S0, dS, steps, variant, phis, rng_seed, L):
        self.widths, self.S0, self.dS, self.steps, self.variant, self.L = widths, S0, dS, steps, variant, L
        self.chunks = partition_blocks(widths)
        self.B, self.n = len(widths), int(sum(widths))
        self.cands = overlap_candidates(widths, variant["mechanism"], v=variant.get("v", 1), shift=str(variant.get("shift", "half")))
        self.specs = self.chunks + [c[1] for c in self.cands]
        self.Btot = len(self.specs)
        self.phis, self.T = phis, len(phis)
        self.refine = int(cfg["likelihood_refine"].get(str(self.n), 3))
        self.G = 1 << (self.n + self.refine)
        rng = generator(rng_seed, 1)
        self.streams = make_streams(phis, self.specs, L, rng)  # chunks first: identical across variants
        jrng = generator(rng_seed, 2)
        self.jit = [jrng.random((self.T, s.M)) * 0.5 for s in self.specs]
        self.u = np.array([s.u_queries_per_shot for s in self.specs])
        self.eps = {"awqpe": 0.9, "awqpe_eps_safe": safe_epsilon(widths)}
        self.ext_capable = variant["mechanism"] == "ext"
        self.decoders = ["awqpe", "awqpe_eps_safe"] + (["awqpe_ext", "awqpe_eps_safe_ext"] if self.ext_capable else []) + ["likelihood"]

    def decode(self, counts, shots, loglik):
        n, B = self.n, self.B
        out = {}
        for name, eps in self.eps.items():
            out[name] = awqpe_vectorised(counts[:B], self.widths, eps, jitter=self.jit[:B])["estimate"] / 2.0**n
            if self.ext_capable:
                ext = {j: (counts[B + i], ve, shots[:, B + i] > 0, self.jit[B + i]) for i, (j, _, ve) in enumerate(self.cands)}
                out[name + "_ext"] = awqpe_ext_decode(counts[:B], self.widths, eps, self.jit[:B], ext) / 2.0**n
        out["likelihood"] = np.argmax(loglik, axis=1) / loglik.shape[1]
        return out

    def score(self, phi_hat):
        err = circular_error(phi_hat, self.phis)
        return err, err <= 2.0**-self.n + 1e-15


def lowerhalf_candidates(ctx, counts) -> np.ndarray:
    return lowerhalf_mask(counts[:ctx.B], ctx.widths, ctx.jit[:ctx.B], ctx.cands)


def boundary_states(ctx, chunk_counts, eps):
    """Per internal boundary j (T, B-1): Algorithm-1 ambiguity flag of the upper chunk j, and whether the
    corrected lower part (chunks j+1..B) is exactly 10..0 -- both from observed counts only."""
    r = awqpe_vectorised(chunk_counts, ctx.widths, eps, jitter=ctx.jit[:ctx.B])
    corr, w = r["corrected"], ctx.widths
    B = ctx.B
    lh = np.zeros((ctx.T, B - 1), dtype=bool)
    for j in range(B - 1):
        m = corr[:, j + 1] == (1 << (w[j + 1] - 1))
        for k in range(j + 2, B):
            m &= corr[:, k] == 0
        lh[:, j] = m
    return r["flags"][:, : B - 1], lh


def baseline_path(ctx, policy, targets, seed, arm, cap_factor=6):
    """Chunk-only allocation (uniform or P5 eig), decoded after EVERY batch, run until every trial's
    cumulative U-queries reach max(targets) (or a cap). Decisions never depend on the target, so one
    path serves every U-matched comparison (docs/DECISIONS.md D-022).

    Returns U (K+1, T), shots (K+1, T) and {decoder: phi_hat (K+1, T)}; row k is the state after k batches.
    """
    T, B, dS = ctx.T, ctx.B, ctx.dS
    shots = np.zeros((T, ctx.Btot), dtype=np.int64)
    shots[:, :B] = ctx.S0
    counts = [batch_counts(ctx.streams[b], np.zeros(T, dtype=np.int64), ctx.S0, ctx.specs[b].M) if b < B
              else np.zeros((T, ctx.specs[b].M), dtype=np.int64) for b in range(ctx.Btot)]
    gp = GridPosterior(ctx.n, T, ctx.refine)
    for b in range(B):
        gp.add_counts(ctx.specs[b], counts[b])
    u = shots @ ctx.u
    prng = generator(seed, stable_int(arm), 11)
    U, S, P = [u.copy()], [shots.sum(axis=1)], {d: [v] for d, v in ctx.decode(counts, shots, gp.loglik).items()}
    goal = float(np.max(targets))
    for step in range(cap_factor * (ctx.steps + 2)):
        if (u >= goal).all():
            break
        E = np.zeros((T, ctx.Btot), dtype=bool)
        E[:, :B] = shots[:, :B] + dS <= ctx.L
        if policy == "uniform":
            chosen = np.full(T, step % B)
        else:
            st = AllocationState(ctx.specs, ctx.n, counts, shots, gp.loglik, ctx.G, 0, prng, step)
            sig, _ = eig_cell_support(st)
            sig = np.where(E, sig, -np.inf)
            chosen = np.argmax(sig + prng.random(sig.shape) * 1e-12, axis=1)
        for b in np.unique(chosen):
            m = (chosen == b) & E[:, b]
            d = batch_counts(ctx.streams[b], np.minimum(shots[:, b], ctx.L - dS), dS, ctx.specs[b].M)
            d[~m] = 0
            counts[b] += d
            shots[m, b] += dS
            gp.add_counts(ctx.specs[b], d)
            u[m] += dS * ctx.u[b]
        U.append(u.copy())
        S.append(shots.sum(axis=1))
        for d, v in ctx.decode(counts, shots, gp.loglik).items():
            P[d].append(v)
    return np.array(U), np.array(S), {d: np.array(v) for d, v in P.items()}


def match_to_target(U, target):
    """Per trial: index of the last path point with U <= target (lower) and the first with U >= target (upper; -1 if none)."""
    K1, T = U.shape
    le = U <= target[None, :]
    lower = np.where(le.any(axis=0), K1 - 1 - np.argmax(le[::-1], axis=0), 0)
    ge = U >= target[None, :]
    upper = np.where(ge.any(axis=0), np.argmax(ge, axis=0), -1)
    return lower, upper


def run_arm(ctx, segments, A_max, uq_target=None, seed=0, record_actions=False, arm=""):
    T, B, Btot, dS = ctx.T, ctx.B, ctx.Btot, ctx.dS
    shots = np.zeros((T, Btot), dtype=np.int64)
    shots[:, :B] = ctx.S0
    counts = [batch_counts(ctx.streams[b], np.zeros(T, dtype=np.int64), ctx.S0, ctx.specs[b].M) if b < B
              else np.zeros((T, ctx.specs[b].M), dtype=np.int64) for b in range(Btot)]
    gp = GridPosterior(ctx.n, T, ctx.refine)
    for b in range(B):
        gp.add_counts(ctx.specs[b], counts[b])
    u_used = shots @ ctx.u
    n_ovl = np.zeros(T, dtype=np.int64)
    uni_ptr = np.zeros(T, dtype=np.int64)
    prng = generator(seed, stable_int(arm), 11)
    actions = []
    step = 0
    for policy, eligible, length in segments:
        it = 0
        while True:
            if length == "uq":
                active = u_used < uq_target
                if not active.any() or it > 50 * (ctx.steps + A_max):
                    break
            else:
                if it >= length:
                    break
                active = np.ones(T, dtype=bool)
            E = np.zeros((T, Btot), dtype=bool)
            if eligible in ("chunks", "all"):
                E[:, :B] = True
            if eligible in ("cands", "all"):
                E[:, B:] = (n_ovl < A_max)[:, None]
            E &= shots + dS <= ctx.L
            if policy == "uniform":
                chosen = uni_ptr % B
            elif policy in ("eig", "lowerhalf"):
                st = AllocationState(ctx.specs, ctx.n, counts, shots, gp.loglik, ctx.G, 0, prng, step)
                sig, _ = eig_cell_support(st)
                if policy == "lowerhalf":
                    bonus = np.zeros((T, Btot))
                    bonus[:, B:] = lowerhalf_candidates(ctx, counts) * 1e6  # structural trigger first; eig breaks ties
                    sig = sig + bonus
                sig = np.where(E, sig, -np.inf)
                chosen = np.argmax(sig + prng.random(sig.shape) * 1e-12, axis=1)
            elif policy.startswith("oracle:"):
                dec = policy.split(":", 1)[1]
                base = ctx.decode(counts, shots, gp.loglik)
                e0, t0 = ctx.score(base[dec])
                gt, ge = np.full((T, Btot), -9.0), np.full((T, Btot), -1e9)
                for b in range(Btot):
                    if not E[:, b].any():
                        continue
                    d = batch_counts(ctx.streams[b], np.minimum(shots[:, b], ctx.L - dS), dS, ctx.specs[b].M)
                    cb = [c + (d if i == b else 0) for i, c in enumerate(counts)]
                    sb = shots.copy()
                    sb[:, b] += dS
                    tmp = GridPosterior(ctx.n, 1, ctx.refine)
                    tmp.loglik = gp.loglik.copy()
                    tmp.add_counts(ctx.specs[b], d)
                    e1, t1 = ctx.score(ctx.decode(cb, sb, tmp.loglik)[dec])
                    gt[:, b] = np.where(E[:, b], t1.astype(int) - t0.astype(int), -9)
                    ge[:, b] = np.where(E[:, b], (e0 - e1) * 2.0**ctx.n, -1e9)
                chosen = greedy_oracle_choice(gt, ge, prng)
            else:
                raise ValueError(policy)
            ok = active & E[np.arange(T), chosen]
            if not ok.any():
                break
            before = ctx.decode(counts, shots, gp.loglik) if record_actions and (chosen[ok] >= B).any() else None
            pre_counts = [c.copy() for c in counts[:B]] if before is not None else None
            for b in np.unique(chosen[ok]):
                m = ok & (chosen == b)
                d = batch_counts(ctx.streams[b], np.minimum(shots[:, b], ctx.L - dS), dS, ctx.specs[b].M)
                d[~m] = 0
                counts[b] += d
                shots[m, b] += dS
                gp.add_counts(ctx.specs[b], d)
                u_used[m] += dS * ctx.u[b]
                if b >= B:
                    n_ovl[m] += 1
            if policy == "uniform":
                uni_ptr[ok] += 1
            if before is not None:
                after = ctx.decode(counts, shots, gp.loglik)
                ov = ok & (chosen >= B)
                bnd = np.array([ctx.cands[c - B][0] for c in chosen[ov]], dtype=np.int64)
                rec = {"trial": np.flatnonzero(ov), "step": step, "boundary": bnd + 1}
                # ambiguity state of the targeted boundary at action time (pre-batch counts), per epsilon
                for tag, eps in (("09", 0.9), ("safe", ctx.eps["awqpe_eps_safe"])):
                    fl, lh = boundary_states(ctx, pre_counts, eps)
                    rec[f"flag{tag}"] = fl[np.flatnonzero(ov), bnd]
                    rec[f"lowerhalf{tag}"] = lh[np.flatnonzero(ov), bnd]
                for dec in ctx.decoders:
                    e0, _ = ctx.score(before[dec])
                    e1, _ = ctx.score(after[dec])
                    rec[f"{dec}__gain_err_lsb"] = ((e0 - e1) * 2.0**ctx.n)[ov]
                actions.append(pd.DataFrame(rec))
            it += 1
            step += 1
    final = ctx.decode(counts, shots, gp.loglik)
    return final, shots, u_used, n_ovl, (pd.concat(actions, ignore_index=True) if actions else pd.DataFrame())


def overlap_shard(cfg: dict, spec: dict) -> pd.DataFrame:
    widths, S0, r, R = spec["widths"], int(spec["S0"]), int(spec["r"]), int(cfg["replicates"])
    dS = S0 if spec["dS"] == "S0" else int(spec["dS"])
    B, n = len(widths), int(sum(widths))
    steps = (r - 1) * S0 * B // dS
    phases = pd.DataFrame(spec["phases"])
    meta = phases.loc[phases.index.repeat(R)].reset_index(drop=True)
    meta["replicate_id"] = np.tile(np.arange(R), len(phases))
    phis = meta["phi"].to_numpy()
    seed = int(generator(cfg["master_seed"], stable_int(cfg["experiment_id"]), stable_int("-".join(map(str, widths))), S0, r,
                         stable_int(str(spec["dS"])), spec["chunk_index"]).integers(0, 2**62))
    A_all = max(int(v["A_max"]) for v in cfg["variants"])
    L = S0 + int(cfg.get("stream_multiplier", 6)) * (steps + A_all) * dS  # identical for every variant
    probs = [block_probabilities(phis, s.offset, s.width) for s in partition_blocks(widths)]
    limit_ok = {"likelihood": np.ones(len(phis), bool)}
    for name, eps in (("awqpe", 0.9), ("awqpe_eps_safe", safe_epsilon(widths))):
        limit_ok[name] = circular_error(awqpe_vectorised(probs, widths, eps)["estimate"] / 2.0**n, phis) <= 2.0**-n + 1e-15
        limit_ok[name + "_ext"] = limit_ok[name]
    rows, acts = [], []
    for variant in cfg["variants"]:
        ctx = Ctx(cfg, widths, S0, dS, steps, variant, phis, seed, L)
        A, trig, vlabel = int(variant["A_max"]), variant.get("trigger", "eig"), variant["label"]
        results = {
            "B1_uniform": run_arm(ctx, [("uniform", "chunks", steps)], A, seed=seed, arm="B1"),
            "B2_p5_eig": run_arm(ctx, [("eig", "chunks", steps)], A, seed=seed, arm="B2"),
            "B3_overlap_alone": run_arm(ctx, [(trig, "cands", A), ("uniform", "chunks", steps - A)], A, seed=seed, record_actions=True, arm="B3"),
            "B4_eig_plus_overlap": run_arm(ctx, [(trig, "all", steps)], A, seed=seed, record_actions=True, arm="B4"),
        }
        # Equal-U-query baselines (D-022): bracket each overlap arm's per-trial U cost on the baseline's own
        # path (same streams, same decisions as B1/B2, decoded after every batch).
        matched = {}
        for base_pol, base_arm, ovl_arm, tag in (("uniform", "B1", "B3_overlap_alone", "B1m"), ("eig", "B2", "B4_eig_plus_overlap", "B2m")):
            target = results[ovl_arm][2].astype(float)
            U, S, P = baseline_path(ctx, base_pol, target, seed, base_arm)
            lo, hi = match_to_target(U, target)
            cols = np.arange(len(phis))
            for side, idx in (("lower", lo), ("upper", hi)):
                valid = idx >= 0
                k = np.where(valid, idx, 0)
                matched[f"{tag}_{side}_{ovl_arm}"] = {
                    "phi": {d: np.where(valid, P[d][k, cols], np.nan) for d in P},
                    "u": np.where(valid, U[k, cols], np.nan), "shots": np.where(valid, S[k, cols], np.nan),
                    "u_target": target, "valid": valid}
        results["R_eig_shots"] = run_arm(ctx, [("eig", "chunks", steps), ("eig", "chunks", A)], A, seed=seed, arm="R1")
        results["R_trigger_overlap"] = run_arm(ctx, [("eig", "chunks", steps), (trig, "cands", A)], A, seed=seed, record_actions=True, arm="R2")
        for dec in cfg["rescue_oracle_decoders"]:
            if dec in ctx.decoders:
                results[f"R_oracle_shots[{dec}]"] = run_arm(ctx, [("eig", "chunks", steps), (f"oracle:{dec}", "chunks", A)], A, seed=seed, arm="R3" + dec)
                results[f"R_oracle_overlap[{dec}]"] = run_arm(ctx, [("eig", "chunks", steps), (f"oracle:{dec}", "cands", A)], A, seed=seed, arm="R4" + dec)
        for name, (final, shots, u_used, n_ovl, actions) in results.items():
            mask = np.zeros(len(phis), dtype=np.int64)
            for i, (j, _, _) in enumerate(ctx.cands):
                mask |= (shots[:, B + i] > 0).astype(np.int64) << j
            for dec in ctx.decoders:
                if name.startswith("R_oracle") and not name.endswith(f"[{dec}]"):
                    continue
                err, tol = ctx.score(final[dec])
                rows.append(pd.DataFrame({"trial": np.arange(len(phis)), "arm": name, "variant": vlabel, "mechanism": variant["mechanism"],
                                          "decoder": dec, "error": err, "tol": tol, "total_shots": shots.sum(axis=1), "u_queries": u_used,
                                          "chunk_shots": shots[:, :B].sum(axis=1), "overlap_shots": shots[:, B:].sum(axis=1),
                                          "n_overlap_actions": n_ovl, "overlap_boundary_mask": mask, "limit_correct": limit_ok[dec]}))
            if len(actions):
                acts.append(actions.assign(arm=name, variant=vlabel))
        for name, mres in matched.items():
            for dec in ctx.decoders:
                err, tol = ctx.score(np.nan_to_num(mres["phi"][dec]))
                rows.append(pd.DataFrame({"trial": np.arange(len(phis)), "arm": name, "variant": vlabel, "mechanism": variant["mechanism"],
                                          "decoder": dec, "error": np.where(mres["valid"], err, np.nan),
                                          "tol": np.where(mres["valid"], tol.astype(float), np.nan),
                                          "total_shots": mres["shots"], "u_queries": mres["u"], "u_target": mres["u_target"],
                                          "u_diff_vs_target": mres["u"] - mres["u_target"], "match_valid": mres["valid"],
                                          "chunk_shots": mres["shots"], "overlap_shots": 0, "n_overlap_actions": 0,
                                          "overlap_boundary_mask": 0, "limit_correct": limit_ok[dec]}))
    out = pd.concat(rows, ignore_index=True).join(meta, on="trial")
    out["tol"] = out["tol"].astype(float)  # NaN = no valid U-match (upper bracket beyond the path cap)
    out["match_valid"] = out["match_valid"].astype("boolean") if "match_valid" in out else pd.NA
    out["row_type"] = "final"
    if acts:
        a = pd.concat(acts, ignore_index=True).join(meta[["phase_id", "stratum", "replicate_id"]], on="trial")
        a["row_type"] = "action"
        out = pd.concat([out, a], ignore_index=True, sort=False)
    out["widths"], out["n"], out["B"], out["S0"], out["r"], out["dS"], out["split"] = "-".join(map(str, widths)), n, B, S0, r, dS, spec["split"]
    return out


def shards(cfg):
    out = []
    for widths in cfg["partitions"]:
        table = make_phase_table(widths, int(cfg["phases_per_stratum"]), cfg["split"], int(cfg["master_seed"]), strata=tuple(cfg["strata"]))
        records = table[["phase_id", "stratum", "phi", "final_residual", "boundary_hardness"]].to_dict("records")
        size = int(cfg.get("phase_chunk", 16))
        for S0 in cfg["S0"]:
            for r in cfg["budget_multipliers"]:
                for dS in cfg["dS"]:
                    for ci, a in enumerate(range(0, len(records), size)):
                        out.append({"shard_id": f"w{'-'.join(map(str, widths))}_S{S0}_r{r}_d{dS}_c{ci}", "widths": widths, "S0": S0, "r": r,
                                    "dS": dS, "chunk_index": ci, "split": cfg["split"], "phases": records[a:a + size]})
    return out


def main() -> None:
    args = common_arguments(__doc__, "p6_overlap_pilot.yaml").parse_args()
    cfg = resolve_config(args)
    run_sharded(args, cfg, shards(cfg), overlap_shard, extra_manifest={"status": cfg.get("status", "PILOT / DEV")})


if __name__ == "__main__":
    main()
