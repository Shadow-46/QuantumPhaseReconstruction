"""Phase 10a: Qiskit Aer shot-sampling validation of the Dirichlet-kernel model (ideal, no noise).

Part A (distribution): for block circuits U = P(2 pi phi) on |1> (the paper's Fig. 1
construction), Aer shot counts vs the kernel probabilities: chi-square goodness of fit
(bins pooled to expected >= 5), total-variation distance, and top-1 agreement with the
kernel's most probable outcome.

Part B (end-to-end): AWQPE (eps 0.9 and eps_safe) and D2 decoded from Aer counts vs from
kernel-sampled counts, same phases, independent shots; success rates compared per cell
(two-proportion z) -- the claim tested is statistical equivalence of the two samplers.

    python -m research.awqpe.run_qiskit_validation [--pilot] [--max-workers N]
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import chi2

from research.awqpe.baseline.awqpe import awqpe_vectorised
from research.awqpe.blocks.geometry import partition_blocks
from research.awqpe.circuits.qiskit_blocks import phase_gate_block
from research.awqpe.decode.likelihood import likelihood_decode
from research.awqpe.evaluation.metrics import circular_error
from research.awqpe.model.kernel import block_probabilities
from research.awqpe.phases import make_phase_table
from research.awqpe.runner.core import common_arguments, resolve_config, run_sharded
from research.awqpe.sim.oracle_sim import batch_block_counts
from research.awqpe.sim.seeds import generator, stable_int
from research.awqpe.verification.theory import safe_epsilon


def aer_counts(phi: float, offset: int, width: int, shots: int, seed: int) -> np.ndarray:
    from qiskit import transpile
    from qiskit_aer import AerSimulator

    sim = AerSimulator(seed_simulator=seed)
    qc = transpile(phase_gate_block(phi, offset, width), sim)
    res = sim.run(qc, shots=shots).result().get_counts()
    out = np.zeros(1 << width, dtype=np.int64)
    for key, c in res.items():
        out[int(key.replace(" ", ""), 2)] += c
    return out


def chisq_pooled(obs: np.ndarray, p: np.ndarray) -> tuple[float, int]:
    exp = p * obs.sum()
    order = np.argsort(exp)
    o, e = obs[order].astype(float), exp[order]
    bo, be, acc_o, acc_e = [], [], 0.0, 0.0
    for oi, ei in zip(o, e):
        acc_o += oi
        acc_e += ei
        if acc_e >= 5:
            bo.append(acc_o)
            be.append(acc_e)
            acc_o = acc_e = 0.0
    if acc_e > 0 and be:
        bo[-1] += acc_o
        be[-1] += acc_e
    bo, be = np.array(bo), np.array(be)
    df = len(be) - 1
    if df < 1:
        return 1.0, 0
    stat = float(((bo - be) ** 2 / be).sum())
    return float(chi2.sf(stat, df)), df


def dist_shard(cfg, spec):
    rows = []
    for rec in spec["phases"]:
        for m in cfg["widths"]:
            for k in cfg["offsets"]:
                seed = int(generator(cfg["master_seed"], stable_int(rec["phase_id"]), m, k).integers(0, 2**31))
                obs = aer_counts(rec["phi"], k, m, int(cfg["dist_shots"]), seed)
                p = block_probabilities(rec["phi"], k, m)
                pv, df = chisq_pooled(obs, p)
                rows.append({"phase_id": rec["phase_id"], "stratum": rec["stratum"], "phi": rec["phi"], "width": m, "offset": k,
                             "shots": int(cfg["dist_shots"]), "chisq_p": pv, "df": df, "tv": 0.5 * float(np.abs(obs / obs.sum() - p).sum()),
                             "top1_aer": int(np.argmax(obs)), "top1_kernel": int(np.argmax(p)),
                             "kernel_top2_ratio": float(np.sort(p)[-2] / np.sort(p)[-1])})
    df_ = pd.DataFrame(rows)
    df_["part"] = "distribution"
    return df_


def e2e_shard(cfg, spec):
    widths, shots, R = spec["widths"], int(spec["shots"]), int(cfg["e2e_replicates"])
    n = sum(widths)
    specs = partition_blocks(widths)
    ph = pd.DataFrame(spec["phases"])
    meta = ph.loc[ph.index.repeat(R)].reset_index(drop=True)
    meta["replicate_id"] = np.tile(np.arange(R), len(ph))
    phis = meta.phi.to_numpy()
    rng = generator(cfg["master_seed"], stable_int("-".join(map(str, widths))), shots, spec["chunk_index"])
    aer = [np.stack([aer_counts(float(phi), s.offset, s.width, shots, int(rng.integers(0, 2**31))) for phi in phis]) for s in specs]
    ker = [batch_block_counts(phis, s, shots, rng) for s in specs]
    frames = []
    for sampler, counts in (("aer", aer), ("kernel", ker)):
        jit = [rng.random(c.shape) * 0.5 for c in counts]
        est = {"awqpe": awqpe_vectorised(counts, widths, 0.9, jitter=jit)["estimate"] / 2.0**n,
               "awqpe_eps_safe": awqpe_vectorised(counts, widths, safe_epsilon(widths), jitter=jit)["estimate"] / 2.0**n,
               "likelihood": likelihood_decode(specs, counts, n, 3)}
        for dec, e in est.items():
            d = meta.copy()
            d["sampler"], d["decoder"] = sampler, dec
            d["error"] = circular_error(e, phis)
            d["tol"] = d["error"] <= 2.0**-n + 1e-15
            frames.append(d)
    out = pd.concat(frames, ignore_index=True)
    out["widths"], out["n"], out["shots_per_block"], out["part"] = "-".join(map(str, widths)), n, shots, "end_to_end"
    return out


def qv_shard(cfg, spec):
    return dist_shard(cfg, spec) if spec["part"] == "distribution" else e2e_shard(cfg, spec)


def shards(cfg):
    out = []
    tab = make_phase_table([4, 4], int(cfg["phases_per_stratum"]), "dev", int(cfg["master_seed"]))
    recs = tab[["phase_id", "stratum", "phi"]].to_dict("records")
    for ci, a in enumerate(range(0, len(recs), 10)):
        out.append({"shard_id": f"dist_c{ci}", "part": "distribution", "phases": recs[a:a + 10]})
    for widths in cfg["e2e_partitions"]:
        t = make_phase_table(widths, int(cfg["e2e_phases_per_stratum"]), "dev", int(cfg["master_seed"]))
        r2 = t[["phase_id", "stratum", "phi"]].to_dict("records")
        for shots in cfg["e2e_shots"]:
            for ci, a in enumerate(range(0, len(r2), 10)):
                out.append({"shard_id": f"e2e_w{'-'.join(map(str, widths))}_s{shots}_c{ci}", "part": "e2e", "widths": widths,
                            "shots": shots, "chunk_index": ci, "phases": r2[a:a + 10]})
    return out


def main() -> None:
    args = common_arguments(__doc__, "p10_qiskit_validation.yaml").parse_args()
    cfg = resolve_config(args)
    run_sharded(args, cfg, shards(cfg), qv_shard, extra_manifest={"status": cfg.get("status", "VALIDATION / DEV")})


if __name__ == "__main__":
    main()
