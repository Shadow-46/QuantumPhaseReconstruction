"""Phase 10b: Qiskit Aer circuit-level noise vs the analytic channels of model/noise.py (D-034).

Part R (readout): Aer ReadoutError (symmetric p) on every measured control vs apply_readout.
          Chi-square goodness of fit of Aer counts to the analytic noisy distribution.
Part J (jitter): control-phase jitter: Aer exact probabilities averaged over K Gaussian offsets
          xi ~ N(0, sigma^2) (block run at phi + xi / 2^k), compared bin by bin (Monte Carlo z) with the
          analytic Fejer damping (jitter_mode "fixed").
Part D (depolarisation): U applied by REPETITION (2^(k+p) labelled controlled-phase gates from
          control p, q = 2^k (2^m - 1) in total), each followed by a two-qubit depolarising error
          p_g; H and inverse QFT noiseless. Exact density-matrix probabilities. The analytic law
          p -> (1 - lam) p + lam / M is fitted per block (least squares lam_hat); reported: the
          residual TV of the best global-depolarising fit, and how lam_hat scales with q
          (the analytic "per_query" law predicts lam = 1 - exp(-gamma q)).

    python -m research.awqpe.run_aer_noise [--pilot] [--max-workers N]
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import QFTGate, UnitaryGate

from research.awqpe.circuits.qiskit_blocks import phase_gate_block
from research.awqpe.model.noise import NoiseSpec, noisy_block_probabilities
from research.awqpe.phases import make_phase_table
from research.awqpe.run_qiskit_validation import chisq_pooled
from research.awqpe.runner.core import common_arguments, resolve_config, run_sharded
from research.awqpe.sim.seeds import generator, stable_int

QUERY_LABEL = "uq"


def _counts(res, width):
    out = np.zeros(1 << width, dtype=np.int64)
    for key, c in res.items():
        out[int(key.replace(" ", ""), 2)] += c
    return out


def readout_row(phi, k, m, p, shots, seed):
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel, ReadoutError

    nm = NoiseModel()
    nm.add_all_qubit_readout_error(ReadoutError([[1 - p, p], [p, 1 - p]]))
    sim = AerSimulator(noise_model=nm, seed_simulator=seed)
    obs = _counts(sim.run(transpile(phase_gate_block(phi, k, m), sim), shots=shots).result().get_counts(), m)
    ref = noisy_block_probabilities(phi, k, m, NoiseSpec(readout_p01=p, readout_p10=p))
    return obs, ref


def jitter_row(phi, k, m, sigma, K, seed):
    """Aer exact (statevector) block probabilities averaged over K Gaussian phase offsets vs the
    analytic damping law. Pooling counts over a finite set of offsets would test the empirical
    K-point mixture rather than the Gaussian one, so exact per-draw probabilities are averaged and
    compared bin by bin with their Monte Carlo standard errors."""
    from qiskit_aer import AerSimulator

    rng = np.random.default_rng(seed)
    sim = AerSimulator(method="statevector")
    circs = []
    for xi in rng.normal(0.0, sigma, K):
        qc = phase_gate_block(phi + xi / 2.0**k, k, m, measure=False)
        qc.save_probabilities(list(range(m)))
        circs.append(qc)
    res = sim.run(transpile(circs, sim)).result()
    P = np.array([res.data(i)["probabilities"] for i in range(K)], dtype=float)
    ref = noisy_block_probabilities(phi, k, m, NoiseSpec(jitter_sigma=sigma, jitter_mode="fixed"))
    se = P.std(axis=0, ddof=1) / np.sqrt(K)
    z = (P.mean(axis=0) - ref) / np.maximum(se, 1e-12)
    return float(np.abs(z).max()), 0.5 * float(np.abs(P.mean(axis=0) - ref).sum())


def repeated_block(phi, k, m):
    """Block circuit with U = P(2 pi phi) applied 2^(k+p) times from control p (labelled gates)."""
    cp = UnitaryGate(np.diag([1, 1, 1, np.exp(2j * np.pi * phi)]), label=QUERY_LABEL)
    qc = QuantumCircuit(m + 1)
    qc.x(m)
    qc.h(range(m))
    for p in range(m):
        for _ in range(1 << (k + p)):
            qc.append(cp, [p, m])
    qc.append(QFTGate(m).inverse(), list(range(m)))
    return qc


def depol_row(phi, k, m, p_g):
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel, depolarizing_error

    nm = NoiseModel(basis_gates=["unitary", "h", "x", "cp", "swap", "u", "cx"])
    for p in range(m):
        nm.add_quantum_error(depolarizing_error(p_g, 2), QUERY_LABEL, [p, m])
    sim = AerSimulator(method="density_matrix", noise_model=nm)
    qc = repeated_block(phi, k, m)
    qc.save_probabilities(list(range(m)))
    qc = transpile(qc, sim, optimization_level=0)
    probs = np.asarray(sim.run(qc).result().data()["probabilities"], dtype=float)
    ideal = noisy_block_probabilities(phi, k, m, None)
    M = 1 << m
    u = np.full(M, 1.0 / M)
    a = u - ideal  # probs ~ ideal + lam * (u - ideal)
    lam = float(np.clip(np.dot(probs - ideal, a) / max(np.dot(a, a), 1e-300), 0.0, 1.0))
    fit = (1 - lam) * ideal + lam * u
    return probs, ideal, lam, 0.5 * float(np.abs(probs - fit).sum()), 0.5 * float(np.abs(probs - ideal).sum())


def aer_noise_shard(cfg, spec):
    rows = []
    for rec in spec["phases"]:
        phi = float(rec["phi"])
        seed = int(generator(cfg["master_seed"], stable_int(rec["phase_id"]), stable_int(spec["part"])).integers(0, 2**31))
        new = []
        if spec["part"] == "readout":
            for (k, m) in cfg["rj_blocks"]:
                for p in cfg["readout_p"]:
                    obs, ref = readout_row(phi, k, m, p, int(cfg["shots"]), seed)
                    pv, df = chisq_pooled(obs, ref)
                    new.append({"param": p, "offset": k, "width": m, "chisq_p": pv, "df": df, "tv": 0.5 * float(np.abs(obs / obs.sum() - ref).sum())})
        elif spec["part"] == "jitter":
            for (k, m) in cfg["rj_blocks"]:
                for s in cfg["jitter_sigma"]:
                    zmax, tv = jitter_row(phi, k, m, s, int(cfg["jitter_draws"]), seed)
                    new.append({"param": s, "offset": k, "width": m, "max_abs_z": zmax, "tv": tv})
        else:
            for (k, m) in cfg["depol_blocks"]:
                for pg in cfg["depol_pg"]:
                    _, _, lam, fit_tv, raw_tv = depol_row(phi, k, m, pg)
                    new.append({"param": pg, "offset": k, "width": m, "q": (1 << k) * ((1 << m) - 1), "lambda_hat": lam, "fit_tv": fit_tv, "noise_tv": raw_tv})
        for r in new:
            r["phase_id"], r["stratum"], r["phi"] = rec["phase_id"], rec["stratum"], phi
        rows += new
    out = pd.DataFrame(rows)
    out["part"] = spec["part"]
    return out


def shards(cfg):
    t = make_phase_table([4, 4], int(cfg["phases_per_stratum"]), "dev", int(cfg["master_seed"]), strata=("S4_uniform", "S1_final_half"))
    recs = t[["phase_id", "stratum", "phi"]].to_dict("records")
    return [{"shard_id": f"{part}_c{ci}", "part": part, "phases": recs[a:a + 2]} for part in ("readout", "jitter", "depol")
            for ci, a in enumerate(range(0, len(recs), 2))]


def main() -> None:
    args = common_arguments(__doc__, "p10b_aer_noise.yaml").parse_args()
    cfg = resolve_config(args)
    run_sharded(args, cfg, shards(cfg), aer_noise_shard, extra_manifest={"status": cfg.get("status", "VALIDATION / DEV")})


if __name__ == "__main__":
    main()
