"""Phase 3: initial-state / eigenphase sensitivity (PRELIMINARY, dev phases).

Separates four things that Paper A conflated:
  state preparation  exact eigenstate vs superposition (mixture of eigencomponents)
  eigenphase         which phase(s) the input carries
  grid location      dyadic at n bits / near-grid / boundary-hard (recorded per component)
  ambiguity          exact top-two ratio of each block's actual (mixture) distribution

Families (all states come from explicit unitaries and are verified):
  A phase_gate       U = P(2 pi phi), |1>                          (eigenstate)
  B diagonal         3-qubit diagonal U, basis states                (eigenstates)
  C conjugated       V D V^dagger with the same D as B, V|j>         (eigenstates)
  D modmul_eigen     U_a with |u_s>, (15,2) r=4, (21,2) r=6, (33,2) r=10 (eigenstates)
  E modmul_one       U_a with computational |1>: the Paper A input  (uniform mixture over s/r)
  F controlled_mix   dominant phi0 with weight w0 in {1, .9, .75, .5} (mixtures)

Every case records ||U psi - e^{2 pi i phi} psi|| (eigenstates) and the max
|statevector - mixture-of-kernels| of real Qiskit block circuits that
prepare psi and apply controlled U^(2^(k+p)).

    python -m research.awqpe.run_state_sensitivity [--pilot] [--max-workers N]
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from research.awqpe.baseline.awqpe import awqpe_vectorised
from research.awqpe.blocks.geometry import BlockSpec, partition_blocks
from research.awqpe.circuits import unitaries as UN
from research.awqpe.circuits.qiskit_blocks import block_statevector_probabilities, general_unitary_block
from research.awqpe.decode.likelihood import likelihood_decode
from research.awqpe.evaluation.metrics import boundary_hardness, circular_error, final_residual
from research.awqpe.model.kernel import block_delta, block_probabilities
from research.awqpe.phases import make_phase_table
from research.awqpe.runner.core import common_arguments, resolve_config, run_sharded
from research.awqpe.sim.oracle_sim import batch_block_counts
from research.awqpe.sim.seeds import generator, stable_int
from research.awqpe.verification.theory import safe_epsilon


def mixture_block_probs(phases, weights, spec) -> np.ndarray:
    return np.asarray(weights) @ block_probabilities(np.asarray(phases), spec.offset, spec.width)


def qiskit_check(inp: UN.SpectralInput, checks) -> float:
    """max |Qiskit statevector - mixture of kernels| over the given (offset, width) blocks."""
    worst = 0.0
    for k, m in checks:
        sv = block_statevector_probabilities(general_unitary_block(inp.unitary, inp.state, k, m, measure=False), m)
        worst = max(worst, float(np.abs(sv - mixture_block_probs(inp.phases, inp.weights, BlockSpec(k, m))).max()))
    return worst


def build_cases(cfg) -> list[dict]:
    checks = [tuple(c) for c in cfg["qiskit_checks"]]
    rng = generator(cfg["master_seed"], stable_int("p3_cases"))
    cases = []

    def add(family, inp: UN.SpectralInput, instance="", dominant=None):
        order = np.argsort(-inp.weights, kind="stable")
        dom = float(inp.phases[order[0]]) if dominant is None else float(dominant)
        cases.append({
            "family": family, "label": inp.label, "instance": instance,
            "state_prep": "eigenstate" if inp.is_eigenstate else "mixture",
            "phases": [float(x) for x in inp.phases], "weights": [float(x) for x in inp.weights],
            "dominant_phase": dom, "dominant_weight": float(inp.weights[order[0]]),
            "eigen_residual": inp.eigen_residual,
            "qiskit_max_abs_diff": qiskit_check(inp, checks),
        })

    table = make_phase_table(cfg["phase_partition"], int(cfg["phases_per_stratum"]), "dev", int(cfg["master_seed"]),
                             strata=("S0_dyadic", "S1_final_half", "S2_boundary", "S3_near_grid", "S4_uniform"))
    for phi in list(table.phi) + [0.25, 1 / 6, 0.1]:
        add("A_phase_gate", UN.phase_gate_input(float(phi)))
    diag_phases = cfg["diagonal_phases"]
    for j in range(len(diag_phases)):
        add("B_diagonal", UN.diagonal_input(diag_phases, j))
        add("C_conjugated", UN.conjugated_input(diag_phases, j, rng))
    for a, N in cfg["modmul_instances"]:
        r = UN.multiplicative_order(a, N)
        for s in range(r):
            add("D_modmul_eigen", UN.modmul_eigenstate_input(a, N, s), instance=f"({N},{a}) r={r}")
        add("E_modmul_one", UN.modmul_one_input(a, N), instance=f"({N},{a}) r={r}", dominant=1.0 / r)
    for phi0 in cfg["mixture_phi0"]:
        for w0 in cfg["mixture_w0"]:
            others_random = [float(x) for x in np.round(rng.random(3), 12)]
            others_shor = [float(np.mod(phi0 + j / 6, 1.0)) for j in range(1, 6)]
            for kind, others in (("random3", others_random), ("shor_like6", others_shor)):
                inp = UN.controlled_mixture_input(float(phi0), float(w0), others)
                add("F_controlled_mix", inp, instance=f"{kind} w0={w0}", dominant=phi0)
    for i, c in enumerate(cases):
        c["case_id"] = f"{c['family']}:{i}"
    return cases


def case_features(case, widths) -> dict:
    n = int(sum(widths))
    phi = case["dominant_phase"]
    ratios = []
    for s in partition_blocks(widths):
        p = np.sort(mixture_block_probs(case["phases"], case["weights"], s))[::-1]
        ratios.append(p[1] / p[0])
    return {
        "n_components": len(case["phases"]),
        "is_dyadic_n": bool(abs(phi * 2**n - round(phi * 2**n)) < 1e-9),
        "final_residual": float(final_residual(phi, n)),
        "boundary_hardness": float(boundary_hardness(phi, widths)),
        "true_min_top2_ratio": float(min(ratios)), "true_max_top2_ratio": float(max(ratios)),
    }


def state_shard(cfg: dict, spec: dict) -> pd.DataFrame:
    widths, shots, R = spec["widths"], int(spec["shots"]), int(cfg["replicates"])
    n = int(sum(widths))
    specs = partition_blocks(widths)
    eps_list = [("awqpe", 0.9), ("awqpe_eps_safe", safe_epsilon(widths))]
    frames = []
    for case in spec["cases"]:
        rng = generator(cfg["master_seed"], stable_int(case["case_id"]), stable_int("-".join(map(str, widths))), shots)
        J = len(case["phases"])
        phases = np.tile(np.asarray(case["phases"]), (R, 1))
        weights = np.tile(np.asarray(case["weights"]), (R, 1))
        counts = [batch_block_counts(phases if J > 1 else phases[:, 0], s, shots, rng, weights=weights if J > 1 else None) for s in specs]
        jitter = [rng.random(c.shape) * 0.5 for c in counts]
        results = []
        for name, eps in eps_list:
            r = awqpe_vectorised(counts, widths, eps, jitter=jitter)
            results.append((name, eps, r["estimate"] / 2.0**n, r))
        results.append(("likelihood", np.nan, likelihood_decode(specs, counts, n, int(cfg.get("likelihood_refine", 4))), None))
        comps = np.asarray(case["phases"])[np.asarray(case["weights"]) >= 0.05]
        t1_dom = [int(np.mod(np.floor(block_delta(case["dominant_phase"], s.offset) * s.M + 0.5), s.M)) for s in specs]
        feats = case_features(case, widths)
        for name, eps, phi_hat, r in results:
            err_dom = circular_error(phi_hat, case["dominant_phase"])
            err_any = np.min(circular_error(phi_hat[:, None], comps[None, :]), axis=1)
            df = pd.DataFrame({
                "replicate_id": np.arange(R), "decoder": name, "epsilon": eps, "phi_hat": phi_hat,
                "error_dominant": err_dom, "error_any": err_any,
                "tol_dominant": err_dom <= 2.0**-n + 1e-15, "tol_any": err_any <= 2.0**-n + 1e-15,
            })
            if r is not None:
                df["flags_mask"] = (r["flags"] * (1 << np.arange(len(widths)))).sum(axis=1)
                df["blocks_matching_dominant"] = (r["top1"] == np.asarray(t1_dom)[None, :]).sum(axis=1)
            else:
                df["flags_mask"] = -1
                df["blocks_matching_dominant"] = -1
            for key in ("case_id", "family", "label", "instance", "state_prep", "dominant_phase", "dominant_weight", "eigen_residual", "qiskit_max_abs_diff"):
                df[key] = case[key]
            for key, val in feats.items():
                df[key] = val
            df["widths"], df["n"], df["shots_per_block"] = "-".join(map(str, widths)), n, shots
            frames.append(df)
    return pd.concat(frames, ignore_index=True)


def main() -> None:
    args = common_arguments(__doc__, "p3_state_sensitivity.yaml").parse_args()
    cfg = resolve_config(args)
    cases = build_cases(cfg)
    size = int(cfg.get("cases_per_shard", 20))
    shards = []
    for widths in cfg["partitions"]:
        for shots in cfg["shots_per_block"]:
            for ci, a in enumerate(range(0, len(cases), size)):
                shards.append({"shard_id": f"w{'-'.join(map(str, widths))}_s{shots}_c{ci}", "widths": widths, "shots": int(shots), "cases": cases[a:a + size]})
    summary = [{k: v for k, v in c.items() if k not in ("phases", "weights")} | {"n_components": len(c["phases"])} for c in cases]
    run_sharded(args, cfg, shards, state_shard, extra_manifest={"status": "PRELIMINARY (dev phases)", "case_verification": summary})


if __name__ == "__main__":
    main()
