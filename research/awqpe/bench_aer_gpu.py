"""Qiskit Aer CPU vs GPU benchmark on this project's actual AWQPE circuit workloads (infrastructure only).

Runs the circuits of P10a/P10b exactly as those phases build them (phase_gate_block; the
P10b repeated-query block with the same labelled two-qubit depolarising noise model) on
AerSimulator(device="CPU") and AerSimulator(device="GPU"), with identical circuits, shots and
seeds. It produces timings and a CPU-vs-GPU correctness check; it changes no scientific
result, parameter or frozen artefact.

Workloads (as used by the phases):
  p10a_e2e    one run() call per block circuit (as run_qiskit_validation.aer_counts), shots/block
  p10a_dist   one run() call per block circuit, 4096 shots, widths 2-5, offsets 0/2/5
  p10b_jitter one batched run() of K statevector circuits with save_probabilities (as jitter_row)
  p10b_depol  one run() per repeated-query density-matrix circuit with noise (as depol_row)

Correctness (CPU vs GPU): exact-probability workloads compare max |p_cpu - p_gpu|; counts
workloads compare both devices' counts with the kernel (pooled chi-square p-values) and with
each other (TV). Timing: build (circuit construction), transpile, execute (run+result), total;
each measured `repeat` times after one warm-up call per device (warm-up reported separately).

Run from WSL:
    PYTHONPATH=<repo> python -m research.awqpe.bench_aer_gpu --pilot --out bench.json
"""

from __future__ import annotations

import argparse
import json
import platform
import resource
import subprocess
import threading
import time

import numpy as np
from qiskit import transpile

from research.awqpe.circuits.qiskit_blocks import phase_gate_block
from research.awqpe.model.kernel import block_probabilities
from research.awqpe.run_aer_noise import QUERY_LABEL, repeated_block
from research.awqpe.run_qiskit_validation import chisq_pooled


class GpuMemPoller:
    """Peak GPU memory used (MiB, whole device) sampled with nvidia-smi while a block runs."""

    def __init__(self, period=0.05):
        self.period, self.peak, self._stop = period, 0, threading.Event()

    def _read(self):
        try:
            out = subprocess.run(["nvidia-smi", "--query-gpu=memory.used", "--format=csv,noheader,nounits"],
                                 capture_output=True, text=True, timeout=5).stdout
            return int(out.strip().splitlines()[0])
        except Exception:
            return -1

    def __enter__(self):
        self.base = self._read()
        self._t = threading.Thread(target=self._loop, daemon=True)
        self._t.start()
        return self

    def _loop(self):
        while not self._stop.is_set():
            self.peak = max(self.peak, self._read())
            time.sleep(self.period)

    def __exit__(self, *a):
        self._stop.set()
        self._t.join()


def _counts(res, i, m):
    out = np.zeros(1 << m, dtype=np.int64)
    for key, c in res.get_counts(i).items():
        out[int(key.replace(" ", ""), 2)] += c
    return out


def sim_for(device, method="automatic", noise_model=None, seed=None):
    from qiskit_aer import AerSimulator

    kw = {"device": device, "method": method}
    if noise_model is not None:
        kw["noise_model"] = noise_model
    if seed is not None:
        kw["seed_simulator"] = seed
    return AerSimulator(**kw)


def depol_noise_model(m, p_g):
    """Same noise model as run_aer_noise.depol_row (mirrored here only to pass `device`)."""
    from qiskit_aer.noise import NoiseModel, depolarizing_error

    nm = NoiseModel(basis_gates=["unitary", "h", "x", "cp", "swap", "u", "cx"])
    for p in range(m):
        nm.add_quantum_error(depolarizing_error(p_g, 2), QUERY_LABEL, [p, m])
    return nm


def phases(N, seed=20261003):
    return np.random.default_rng(seed).random(N)


def w_counts_per_call(device, blocks, shots, seed0):
    """P10a pattern: one AerSimulator + transpile + run per circuit. blocks = [(phi, k, m)]."""
    t_build = t_tr = t_ex = 0.0
    counts = []
    for i, (phi, k, m) in enumerate(blocks):
        a = time.perf_counter()
        qc = phase_gate_block(phi, k, m)
        b = time.perf_counter()
        sim = sim_for(device, seed=seed0 + i)
        tq = transpile(qc, sim)
        c = time.perf_counter()
        res = sim.run(tq, shots=shots).result()
        d = time.perf_counter()
        counts.append(_counts(res, 0, m))
        t_build, t_tr, t_ex = t_build + b - a, t_tr + c - b, t_ex + d - c
    return {"build_s": t_build, "transpile_s": t_tr, "exec_s": t_ex}, counts


def w_jitter_batch(device, phi, k, m, K, seed):
    """P10b jitter pattern: one batched statevector run of K circuits with save_probabilities."""
    rng = np.random.default_rng(seed)
    a = time.perf_counter()
    circs = []
    for xi in rng.normal(0.0, 0.04, K):
        qc = phase_gate_block(phi + xi / 2.0**k, k, m, measure=False)
        qc.save_probabilities(list(range(m)))
        circs.append(qc)
    b = time.perf_counter()
    sim = sim_for(device, method="statevector")
    tq = transpile(circs, sim)
    c = time.perf_counter()
    res = sim.run(tq).result()
    d = time.perf_counter()
    P = np.array([res.data(i)["probabilities"] for i in range(K)], dtype=float)
    return {"build_s": b - a, "transpile_s": c - b, "exec_s": d - c}, P


def w_depol(device, cases, p_g):
    """P10b depolarisation pattern: one density-matrix run per repeated-query circuit."""
    t_build = t_tr = t_ex = 0.0
    probs = []
    for phi, k, m in cases:
        a = time.perf_counter()
        qc = repeated_block(phi, k, m)
        qc.save_probabilities(list(range(m)))
        b = time.perf_counter()
        sim = sim_for(device, method="density_matrix", noise_model=depol_noise_model(m, p_g))
        tq = transpile(qc, sim, optimization_level=0)
        c = time.perf_counter()
        res = sim.run(tq).result()
        d = time.perf_counter()
        probs.append(np.asarray(res.data()["probabilities"], dtype=float))
        t_build, t_tr, t_ex = t_build + b - a, t_tr + c - b, t_ex + d - c
    return {"build_s": t_build, "transpile_s": t_tr, "exec_s": t_ex}, probs


def workloads(pilot):
    N = 4 if pilot else 24
    ph = phases(N)
    wl = []
    for widths in ([3, 2, 3], [4, 4], [4, 4, 4, 4]):
        ks = np.cumsum([0] + widths[:-1])
        blocks = [(float(p), int(k), int(m)) for p in ph for k, m in zip(ks, widths)]
        for shots in ((16,) if pilot else (16, 1024)):
            wl.append(dict(name="p10a_e2e", n=sum(widths), widths="-".join(map(str, widths)), width=max(widths),
                           shots=shots, n_circuits=len(blocks), kind="counts", blocks=blocks))
    dist = [(float(p), k, m) for p in ph[: max(2, N // 4)] for m in (2, 3, 4, 5) for k in (0, 2, 5)]
    wl.append(dict(name="p10a_dist", n=None, widths="2..5", width=5, shots=4096, n_circuits=len(dist), kind="counts", blocks=dist))
    K = 200 if pilot else 2000
    wl.append(dict(name="p10b_jitter", n=None, widths="4", width=4, shots=0, n_circuits=K, kind="jitter", phi=float(ph[0]), k=2, m=4, K=K))
    cases = [(float(ph[0]), k, m) for k in range(4) for m in (2, 3, 4)]
    wl.append(dict(name="p10b_depol", n=None, widths="2..4", width=4, shots=0, n_circuits=len(cases), kind="depol", cases=cases, p_g=0.01))
    return wl


def run_one(w, device):
    rss0 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    with GpuMemPoller() as pm:
        t0 = time.perf_counter()
        if w["kind"] == "counts":
            t, out = w_counts_per_call(device, w["blocks"], w["shots"], 1000)
        elif w["kind"] == "jitter":
            t, out = w_jitter_batch(device, w["phi"], w["k"], w["m"], w["K"], 7)
        else:
            t, out = w_depol(device, w["cases"], w["p_g"])
        t["total_s"] = time.perf_counter() - t0
    t["rss_peak_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    t["rss_start_mb"] = rss0
    t["gpu_mem_peak_mib"], t["gpu_mem_base_mib"] = pm.peak, pm.base
    return t, out


def correctness(w, out_cpu, out_gpu):
    if w["kind"] == "counts":
        pc, pg, tv = [], [], []
        for (phi, k, m), a, b in zip(w["blocks"], out_cpu, out_gpu):
            ref = block_probabilities(phi, k, m)
            for arr, lst in ((a, pc), (b, pg)):
                p, df = chisq_pooled(arr, ref)
                if df > 0:
                    lst.append(p)
            tv.append(0.5 * float(np.abs(a / a.sum() - b / b.sum()).sum()))
        return {"frac_p_lt_005_cpu_vs_kernel": float(np.mean(np.array(pc) < 0.05)) if pc else None,
                "frac_p_lt_005_gpu_vs_kernel": float(np.mean(np.array(pg) < 0.05)) if pg else None,
                "tested_blocks": len(pg), "mean_tv_cpu_gpu": float(np.mean(tv)), "identical_counts": bool(all((a == b).all() for a, b in zip(out_cpu, out_gpu)))}
    if w["kind"] == "jitter":
        return {"max_abs_prob_diff": float(np.abs(out_cpu - out_gpu).max())}
    return {"max_abs_prob_diff": float(max(np.abs(a - b).max() for a, b in zip(out_cpu, out_gpu)))}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pilot", action="store_true")
    ap.add_argument("--repeat", type=int, default=3)
    ap.add_argument("--devices", default="CPU,GPU")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    import qiskit
    import qiskit_aer

    devices = args.devices.split(",")
    env = {"python": platform.python_version(), "qiskit": qiskit.__version__, "qiskit_aer": qiskit_aer.__version__,
           "available_devices": list(qiskit_aer.AerSimulator().available_devices()), "platform": platform.platform()}
    warm = {}
    for dev in devices:  # warm-up: first call pays library / CUDA context initialisation
        t0 = time.perf_counter()
        w_counts_per_call(dev, [(0.3137, 0, 3)], 16, 1)
        warm[dev] = time.perf_counter() - t0
    records, checks = [], []
    for w in workloads(args.pilot):
        outs = {}
        for rep in range(args.repeat):
            for dev in devices:
                t, out = run_one(w, dev)
                outs[dev] = out
                rec = {k: w[k] for k in ("name", "n", "widths", "width", "shots", "n_circuits")}
                rec.update(device=dev, repeat=rep, **t)
                records.append(rec)
                print(json.dumps(rec), flush=True)
        if set(devices) >= {"CPU", "GPU"}:
            c = correctness(w, outs["CPU"], outs["GPU"])
            c.update(name=w["name"], widths=w["widths"], shots=w["shots"])
            checks.append(c)
            print("CHECK", json.dumps(c), flush=True)
    json.dump({"env": env, "warmup_s": warm, "pilot": args.pilot, "records": records, "checks": checks}, open(args.out, "w"), indent=1)


if __name__ == "__main__":
    main()
