# Qiskit Aer GPU: installation attempt and CPU-vs-GPU benchmark (2026-09-28)

## Environment inspected (existing WSL `.venv`)

| Item | Value |
|---|---|
| Distribution | WSL2 Ubuntu 24.04.4, kernel 6.6.87.2-microsoft-standard-WSL2 |
| `.venv` location | `QuantumPhaseReconstruction/.venv` |
| Python | 3.12.3 |
| Qiskit stack | qiskit 2.5.0, qiskit-aer 0.17.2 (CPU build), numpy 2.5.1, scipy 1.18.0 |
| GPU | NVIDIA GeForce RTX 5050 Laptop, 8 GB (Blackwell, sm_120) |
| Driver | 595.97, CUDA 13.2 capable. Visible in WSL via `/usr/lib/wsl/lib/libcuda.so.1` |
| Toolchain | g++ 13.3, make, git and Python headers present. No system CUDA toolkit (nvcc), no cmake. |
| Privileges | sudo requires a password, so system packages cannot be installed unattended. |

A `pip freeze` of the `.venv` was taken before any change.

## What was tried, and the outcome

1. **`qiskit-aer-gpu` from PyPI: incompatible with Qiskit 2.5.0.**
   - PyPI's newest GPU wheel is 0.15.1, and no 0.16 or 0.17 GPU wheel exists.
   - It was installed only in a throwaway venv, never in `.venv`, together with the CUDA-12 runtime wheels (see the table below).
   - It fails at import: `ImportError: cannot import name 'convert_to_target' from 'qiskit.providers'`. Aer 0.15 predates Qiskit 2.x, which removed that function.
   - Using it would mean downgrading Qiskit, so it was rejected.
2. **Building qiskit-aer 0.17.2 (the installed version) from its PyPI sdist with the CUDA Thrust backend: blocked by the toolchain.**
   - The build was targeted at sm_120, in scratch under `~/aer_gpu_scratch`.
   - The pip wheel `nvidia-cuda-nvcc-cu12` 12.9.86 contains only `ptxas`, `libnvvm` and headers. It has **no `nvcc` driver** and no `cicc` or `fatbinary`, so CUDA C++ cannot be compiled from pip packages alone.
   - A build needs a full CUDA toolkit ≥ 12.8, the first version with sm_120 support. It can come from:
     - (a) apt `cuda-toolkit-12-9`, which needs sudo;
     - (b) NVIDIA's runfile installed in user space (`--toolkit --toolkitpath=$HOME/cuda-12.9`, about 4 GB download);
     - (c) conda-forge `cuda-nvcc`/`cuda-toolkit` in a user-space micromamba env.
   - None of these is a package in the existing `.venv`, so this goes beyond the installation that was authorised. **Stopped for a decision.**

**Net result.** **Nothing was installed or changed in the existing `.venv`.** It still has qiskit 2.5.0 and qiskit-aer 0.17.2 (CPU). GPU verification (import, `available_devices()` containing "GPU", a minimal GPU run) could not be completed.

**Scratch artefacts.** These are outside the repository and outside `.venv`; remove them with `rm -rf ~/aer_gpu_scratch`, about 2.5 GB. Packages installed there:

| Package | Version |
|---|---|
| qiskit-aer-gpu | 0.15.1 |
| nvidia-cuda-runtime-cu12 | 12.9.79 |
| nvidia-cuda-nvrtc-cu12 | 12.9.86 |
| nvidia-cublas-cu12 | 12.9.2.10 |
| nvidia-cusolver-cu12 | 11.7.5.82 |
| nvidia-cusparse-cu12 | 12.5.10.65 |
| nvidia-nvjitlink-cu12 | 12.9.86 |
| custatevec-cu12 | 1.15.0 |
| nvidia-cuda-nvcc-cu12 | 12.9.86 |
| nvidia-cuda-cccl-cu12 | (latest) |
| cmake | <4 |
| ninja, pybind11, scikit-build | (latest) |

**Operational note.** WSL's `/tmp` is tmpfs and is wiped when the WSL VM idles out. Scratch work must live under `$HOME`.

## Benchmark tool and CPU pilot (the GPU arm is pending)

**Tool.** `research/awqpe/bench_aer_gpu.py` runs the actual P10a/P10b circuits, as those phases build and call them:
- `phase_gate_block`, with one run() per block;
- the P10b batched statevector jitter run;
- the P10b repeated-query density-matrix run with labelled two-qubit depolarising noise.

**What it measures.** Identical circuits, shots and seeds on `device="CPU"` and `device="GPU"`. It times construction, transpilation, execution and total separately, and records RSS and whole-device GPU memory (sampled with nvidia-smi). It also runs a CPU-vs-GPU correctness check: exact probabilities for the statevector and density-matrix workloads, and chi-square against the kernel plus CPU/GPU total-variation distance for counts. It changes no scientific result.

**CPU-only pilot.** Run in the existing `.venv` (Aer 0.17.2 CPU), 2 repeats, after a warm-up of 1.13 s. Output: `research/results/bench_aer/bench_cpu_pilot_20260928.json`. Mean seconds per workload:

| Workload | n / widths | Shots | Circuits | Build | Transpile | Execute | Total | Execute share | Max total speedup if execution cost 0 |
|---|---|---|---|---|---|---|---|---|---|
| p10a_e2e | 8, [3,2,3] | 16 | 12 | 0.006 | 0.570 | 0.028 | 0.605 | 5% | 1.05x |
| p10a_e2e | 8, [4,4] | 16 | 8 | 0.005 | 0.391 | 0.029 | 0.425 | 7% | 1.07x |
| p10a_e2e | 16, [4,4,4,4] | 16 | 16 | 0.008 | 0.770 | 0.058 | 0.837 | 7% | 1.07x |
| p10a_dist | widths 2-5 | 4096 | 24 | 0.014 | 1.151 | 0.276 | 1.442 | 19% | 1.24x |
| p10b_jitter | 4 | exact | 200 | 0.050 | 1.483 | 0.049 | 1.583 | 3% | 1.03x |
| p10b_depol | 2-4 (density matrix) | exact | 12 | 0.009 | 0.139 | 0.050 | 0.198 | 25% | 1.34x |

**What the CPU pilot already establishes.**
- For this project's actual Aer workload, simulation is a small fraction of wall time: 3-25%. Transpilation dominates.
- Even an infinitely fast GPU could reduce total time by at most 1.03-1.34x (Amdahl's bound).
- GPU execution of 3-6-qubit circuits also adds per-call launch and transfer overhead, so a net slowdown is plausible.
- The GPU timings themselves remain unmeasured until a GPU-enabled Aer exists.
- **All scientific phases' Monte Carlo work used the NumPy Dirichlet-kernel simulator, not Aer.** Aer GPU therefore cannot speed up P2-P9. It could only affect the P10 validation runs.
