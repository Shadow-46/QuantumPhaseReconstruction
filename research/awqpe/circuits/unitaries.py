"""Unitaries, verified eigenstates, and spectral decomposition of input states.

A state is only called an eigenstate after ||U|u> - e^{2 pi i phi}|u>|| is
computed and reported. Any input state |psi> is decomposed as
sum_j c_j |u_j>; AWQPE blocks then see the mixture sum_j |c_j|^2 p(y | phi_j)
because every shot collapses onto one eigencomponent independently.

Cases used by Phase 3:
  phase_gate        U = diag(1, e^{2 pi i phi}); |1> is an exact eigenstate.
  diagonal          U = diag(e^{2 pi i phi_j}) on q qubits; each |j> is exact.
  conjugated        U = V D V^dagger with Haar-random V; eigenvectors V|j>.
  modmul            U_a|y> = |a y mod N> (the Paper A unitary, imported from
                    the frozen Circuits/modular_multiplication.py);
                    |u_s> = r^-1/2 sum_t e^{-2 pi i s t / r} |a^t mod N>.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from math import gcd
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parents[3]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))


@dataclass(frozen=True)
class SpectralInput:
    """What a block sees: eigenphases with weights, plus provenance."""

    label: str
    unitary: np.ndarray
    state: np.ndarray
    phases: np.ndarray  # eigenphases in [0,1) with non-negligible weight
    weights: np.ndarray  # |<u_j|psi>|^2, summing to 1
    eigen_residual: float  # ||U psi - e^{2 pi i phi} psi|| if psi is claimed eigenstate, else NaN
    is_eigenstate: bool


def eigen_residual(U: np.ndarray, v: np.ndarray, phi: float) -> float:
    return float(np.linalg.norm(U @ v - np.exp(2j * np.pi * phi) * v))


def spectral_decomposition(U: np.ndarray, psi: np.ndarray, tol: float = 1e-10) -> tuple[np.ndarray, np.ndarray]:
    """Eigenphases of U carrying weight in psi, merged when degenerate."""
    from scipy.linalg import schur

    psi = psi / np.linalg.norm(psi)
    # U is unitary (normal): its complex Schur form is diagonal with orthonormal
    # Schur vectors, which stay valid eigenvectors even for degenerate eigenvalues.
    T, Z = schur(U.astype(complex), output="complex")
    phases = np.mod(np.angle(np.diag(T)) / (2 * np.pi), 1.0)
    amps = np.abs(Z.conj().T @ psi) ** 2
    keys = np.mod(np.round(phases, 10), 1.0)
    uniq = np.unique(keys)
    w = np.array([amps[keys == u].sum() for u in uniq])
    keep = w > tol
    return uniq[keep], w[keep] / w[keep].sum()


def _make(label: str, U: np.ndarray, psi: np.ndarray, claimed_phase: float | None) -> SpectralInput:
    phases, weights = spectral_decomposition(U, psi)
    res = eigen_residual(U, psi / np.linalg.norm(psi), claimed_phase) if claimed_phase is not None else float("nan")
    return SpectralInput(label, U, psi, phases, weights, res, bool(claimed_phase is not None and res < 1e-9))


def phase_gate_input(phi: float) -> SpectralInput:
    U = np.diag([1.0, np.exp(2j * np.pi * phi)])
    return _make("phase_gate|1>", U, np.array([0, 1], dtype=complex), phi)


def diagonal_input(phases, index: int) -> SpectralInput:
    phases = np.asarray(phases, dtype=float)
    U = np.diag(np.exp(2j * np.pi * phases))
    psi = np.zeros(len(phases), dtype=complex)
    psi[index] = 1.0
    return _make(f"diagonal|{index}>", U, psi, float(phases[index]))


def haar_unitary(dim: int, rng: np.random.Generator) -> np.ndarray:
    z = (rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))) / np.sqrt(2)
    q, r = np.linalg.qr(z)
    return q * (np.diag(r) / np.abs(np.diag(r)))


def conjugated_input(phases, index: int, rng: np.random.Generator) -> SpectralInput:
    phases = np.asarray(phases, dtype=float)
    V = haar_unitary(len(phases), rng)
    U = V @ np.diag(np.exp(2j * np.pi * phases)) @ V.conj().T
    return _make(f"VDV+|v{index}>", U, V[:, index].copy(), float(phases[index]))


def modmul_matrix(a: int, N: int) -> np.ndarray:
    from Circuits.modular_multiplication import ModularMultiplicationOperator  # frozen Paper A module, read-only

    return ModularMultiplicationOperator(a, N).matrix()


def multiplicative_order(a: int, N: int) -> int:
    if gcd(a, N) != 1:
        raise ValueError("a must be coprime to N.")
    r, x = 1, a % N
    while x != 1:
        x = (x * a) % N
        r += 1
    return r


def modmul_eigenstate_input(a: int, N: int, s: int) -> SpectralInput:
    U = modmul_matrix(a, N)
    r = multiplicative_order(a, N)
    psi = np.zeros(U.shape[0], dtype=complex)
    for t in range(r):
        psi[pow(a, t, N)] += np.exp(-2j * np.pi * s * t / r) / np.sqrt(r)
    return _make(f"U_{a}mod{N}|u_{s}>", U, psi, s / r)


def modmul_one_input(a: int, N: int) -> SpectralInput:
    """The Paper A input: computational |1> on the work register (NOT an eigenstate)."""
    U = modmul_matrix(a, N)
    psi = np.zeros(U.shape[0], dtype=complex)
    psi[1] = 1.0
    return _make(f"U_{a}mod{N}|1>", U, psi, None)


def controlled_mixture_input(phi0: float, w0: float, others) -> SpectralInput:
    """Diagonal U with target component phi0 of weight w0 and the rest spread equally over `others`."""
    others = list(others)
    phases = np.array([phi0] + others)
    amps = np.sqrt(np.array([w0] + [(1 - w0) / len(others)] * len(others))) if others else np.array([1.0])
    # Pad to a qubit register (power-of-two dimension) with zero-amplitude eigenvalues.
    dim = 1 << max(1, int(np.ceil(np.log2(len(phases)))))
    phases = np.concatenate([phases, np.zeros(dim - len(phases))])
    amps = np.concatenate([amps, np.zeros(dim - len(amps))])
    U = np.diag(np.exp(2j * np.pi * phases))
    return _make(f"mixture(w0={w0})", U, amps.astype(complex), phi0 if w0 == 1.0 else None)
