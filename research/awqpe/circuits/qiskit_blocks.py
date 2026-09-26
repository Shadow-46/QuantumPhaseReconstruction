"""Qiskit circuits for single AWQPE blocks (Algorithm 1, lines 5-12).

A block with offset k and width m allocates m control qubits plus the target
register, prepares the target eigenstate, applies H to every control, applies
controlled U^(2^(k+p)) from control p (p = 0..m-1), applies the inverse QFT to
the controls and measures them. Control p is classical bit p, so the integer
value of a Qiskit counts key equals the kernel outcome y.

Unitaries and eigenstates are supplied by research.awqpe.circuits.unitaries,
which verifies U|u> = exp(2 pi i phi)|u> before a circuit is built.
"""

from __future__ import annotations

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import CPhaseGate, QFTGate


def phase_gate_block(phi: float, offset: int, width: int, measure: bool = True) -> QuantumCircuit:
    """Block circuit for U = P(2 pi phi) on one target prepared in |1>.

    |1> is an exact eigenstate of the phase gate with eigenphase phi; this is
    the construction drawn in the paper's Fig. 1 (X on the target qubit).
    Controlled powers are applied as a single CPhase(2 pi frac(2^(k+p) phi)),
    i.e. U^(2^(k+p)) synthesised exactly rather than by repetition; circuit
    depth for noise studies is modelled separately (see model/noise.py).
    """
    controls = list(range(width))
    target = width
    qc = QuantumCircuit(width + 1, width if measure else 0, name=f"AWQPE_k{offset}_m{width}")
    qc.x(target)
    qc.h(controls)
    for p in controls:
        angle = 2.0 * np.pi * float(np.mod(np.ldexp(phi, offset + p), 1.0))
        qc.append(CPhaseGate(angle), [p, target])
    qc.append(QFTGate(width).inverse(), controls)
    if measure:
        qc.measure(controls, controls)
    return qc


def general_unitary_block(U: np.ndarray, psi: np.ndarray, offset: int, width: int, measure: bool = True) -> QuantumCircuit:
    """Block circuit for an arbitrary unitary U on target register prepared in psi.

    Controlled powers U^(2^(k+p)) are exact matrix powers wrapped as
    UnitaryGate(...).control(1); psi is loaded with StatePreparation. Used to
    check at circuit level that non-eigenstate inputs produce the eigenphase
    mixture predicted by circuits/unitaries.spectral_decomposition.
    """
    from qiskit.circuit.library import StatePreparation, UnitaryGate

    nt = int(np.log2(U.shape[0]))
    controls = list(range(width))
    targets = list(range(width, width + nt))
    qc = QuantumCircuit(width + nt, width if measure else 0)
    qc.append(StatePreparation(psi / np.linalg.norm(psi)), targets)
    qc.h(controls)
    dim = U.shape[0]
    P0, P1 = np.diag([1.0, 0.0]), np.diag([0.0, 1.0])
    for p in controls:
        Up = np.linalg.matrix_power(U, 1 << (offset + p))
        # Controlled-U^(2^(k+p)) as an explicit matrix on qargs [control, *targets]
        # (Qiskit little-endian: the control is the least significant qarg). This
        # avoids gate synthesis of large controlled unitaries; it is exact.
        cu = np.kron(np.eye(dim), P0) + np.kron(Up, P1)
        qc.append(UnitaryGate(cu, check_input=False), [p, *targets])
    qc.append(QFTGate(width).inverse(), controls)
    if measure:
        qc.measure(controls, controls)
    return qc


def block_statevector_probabilities(qc: QuantumCircuit, width: int) -> np.ndarray:
    """Exact control-register outcome distribution p[y] of an unmeasured block."""
    from qiskit.quantum_info import Statevector

    probs = Statevector.from_instruction(qc).probabilities(qargs=list(range(width)))
    return np.asarray(probs, dtype=float)
