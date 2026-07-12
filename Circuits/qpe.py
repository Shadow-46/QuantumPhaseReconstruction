"""
Purpose
    Build and validate standard quantum phase estimation circuits.
Theory
    QPE prepares a uniform phase register, applies controlled powers of a
    unitary to an eigenstate, and uses IQFT to convert phase into measured bits.
Inputs
    Phase precision, a unitary gate, work-register size, and eigenstate bits.
Outputs
    Qiskit circuits plus lightweight exact-phase verification utilities.
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from collections import Counter
from fractions import Fraction
from typing import Mapping

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit import Gate
from qiskit.circuit.library import PhaseGate
from qiskit.quantum_info import Statevector

from Circuits.iqft import inverse_qft_circuit


def bitstring_to_fraction(bitstring: str) -> Fraction:
    """Convert a measured QPE bitstring into an estimated phase fraction."""
    return Fraction(int(bitstring, 2), 2 ** len(bitstring))


def most_likely_bitstring(counts: Mapping[str, int]) -> str:
    """Return the most frequently observed bitstring."""
    if not counts:
        raise ValueError("counts must not be empty.")
    return max(counts.items(), key=lambda item: item[1])[0].replace(" ", "")


def build_qpe_circuit(unitary: Gate, num_phase_qubits: int, num_work_qubits: int, eigenstate: int = 1, measure: bool = True) -> QuantumCircuit:
    """Build a standard QPE circuit for a supplied unitary and eigenstate basis value."""
    if num_phase_qubits < 1 or num_work_qubits < 1:
        raise ValueError("phase and work registers must contain at least one qubit.")
    if not 0 <= eigenstate < 2**num_work_qubits:
        raise ValueError("eigenstate basis value does not fit in the work register.")
    circuit = QuantumCircuit(num_phase_qubits + num_work_qubits, num_phase_qubits if measure else 0, name="QPE")
    phase = list(range(num_phase_qubits))
    work = list(range(num_phase_qubits, num_phase_qubits + num_work_qubits))
    for qubit in phase:
        circuit.h(qubit)
    for index, bit in enumerate(format(eigenstate, f"0{num_work_qubits}b")[::-1]):
        if bit == "1":
            circuit.x(work[index])
    for control in range(num_phase_qubits):
        powered = unitary.power(2**control).control(1)
        circuit.append(powered, [phase[control], *work])
    circuit.append(inverse_qft_circuit(num_phase_qubits).to_gate(label="IQFT"), phase)
    if measure:
        circuit.measure(phase, range(num_phase_qubits))
    return circuit


def exact_phase_counts(phase: Fraction, num_phase_qubits: int) -> dict[str, int]:
    """Return exact statevector counts for a one-qubit phase eigenstate."""
    gate = PhaseGate(2 * np.pi * float(phase))
    circuit = build_qpe_circuit(gate, num_phase_qubits, 1, eigenstate=1, measure=False)
    state = Statevector.from_instruction(circuit)
    probs = state.probabilities_dict(qargs=list(range(num_phase_qubits)))
    scaled = Counter()
    for raw, probability in probs.items():
        bitstring = raw
        if probability > 1e-9:
            scaled[bitstring] += int(round(probability * 1_000_000))
    return dict(scaled)


def verify_known_phases(num_phase_qubits: int = 4) -> bool:
    """Verify exact QPE recovery for phases 1/2, 1/4, and 3/8."""
    for phase in (Fraction(1, 2), Fraction(1, 4), Fraction(3, 8)):
        counts = exact_phase_counts(phase, num_phase_qubits)
        estimate = bitstring_to_fraction(most_likely_bitstring(counts))
        if estimate != phase:
            return False
    return True


def self_test() -> None:
    """Run exact phase checks used by repository smoke tests."""
    assert verify_known_phases()


if __name__ == "__main__":
    self_test()
    print("qpe self-test passed")
