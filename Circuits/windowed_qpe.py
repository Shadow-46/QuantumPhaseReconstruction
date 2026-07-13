"""
Purpose
    Build block-local windowed phase-estimation circuits and exact window
    samples.
Theory
    Windowed QPE (arXiv:2507.22460, arXiv:2509.05010) partitions the phase
    register into small independent blocks. Each block allocates only
    `width` phase qubits and applies controlled unitaries U^(2^(offset+j))
    for j = 0..width-1, where offset is the block's exponent offset (the
    window's start position). This differs from measuring a subset of a
    single full-precision QPE register: the whole point of windowing is that
    no single circuit ever materializes more than `width` phase qubits, so
    qubit count and depth are decoupled from the total precision n.
Inputs
    A ModularMultiplicationOperator, work-register eigenstate, and window
    bounds.
Outputs
    Qiskit circuits, Aer-sampled counts, and deterministic per-window
    verification counts.
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from dataclasses import dataclass
from fractions import Fraction

from qiskit import QuantumCircuit

from Circuits.iqft import inverse_qft_circuit
from Circuits.modular_multiplication import ModularMultiplicationOperator
from Reconstruction.candidate_generation import exact_window_counts


@dataclass(frozen=True)
class WindowSpec:
    """A contiguous measured phase-register window."""
    start: int
    width: int
    total_precision: int

    @property
    def end(self) -> int:
        """Return the exclusive end index of this window."""
        return self.start + self.width

    def validate(self) -> None:
        """Validate the window geometry."""
        if self.start < 0 or self.width < 1 or self.end > self.total_precision:
            raise ValueError("invalid window geometry.")


def make_overlapping_windows(total_precision: int, window_size: int, overlap: int) -> list[WindowSpec]:
    """Create a left-to-right covering of overlapping phase-bit windows."""
    if total_precision < 1 or window_size < 1:
        raise ValueError("precision and window size must be positive.")
    if overlap < 0 or overlap >= window_size:
        raise ValueError("overlap must satisfy 0 <= overlap < window_size.")
    step = window_size - overlap
    starts = list(range(0, total_precision, step))
    specs: list[WindowSpec] = []
    for start in starts:
        width = min(window_size, total_precision - start)
        if width > 0:
            specs.append(WindowSpec(start, width, total_precision))
        if start + width == total_precision:
            break
    return specs


def build_windowed_qpe_circuit(op: ModularMultiplicationOperator, eigenstate: int, spec: WindowSpec) -> QuantumCircuit:
    """Build a shallow, block-local QPE circuit for a single window.

    Only `spec.width` phase qubits and `op.qubits` work qubits are ever
    allocated: physical phase qubit j (0-indexed) is controlled by
    U^(2^(spec.start + j)), so the block's own IQFT and measurement directly
    yield the best `spec.width`-bit approximation of bits
    [spec.start, spec.end) of the full-precision phase, matching Algorithm 2
    of arXiv:2509.05010.
    """
    spec.validate()
    if not 0 <= eigenstate < op.dimension:
        raise ValueError("eigenstate basis value does not fit in the work register.")
    work_qubits = op.qubits
    circuit = QuantumCircuit(spec.width + work_qubits, spec.width, name=f"WQPE_{spec.start}_{spec.width}")
    phase = list(range(spec.width))
    work = list(range(spec.width, spec.width + work_qubits))
    for qubit in phase:
        circuit.h(qubit)
    for index, bit in enumerate(format(eigenstate, f"0{work_qubits}b")[::-1]):
        if bit == "1":
            circuit.x(work[index])
    for local_index in range(spec.width):
        gate = op.controlled_power_gate(2 ** (spec.start + local_index))
        circuit.append(gate, [phase[local_index], *work])
    circuit.append(inverse_qft_circuit(spec.width).to_gate(label="IQFT"), phase)
    circuit.measure(phase, range(spec.width))
    return circuit


def run_windowed_qpe_block(op: ModularMultiplicationOperator, eigenstate: int, spec: WindowSpec, shots: int, method: str = "automatic", seed: int | None = None, noise_model=None) -> dict[str, int]:
    """Build one block-local windowed QPE circuit and sample it on Aer."""
    from Simulation.backend import run_counts

    circuit = build_windowed_qpe_circuit(op, eigenstate, spec)
    return run_counts(circuit, shots=shots, method=method, seed=seed, noise_model=noise_model)


def deterministic_window_samples(phase: Fraction, specs: list[WindowSpec]) -> list[dict[str, int]]:
    """Return exact local counts for a phase exactly representable at the requested precision."""
    if not specs:
        raise ValueError("at least one window is required.")
    precision = specs[0].total_precision
    phase_value = int(phase * (2**precision))
    return [exact_window_counts(phase_value, precision, spec.start, spec.width) for spec in specs]


def _phase_register_probabilities(circuit: QuantumCircuit, num_phase_qubits: int) -> dict[str, float]:
    """Return exact phase-register marginal probabilities via statevector."""
    from qiskit.quantum_info import Statevector

    unmeasured = circuit.copy()
    unmeasured.remove_final_measurements(inplace=True)
    state = Statevector.from_instruction(unmeasured)
    return state.probabilities_dict(qargs=list(range(num_phase_qubits)))


def self_test() -> None:
    """Verify every deterministic window of the phase 1/4 at six-bit precision
    and cross-check the block-local circuit's exact distribution against a
    full-register standard QPE run with the same exponent offset."""
    specs = make_overlapping_windows(6, 3, 1)
    samples = deterministic_window_samples(Fraction(1, 4), specs)
    assert samples == [{"010": 1}, {"000": 1}, {"00": 1}]

    # N=15, a=2 has order 4; the standard Shor initialization |1> is an equal
    # superposition of the four eigenstates with phases k/4, k=0..3. A block
    # circuit with exponent offset 0 must therefore reproduce exactly the
    # same distribution as an ordinary `width`-qubit standard QPE run (which
    # allocates the same number of phase qubits with exponents 2^0..2^(w-1)).
    from Circuits.qpe import build_qpe_circuit

    op = ModularMultiplicationOperator(2, 15)
    zero_offset_spec = WindowSpec(start=0, width=2, total_precision=4)
    windowed_circuit = build_windowed_qpe_circuit(op, eigenstate=1, spec=zero_offset_spec)
    assert windowed_circuit.num_qubits == zero_offset_spec.width + op.qubits
    standard_circuit = build_qpe_circuit(op.gate(), zero_offset_spec.width, op.qubits, eigenstate=1, measure=False)
    windowed_probs = _phase_register_probabilities(windowed_circuit, zero_offset_spec.width)
    standard_probs = _phase_register_probabilities(standard_circuit, zero_offset_spec.width)
    assert windowed_probs == standard_probs

    # A block circuit further from the MSB (exponent offset > 0) must stay
    # exactly `width` phase qubits wide, never growing with total_precision:
    # this is the qubit-count reduction the windowed formulation promises.
    far_spec = WindowSpec(start=6, width=2, total_precision=10)
    far_circuit = build_windowed_qpe_circuit(op, eigenstate=1, spec=far_spec)
    assert far_circuit.num_qubits == far_spec.width + op.qubits


if __name__ == "__main__":
    self_test()
    print("windowed_qpe self-test passed")
