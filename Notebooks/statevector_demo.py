"""
Purpose
    Provide a scriptable statevector demonstration for exact QPE phases.
Theory
    Notebook-adjacent demos should be executable Python so they participate in
    repository validation while showing phase peaks without running notebooks.
Inputs
    Built-in exact phases 1/2, 1/4, and 3/8.
Outputs
    Printed exact QPE count dictionaries.
Author
    Sanjay
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from fractions import Fraction

from Circuits.qpe import exact_phase_counts


def main() -> int:
    """Print exact statevector-derived phase-estimation peaks."""
    for phase in (Fraction(1, 2), Fraction(1, 4), Fraction(3, 8)):
        print(phase, exact_phase_counts(phase, 4))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
