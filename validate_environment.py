"""
Purpose
    Validate the Python environment and run repository self-tests.
Theory
    Research runs should fail early when a virtual environment contains a
    corrupted binary package or a dependency version outside the supported
    range. This script checks imports, versions, and executable module tests.
Inputs
    The active Python interpreter and installed packages in the current venv.
Outputs
    A validation report and a nonzero exit code on dependency or self-test
    failure.
Author
    Sanjay
"""

from __future__ import annotations

from dataclasses import dataclass
from importlib import import_module
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


@dataclass(frozen=True)
class Requirement:
    """A dependency with an inclusive lower bound and optional exclusive upper bound."""

    distribution: str
    module: str
    minimum: tuple[int, ...]
    maximum_exclusive: tuple[int, ...] | None = None


REQUIREMENTS = (
    Requirement("qiskit", "qiskit", (2, 5)),
    Requirement("qiskit-aer", "qiskit_aer", (0, 15)),
    Requirement("numpy", "numpy", (2, 0), (3, 0)),
    Requirement("scipy", "scipy", (1, 14)),
    Requirement("matplotlib", "matplotlib", (3, 9)),
    Requirement("networkx", "networkx", (3, 3)),
    Requirement("pandas", "pandas", (2, 2), (3, 0)),
)


SELF_TEST_MODULES = (
    "Circuits.modular_multiplication",
    "Circuits.iqft",
    "Circuits.qpe",
    "Circuits.windowed_qpe",
    "Reconstruction.candidate_generation",
    "Reconstruction.carry",
    "Reconstruction.stitching",
    "Reconstruction.continued_fraction",
    "Algorithms.standard_shor",
    "Algorithms.paper_algorithm",
    "Algorithms.adaptive_algorithm",
    "Evaluation.metrics",
    "Evaluation.experiments",
    "Evaluation.plots",
    "Simulation.backend",
    "Simulation.noise",
)


def parse_version(value: str) -> tuple[int, ...]:
    """Parse the numeric prefix of a Python package version string."""

    parts: list[int] = []
    for token in value.replace("-", ".").split("."):
        digits = "".join(ch for ch in token if ch.isdigit())
        if not digits:
            break
        parts.append(int(digits))
    return tuple(parts)


def version_satisfies(installed: tuple[int, ...], requirement: Requirement) -> bool:
    """Return True when an installed version satisfies a Requirement."""

    lower_ok = installed >= requirement.minimum
    upper_ok = requirement.maximum_exclusive is None or installed < requirement.maximum_exclusive
    return lower_ok and upper_ok


def validate_dependencies() -> None:
    """Import and version-check every required dependency."""

    print(f"python {sys.version.split()[0]} at {sys.executable}")
    for requirement in REQUIREMENTS:
        try:
            installed_text = version(requirement.distribution)
        except PackageNotFoundError as exc:
            raise RuntimeError(f"missing dependency: {requirement.distribution}") from exc
        installed = parse_version(installed_text)
        if not version_satisfies(installed, requirement):
            upper = f", <{'.'.join(map(str, requirement.maximum_exclusive))}" if requirement.maximum_exclusive else ""
            minimum = ".".join(map(str, requirement.minimum))
            raise RuntimeError(
                f"{requirement.distribution} {installed_text} is outside supported range >={minimum}{upper}"
            )
        import_module(requirement.module)
        print(f"{requirement.distribution} {installed_text}")


def run_self_tests() -> None:
    """Run all repository module self-tests."""

    for module_name in SELF_TEST_MODULES:
        module = import_module(module_name)
        module.self_test()
        print(f"{module_name}: ok")


def main() -> int:
    """Validate dependencies and execute the project self-test suite."""

    validate_dependencies()
    run_self_tests()
    print("environment validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
