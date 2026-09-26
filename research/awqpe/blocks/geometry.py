"""Block geometry and quantum-cost accounting for AWQPE.

AWQPE partitions n phase bits into non-overlapping chunks [m_1, ..., m_B]
with m_i > 1 (Algorithm 1 input). Block i has offset k_i = sum_{j<i} m_j and
applies U^(2^(k_i+p)), p = 0..m_i-1, so it measures the m_i most significant
bits of frac(2^(k_i) phi) (bits k_i+1 .. k_i+m_i of phi, MSB first).

Overlap is NOT part of AWQPE. Extra blocks created by overlap policies are
represented by the same BlockSpec (any offset/width) and tagged with a role so
that the faithful baseline never sees them.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, order=True)
class BlockSpec:
    """One independently executed block circuit."""

    offset: int
    width: int
    role: str = "chunk"  # "chunk" (AWQPE partition) | "ext" | "bridge" (overlap extensions)

    def __post_init__(self) -> None:
        if self.offset < 0 or self.width < 1:
            raise ValueError("invalid block geometry.")
        if self.role not in {"chunk", "ext", "bridge"}:
            raise ValueError(f"unknown block role {self.role!r}.")

    @property
    def M(self) -> int:
        return 1 << self.width

    @property
    def u_queries_per_shot(self) -> int:
        """Controlled-U applications per execution: 2^k (2^m - 1)."""
        return (1 << self.offset) * ((1 << self.width) - 1)

    @property
    def max_power(self) -> int:
        """Largest controlled power 2^(k+m-1) in the circuit (depth driver)."""
        return 1 << (self.offset + self.width - 1)

    @property
    def key(self) -> tuple[int, int]:
        """Integer key used for seed derivation (geometry only, role-independent)."""
        return (self.offset, self.width)


def partition_blocks(widths, allow_width_one: bool = False) -> list[BlockSpec]:
    """AWQPE partition [m_1..m_B] -> chunk BlockSpecs, MSB block first."""
    widths = [int(m) for m in widths]
    if not widths:
        raise ValueError("at least one block is required.")
    if not allow_width_one and any(m < 2 for m in widths):
        raise ValueError("AWQPE requires every block width m_i > 1.")
    specs, k = [], 0
    for m in widths:
        specs.append(BlockSpec(k, m))
        k += m
    return specs


def extension_block(chunk: BlockSpec, extra: int) -> BlockSpec:
    """Overlap form O-ext: same offset, width m+v (overlaps the next chunk's top v bits)."""
    if extra < 1:
        raise ValueError("extension must add at least one bit.")
    return BlockSpec(chunk.offset, chunk.width + extra, role="ext")


def bridge_block(upper: BlockSpec, lower: BlockSpec, width: int, shift: int) -> BlockSpec:
    """Overlap form O-bridge: a block straddling the boundary between two chunks.

    Its offset is (boundary - shift), so its top `shift` bits overlap the upper
    chunk's least significant bits and the rest overlap the lower chunk.
    """
    boundary = upper.offset + upper.width
    if lower.offset != boundary:
        raise ValueError("bridge requires adjacent chunks.")
    if not 1 <= shift < width:
        raise ValueError("bridge shift must satisfy 1 <= shift < width.")
    if shift > upper.width:
        raise ValueError("bridge cannot reach above the upper chunk.")
    return BlockSpec(boundary - shift, width, role="bridge")
