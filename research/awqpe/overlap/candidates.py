"""Overlap candidate blocks for Phase 6 (decision side: geometry only, no truth).

For every internal boundary j (between chunk j and chunk j+1, 0-indexed):

  O-ext    BlockSpec(k_j, m_j + v, role="ext"): chunk j re-measured with v extra
           controls toward the LSB, so its lowest v bits overlap chunk j+1's top
           v bits (v capped at m_{j+1}).
  O-bridge BlockSpec(K_j - s, w, role="bridge"): a width-w block whose top s bits
           are chunk j's LSBs and whose remaining w - s bits are chunk j+1's MSBs
           (K_j = boundary position). w defaults to m_j; s is 1 or ceil(w/2).

The two mechanisms are never mixed within one experimental arm.
"""

from __future__ import annotations

from math import ceil

from research.awqpe.blocks.geometry import BlockSpec, bridge_block, extension_block, partition_blocks


def overlap_candidates(widths, mechanism: str, v: int = 1, shift: str = "half", bridge_width: int | None = None) -> list[tuple[int, BlockSpec, int]]:
    """List of (boundary j, spec, v_effective) for every internal boundary."""
    chunks = partition_blocks(widths)
    out = []
    for j in range(len(chunks) - 1):
        up, lo = chunks[j], chunks[j + 1]
        if mechanism == "ext":
            v_eff = min(int(v), lo.width)
            out.append((j, extension_block(up, v_eff), v_eff))
        elif mechanism == "bridge":
            w = int(bridge_width or up.width)
            s = 1 if shift == "1" else min(int(ceil(w / 2)), up.width)
            s = max(1, min(s, w - 1))
            out.append((j, bridge_block(up, lo, w, s), 0))
        else:
            raise ValueError(f"unknown overlap mechanism {mechanism!r}")
    return out
