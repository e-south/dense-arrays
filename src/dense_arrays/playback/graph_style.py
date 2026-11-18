"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/playback/graph_style.py

Renderer-neutral arrow dimensions relative to graph nucleotide type.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import asdict, dataclass, fields


@dataclass(frozen=True, slots=True)
class GraphArrowDimensions:
    """Resolved physical geometry and opacity, independent of a drawing backend."""

    shaft_width_pt: float
    head_length_pt: float
    head_width_pt: float
    opacity: float


@dataclass(frozen=True, slots=True)
class GraphArrowStyle:
    """Size shafts and filled heads in multiples of the graph font's point size.

    Emphasis changes color, not geometry. Heads have no stroked outline.
    """

    shaft_width_em: float = 0.07
    head_length_em: float = 0.53
    head_width_em: float = 0.385
    opacity: float = 0.55

    def __post_init__(self) -> None:
        """Reject ambiguous units and invalid dimensions before rendering."""
        for item in fields(self):
            value = getattr(self, item.name)
            if isinstance(value, bool) or not isinstance(value, (float, int)):
                msg = f"graph_arrows.{item.name} must be a positive finite number"
                raise TypeError(msg)
            if not math.isfinite(value) or value <= 0:
                msg = f"graph_arrows.{item.name} must be a positive finite number"
                raise ValueError(msg)

        if self.opacity > 1:
            msg = "graph_arrows.opacity must be at most 1"
            raise ValueError(msg)

    def resolve(self, font_size_pt: float) -> GraphArrowDimensions:
        """Resolve physical arrow dimensions once from the graph font size."""
        if (
            isinstance(font_size_pt, bool)
            or not isinstance(font_size_pt, (float, int))
            or not math.isfinite(font_size_pt)
            or font_size_pt <= 0
        ):
            msg = "graph font size must be a positive finite number"
            raise ValueError(msg)
        return GraphArrowDimensions(
            self.shaft_width_em * font_size_pt,
            self.head_length_em * font_size_pt,
            self.head_width_em * font_size_pt,
            self.opacity,
        )

    @classmethod
    def from_mapping(cls, raw: Mapping[str, object]) -> GraphArrowStyle:
        """Admit the same fields for direct API and JSON CLI use."""
        if not isinstance(raw, Mapping):
            msg = "graph_arrows must be a mapping"
            raise TypeError(msg)
        unknown = set(raw) - {item.name for item in fields(cls)}
        if unknown:
            msg = f"unknown graph_arrows fields: {sorted(unknown)}"
            raise ValueError(msg)
        return cls(**raw)

    def to_mapping(self) -> dict[str, float]:
        """Return detached, JSON-compatible settings."""
        return asdict(self)
