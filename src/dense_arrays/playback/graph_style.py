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

from .colors import validate_color


@dataclass(frozen=True, slots=True)
class GraphArrowDimensions:
    """Resolved physical geometry and opacity, independent of a drawing backend."""

    shaft_width_pt: float
    head_length_pt: float
    head_width_pt: float
    opacity: float
    head_shape: str


@dataclass(frozen=True, slots=True)
class GraphArrowStyle:
    """Size shafts and filled heads in multiples of the graph font's point size.

    Emphasis changes color, not geometry. Heads have no stroked outline.
    """

    shaft_width_em: float = 0.07
    head_length_em: float = 0.43
    head_width_em: float = 0.5
    opacity: float = 1.0
    head_shape: str = "triangle"
    resting_color: str = "#9A9A9A"
    completed_color: str = "#505050"
    active_color: str = "#202020"

    def __post_init__(self) -> None:
        """Reject ambiguous units and invalid dimensions before rendering."""
        for name in ("shaft_width_em", "head_length_em", "head_width_em", "opacity"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (float, int)):
                msg = f"graph_arrows.{name} must be a positive finite number"
                raise TypeError(msg)
            if (
                not math.isfinite(value)
                or value < 0
                or (name != "opacity" and value == 0)
            ):
                msg = f"graph_arrows.{name} must be a positive finite number"
                raise ValueError(msg)

        if self.head_shape not in ("triangle", "notched"):
            msg = "graph_arrows.head_shape must be triangle or notched"
            raise ValueError(msg)
        for name in ("resting_color", "completed_color", "active_color"):
            try:
                validate_color(getattr(self, name))
            except ValueError as error:
                msg = f"graph_arrows.{name}: {error}"
                raise ValueError(msg) from error
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
            self.head_shape,
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

    def to_mapping(self) -> dict[str, object]:
        """Return detached, JSON-compatible settings."""
        return asdict(self)
