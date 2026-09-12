"""Generic presentation choices and explicit color validation.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .models import PlaybackStep

DEFAULT_STEP_COLORS = (
    "#67BFA5",
    "#D883A4",
    "#7BA4D9",
    "#C08A56",
    "#5DA79F",
    "#D1B06C",
    "#74C0CB",
    "#86A5D8",
)
_MAX_GRAPH_FRACTION = 0.5
_PROFILES = {"categorical", "uniform", "constraints"}


def validate_color(color: str) -> None:
    """Require a portable opaque RGB color for raster output."""
    if not isinstance(color, str) or re.fullmatch(r"#[0-9a-fA-F]{6}", color) is None:
        msg = "colors must use #RRGGBB notation"
        raise ValueError(msg)


@dataclass(frozen=True, slots=True)
class LegendEntry:
    """One caller-authored legend item."""

    key: str
    label: str
    color: str

    def __post_init__(self) -> None:
        """Validate portable presentation values without coercion."""
        if not isinstance(self.key, str) or not self.key.strip():
            msg = "legend key must be nonempty text"
            raise ValueError(msg)
        if not isinstance(self.label, str) or not self.label.strip():
            msg = "legend label must be nonempty text"
            raise ValueError(msg)
        validate_color(self.color)


@dataclass(frozen=True, slots=True)
class PlaybackPresentation:
    """Presentation controls for raster scenes.

    Required authority, order, and failure evidence is always visible.
    ``show_authority_notice`` additionally displays the plan's detailed notices.
    ``full`` shows context relations; ``reduced`` shows traversal relations only.
    ``none`` hides the graph. Distance brackets refer only to declared constraints;
    ``always`` also states when no constraint was declared.
    """

    color_profile: str = "categorical"
    legend_entries: tuple[LegendEntry, ...] = ()
    graph_detail: str = "full"
    graph_fraction: float = 0.35
    show_edge_costs: bool = True
    show_authority_notice: bool = False
    show_distance_bracket: str = "when_declared"

    def __post_init__(self) -> None:
        """Validate portable presentation values without coercion."""
        if self.color_profile not in _PROFILES:
            msg = (
                "color_profile must be categorical, uniform, or constraints; "
                "use explicit color_overrides for caller palettes"
            )
            raise ValueError(msg)
        if self.graph_detail not in {"full", "reduced", "none"}:
            msg = "graph_detail must be full, reduced, or none"
            raise ValueError(msg)
        if isinstance(self.graph_fraction, bool) or not isinstance(
            self.graph_fraction, (int, float)
        ):
            msg = "graph_fraction must be numeric"
            raise TypeError(msg)
        if (
            not math.isfinite(self.graph_fraction)
            or not 0 <= self.graph_fraction <= _MAX_GRAPH_FRACTION
        ):
            msg = "graph_fraction must be finite and between 0.0 and 0.5"
            raise ValueError(msg)
        if (self.graph_detail == "none") != (self.graph_fraction == 0):
            msg = "graph_fraction must be zero exactly when graph_detail is none"
            raise ValueError(msg)
        for name in ("show_edge_costs", "show_authority_notice"):
            if not isinstance(getattr(self, name), bool):
                msg = f"{name} must be a boolean"
                raise TypeError(msg)
        if self.show_distance_bracket not in {"never", "when_declared", "always"}:
            msg = "show_distance_bracket must be never, when_declared, or always"
            raise ValueError(msg)
        object.__setattr__(self, "legend_entries", _legend_entries(self.legend_entries))


def _legend_entries(values: tuple[LegendEntry, ...]) -> tuple[LegendEntry, ...]:
    entries = tuple(values)
    if any(not isinstance(entry, LegendEntry) for entry in entries):
        msg = "legend_entries must contain LegendEntry values"
        raise TypeError(msg)
    if len({entry.key for entry in entries}) != len(entries):
        msg = "legend keys must be unique"
        raise ValueError(msg)
    return entries


def step_color(step: PlaybackStep, index: int, profile: str = "categorical") -> str:
    """Resolve a generic palette from declared kind and placement order."""
    if profile not in _PROFILES:
        msg = "unknown color profile; use explicit caller colors"
        raise ValueError(msg)
    if profile == "uniform":
        return "#667673"
    if profile == "constraints":
        return "#63558D" if step.placement_kind == "fixed_element" else "#687774"
    return DEFAULT_STEP_COLORS[index % len(DEFAULT_STEP_COLORS)]


def constraint_relation_color(profile: str) -> str:
    """Resolve the generic constraint relation color."""
    if profile not in _PROFILES:
        msg = "unknown color profile"
        raise ValueError(msg)
    return "#63558D" if profile == "constraints" else "#687876"
