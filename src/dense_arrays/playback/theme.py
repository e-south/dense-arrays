"""Semantic presentation profiles shared by graph and duplex playback."""

from __future__ import annotations

from dataclasses import dataclass

from .models import PlaybackStep

_DEFAULT_COLORS = (
    "#67BFA5",
    "#D883A4",
    "#7BA4D9",
    "#C08A56",
    "#5DA79F",
    "#D1B06C",
    "#74C0CB",
    "#86A5D8",
)
_SLATE = "#667673"
_CONSTRAINT_SLATE = "#687774"
_FIXED_PURPLE = "#63558D"
_SECG_COLORS = {
    "baer": "#9C572B",
    "cpxr": "#267C73",
    "lexa": "#A34770",
}


@dataclass(frozen=True, slots=True)
class LegendEntry:
    key: str
    label: str
    color: str


@dataclass(frozen=True, slots=True)
class PlaybackPresentation:
    color_profile: str = "categorical"
    legend_entries: tuple[LegendEntry, ...] = ()
    graph_detail: str = "full"
    graph_fraction: float = 0.35
    show_edge_costs: bool = True
    show_authority_notice: bool = False
    show_distance_bracket: str = "when_declared"

    def __post_init__(self) -> None:
        if self.graph_detail not in {"full", "reduced", "inset", "none"}:
            raise ValueError(
                "graph_detail must be 'full', 'reduced', 'inset', or 'none'"
            )
        if not 0.0 <= float(self.graph_fraction) <= 0.5:
            raise ValueError("graph_fraction must be between 0.0 and 0.5")
        if self.graph_detail == "none" and float(self.graph_fraction) != 0.0:
            raise ValueError("graph_fraction must be 0.0 when graph_detail is 'none'")
        if self.graph_detail != "none" and float(self.graph_fraction) <= 0.0:
            raise ValueError("graph_fraction must be positive when a graph is shown")
        if self.show_distance_bracket not in {"never", "when_declared", "always"}:
            raise ValueError(
                "show_distance_bracket must be 'never', 'when_declared', or 'always'"
            )


def step_category(step: PlaybackStep) -> str:
    label = str(step.label or "")
    identity = f"{label} {step.placement_id}".lower()
    if step.placement_kind == "fixed_element":
        if "downstream" in identity or "-10" in identity:
            return "fixed_downstream"
        return "fixed_upstream"
    if "baer" in identity or "bayr" in identity:
        return "baer"
    if "cpxr" in identity:
        return "cpxr"
    if "lexa" in identity:
        return "lexa"
    return "binding_site"


def step_color(step: PlaybackStep, index: int, profile: str = "categorical") -> str:
    category = step_category(step)
    if profile == "uniform":
        return _SLATE
    if profile == "constraints":
        return _FIXED_PURPLE if category.startswith("fixed_") else _CONSTRAINT_SLATE
    if profile == "secg":
        if category in _SECG_COLORS:
            return _SECG_COLORS[category]
        return _FIXED_PURPLE if category.startswith("fixed_") else _SLATE
    if category == "fixed_downstream":
        return "#C886D1"
    if category == "fixed_upstream":
        return "#7D86D1"
    return _DEFAULT_COLORS[index % len(_DEFAULT_COLORS)]


def constraint_relation_color(profile: str) -> str:
    return _FIXED_PURPLE if profile in {"constraints", "secg"} else "#687876"


def step_text_color(
    step: PlaybackStep, index: int, profile: str = "categorical"
) -> str:
    """Return the publication text color for a step fill."""
    del step, index, profile
    return "#FFFFFF"


def legend_entries_for_profile(profile: str) -> tuple[LegendEntry, ...]:
    if profile != "secg":
        return ()
    return (
        LegendEntry("baer", "BaeR binding site", _SECG_COLORS["baer"]),
        LegendEntry("cpxr", "CpxR binding site", _SECG_COLORS["cpxr"]),
        LegendEntry("lexa", "LexA binding site", _SECG_COLORS["lexa"]),
        LegendEntry("fixed", "Fixed promoter element", _FIXED_PURPLE),
    )
