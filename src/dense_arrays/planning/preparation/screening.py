"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/preparation/screening.py

Resolve optional screens and verify their saved observation bindings.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING

from dense_arrays._record_validation import object_fields, required_text
from dense_arrays.parts.screening import PWMExclusion

from .motifs import MotifSource, resolve_motif

if TYPE_CHECKING:
    from pathlib import Path

    from dense_arrays.parts import PreparationSpec
    from dense_arrays.parts.candidates import Candidate


@dataclass(frozen=True)
class BoundExclusion:
    """A named screen with ordered, immutable motif/scorer bindings."""

    rule_id: str
    motifs: tuple[MotifSource, ...]

    def __post_init__(self) -> None:
        """Reject repeated motif identities and incomplete bindings."""
        required_text(self.rule_id, field_name="bound screen rule_id")
        if (
            not isinstance(self.motifs, (tuple, list))
            or not self.motifs
            or any(not isinstance(x, MotifSource) for x in self.motifs)
        ):
            msg = "bound screen requires ordered motif sources"
            raise TypeError(msg)
        object.__setattr__(self, "motifs", tuple(self.motifs))
        if len({x.input.motif.motif_id for x in self.motifs}) != len(self.motifs):
            msg = "PWM exclusion motif identities must be unique within each rule"
            raise ValueError(msg)

    def to_dict(self, *, base: Path | None = None) -> dict[str, object]:
        """Save the rule binding separately from its declared threshold policy."""
        return {
            "rule_id": self.rule_id,
            "motifs": [x.to_dict(base=base) for x in self.motifs],
        }

    @classmethod
    def from_dict(cls, value: object, *, base: Path | None = None) -> BoundExclusion:
        """Restore complete bindings without invoking a tool or opening an input."""
        data = object_fields(value, {"rule_id", "motifs"}, "bound screen")
        if set(data) != {"rule_id", "motifs"} or not isinstance(data["motifs"], list):
            msg = "incomplete bound screen"
            raise ValueError(msg)
        return cls(
            data["rule_id"],
            tuple(MotifSource.from_dict(x, base=base) for x in data["motifs"]),
        )


def resolve_screens(
    request: PreparationSpec,
) -> tuple[PreparationSpec, tuple[BoundExclusion, ...]]:
    """Bind every declared motif before candidate generation or output ownership."""
    declarations, bindings = [], []
    for rule in request.screening:
        if not isinstance(rule, PWMExclusion):
            declarations.append(rule)
            continue
        resolved = [resolve_motif(motif, rule.scoring) for motif in rule.motifs]
        bindings.append(BoundExclusion(rule.id, tuple(item[1] for item in resolved)))
        declarations.append(
            replace(
                rule,
                motifs=tuple(item[0] for item in resolved),
                scoring=resolved[0][1].scoring.settings,
            )
        )
    return replace(request, screening=tuple(declarations)), tuple(bindings)


def validate_screens(
    request: PreparationSpec, screens: tuple[BoundExclusion, ...]
) -> None:
    """Check declared identities, frozen policies and per-call scoring bounds."""
    declared = [rule for rule in request.screening if isinstance(rule, PWMExclusion)]
    if len(declared) != len(screens):
        msg = "sampled screen bindings disagree with the request"
        raise ValueError(msg)
    for rule, bound in zip(declared, screens, strict=True):
        if rule.id != bound.rule_id or len(rule.motifs) != len(bound.motifs):
            msg = "sampled screen identities disagree with the request"
            raise ValueError(msg)
        for source, motif in zip(rule.motifs, bound.motifs, strict=True):
            if (
                source.path != motif.input.path
                or rule.scoring != motif.scoring.settings
                or source.window != motif.window
                or source.format != motif.input.format
                or (
                    source.motif_ids
                    and source.motif_ids != (motif.input.motif.motif_id,)
                )
            ):
                msg = "screen source/scorer binding disagrees with the request"
                raise ValueError(msg)
            batch = min(request.budget.candidates, request.budget.batch_size)
            if (
                motif.window_bound(request.sampling.maximum_length, batch)
                > rule.scoring.limits.windows
            ):
                msg = "screen batch exceeds scoring window limit, including calibration"
                raise ValueError(msg)


def verify_observations(
    candidate: Candidate, screens: tuple[BoundExclusion, ...]
) -> None:
    """Require complete observations, or an ordered prefix for a failed batch."""
    expected = [
        (screen.rule_id, motif.scoring) for screen in screens for motif in screen.motifs
    ]
    observed = candidate.screening
    if len(observed) > len(expected) or (
        not candidate.error and len(observed) != len(expected)
    ):
        msg = "candidate has missing or excess screen observations"
        raise ValueError(msg)
    for item, (rule_id, scoring) in zip(observed, expected, strict=False):
        if item.rule_id != rule_id or item.binding_id != scoring.binding_id:
            msg = "screen observation disagrees with its bound motif/scorer"
            raise ValueError(msg)
        scoring.verify_hit(item.hit)


def screening_preview(
    request: PreparationSpec, screens: tuple[BoundExclusion, ...]
) -> dict[str, object]:
    """Expose total motif count and full-budget oriented-window work."""
    if not screens:
        return {}
    count = request.budget.candidates
    calls = (count + request.budget.batch_size - 1) // request.budget.batch_size
    windows = 0
    for screen in screens:
        for motif in screen.motifs:
            strands = 1 if motif.scoring.settings.strands == "single" else 2
            windows += (
                motif.window_bound(request.sampling.maximum_length, count)
                + (calls - 1) * strands
            )
    return {
        "screening_motifs": sum(len(s.motifs) for s in screens),
        "screening_window_bound": windows,
        **(
            {
                "screening_windows": [
                    {
                        "rule_id": screen.rule_id,
                        "motif_id": motif.motif.motif_id,
                        "window": motif.selection.to_dict(),
                    }
                    for screen in screens
                    for motif in screen.motifs
                    if motif.selection
                ]
            }
            if any(m.window for s in screens for m in s.motifs)
            else {}
        ),
    }
