"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/screening/models.py

Declared optional PWM screens and persisted motif-hit observations.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dataclasses import dataclass

from dense_arrays._record_validation import digest, object_fields, required_text
from dense_arrays.parts.motifs.artifacts import PWMArtifact
from dense_arrays.parts.scoring import FimoHit, FimoScoring
from dense_arrays.parts.scoring.configuration import finite
from dense_arrays.sequence import reverse_complement


@dataclass(frozen=True)
class PWMExclusion:
    """Reject a qualifying hit or a named score strictly above a threshold."""

    id: str
    motifs: tuple[PWMArtifact, ...]
    scoring: FimoScoring
    reject: str = "any_hit"
    score_field: str | None = None
    threshold: float | None = None

    def __post_init__(self) -> None:
        """Require explicit scoring and nonempty, ordered motif sources."""
        required_text(self.id, field_name="PWM exclusion id")
        if (
            not isinstance(self.motifs, (tuple, list))
            or not self.motifs
            or any(not isinstance(source, PWMArtifact) for source in self.motifs)
        ):
            msg = "PWM exclusion motifs require a nonempty array of PWMArtifact"
            raise TypeError(msg)
        object.__setattr__(self, "motifs", tuple(self.motifs))
        if not isinstance(self.scoring, FimoScoring):
            msg = "PWM exclusion requires FimoScoring"
            raise TypeError(msg)
        if self.reject == "any_hit":
            if self.score_field is not None or self.threshold is not None:
                msg = "any_hit exclusion cannot include a score field or threshold"
                raise ValueError(msg)
        elif self.reject == "score_above":
            if self.score_field not in {"raw", "fraction_of_max"}:
                msg = "score_above requires raw or fraction_of_max"
                raise ValueError(msg)
            object.__setattr__(
                self, "threshold", finite(self.threshold, "exclusion threshold")
            )
        else:
            msg = "PWM exclusion reject must be any_hit or score_above"
            raise ValueError(msg)

    def rejects(self, hit: FimoHit | None) -> bool:
        """Interpret a recorded qualifying hit without rescanning DNA."""
        if hit is None:
            return False
        if self.reject == "any_hit":
            return True
        score = hit.raw if self.score_field == "raw" else hit.fraction_of_max
        if score is None:
            msg = "fraction_of_max exclusion requires a positive theoretical maximum"
            raise ValueError(msg)
        return score > self.threshold


@dataclass(frozen=True)
class ScreenObservation:
    """One bound motif's best qualifying hit, including an explicit absent hit."""

    rule_id: str
    binding_id: str
    hit: FimoHit | None

    def __post_init__(self) -> None:
        """Bind valid labeled evidence to its rule and scoring model."""
        required_text(self.rule_id, field_name="screen rule_id")
        digest(self.binding_id, field_name="screen binding_id")
        if self.hit is not None and not isinstance(self.hit, FimoHit):
            msg = "screen hit requires FimoHit or None"
            raise TypeError(msg)

    def verify_sequence(self, sequence: str) -> None:
        """Check orientation and interval against the recorded full candidate."""
        if self.hit is not None:
            hit = self.hit
            core = sequence[hit.start : hit.end]
            if hit.strand == "reverse":
                core = reverse_complement(core)
            if core != hit.core:
                msg = "screen hit geometry disagrees with candidate sequence"
                raise ValueError(msg)

    def to_dict(self) -> dict[str, object]:
        """Keep absent hits distinct from missing observations."""
        return {
            "rule_id": self.rule_id,
            "binding_id": self.binding_id,
            "hit": self.hit.to_dict() if self.hit is not None else None,
        }

    @classmethod
    def from_dict(cls, value: object) -> "ScreenObservation":
        """Reject incomplete observations and validate labeled score evidence."""
        keys = {"rule_id", "binding_id", "hit"}
        data = object_fields(value, keys, "screen observation")
        if set(data) != keys:
            msg = "incomplete screen observation"
            raise ValueError(msg)
        data["hit"] = (
            FimoHit.from_dict(data["hit"]) if data["hit"] is not None else None
        )
        return cls(**data)
