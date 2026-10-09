"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/candidates.py

Recorded candidate observations and explicit retention decisions.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass

from dense_arrays._record_validation import integer, object_fields, required_text
from dense_arrays.parts.models import Part
from dense_arrays.parts.retention import MMRDecision
from dense_arrays.parts.scoring import FimoHit
from dense_arrays.parts.screening import ScreenObservation
from dense_arrays.parts.serialization import part_from_dict, part_to_dict


@dataclass(frozen=True)
class Candidate:
    """One mined sequence with frozen screening and representative evidence."""

    index: int
    part: Part
    reasons: tuple[str, ...] = ()
    representative: int | None = None
    rank: int | None = None
    retained: bool = False
    error: str | None = None
    selection: MMRDecision | None = None
    screening: tuple[ScreenObservation, ...] = ()
    recipe_id: str | None = None
    recipe_index: int | None = None
    score_band: int | None = None

    def __post_init__(self) -> None:
        """Reject contradictory decisions and malformed source/core geometry."""
        integer(self.index, field_name="candidate.index", minimum=1)
        self._validate_origin()
        if not isinstance(self.part, Part):
            msg = "candidate requires Part"
            raise TypeError(msg)
        self._validate_screening()
        if not isinstance(self.reasons, (tuple, list)):
            msg = "candidate reasons must be an ordered array"
            raise TypeError(msg)
        for reason in self.reasons:
            required_text(reason, field_name="candidate.reason")
        object.__setattr__(self, "reasons", tuple(self.reasons))
        if len(set(self.reasons)) != len(self.reasons):
            msg = "candidate reasons must be unique"
            raise ValueError(msg)
        for name in ("representative", "rank"):
            if getattr(self, name) is not None:
                integer(getattr(self, name), field_name=f"candidate.{name}", minimum=1)
        if not isinstance(self.retained, bool):
            msg = "candidate retained flag must be boolean"
            raise TypeError(msg)
        if self.error is not None:
            required_text(self.error, field_name="candidate.error")
        self._validate_score_and_selection()
        self._validate_score_band()

    def _validate_score_band(self) -> None:
        if self.score_band is not None:
            integer(self.score_band, field_name="candidate.score_band", minimum=1)
            if self.representative != self.index or self.score is None:
                msg = "score bands require a scored eligible representative"
                raise ValueError(msg)

    def _validate_origin(self) -> None:
        """Keep recipe-local coordinates complete and separate from pool order."""
        if (self.recipe_id is None) != (self.recipe_index is None):
            msg = "candidate recipe ID and local index must be declared together"
            raise ValueError(msg)
        if self.recipe_id is not None:
            required_text(self.recipe_id, field_name="candidate.recipe_id")
            integer(self.recipe_index, field_name="candidate.recipe_index", minimum=1)

    def _validate_screening(self) -> None:
        if not isinstance(self.screening, (list, tuple)) or any(
            not isinstance(observation, ScreenObservation)
            for observation in self.screening
        ):
            msg = "candidate screening requires ordered ScreenObservation records"
            raise TypeError(msg)
        object.__setattr__(self, "screening", tuple(self.screening))
        for observation in self.screening:
            observation.verify_sequence(self.part.sequence)

    def _validate_score_and_selection(self) -> None:
        if self.selection is not None:
            if (
                not isinstance(self.selection, MMRDecision)
                or self.representative != self.index
            ):
                msg = "MMR decisions require an eligible representative"
                raise ValueError(msg)
            if (self.selection.utility is not None) != self.retained:
                msg = "MMR utility belongs to retained candidates only"
                raise ValueError(msg)
        score = self.score
        if score is not None and (
            score.core != self.part.core_sequence
            or (score.start, score.end, score.strand)
            != (self.part.core_start, self.part.core_end, self.part.core_orientation)
        ):
            msg = "candidate score and part core geometry disagree"
            raise ValueError(msg)
        if (self.error or self.reasons) and any(
            x is not None for x in (self.representative, self.rank)
        ):
            msg = "rejected candidates cannot carry selection decisions"
            raise ValueError(msg)
        if self.retained and (self.rank is None or self.representative != self.index):
            msg = "retained candidates require their own representative and rank"
            raise ValueError(msg)

    @property
    def score(self) -> FimoHit | None:
        """Decode the labeled score evidence recorded with the candidate part."""
        value = self.part.metadata.get("score")
        return None if value is None else FimoHit.from_dict(value)

    def to_dict(self) -> dict[str, object]:
        """Publish all rejection reasons and the retained representative mapping."""
        return {
            "schema": "dense_arrays.preparation_candidate.v1",
            "index": self.index,
            "part": part_to_dict(self.part),
            "reasons": list(self.reasons),
            "representative": self.representative,
            "rank": self.rank,
            "retained": self.retained,
            "error": self.error,
            **({"score_band": self.score_band} if self.score_band is not None else {}),
            **(
                {"recipe_id": self.recipe_id, "recipe_index": self.recipe_index}
                if self.recipe_id is not None
                else {}
            ),
            **(
                {"screening": [item.to_dict() for item in self.screening]}
                if self.screening
                else {}
            ),
            **(
                {"selection": self.selection.to_dict()}
                if self.selection is not None
                else {}
            ),
        }

    @classmethod
    def from_dict(cls, value: object) -> Candidate:
        """Validate complete recorded candidate evidence."""
        keys = {
            "schema",
            "index",
            "part",
            "reasons",
            "representative",
            "rank",
            "retained",
            "error",
        }
        data = object_fields(
            value,
            keys
            | {"selection", "screening", "recipe_id", "recipe_index", "score_band"},
            "preparation candidate",
        )
        recipe_id = data.pop("recipe_id", None)
        score_band = data.pop("score_band", None)
        recipe_index = data.pop("recipe_index", None)
        screening = data.pop("screening", [])
        if not isinstance(screening, list):
            msg = "candidate screening must be an ordered array"
            raise TypeError(msg)
        selection = data.pop("selection", None)
        if (
            set(data) != keys
            or data.pop("schema") != "dense_arrays.preparation_candidate.v1"
        ):
            msg = "unsupported or incomplete preparation candidate"
            raise ValueError(msg)
        data["part"] = part_from_dict(data["part"])
        data["screening"] = tuple(
            ScreenObservation.from_dict(item) for item in screening
        )
        if selection is not None:
            data["selection"] = MMRDecision.from_dict(selection)
        return cls(
            **data,
            recipe_id=recipe_id,
            recipe_index=recipe_index,
            score_band=score_band,
        )
