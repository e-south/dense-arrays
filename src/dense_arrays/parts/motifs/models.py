"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/motifs/models.py

Immutable probability and score matrices with independent model identities.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import cached_property
from numbers import Real
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    immutable_json_mapping,
    mutable_json,
    object_fields,
    required_text,
    semantic_digest,
)

if TYPE_CHECKING:
    from collections.abc import Mapping

BASES = "ACGT"
PROBABILITY_TOLERANCE = 1e-3
UNIT_ROUNDING_TOLERANCE = 1e-15
MOTIF_SCHEMA = "dense_arrays.motif.v1"
MODEL_SCHEMA = "dense_arrays.motif_model.v1"


def numeric_row(
    value: object, *, probability: bool, positive: bool = False
) -> tuple[float, ...]:
    """Validate exact A/C/G/T columns; normalize only near-unit probabilities."""
    if not isinstance(value, (tuple, list)) or len(value) != len(BASES):
        msg = "motif rows require four ordered ACGT values"
        raise ValueError(msg)
    if any(
        isinstance(v, bool) or not isinstance(v, Real) or not math.isfinite(v)
        for v in value
    ):
        msg = "motif matrix values must be finite numbers"
        raise ValueError(msg)
    row = tuple(float(v) for v in value)
    if probability:
        if any(v < 0 or (positive and v == 0) for v in row):
            msg = "motif probabilities must be nonnegative; background must be positive"
            raise ValueError(msg)
        try:
            total = math.fsum(row)
        except OverflowError as err:
            msg = "motif probability rows must sum to one"
            raise ValueError(msg) from err
        if not math.isfinite(total) or abs(total - 1) > PROBABILITY_TOLERANCE:
            msg = "motif probability rows must sum to one within 0.001"
            raise ValueError(msg)
        if abs(total - 1) > UNIT_ROUNDING_TOLERANCE:
            row = tuple(v / total for v in row)
    return row


def _matrix(value: object, *, probability: bool) -> tuple[tuple[float, ...], ...]:
    if not isinstance(value, (list, tuple)) or not value:
        msg = "motif matrices require nonempty ordered rows"
        raise ValueError(msg)
    return tuple(numeric_row(row, probability=probability) for row in value)


@dataclass(frozen=True)
class Motif:
    """A labeled probability matrix with an optional declared score matrix.

    ``log_odds=None`` records unavailable supplied scores. Scoring such a model
    requires an explicit backend; probabilities alone imply no smoothing rule.
    """

    motif_id: str
    probabilities: tuple[tuple[float, ...], ...]
    background: tuple[float, ...]
    log_odds: tuple[tuple[float, ...], ...] | None
    producer: str
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Freeze matrix data without inferring a log base or smoothing rule."""
        required_text(self.motif_id, field_name="motif_id")
        required_text(self.producer, field_name="motif.producer")
        object.__setattr__(
            self, "probabilities", _matrix(self.probabilities, probability=True)
        )
        object.__setattr__(
            self,
            "background",
            numeric_row(self.background, probability=True, positive=True),
        )
        if self.log_odds is not None:
            object.__setattr__(
                self, "log_odds", _matrix(self.log_odds, probability=False)
            )
        if self.log_odds is not None and len(self.log_odds) != self.width:
            msg = "log-odds rows must match the probability matrix width"
            raise ValueError(msg)
        object.__setattr__(self, "metadata", immutable_json_mapping(self.metadata))

    @property
    def width(self) -> int:
        """Return the number of aligned motif positions."""
        return len(self.probabilities)

    def model(self) -> dict[str, object]:
        """Separate computational identity from labels and producer annotations."""
        return {
            "schema": MODEL_SCHEMA
            if self.log_odds is not None
            else "dense_arrays.motif_model.v2",
            "alphabet": BASES,
            "probabilities": [list(row) for row in self.probabilities],
            "background": list(self.background),
            "log_odds": [list(row) for row in self.log_odds]
            if self.log_odds is not None
            else None,
            "normalization": "near_unit_0.001.v1",
            "score_units": "declared_log_odds" if self.log_odds is not None else None,
        }

    @cached_property
    def model_id(self) -> str:
        """Fingerprint the exact matrices, normalization and score interpretation."""
        return semantic_digest(self.model())

    def to_dict(self) -> dict[str, object]:
        """Serialize normalized matrices and immutable source annotations."""
        return {
            **self.model(),
            "schema": MOTIF_SCHEMA
            if self.log_odds is not None
            else "dense_arrays.motif.v2",
            "model_id": self.model_id,
            "motif_id": self.motif_id,
            "producer": self.producer,
            "metadata": mutable_json(self.metadata),
        }

    @classmethod
    def from_dict(cls, value: object) -> Motif:
        """Read normalized motif evidence and reject changed interpretations."""
        keys = {
            "schema",
            "model_id",
            "motif_id",
            "producer",
            "metadata",
            "alphabet",
            "probabilities",
            "background",
            "log_odds",
            "normalization",
            "score_units",
        }
        data = object_fields(value, keys, "motif")
        if set(data) != keys or data["schema"] not in {
            MOTIF_SCHEMA,
            "dense_arrays.motif.v2",
        }:
            msg = "unsupported or incomplete motif record"
            raise ValueError(msg)
        result = cls(
            data["motif_id"],
            data["probabilities"],
            data["background"],
            data["log_odds"],
            data["producer"],
            data["metadata"],
        )
        if result.to_dict() != data:
            msg = "motif model identity or normalized fields mismatch"
            raise ValueError(msg)
        return result
