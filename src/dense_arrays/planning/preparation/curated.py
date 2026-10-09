"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/preparation/curated.py

Preview and freeze curated preparation without publishing or running tools.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    canonical_json,
    mutable_json,
    object_fields,
    semantic_digest,
)
from dense_arrays.artifacts.publication import write_new
from dense_arrays.parts import Part, PreparationSpec
from dense_arrays.parts.ingestion import validate_parts
from dense_arrays.parts.provenance import ImportReport
from dense_arrays.parts.serialization import (
    part_from_dict,
    part_to_dict,
)
from dense_arrays.planning.resolution import InputBinding

if TYPE_CHECKING:
    from collections.abc import Mapping

from .requests import preparation_from_dict, preparation_to_dict

PREPARATION_PLAN_SCHEMA = "dense_arrays.preparation_plan.v1"
PREPARATION_POLICIES = MappingProxyType(
    {
        "import": "strict_curated.v1",
        "retention": "ordered_filter.v1",
        "duplicates": "retain.v1",
    }
)


@dataclass(frozen=True, repr=False)
class CuratedPreparation:
    """An immutable normalized snapshot with exact curated retention counts."""

    request: PreparationSpec
    parts: tuple[Part, ...]
    input: InputBinding
    import_report: ImportReport
    retained_indices: tuple[int, ...] = field(init=False)
    plan_id: str = field(init=False)

    def __post_init__(self) -> None:
        """Resolve deterministic retention and content identity once."""
        if not isinstance(self.request, PreparationSpec) or not isinstance(
            self.input, InputBinding
        ):
            msg = "preparation plans require typed requests and input bindings"
            raise TypeError(msg)
        object.__setattr__(self, "parts", validate_parts(self.parts))
        if not isinstance(
            self.import_report, ImportReport
        ) or self.import_report.rows != len(self.parts):
            msg = "import report must describe the preparation source"
            raise ValueError(msg)
        selected = self.request.retain.select
        if selected is not None:
            selected.validate(self.parts)
        indices = tuple(
            i
            for i, part in enumerate(self.parts)
            if selected is None or selected.matches(part)
        )
        object.__setattr__(self, "retained_indices", indices)
        request = preparation_to_dict(self.request)
        request["source"].pop("table")
        identity = {
            "schema": PREPARATION_PLAN_SCHEMA,
            "request": request,
            "parts": [part_to_dict(p) for p in self.parts],
            "input_digest": self.input.sha256,
            "import_report": self.import_report.to_dict(),
            "policies": dict(PREPARATION_POLICIES),
        }
        object.__setattr__(self, "plan_id", semantic_digest(identity))

    def __repr__(self) -> str:
        """Show bounded identity and counts without listing every part."""
        return (
            f"CuratedPreparation({self.plan_id[:12]}, source={len(self.parts)}, "
            f"retained={len(self.retained_indices)})"
        )

    @property
    def preview(self) -> Mapping[str, object]:
        """Report deterministic curated counts without a scoring or mining claim."""
        return MappingProxyType(
            {
                "source_parts": len(self.parts),
                "retained_parts": len(self.retained_indices),
                "candidate_budget": 0,
                "required_tools": (),
                "retained_count_status": "exact",
                "screening_stages": ("strict_import", "part_filter"),
            }
        )

    def to_dict(self, *, base: Path | None = None) -> dict[str, object]:
        """Serialize bound input, complete provenance and every resolved policy."""
        return {
            "schema": PREPARATION_PLAN_SCHEMA,
            "plan_id": self.plan_id,
            "request": preparation_to_dict(self.request, base=base),
            "parts": [part_to_dict(p) for p in self.parts],
            "input": self.input.to_dict(base=base),
            "import_report": self.import_report.to_dict(),
            "policies": dict(PREPARATION_POLICIES),
            "retained_indices": list(self.retained_indices),
            "preview": mutable_json(dict(self.preview)),
        }

    def verify_inputs(self) -> None:
        """Reject changed source bytes; execution consumes this bound snapshot."""
        if (
            hashlib.sha256(self.input.path.read_bytes()).hexdigest()
            != self.input.sha256
        ):
            msg = f"input changed since planning: {self.input.path}; create a new plan"
            raise ValueError(msg)

    def write(self, path: str | Path) -> None:
        """Create a plan file atomically without replacing an existing destination."""
        path = Path(path).absolute()
        write_new(path, canonical_json(self.to_dict(base=path.parent)) + "\n")

    @classmethod
    def from_dict(
        cls, value: object, *, base: Path | None = None
    ) -> CuratedPreparation:
        """Reject incomplete, altered or unsupported resolved plans."""
        keys = {
            "schema",
            "plan_id",
            "request",
            "parts",
            "input",
            "import_report",
            "policies",
            "retained_indices",
            "preview",
        }
        data = object_fields(value, keys, "preparation plan")
        if (
            set(data) != keys
            or data["schema"] != PREPARATION_PLAN_SCHEMA
            or data["policies"] != PREPARATION_POLICIES
        ):
            msg = "unsupported or incomplete preparation plan schema or policies"
            raise ValueError(msg)
        plan = cls(
            preparation_from_dict(data["request"], base=base),
            tuple(part_from_dict(p) for p in data["parts"]),
            InputBinding.from_dict(data["input"], base=base),
            ImportReport.from_dict(data["import_report"]),
        )
        normalized = mutable_json(data)
        normalized["request"]["source"]["table"] = str(plan.request.source.table)
        normalized["input"]["path"] = str(plan.input.path)
        if plan.to_dict() != normalized:
            msg = (
                "preparation plan digest or normalized fields mismatch; "
                "defaults are not reapplied"
            )
            raise ValueError(msg)
        return plan
