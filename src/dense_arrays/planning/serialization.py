"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/serialization.py

Versioned, JSON-compatible design requests with strict nested fields.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import os
from collections.abc import Mapping
from dataclasses import fields, replace
from types import MappingProxyType
from typing import TYPE_CHECKING

from dense_arrays._record_validation import mutable_json, object_fields
from dense_arrays.parts import (
    BoundParts,
    Normalization,
    Part,
    PartFilter,
    PartSelector,
    PartTable,
    PoolHandle,
    PoolSource,
)
from dense_arrays.parts.serialization import table_to_dict
from dense_arrays.planning.batches import BatchSchedule, CandidateBatch, Resampling
from dense_arrays.planning.libraries import LibraryExclusion
from dense_arrays.planning.lineage import Lineage
from dense_arrays.planning.models import (
    Assembly,
    DesignSpec,
    Length,
    Limits,
    Padding,
    Target,
)
from dense_arrays.planning.requirements import (
    GC,
    Avoid,
    Fixed,
    GroupCoverage,
    Occurrences,
    Requirement,
    Spacing,
    StartWindow,
)

if TYPE_CHECKING:
    from pathlib import Path

DESIGN_SCHEMA = "dense_arrays.design.v1"
PLAN_SCHEMA = "dense_arrays.generation_plan.v1"
POLICIES = MappingProxyType(
    {
        "packing": "oriented_path.v1",
        "objective": "selected_occurrences.v1",
        "proof": "optimal_required.v1",
        "enumeration": "ordered_path_exclusion.v1",
        "uniqueness": "exact_sequence_per_cell.v1",
        "randomness": "logical_work_sha256.v1",
        "backend": "CBC",
        "assembly": "single_side_uniform_shake256.v1",
    }
)


def _record(value: object) -> dict[str, object]:
    return {f.name: mutable_json(getattr(value, f.name)) for f in fields(value)}


def policies_for(request: DesignSpec) -> dict[str, str]:
    """Bind only requested additions; default plan policy identities stay stable."""
    return {
        **POLICIES,
        **(
            {
                "objective": "greedy_occurrences_then_length.v1",
                "proof": "heuristic_unproven.v1",
                "enumeration": "one_greedy_candidate_per_batch.v1",
                "backend": "greedy_multistart.v1",
            }
            if request.search == "greedy"
            else {}
        ),
        **(
            {"packing_preference": "underused_parts.v1"}
            if request.packing_preference
            else {}
        ),
    }


def requirement_to_dict(requirement: Requirement) -> dict[str, object]:
    """Encode one rule without serializing its unrelated part collection."""
    value = _record(requirement)
    if isinstance(requirement, Occurrences):
        value.update(kind="occurrences", select=_record(requirement.select))
    elif isinstance(requirement, GroupCoverage):
        value["kind"] = "group_coverage"
    elif isinstance(requirement, Fixed):
        value.update(
            kind="fixed",
            start=None if requirement.start is None else _record(requirement.start),
        )
    elif isinstance(requirement, Spacing):
        value["kind"] = "spacing"
    elif isinstance(requirement, Avoid):
        value["kind"] = "avoid"
    elif isinstance(requirement, GC):
        value["kind"] = "gc"
    else:
        msg = f"unsupported requirement type: {type(requirement).__name__}"
        raise TypeError(msg)
    return value


def request_to_dict(
    request: DesignSpec, *, base: Path | None = None
) -> dict[str, object]:
    """Serialize explicit sources or resolved parts without reading external inputs."""
    requirements = [requirement_to_dict(r) for r in request.requirements]
    return {
        "schema": DESIGN_SCHEMA,
        "parts": parts_to_dict(request.parts, base=base),
        "length": _record(request.length),
        "requirements": requirements,
        "target": _record(request.target),
        "seed": request.seed,
        "limits": request.limits.to_dict(),
        "strands": request.strands,
        **({"search": request.search} if request.search != "exact" else {}),
        **(
            {"packing_preference": request.packing_preference}
            if request.packing_preference
            else {}
        ),
        "assembly": None
        if request.assembly is None
        else {
            "padding": None
            if request.assembly.padding is None
            else _record(request.assembly.padding)
        },
        **(
            {"resampling": request.resampling.to_dict()}
            if request.resampling is not None
            else {}
        ),
        **({"batch": request.batch.to_dict()} if request.batch is not None else {}),
        **(
            {"schedule": request.schedule.to_dict()}
            if request.schedule is not None
            else {}
        ),
        **(
            {"lineage": request.lineage.to_dict()}
            if request.lineage is not None
            else {}
        ),
        **(
            {"exclude": request.exclude.to_dict(base=base)}
            if request.exclude is not None
            else {}
        ),
    }


def parts_to_dict(
    source: tuple[Part, ...] | PartTable | PoolSource | BoundParts, *, base: Path | None
) -> dict[str, object] | list[dict[str, object]]:
    """Keep locators relative to their document and preserve source predicates."""
    if isinstance(source, BoundParts):
        return source.to_dict(base=base)
    if isinstance(source, PartTable):
        value = table_to_dict(source)
        if base is not None:
            value["table"] = os.path.relpath(source.table.absolute(), base)
        return value
    if isinstance(source, PoolSource):
        if isinstance(source.pool, PoolHandle):
            msg = (
                "a PoolHandle identity cannot be encoded as a path-only request; "
                "export the resolved plan to preserve its bound evidence"
            )
            raise TypeError(msg)
        selected = None if source.select is None else source.select.to_dict()
        if selected is not None:
            selected.pop("schema")
        return {
            "pool": str(source.path)
            if base is None
            else os.path.relpath(source.path.absolute(), base),
            "select": selected,
        }
    return [_record(part) for part in source]


def requirement_from_dict(value: object) -> Requirement:
    """Parse one declared rule through the shared typed requirement constructors."""
    if not isinstance(value, Mapping):
        msg = "requirement must be an object with kind"
        raise TypeError(msg)
    value = dict(value)
    kind = value.pop("kind", None)
    if kind == "occurrences":
        selector = value.get("select")
        value["select"] = PartSelector(
            **object_fields(selector, {"part_ids", "groups"}, "select")
        )
        return Occurrences(**value)
    if kind == "group_coverage":
        return GroupCoverage(**value)
    if kind == "fixed":
        if value.get("start") is not None:
            value["start"] = StartWindow(**value["start"])
        return Fixed(**value)
    if kind == "spacing":
        return Spacing(**value)
    if kind == "avoid":
        return Avoid(**value)
    if kind == "gc":
        return GC(**value)
    msg = f"unsupported requirement kind {kind!r}"
    raise ValueError(msg)


def request_from_dict(value: object, *, base: Path | None = None) -> DesignSpec:
    """Parse declared requests without invoking preparation or optimization."""
    value = object_fields(
        value, {"schema", *(f.name for f in fields(DesignSpec))}, "design"
    )
    if value.pop("schema", None) != DESIGN_SCHEMA:
        msg = "unsupported design schema"
        raise ValueError(msg)
    source = parts_from_dict(value.pop("parts", None), base=base)
    for name, record_type in (
        ("length", Length),
        ("target", Target),
        ("limits", Limits),
    ):
        if name in value:
            value[name] = record_type(**value[name])
    requirements = value.pop("requirements", [])
    for name, record_type in (
        ("schedule", BatchSchedule),
        ("resampling", Resampling),
        ("batch", CandidateBatch),
        ("lineage", Lineage),
    ):
        if value.get(name) is not None:
            value[name] = record_type.from_dict(value[name])
    if value.get("exclude") is not None:
        value["exclude"] = LibraryExclusion.from_dict(value["exclude"], base=base)
    if value.get("assembly") is not None:
        assembly = object_fields(value["assembly"], {"padding"}, "assembly")
        if assembly.get("padding") is not None:
            assembly["padding"] = Padding(**assembly["padding"])
        value["assembly"] = Assembly(**assembly)
    if not isinstance(requirements, list):
        msg = "requirements must be an array"
        raise TypeError(msg)
    return DesignSpec(
        parts=source,
        requirements=tuple(requirement_from_dict(r) for r in requirements),
        **value,
    )


def parts_from_dict(
    source: object, *, base: Path | None
) -> tuple[Part, ...] | PoolSource | PartTable | BoundParts:
    """Resolve explicit source kinds without mixing execution settings into import."""
    if isinstance(source, Mapping) and "schema" in source:
        return BoundParts.from_dict(source, base=base)
    if isinstance(source, list):
        source = tuple(
            Part(**object_fields(item, {f.name for f in fields(Part)}, "part"))
            for item in source
        )
    elif isinstance(source, Mapping) and "pool" in source:
        source = object_fields(source, {"pool", "select"}, "pool source")
        if source.get("select") is not None:
            source["select"] = PartFilter.from_dict(source["select"], declared=False)
        source = PoolSource(**source)
        if base is not None and not source.path.is_absolute():
            source = replace(source, pool=base / source.path)
    elif isinstance(source, Mapping):
        source = dict(source)
        normalization = source.pop("normalization", {})
        source["normalization"] = Normalization(**normalization)
        source = PartTable(**source)
        if base is not None and not source.table.is_absolute():
            source = replace(source, table=base / source.table)
    else:
        msg = "parts must be a typed table source or an array of part records"
        raise TypeError(msg)
    return source
