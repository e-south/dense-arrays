"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/serialization.py

Strict source and part encodings shared by preparation, planning and artifacts.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dataclasses import fields

from dense_arrays._record_validation import mutable_json, object_fields
from dense_arrays.parts.models import Normalization, Part, PartTable


def part_to_dict(part: Part) -> dict[str, object]:
    """Encode the supplied occurrence without introducing a new identity."""
    return {f.name: mutable_json(getattr(part, f.name)) for f in fields(Part)}


def part_from_dict(value: object) -> Part:
    """Read one strict occurrence record, including its metadata namespace."""
    return Part(**object_fields(value, {f.name for f in fields(Part)}, "part"))


def table_to_dict(source: PartTable) -> dict[str, object]:
    """Serialize table settings explicitly; paths remain outside content identity."""
    value = {f.name: mutable_json(getattr(source, f.name)) for f in fields(PartTable)}
    value["table"] = str(source.table)
    value["normalization"] = {
        f.name: getattr(source.normalization, f.name) for f in fields(Normalization)
    }
    if source.sheet is None:
        value.pop("sheet")
    return value


def table_from_dict(value: object) -> PartTable:
    """Reject unknown source options and normalization settings."""
    data = object_fields(value, {f.name for f in fields(PartTable)}, "table")
    data["normalization"] = Normalization(**data.get("normalization", {}))
    return PartTable(**data)
