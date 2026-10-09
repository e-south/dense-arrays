"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/ingestion.py

Read one immutable byte snapshot through the curated part contract.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import TYPE_CHECKING

from dense_arrays._record_validation import integer, records, required_text
from dense_arrays.parts.models import Part, PartTable
from dense_arrays.parts.provenance import ImportReport, Transformation
from dense_arrays.parts.tables import open_table
from dense_arrays.parts.tables.diagnostics import RowDiagnostic, TableImportError

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence


@dataclass(frozen=True)
class ImportedParts:
    """Validated parts and bounded import provenance from one source snapshot."""

    parts: tuple[Part, ...]
    source_digest: str | None = None
    ignored_columns: tuple[str, ...] = ()
    transformations: tuple[tuple[int, str, str, str], ...] = ()
    report: ImportReport | None = None


def validate_parts(parts: Sequence[Part]) -> tuple[Part, ...]:
    """Freeze a nonempty pool with unique supplied occurrence identities."""
    result = records(parts, Part, field_name="parts")
    if not result:
        msg = "parts must contain at least one eligible occurrence"
        raise ValueError(msg)
    ids = [part.part_id for part in result]
    if len(ids) != len(set(ids)):
        msg = "part_id values must be unique; repeated DNA needs distinct IDs"
        raise ValueError(msg)
    return result


def read_parts(source: PartTable | Sequence[Part]) -> ImportedParts:
    """Validate typed parts or parse one captured table byte sequence."""
    if not isinstance(source, PartTable):
        parts = validate_parts(source)
        return ImportedParts(parts, report=ImportReport("inline", len(parts)))
    data = source.table.read_bytes()
    with open_table(source, data) as table:
        columns = _resolve_columns(source, table.header)
        mapped = set(columns.values()) | set(source.metadata_columns.values())
        ignored = set(table.header) - mapped
        parts, changes = _read_rows(source, columns, table.rows(tuple(sorted(mapped))))
    return ImportedParts(
        validate_parts(parts),
        hashlib.sha256(data).hexdigest(),
        tuple(sorted(ignored)),
        tuple(changes),
        ImportReport(
            "table",
            len(parts),
            columns,
            tuple(sorted(ignored)),
            tuple(
                Transformation(row, field, columns[field], before, after)
                for row, field, before, after in changes
            ),
            source.normalization,
            source.id_policy,
        ),
    )


def _read_rows(
    source: PartTable,
    columns: dict[str, str],
    rows: Iterator[dict[str, object]],
) -> tuple[list[Part], list[tuple[int, str, str, str]]]:
    parts = []
    changes = []
    seen = {}
    diagnostics = []
    invalid_rows = number = 0
    locations = {
        **columns,
        **{
            f"metadata.{name}": column
            for name, column in source.metadata_columns.items()
        },
    }
    for number, row in enumerate(rows, 1):
        issue = _duplicate_row(columns, row, number, seen)
        if issue is None:
            try:
                part, changed = _read_row(source, columns, row, number)
            except (ValueError, TypeError) as err:
                issue = RowDiagnostic(number, locations, str(err))
        if issue is not None:
            invalid_rows += 1
            parts.clear()
            changes.clear()
            if len(diagnostics) < TableImportError.sample_limit:
                diagnostics.append(issue)
        elif not invalid_rows:
            parts.append(part)
            changes.extend(changed)
    if invalid_rows:
        raise TableImportError(source.table, number, invalid_rows, diagnostics)
    return parts, changes


def _duplicate_row(
    columns: dict[str, str],
    row: dict[str, object],
    number: int,
    seen: dict[str, int],
) -> RowDiagnostic | None:
    if "part_id" not in columns:
        return None
    identity = row[columns["part_id"]]
    if not isinstance(identity, str) or not identity.strip():
        return None
    previous = seen.setdefault(identity, number)
    if previous == number:
        return None
    return RowDiagnostic(
        number,
        {"part_id": columns["part_id"]},
        f"part_id must be unique; already supplied at row {previous}",
        related_row=previous,
    )


def _resolve_columns(source: PartTable, header: tuple[object, ...]) -> dict[str, str]:
    if (
        not header
        or any(not isinstance(c, str) or not c.strip() for c in header)
        or len(header) != len(set(header))
    ):
        msg = f"{source.table}: expected nonempty, unique header columns"
        raise ValueError(msg)
    names = set(Part.__dataclass_fields__) - {"metadata"}
    columns = {name: name for name in names if name in header}
    columns.update(source.columns)
    required = {"sequence"} | ({"part_id"} if source.id_policy == "provided" else set())
    for name in required:
        columns.setdefault(name, name)
    if source.id_policy == "row" and "part_id" in columns:
        msg = "id_policy=row cannot discard a provided part_id column"
        raise ValueError(msg)
    missing = (set(columns.values()) | set(source.metadata_columns.values())) - set(
        header
    )
    if missing:
        msg = f"{source.table}: missing columns {sorted(missing)}"
        raise ValueError(msg)
    return columns


def _read_row(
    source: PartTable, columns: dict[str, str], row: dict[str, object], number: int
) -> tuple[Part, list[tuple[int, str, str, str]]]:
    fields = {name: row[column] for name, column in columns.items()}
    if source.id_policy == "row":
        fields["part_id"] = f"row:{number}"
    changes = []
    original = required_text(fields["sequence"], field_name="sequence")
    normalized = original
    if source.normalization.trim_outer_whitespace:
        normalized = normalized.strip()
    if source.normalization.uppercase:
        normalized = normalized.upper()
    if normalized != original:
        fields["sequence"] = normalized
        changes.append((number, "sequence", original, normalized))
    for name in ("group", "source", "core_start", "core_end", "core_orientation"):
        if name in fields and fields[name] == "":
            fields[name] = None
    for name in ("core_start", "core_end"):
        if fields.get(name) is not None:
            fields[name] = integer(
                int(fields[name]) if isinstance(fields[name], str) else fields[name],
                field_name=name,
            )
    fields["metadata"] = {
        key: row[column] for key, column in source.metadata_columns.items()
    }
    return Part(**fields), changes
