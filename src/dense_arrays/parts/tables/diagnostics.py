"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/tables/diagnostics.py

Bounded row diagnostics for rejected part-table imports.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    immutable_json_mapping,
    integer,
    records,
    required_text,
)

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence
    from pathlib import Path


@dataclass(frozen=True)
class RowDiagnostic:
    """One failed row with its field mappings and optional conflicting row."""

    row: int
    columns: Mapping[str, str]
    message: str
    related_row: int | None = None

    def __post_init__(self) -> None:
        """Retain one-based locations and immutable field-to-column mappings."""
        integer(self.row, field_name="diagnostic.row", minimum=1)
        required_text(self.message, field_name="diagnostic.message")
        object.__setattr__(self, "columns", immutable_json_mapping(self.columns))
        for name, column in self.columns.items():
            required_text(name, field_name="diagnostic.field")
            required_text(column, field_name="diagnostic.column")
        if self.related_row is not None:
            integer(self.related_row, field_name="diagnostic.related_row", minimum=1)
            if self.related_row >= self.row:
                msg = "a duplicate ID must refer to an earlier source row"
                raise ValueError(msg)

    def to_dict(self) -> dict[str, object]:
        """Return row evidence without retaining the source row's values."""
        return {
            "row": self.row,
            "columns": dict(self.columns),
            "message": self.message,
            "related_row": self.related_row,
        }


class TableImportError(ValueError):
    """Reject a table with exact invalid-row counts and the first 20 diagnostics.

    Each invalid row contributes once, using its first validation failure.
    Structural parser failures have no reliable total and raise immediately.
    """

    sample_limit = 20

    def __init__(
        self,
        source: Path,
        rows: int,
        invalid_rows: int,
        diagnostics: Sequence[RowDiagnostic],
    ) -> None:
        """Bind complete row counts to a bounded, ordered diagnostic sample."""
        self.source = source
        self.rows = integer(rows, field_name="import.rows", minimum=1)
        self.invalid_rows = integer(
            invalid_rows, field_name="import.invalid_rows", minimum=1
        )
        self.diagnostics = records(
            diagnostics, RowDiagnostic, field_name="import.diagnostics"
        )
        locations = [item.row for item in self.diagnostics]
        if (
            invalid_rows > rows
            or len(self.diagnostics) != min(invalid_rows, self.sample_limit)
            or locations != sorted(set(locations))
            or any(row > rows for row in locations)
        ):
            msg = "import diagnostic counts and row locations must agree"
            raise ValueError(msg)
        label = "row" if invalid_rows == 1 else "rows"
        lines = [
            (
                f"{source}: {invalid_rows} invalid {label} among {rows}; "
                f"showing {len(self.diagnostics)}. No parts imported."
            )
        ]
        for item in self.diagnostics:
            mapping = ", ".join(
                f"{column}→{name}" for name, column in sorted(item.columns.items())
            )
            lines.append(f"row {item.row}, columns {mapping}: {item.message}")
        super().__init__("\n".join(lines))

    def to_dict(self) -> dict[str, object]:
        """Expose the same counted report to Python and CLI callers."""
        return {
            "source": str(self.source),
            "rows": self.rows,
            "invalid_rows": self.invalid_rows,
            "sample_limit": self.sample_limit,
            "omitted_rows": self.invalid_rows - len(self.diagnostics),
            "diagnostics": [item.to_dict() for item in self.diagnostics],
        }
