"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/receipts.py

Portable publication receipts shared by data exports and visual artifacts.

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
    mutable_json,
    required_text,
)

if TYPE_CHECKING:
    from collections.abc import Mapping


@dataclass(frozen=True)
class ExportReceipt:
    """Published destination and exact source identities, without owning open files."""

    destination: str
    format: str
    view: str
    records: int
    sources: tuple[Mapping[str, object], ...]
    design_refs: tuple[str, ...]
    files: tuple[Mapping[str, object], ...]
    selection: Mapping[str, object] | None = None

    def __post_init__(self) -> None:
        """Detach publication evidence from mutable caller collections."""
        for name in ("destination", "format", "view"):
            required_text(getattr(self, name), field_name=name)
        integer(self.records, field_name="records", minimum=0)
        object.__setattr__(
            self, "sources", tuple(immutable_json_mapping(s) for s in self.sources)
        )
        object.__setattr__(
            self, "files", tuple(immutable_json_mapping(f) for f in self.files)
        )
        object.__setattr__(self, "design_refs", tuple(self.design_refs))
        if self.selection is not None:
            object.__setattr__(
                self, "selection", immutable_json_mapping(self.selection)
            )

    def to_dict(self) -> dict[str, object]:
        """Encode successful publication; failures do not return a success receipt."""
        return {
            "schema": "dense_arrays.export_receipt.v1",
            "destination": self.destination,
            "format": self.format,
            "view": self.view,
            "records": self.records,
            "sources": mutable_json(self.sources),
            "design_refs": list(self.design_refs),
            "files": mutable_json(self.files),
            **(
                {"selection": mutable_json(self.selection)}
                if self.selection is not None
                else {}
            ),
        }
