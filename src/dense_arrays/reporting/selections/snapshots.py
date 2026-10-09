"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/selections/snapshots.py

Frozen selection membership, source revisions and allocation evidence.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import os
from collections import Counter
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path

from dense_arrays._record_validation import (
    digest,
    immutable_json_mapping,
    integer,
    mutable_json,
    object_fields,
    records,
    required_text,
    semantic_digest,
)
from dense_arrays.artifacts.reading import ReadCost, ReadLimits

from .requests import LibrarySelection, resolve_quotas

_CELL_COMPONENTS = 2
_DESIGN_COMPONENTS = 3

SNAPSHOT_SCHEMA = "dense_arrays.selection-snapshot.v1"


@dataclass(frozen=True)
class SelectionSource:
    """A committed source identity; its locator is outside semantic identity."""

    kind: str
    source_id: str
    revision: int
    manifest_digest: str
    path: Path
    cells: tuple[str, ...]

    def __post_init__(self) -> None:
        """Keep source identity and namespaces explicit even for empty populations."""
        if self.kind not in {"run", "bundle"}:
            msg = "selection sources must be runs or bundles"
            raise ValueError(msg)
        required_text(self.source_id, field_name="source_id")
        integer(self.revision, field_name="revision", minimum=0)
        digest(self.manifest_digest, field_name="manifest_digest")
        object.__setattr__(self, "path", Path(self.path).absolute())
        cells = records(self.cells, str, field_name="source cells")
        if len(set(cells)) != len(cells) or any(
            len(c.split("/")) != _CELL_COMPONENTS or not all(c.split("/"))
            for c in cells
        ):
            msg = "source cells must be distinct full run/cell references"
            raise ValueError(msg)
        object.__setattr__(self, "cells", cells)

    def content(self) -> dict[str, object]:
        """Encode the pinned evidence independently of its location."""
        return {
            "kind": self.kind,
            "source_id": self.source_id,
            "revision": self.revision,
            "manifest_digest": self.manifest_digest,
            "cells": list(self.cells),
        }

    def to_dict(self, *, base: Path | None = None) -> dict[str, object]:
        """Locate source evidence relative to a saved snapshot when possible."""
        return {
            **self.content(),
            "path": str(self.path)
            if base is None
            else os.path.relpath(self.path, base),
        }

    @classmethod
    def from_dict(cls, value: object, *, base: Path | None = None) -> SelectionSource:
        """Read an explicit locator without substituting another revision."""
        data = object_fields(
            value,
            {"kind", "source_id", "revision", "manifest_digest", "path", "cells"},
            "selection source",
        )
        if set(data) != {
            "kind",
            "source_id",
            "revision",
            "manifest_digest",
            "path",
            "cells",
        }:
            msg = "selection source is missing required fields"
            raise ValueError(msg)
        path = Path(required_text(data.pop("path"), field_name="source path"))
        return cls(
            **data, path=path if path.is_absolute() else (base or Path.cwd()) / path
        )


@dataclass(frozen=True)
class SelectedDesign:
    """Full identity and exact record content selected at a pinned revision."""

    reference: str
    record_digest: str

    def __post_init__(self) -> None:
        """Reject local-only references and malformed record digests."""
        required_text(self.reference, field_name="design reference")
        if len(self.reference.split("/")) != _DESIGN_COMPONENTS or not all(
            self.reference.split("/")
        ):
            msg = "selected designs require full run/cell/design references"
            raise ValueError(msg)
        digest(self.record_digest, field_name="record_digest")

    def to_dict(self) -> dict[str, str]:
        """Serialize a membership record without copying its source design."""
        return {"reference": self.reference, "record_digest": self.record_digest}


class SelectionShortfall(ValueError):  # noqa: N818 - public domain outcome
    """A valid allocation could not be filled from its eligible source population."""

    def __init__(self, counts: Mapping[str, Mapping[str, int]]) -> None:
        """Preserve exact deficits rather than silently redistributing quotas."""
        self.counts = immutable_json_mapping(counts)
        requested = sum(c["requested"] for c in counts.values())
        available = sum(c["available"] for c in counts.values())
        selected = sum(c["selected"] for c in counts.values())
        super().__init__(
            f"selection shortfall: requested {requested}, available {available}, "
            f"selected {selected}; inspect cell counts or explicitly allow_partial"
        )


@dataclass(frozen=True, repr=False)
class SelectionSnapshot:
    """An inspectable selection with lazy reference access and bounded repr."""

    sources: tuple[SelectionSource, ...]
    request: LibrarySelection
    members: tuple[SelectedDesign, ...]
    counts: Mapping[str, Mapping[str, int]]
    algorithm: str
    read_limits: ReadLimits = field(default_factory=ReadLimits, compare=False)

    def __post_init__(self) -> None:
        """Reconcile membership, allocation and declared namespaces before reuse."""
        object.__setattr__(
            self,
            "sources",
            records(self.sources, SelectionSource, field_name="sources"),
        )
        object.__setattr__(
            self, "members", records(self.members, SelectedDesign, field_name="members")
        )
        if not self.sources or not isinstance(self.request, LibrarySelection):
            msg = "selection snapshots require source bindings and LibrarySelection"
            raise ValueError(msg)
        expected = (
            "sha256_priority.v1"
            if self.request.take and self.request.take.policy == "random"
            else "source_order.v1"
        )
        if self.algorithm != expected:
            msg = "selection algorithm does not match its request"
            raise ValueError(msg)
        if not isinstance(self.read_limits, ReadLimits):
            msg = "selection snapshots require ReadLimits"
            raise TypeError(msg)
        _validate_counts(self)
        object.__setattr__(self, "counts", immutable_json_mapping(self.counts))

    @property
    def requested(self) -> int:
        """Requested records in the declared total or cell allocations."""
        return sum(c["requested"] for c in self.counts.values())

    @property
    def available(self) -> int:
        """Eligible records in the allocation's scope, including zero-quota cells."""
        return sum(c["available"] for c in self.counts.values())

    @property
    def selected(self) -> int:
        """Number of materialized full design identities."""
        return len(self.members)

    @property
    def shortfall(self) -> int:
        """Unfilled quotas; excess availability elsewhere never offsets a deficit."""
        return sum(c["shortfall"] for c in self.counts.values())

    @property
    def status(self) -> str:
        """Qualify membership completeness separately from native run completion."""
        return "partial" if self.shortfall else "complete"

    def references(self) -> Iterator[str]:
        """Iterate saved order without sampling or reopening source artifacts."""
        return (member.reference for member in self.members)

    def summary(self) -> dict[str, object]:
        """Expose allocation evidence without membership expansion or local paths."""
        return {
            "schema": "dense_arrays.selection_summary.v1",
            "snapshot_id": self.snapshot_id,
            "status": self.status,
            "requested": self.requested,
            "available": self.available,
            "selected": self.selected,
            "shortfall": self.shortfall,
            "counts": mutable_json(self.counts),
            "request": self.request.to_dict(),
            "algorithm": self.algorithm,
        }

    def _content(self) -> dict[str, object]:
        return {
            "schema": SNAPSHOT_SCHEMA,
            "sources": [s.content() for s in self.sources],
            "request": self.request.to_dict(),
            "members": [m.to_dict() for m in self.members],
            "counts": mutable_json(self.counts),
            "algorithm": self.algorithm,
        }

    @cached_property
    def snapshot_id(self) -> str:
        """Bind policy, evidence and membership without binding filesystem locations."""
        return semantic_digest(self._content())

    @property
    def cost(self) -> ReadCost:
        """Membership is already materialized; source reuse has a separate read cost."""
        return ReadCost(
            self.snapshot_id, 0, "manifest", "selection", 0, self.read_limits
        )

    def to_dict(self, *, base: Path | None = None) -> dict[str, object]:
        """Serialize explicit membership and source locators as a reusable document."""
        return {
            **self._content(),
            "snapshot_id": self.snapshot_id,
            "sources": [s.to_dict(base=base) for s in self.sources],
        }

    @classmethod
    def from_dict(cls, value: object, *, base: Path | None = None) -> SelectionSnapshot:
        """Check the complete snapshot before exposing saved membership."""
        data = object_fields(
            value,
            {
                "schema",
                "snapshot_id",
                "sources",
                "request",
                "members",
                "counts",
                "algorithm",
            },
            "selection snapshot",
        )
        if set(data) != {
            "schema",
            "snapshot_id",
            "sources",
            "request",
            "members",
            "counts",
            "algorithm",
        }:
            msg = "selection snapshot is missing required fields"
            raise ValueError(msg)
        if data.pop("schema", None) != SNAPSHOT_SCHEMA:
            msg = "unsupported selection-snapshot schema"
            raise ValueError(msg)
        identity = data.pop("snapshot_id", None)
        result = cls(
            sources=tuple(
                SelectionSource.from_dict(s, base=base) for s in data.pop("sources")
            ),
            request=LibrarySelection.from_dict(data.pop("request")),
            members=tuple(
                SelectedDesign(
                    **object_fields(
                        m, {"reference", "record_digest"}, "selected design"
                    )
                )
                for m in data.pop("members")
            ),
            **data,
        )
        if identity != result.snapshot_id:
            msg = "selection snapshot digest mismatch"
            raise ValueError(msg)
        return result

    def __repr__(self) -> str:
        """Show useful allocation totals without expanding membership or locators."""
        return (
            f"SelectionSnapshot(selected={self.selected}, requested={self.requested}, "
            f"available={self.available}, status={self.status}, "
            f"sources={len(self.sources)})"
        )


def _validate_counts(snapshot: SelectionSnapshot) -> None:
    """Validate saved allocation evidence without regenerating a random sample."""
    cells = {c for s in snapshot.sources for c in s.cells}
    _validate_allocation(snapshot.request, snapshot.counts, cells)
    observed = Counter(m.reference.rsplit("/", 1)[0] for m in snapshot.members)
    if (
        not observed.keys() <= cells
        or len(set(snapshot.references())) != snapshot.selected
    ):
        msg = "selection contains unknown cells or duplicate full design references"
        raise ValueError(msg)
    per_cell = (
        snapshot.request.take is not None and snapshot.request.take.per_cell is not None
    )
    for key, value in snapshot.counts.items():
        expected = observed[key] if per_cell else snapshot.selected
        if value["selected"] != expected:
            msg = "selection counts do not reconcile with saved membership"
            raise ValueError(msg)
    if sum(c["selected"] for c in snapshot.counts.values()) != snapshot.selected:
        msg = "selection counts omit selected cells"
        raise ValueError(msg)


def validate_summary(value: object, *, cells: set[str]) -> dict[str, object]:
    """Check allocation metadata retained by exports without requiring source paths."""
    data = object_fields(
        value,
        {
            "schema",
            "snapshot_id",
            "status",
            "requested",
            "available",
            "selected",
            "shortfall",
            "counts",
            "request",
            "algorithm",
        },
        "selection summary",
    )
    if data.get("schema") != "dense_arrays.selection_summary.v1":
        msg = "unsupported selection-summary schema"
        raise ValueError(msg)
    digest(data.get("snapshot_id"), field_name="snapshot_id")
    request = LibrarySelection.from_dict(data["request"])
    _validate_allocation(request, data["counts"], cells)
    expected = (
        "sha256_priority.v1"
        if request.take and request.take.policy == "random"
        else "source_order.v1"
    )
    if data["algorithm"] != expected:
        msg = "selection summary algorithm does not match request"
        raise ValueError(msg)
    for name in ("requested", "available", "selected", "shortfall"):
        integer(data.get(name), field_name=name, minimum=0)
        if data[name] != sum(c[name] for c in data["counts"].values()):
            msg = "selection summary totals do not reconcile"
            raise ValueError(msg)
    if (
        data["status"] != ("partial" if data["shortfall"] else "complete")
        or data["requested"] - data["selected"] != data["shortfall"]
    ):
        msg = "selection summary status or shortfall mismatch"
        raise ValueError(msg)
    return data


def _validate_allocation(
    request: LibrarySelection,
    counts: object,
    cells: set[str],
) -> None:
    """Use one allocation policy for full snapshots and portable summaries."""
    expected_quotas = resolve_quotas(request.take, cells)
    counts = object_fields(counts, set(expected_quotas), "allocation counts")
    if set(counts) != set(expected_quotas):
        msg = "selection counts omit declared allocation quotas"
        raise ValueError(msg)
    for key, value in counts.items():
        row = object_fields(
            value,
            {"requested", "available", "selected", "shortfall"},
            "selection counts",
        )
        for name in ("requested", "available", "selected", "shortfall"):
            integer(row.get(name), field_name=name, minimum=0)
        expected = expected_quotas[key]
        requested = row["available"] if expected is None else expected
        if row["requested"] != requested:
            msg = "selection counts disagree with the requested allocation"
            raise ValueError(msg)
        if (
            row["selected"] != min(requested, row["available"])
            or row["shortfall"] != requested - row["selected"]
        ):
            msg = "selection counts do not reconcile with the declared allocation"
            raise ValueError(msg)
    if any(c["shortfall"] for c in counts.values()) and (
        request.take is None or request.take.shortfall != "allow_partial"
    ):
        raise SelectionShortfall(counts)
