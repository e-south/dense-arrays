"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/matrices/requests.py

Typed substitutions over declared part and requirement identities.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass, field, fields, replace
from types import MappingProxyType
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    integer,
    mutable_json,
    object_fields,
    records,
)
from dense_arrays.parts import BoundParts, Part, PartTable, PoolSource
from dense_arrays.parts.ingestion import validate_parts
from dense_arrays.planning.batches import BatchSchedule, CandidateBatch, Resampling
from dense_arrays.planning.batches.schedules import selection_from_dict
from dense_arrays.planning.libraries import LibraryExclusion, ParentLibrary
from dense_arrays.planning.models import DesignSpec
from dense_arrays.planning.requirements import REQUIREMENT_TYPES, Requirement
from dense_arrays.planning.serialization import (
    parts_from_dict,
    parts_to_dict,
    request_from_dict,
    request_to_dict,
    requirement_from_dict,
    requirement_to_dict,
)

from .allocation import Allocation

if TYPE_CHECKING:
    from pathlib import Path

MATRIX_SCHEMA = "dense_arrays.matrix.v1"


def label(value: object) -> None:
    """Keep cell labels readable and unambiguous in full references."""
    if (
        not isinstance(value, str)
        or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value) is None
    ):
        msg = (
            "matrix labels require letters/digits followed by letters/digits, _, . or -"
        )
        raise ValueError(msg)


@dataclass(frozen=True)
class Variant:
    """Replace named parts/rules or explicitly add rules; empty retains the base."""

    parts: tuple[Part, ...] = ()
    requirements: tuple[Requirement, ...] = ()
    add_requirements: tuple[Requirement, ...] = field(default=(), kw_only=True)

    def __post_init__(self) -> None:
        """Freeze complete substitutions and reject repeated identities."""
        object.__setattr__(
            self, "parts", records(self.parts, Part, field_name="variant.parts")
        )
        for name in ("requirements", "add_requirements"):
            values = getattr(self, name)
            if not isinstance(values, (tuple, list)) or any(
                not isinstance(r, REQUIREMENT_TYPES) for r in values
            ):
                msg = f"variant.{name} must contain typed requirements"
                raise TypeError(msg)
            object.__setattr__(self, name, tuple(values))
        rules = self.requirements + self.add_requirements
        if len({p.part_id for p in self.parts}) != len(self.parts) or len(
            {r.id for r in rules}
        ) != len(rules):
            msg = "variant contains repeated part or requirement identities"
            raise ValueError(msg)

    def to_dict(self) -> dict[str, object]:
        """Serialize complete part records and typed requirement variants."""
        return {
            "parts": [
                {f.name: mutable_json(getattr(part, f.name)) for f in fields(Part)}
                for part in self.parts
            ],
            "requirements": [requirement_to_dict(r) for r in self.requirements],
            **(
                {
                    "add_requirements": [
                        requirement_to_dict(r) for r in self.add_requirements
                    ]
                }
                if self.add_requirements
                else {}
            ),
        }

    @classmethod
    def from_dict(cls, value: object) -> Variant:
        """Parse a compact substitution without reading or changing the base."""
        data = object_fields(
            value, {"parts", "requirements", "add_requirements"}, "variant"
        )
        return cls(
            tuple(Part(**p) for p in data.get("parts", [])),
            tuple(requirement_from_dict(r) for r in data.get("requirements", [])),
            add_requirements=tuple(
                requirement_from_dict(r) for r in data.get("add_requirements", [])
            ),
        )


@dataclass(frozen=True, repr=False)
class MatrixSpec:
    """A base recipe, named alternatives, bounded pairing and explicit allocation."""

    base: DesignSpec
    axes: Mapping[str, Mapping[str, Variant]]
    allocation: Allocation
    max_cells: int
    pairing: str = "cross_product"
    pairs: tuple[Mapping[str, str], ...] = ()
    exclude: Mapping[str, LibraryExclusion] = field(default_factory=dict)
    batches: Mapping[str, CandidateBatch | BatchSchedule | Resampling] = field(
        default_factory=dict
    )
    sources: Mapping[str, PartTable | PoolSource | BoundParts | tuple[Part, ...]] = (
        field(default_factory=dict)
    )

    def __post_init__(self) -> None:
        """Freeze declared order and reject implicit pairing or target ownership."""
        if not isinstance(self.base, DesignSpec) or not isinstance(
            self.allocation, Allocation
        ):
            msg = "matrix requires a DesignSpec base and Allocation"
            raise TypeError(msg)
        if self.base.target.count != 1:
            msg = "matrix allocation owns targets; omit the base target"
            raise ValueError(msg)
        object.__setattr__(self, "exclude", _freeze_exclusions(self.exclude))
        if self.base.batch is not None or self.base.schedule is not None:
            msg = "matrix batches require explicit per-cell declarations"
            raise ValueError(msg)
        if not isinstance(self.batches, Mapping) or any(
            not isinstance(k, str)
            or not isinstance(v, (CandidateBatch, BatchSchedule, Resampling))
            for k, v in self.batches.items()
        ):
            msg = (
                "matrix batches must map cells to "
                "CandidateBatch, BatchSchedule or Resampling"
            )
            raise TypeError(msg)
        object.__setattr__(self, "batches", MappingProxyType(dict(self.batches)))
        if not isinstance(self.sources, Mapping) or any(
            not isinstance(key, str) for key in self.sources
        ):
            msg = "matrix sources must map cell IDs to part sources"
            raise TypeError(msg)
        object.__setattr__(
            self,
            "sources",
            MappingProxyType(
                {
                    key: source
                    if isinstance(source, (PartTable, PoolSource, BoundParts))
                    else validate_parts(source)
                    for key, source in self.sources.items()
                }
            ),
        )
        integer(self.max_cells, field_name="max_cells", minimum=1)
        if self.pairing not in {"cross_product", "zip", "explicit"}:
            msg = "matrix pairing must be cross_product, zip or explicit"
            raise ValueError(msg)
        object.__setattr__(self, "axes", _freeze_axes(self.axes))
        if not isinstance(self.pairs, (tuple, list)) or any(
            not isinstance(pair, Mapping) for pair in self.pairs
        ):
            msg = "matrix pairs must be an array of choice mappings"
            raise TypeError(msg)
        object.__setattr__(
            self, "pairs", tuple(MappingProxyType(dict(pair)) for pair in self.pairs)
        )
        if (self.pairing == "explicit") != bool(self.pairs):
            msg = "explicit pairing requires pairs; other pairing modes prohibit them"
            raise ValueError(msg)

    def with_changes(self, **changes: object) -> MatrixSpec:
        """Return a validated immutable edit without reading sources."""
        return replace(self, **changes)

    def __repr__(self) -> str:
        """Show shape without expanding the variants or source collection."""
        return (
            f"MatrixSpec(axes={len(self.axes)}, pairing={self.pairing!r}, "
            f"max_cells={self.max_cells})"
        )

    def to_dict(self, *, base: Path | None = None) -> dict[str, object]:
        """Use ordered arrays so canonical JSON cannot change expansion order."""
        return {
            "schema": MATRIX_SCHEMA,
            "base": request_to_dict(self.base, base=base),
            "axes": [
                {
                    "name": name,
                    "choices": [
                        {"name": choice, "variant": variant.to_dict()}
                        for choice, variant in options.items()
                    ],
                }
                for name, options in self.axes.items()
            ],
            "allocation": self.allocation.to_dict(),
            "max_cells": self.max_cells,
            "pairing": self.pairing,
            "pairs": [dict(pair) for pair in self.pairs],
            **(
                {
                    "sources": {
                        cell: parts_to_dict(source, base=base)
                        for cell, source in self.sources.items()
                    }
                }
                if self.sources
                else {}
            ),
            **(
                {
                    "batches": {
                        cell: batch.to_dict() for cell, batch in self.batches.items()
                    }
                }
                if self.batches
                else {}
            ),
            **(
                {
                    "exclude": {
                        cell: scope.to_dict() for cell, scope in self.exclude.items()
                    }
                }
                if self.exclude
                else {}
            ),
        }

    @classmethod
    def from_dict(cls, value: object, *, base: Path | None = None) -> MatrixSpec:
        """Read compact named axes or their explicitly ordered native encoding."""
        data = object_fields(
            value,
            {
                "schema",
                "base",
                "axes",
                "allocation",
                "max_cells",
                "pairing",
                "pairs",
                "exclude",
                "batches",
                "sources",
            },
            "matrix",
        )
        if data.pop("schema", None) != MATRIX_SCHEMA:
            msg = "unsupported matrix schema"
            raise ValueError(msg)
        missing = {"base", "axes", "allocation", "max_cells"} - set(data)
        if missing:
            msg = f"matrix is missing required fields: {', '.join(sorted(missing))}"
            raise ValueError(msg)
        data["base"] = request_from_dict(data["base"], base=base)
        data["allocation"] = Allocation(**data["allocation"])
        data["axes"] = _axes_from_dict(data["axes"])
        if "sources" in data:
            if not isinstance(data["sources"], Mapping):
                msg = "matrix sources must map cell IDs to part sources"
                raise TypeError(msg)
            data["sources"] = {
                cell: parts_from_dict(source, base=base)
                for cell, source in data["sources"].items()
            }
        if "batches" in data:
            if not isinstance(data["batches"], Mapping):
                msg = "matrix batches must be a cell mapping"
                raise TypeError(msg)
            data["batches"] = {
                cell: selection_from_dict(batch)
                for cell, batch in data["batches"].items()
            }
        if "exclude" in data:
            if not isinstance(data["exclude"], Mapping):
                msg = "matrix exclude must map cell identities to exclusions"
                raise TypeError(msg)
            data["exclude"] = {
                cell: LibraryExclusion.from_dict(scope, base=base)
                for cell, scope in data["exclude"].items()
            }
        return cls(**data)


def _freeze_axes(value: object) -> Mapping[str, Mapping[str, Variant]]:
    if not isinstance(value, Mapping) or not value:
        msg = "matrix axes must be a nonempty mapping"
        raise TypeError(msg)
    axes = {}
    for name, options in value.items():
        label(name)
        if not isinstance(options, Mapping) or not options:
            msg = f"axis {name!r} requires named variants"
            raise TypeError(msg)
        for choice, variant in options.items():
            label(choice)
            if not isinstance(variant, Variant):
                msg = f"axis {name!r} choice {choice!r} requires Variant"
                raise TypeError(msg)
        axes[name] = MappingProxyType(dict(options))
    return MappingProxyType(axes)


def _named_records(value: object, field: str) -> dict[str, object]:
    result = {}
    for record in records(value, Mapping, field_name=field):
        data = object_fields(record, {"name", field}, f"matrix {field}")
        if {"name", field} - set(data):
            msg = f"matrix {field} record requires name and {field}"
            raise ValueError(msg)
        label(data.get("name"))
        if data["name"] in result:
            msg = f"duplicate matrix {field} name {data['name']!r}"
            raise ValueError(msg)
        result[data["name"]] = data[field]
    return result


def _axes_from_dict(value: object) -> dict[str, dict[str, Variant]]:
    axes = value if isinstance(value, Mapping) else _named_records(value, "choices")
    result = {}
    for name, options in axes.items():
        choices = (
            options
            if isinstance(options, Mapping)
            else _named_records(options, "variant")
        )
        result[name] = {choice: Variant.from_dict(v) for choice, v in choices.items()}
    return result


def _freeze_exclusions(value: object) -> Mapping[str, LibraryExclusion]:
    if not isinstance(value, Mapping):
        msg = "matrix exclude must map cell identities to exclusions"
        raise TypeError(msg)
    for cell, scope in value.items():
        if not isinstance(scope, LibraryExclusion) or not isinstance(
            scope.source, ParentLibrary
        ):
            msg = "matrix exclusions require frozen per-cell accepted libraries"
            raise TypeError(msg)
        if dict(scope.cell_mapping) != {cell: cell}:
            msg = "matrix exclusion mapping must match its declared cell"
            raise ValueError(msg)
    return MappingProxyType(dict(value))
