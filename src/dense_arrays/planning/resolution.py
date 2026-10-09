"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/resolution.py

Resolve curated requests into immutable, input-bound generation plans.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import hashlib
import os
from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING

from dense_arrays._record_validation import canonical_json, digest, object_fields
from dense_arrays.artifacts.publication import write_new
from dense_arrays.parts import BoundParts, Part, PartTable, PoolSource
from dense_arrays.parts.ingestion import read_parts
from dense_arrays.parts.provenance import ImportReport
from dense_arrays.planning.batches.bindings import maximum_offered_parts
from dense_arrays.planning.evidence import PlanEvidence
from dense_arrays.planning.libraries import ParentLibrary
from dense_arrays.planning.models import DesignSpec
from dense_arrays.planning.serialization import (
    PLAN_SCHEMA,
    policies_for,
    request_from_dict,
    request_to_dict,
)

if TYPE_CHECKING:
    from dense_arrays.planning.libraries import ExcludedDesign


@dataclass(frozen=True)
class InputBinding:
    """A source location outside semantic identity and its exact byte digest."""

    path: Path
    sha256: str

    def __post_init__(self) -> None:
        """Detach the path and validate the source fingerprint."""
        object.__setattr__(self, "path", Path(self.path).absolute())
        digest(self.sha256, field_name="input.sha256")

    def verify(self) -> None:
        """Check source bytes against the frozen fingerprint."""
        if hashlib.sha256(self.path.read_bytes()).hexdigest() != self.sha256:
            msg = f"input changed since planning: {self.path}; create a new plan"
            raise ValueError(msg)

    def to_dict(self, *, base: Path | None = None) -> dict[str, str]:
        """Encode a locator independently of the byte identity it binds."""
        return {
            "path": str(self.path)
            if base is None
            else os.path.relpath(self.path, base),
            "sha256": self.sha256,
        }

    @classmethod
    def from_dict(cls, value: object, *, base: Path | None = None) -> InputBinding:
        """Resolve relative input locators against their containing document."""
        data = object_fields(value, {"path", "sha256"}, "input binding")
        if base is not None:
            data["path"] = base / data["path"]
        return cls(**data)


@dataclass(frozen=True, repr=False)
class GenerationPlan:
    """An immutable normalized request with frozen v1 policies and source bindings."""

    request: DesignSpec
    inputs: tuple[InputBinding, ...] = ()
    import_report: ImportReport | None = None
    parent: ParentLibrary | None = None
    embedded_input_digests: tuple[str, ...] = ()
    plan_id: str = field(init=False)
    evidence: PlanEvidence = field(init=False, repr=False)
    _resume_supported: bool = field(
        default=True, repr=False, compare=False, kw_only=True
    )

    def __post_init__(self) -> None:
        """Enforce supported semantics before a runtime can own any destination."""
        if not isinstance(self._resume_supported, bool):
            msg = "preview.resume_supported must be a boolean"
            raise TypeError(msg)
        if not isinstance(self.inputs, (list, tuple)) or any(
            not isinstance(item, InputBinding) for item in self.inputs
        ):
            msg = "inputs must contain InputBinding records"
            raise TypeError(msg)
        object.__setattr__(self, "inputs", tuple(self.inputs))
        if not isinstance(self.embedded_input_digests, (list, tuple)):
            msg = "embedded_input_digests must be an ordered array"
            raise TypeError(msg)
        if self.inputs and self.embedded_input_digests:
            msg = "generation plans cannot mix bound files and embedded origins"
            raise ValueError(msg)
        object.__setattr__(
            self, "embedded_input_digests", tuple(self.embedded_input_digests)
        )
        evidence = PlanEvidence(
            self.request,
            self.embedded_input_digests or tuple(i.sha256 for i in self.inputs),
            self.import_report,
            self.parent,
        )
        object.__setattr__(self, "import_report", evidence.import_report)
        object.__setattr__(self, "plan_id", evidence.plan_id)
        object.__setattr__(self, "evidence", evidence)

    def __repr__(self) -> str:
        """Summarize immutable plan identity without serializing the part collection."""
        return (
            f"GenerationPlan({self.plan_id[:12]}, {len(self.request.parts)} parts, "
            f"target={self.request.target.count})"
        )

    @property
    def collection_id(self) -> str:
        """Identify supplied parts independently of target, seed and run lineage."""
        return self.evidence.collection_id

    def admit_work(self) -> None:
        """Bound the largest active model without changing saved evidence."""
        if self.request.target.count:
            nodes = maximum_offered_parts(self.request) * (
                1 if self.request.strands == "single" else 2
            )
            self.request.limits.admit_model(nodes)

    @property
    def exclusions(self) -> tuple[ExcludedDesign, ...]:
        """Frozen sequence exclusions, independent of additional-target semantics."""
        return self.evidence.exclusions

    @property
    def preview(self) -> Mapping[str, object]:
        """Cheap model dimensions; no adjacency matrix or backend allocation."""
        count = len(self.request.parts)
        schedule = self.request.schedule
        offered = maximum_offered_parts(self.request)
        resampling = self.request.resampling
        nodes = offered * (1 if self.request.strands == "single" else 2)
        return MappingProxyType(
            {
                "parts": count,
                **(
                    {"offered_parts": offered, "batch_id": self.request.batch.batch_id}
                    if self.request.batch is not None
                    else {}
                ),
                **(
                    {
                        "batches": len(schedule.batches),
                        "max_offered_parts": offered,
                        "attempts_per_batch": schedule.attempts_per_batch,
                        **(
                            {"accepted_per_batch": schedule.accepted_per_batch}
                            if schedule.accepted_per_batch is not None
                            else {}
                        ),
                    }
                    if schedule is not None
                    else {}
                ),
                **(
                    {
                        "offered_parts": offered,
                        "max_batches": resampling.max_batches,
                        "attempts_per_batch": resampling.attempts_per_batch,
                        "accepted_per_batch": resampling.accepted_per_batch,
                        "feedback": None
                        if resampling.feedback is None
                        else resampling.feedback.to_dict(),
                    }
                    if resampling is not None
                    else {}
                ),
                "cells": 1,
                "target": self.request.target.count,
                "oriented_nodes": nodes,
                "path_variables": 0
                if self.request.search == "greedy"
                else nodes * (nodes - 1) + 2 * nodes,
                "feasibility": "unknown",
                "proof_scope": None
                if self.request.search == "greedy"
                else "offered_packing_model",
                **({"search": "greedy"} if self.request.search == "greedy" else {}),
                **(
                    {"packing_preference": self.request.packing_preference}
                    if self.request.packing_preference
                    else {}
                ),
                "solver": None if self.request.search == "greedy" else "CBC",
                "resume_supported": self._resume_supported,
                **(
                    {"input_verification": "embedded_parts"}
                    if self.embedded_input_digests
                    else {}
                ),
                **(
                    {
                        "excluded_sequences": len(self.exclusions),
                        "exclusion_scope": {
                            "cell_mapping": dict(self.request.exclude.cell_mapping),
                            "uniqueness": self.request.exclude.uniqueness,
                        },
                    }
                    if self.request.exclude is not None
                    else {}
                ),
                **(
                    {
                        "excluded_sequences": len(self.parent.exclusions),
                        "parent_run_id": self.parent.run_id,
                        "target_scope": "additional",
                    }
                    if self.parent is not None
                    else {}
                ),
            }
        )

    def verify_inputs(self) -> None:
        """Check current bytes; execution uses the already-bound immutable parts."""
        if self.parent is not None:
            # Extension inputs were verified from native parent evidence at planning.
            return
        for item in self.inputs:
            item.verify()

    def to_dict(self, *, base: Path | None = None) -> dict[str, object]:
        """Return a portable logical plan with explicit policy/default versions."""
        return {
            "schema": PLAN_SCHEMA,
            "plan_id": self.plan_id,
            "request": request_to_dict(self.request),
            "policies": policies_for(self.request),
            "inputs": [i.to_dict(base=base) for i in self.inputs],
            "preview": dict(self.preview),
            "import_report": self.import_report.to_dict(),
            **(
                {"embedded_input_digests": list(self.embedded_input_digests)}
                if self.embedded_input_digests
                else {}
            ),
            **({"parent": self.parent.to_dict()} if self.parent is not None else {}),
        }

    def write(self, path: str | Path) -> None:
        """Save a complete plan with relative locators; never replace a destination."""
        path = Path(path).absolute()
        write_new(path, canonical_json(self.to_dict(base=path.parent)) + "\n")

    @classmethod
    def from_dict(cls, value: object, *, base: Path | None = None) -> GenerationPlan:
        """Validate a saved plan without reapplying omitted defaults."""
        keys = {
            "schema",
            "plan_id",
            "request",
            "policies",
            "inputs",
            "preview",
            "import_report",
        }
        value = object_fields(
            value, keys | {"parent", "embedded_input_digests"}, "generation_plan"
        )
        if not keys <= set(value) or value.get("schema") != PLAN_SCHEMA:
            msg = "unsupported or incomplete generation plan schema"
            raise ValueError(msg)
        request = request_from_dict(value["request"])
        if value["policies"] != policies_for(request):
            msg = "unsupported generation plan policy version"
            raise ValueError(msg)
        if request_to_dict(request) != value["request"]:
            msg = (
                "saved plans require every normalized field; defaults are not reapplied"
            )
            raise ValueError(msg)
        if not isinstance(value["inputs"], list):
            msg = "plan inputs must be an array"
            raise TypeError(msg)
        if not isinstance(value["preview"], Mapping):
            msg = "generation plan preview must be an object"
            raise TypeError(msg)
        plan = cls(
            request,
            tuple(InputBinding.from_dict(item, base=base) for item in value["inputs"]),
            ImportReport.from_dict(value["import_report"]),
            ParentLibrary.from_dict(value["parent"]) if "parent" in value else None,
            value.get("embedded_input_digests", ()),
            _resume_supported=value["preview"].get("resume_supported"),
        )
        if plan.plan_id != value["plan_id"]:
            msg = "generation plan digest mismatch"
            raise ValueError(msg)
        if dict(plan.preview) != value["preview"]:
            msg = "generation plan preview does not match its request"
            raise ValueError(msg)
        return plan


def resolve_design(request: DesignSpec) -> GenerationPlan:
    """Read curated inputs once, resolve supported requirements and bind provenance."""
    if not isinstance(request, DesignSpec):
        msg = "plan requires a DesignSpec"
        raise TypeError(msg)
    plan = bind_design(request, resolve_source(request.parts))
    if isinstance(request.parts, BoundParts):
        plan.verify_inputs()
    return plan


def bind_design(request: DesignSpec, source: BoundParts) -> GenerationPlan:
    """Apply frozen source evidence to a recipe without opening its original files."""
    bindings = (
        ()
        if source.locations is None
        else tuple(
            InputBinding(path, fingerprint)
            for path, fingerprint in zip(
                source.locations, source.input_digests, strict=True
            )
        )
    )
    return GenerationPlan(
        replace(request, parts=source.parts),
        bindings,
        source.import_report,
        embedded_input_digests=source.input_digests if source.locations is None else (),
    )


def resolve_source(
    source: tuple[Part, ...] | PartTable | PoolSource | BoundParts,
) -> BoundParts:
    """Import one part source through the same boundary for single and matrix plans."""
    if isinstance(source, BoundParts):
        return source
    if isinstance(source, PoolSource):
        from dense_arrays.artifacts.pools import (  # noqa: PLC0415
            POOL_DATABASE,
            read_pool_source,
        )

        imported = read_pool_source(source)
        binding_path = source.path / POOL_DATABASE
    else:
        imported = read_parts(source)
        binding_path = source.table if isinstance(source, PartTable) else None
    return BoundParts(
        imported.parts,
        imported.report,
        () if binding_path is None else (imported.source_digest,),
        None if binding_path is None else (binding_path,),
    )
