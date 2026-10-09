"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/workflow/inputs.py

Strict declared-schema file dispatch at the application boundary.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

from dense_arrays.artifacts.documents import read_document
from dense_arrays.artifacts.pools import decode_preparation
from dense_arrays.artifacts.run_plans import check_plan_size, decode_plan
from dense_arrays.parts import PartFilter
from dense_arrays.planning import (
    DesignSpec,
    ExtensionSpec,
    GenerationPlan,
    PlanEvidence,
    PreparationPlan,
)
from dense_arrays.planning.evidence import EVIDENCE_SCHEMA
from dense_arrays.planning.extension import EXTENSION_SCHEMA
from dense_arrays.planning.matrices import MatrixPlan, MatrixSpec
from dense_arrays.planning.matrices.requests import MATRIX_SCHEMA
from dense_arrays.planning.matrices.resolution import MATRIX_PLAN_SCHEMA
from dense_arrays.planning.preparation import (
    PREPARATION_PLAN_SCHEMA,
    PREPARE_SCHEMA,
    PREPARE_SET_SCHEMA,
    PREPARE_WINDOWS_SCHEMA,
    SAMPLED_PLAN_SCHEMA,
    SET_PLAN_SCHEMA,
    preparation_from_dict,
)
from dense_arrays.planning.serialization import PLAN_SCHEMA, request_from_dict
from dense_arrays.reporting.design_filters import DesignFilter
from dense_arrays.reporting.filters import AttemptFilter
from dense_arrays.reporting.plans.filters import PLAN_FILTER_SCHEMA, PlanFilter
from dense_arrays.reporting.pools.filters import (
    CANDIDATE_FILTER_SCHEMA,
    CandidateFilter,
)
from dense_arrays.reporting.pools.snapshots import (
    POOL_QUALITY_SCHEMA,
    PoolQualitySnapshot,
)
from dense_arrays.reporting.quality import QualitySnapshot
from dense_arrays.reporting.selections import LibrarySelection, SelectionSnapshot
from dense_arrays.reporting.selections.requests import SELECTION_SCHEMA
from dense_arrays.reporting.selections.snapshots import SNAPSHOT_SCHEMA

if TYPE_CHECKING:
    from pathlib import Path

    from dense_arrays.artifacts.reading import ReadLimits
    from dense_arrays.parts import PreparationSet, PreparationSpec


def read_selection(
    path: Path,
) -> (
    PartFilter
    | CandidateFilter
    | AttemptFilter
    | DesignFilter
    | PlanFilter
    | LibrarySelection
    | SelectionSnapshot
):
    """Read a declared predicate without guessing its meaning from the view."""
    value = read_document(path)
    if isinstance(value, dict) and value.get("schema") == SNAPSHOT_SCHEMA:
        return SelectionSnapshot.from_dict(value, base=path.absolute().parent)
    parsers = {
        CANDIDATE_FILTER_SCHEMA: CandidateFilter,
        PLAN_FILTER_SCHEMA: PlanFilter,
        SELECTION_SCHEMA: LibrarySelection,
        "dense_arrays.design-filter.v1": DesignFilter,
        "dense_arrays.attempt-filter.v1": AttemptFilter,
    }
    cls = (
        parsers.get(value.get("schema"), PartFilter)
        if isinstance(value, dict)
        else PartFilter
    )
    return cls.from_dict(value)


def read_source(
    path: Path,
    *,
    max_identities: int | None = None,
) -> (
    DesignSpec
    | ExtensionSpec
    | GenerationPlan
    | PlanEvidence
    | PreparationSpec
    | PreparationSet
    | PreparationPlan
    | MatrixPlan
    | MatrixSpec
):
    """Read YAML or JSON by declared schema, never by a guessed filename."""
    value = read_document(path)
    if isinstance(value, dict) and value.get("schema") == MATRIX_SCHEMA:
        return MatrixSpec.from_dict(value, base=path.absolute().parent)
    if isinstance(value, dict) and value.get("schema") == EXTENSION_SCHEMA:
        return ExtensionSpec.from_dict(value, base=path.absolute().parent)
    plans = {
        PLAN_SCHEMA,
        MATRIX_PLAN_SCHEMA,
        PREPARATION_PLAN_SCHEMA,
        SAMPLED_PLAN_SCHEMA,
        SET_PLAN_SCHEMA,
    }
    if isinstance(value, dict) and value.get("schema") in plans:
        return (
            decode_plan(value, max_identities, base=path.absolute().parent)
            if value["schema"] in {PLAN_SCHEMA, MATRIX_PLAN_SCHEMA}
            else decode_preparation(
                value, max_identities=max_identities, base=path.absolute().parent
            )
        )
    if isinstance(value, dict) and value.get("schema") == EVIDENCE_SCHEMA:
        check_plan_size(value.get("content", {}), max_identities)
        return PlanEvidence.from_dict(value)
    if isinstance(value, dict) and value.get("schema") in {
        PREPARE_SCHEMA,
        PREPARE_SET_SCHEMA,
        PREPARE_WINDOWS_SCHEMA,
    }:
        return preparation_from_dict(value, base=path.absolute().parent)
    return request_from_dict(value, base=path.absolute().parent)


def read_quality(
    path: Path, limits: ReadLimits | None = None
) -> QualitySnapshot | PoolQualitySnapshot:
    """Load one saved quality document under explicit structural read bounds."""
    value = read_document(path)
    cls = (
        PoolQualitySnapshot
        if isinstance(value, dict) and value.get("schema") == POOL_QUALITY_SCHEMA
        else QualitySnapshot
    )
    return cls.from_dict(value, read_limits=limits)
