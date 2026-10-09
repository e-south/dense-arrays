"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/arrays/exporting.py

Publish supplied-array collections, scalar tables and bounded metadata documents.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import csv
from contextlib import nullcontext
from dataclasses import dataclass, fields
from pathlib import Path
from typing import TYPE_CHECKING, TextIO
from urllib.parse import quote

from dense_arrays._record_validation import canonical_json
from dense_arrays.artifacts.publication import new_text_file
from dense_arrays.artifacts.reading import ReadBudget, ReadCost, ReadLimits
from dense_arrays.artifacts.receipts import ExportReceipt

from .models import ArrayCollection, CollectionSummary
from .projections import ArrayPlacement, ArraySequence
from .publication import publish_collection
from .reading import CollectionView
from .serialization import read_input, write_input
from .storage import file_digest, is_collection, read_parts, reader

if TYPE_CHECKING:
    from collections.abc import Iterator

    from dense_arrays.realized import RealizedArray


@dataclass(frozen=True)
class SuppliedExport:
    """A once-consumed supplied-array source and its explicit work limits."""

    source: ArrayCollection
    read_limits: ReadLimits
    view: str = "arrays"

    @property
    def cost(self) -> ReadCost:
        """Expose an unknown stream length instead of consuming it for a preview."""
        return ReadCost("supplied-arrays", 0, "scan", "arrays", None, self.read_limits)


def handles_export(artifact: object) -> bool:
    """Recognize collection objects and explicitly named line-delimited inputs."""
    return (
        is_collection(artifact)
        or isinstance(artifact, (CollectionView, CollectionSummary))
        or (isinstance(artifact, (str, Path)) and Path(artifact).suffix == ".jsonl")
    )


def resolve_export(  # noqa: PLR0913 - paired Python/CLI query options
    artifact: object,
    *,
    view: str | None,
    format_name: str,
    all_rows: bool,
    selected: object,
    read_limits: ReadLimits | None,
    limit: int | None,
    compare: object,
) -> SuppliedExport | CollectionView | CollectionSummary:
    """Fail unsupported combinations before scanning or creating a destination."""
    from dense_arrays.workflow.operations import inspect  # noqa: PLC0415

    limits = _export_limits(read_limits, compare, limit)
    view = view or (
        artifact.view
        if isinstance(artifact, CollectionView)
        else "summary"
        if isinstance(artifact, CollectionSummary)
        else "arrays"
    )
    if view == "summary":
        if format_name != "json" or all_rows is not False or selected is not None:
            msg = (
                "collection summary exports require JSON "
                "without record-selection options"
            )
            raise ValueError(msg)
        return (
            artifact
            if isinstance(artifact, CollectionSummary)
            else inspect(artifact, read_limits=limits)
        )
    if all_rows is not True:
        msg = "record export requires all=True to declare its complete scope"
        raise ValueError(msg)
    supported = {
        "arrays": {"json", "jsonl", "bundle"},
        "parts": {"json"},
        "sequences": {"json", "csv", "tsv", "fasta"},
        "placements": {"json", "csv", "tsv"},
    }
    if format_name not in supported.get(view, set()):
        msg = "unsupported array collection view/format combination"
        raise ValueError(msg)
    if isinstance(artifact, (str, Path)) and Path(artifact).suffix == ".jsonl":
        artifact = read_input(Path(artifact), limits)
    if isinstance(artifact, ArrayCollection):
        if (
            view != "arrays"
            or format_name not in {"bundle", "jsonl"}
            or selected is not None
        ):
            msg = "publish supplied arrays as bundle or jsonl before querying them"
            raise ValueError(msg)
        return SuppliedExport(artifact, limits)
    if isinstance(artifact, CollectionView):
        if (
            artifact.view != view
            or artifact.limit is not None
            or artifact.after is not None
            or selected is not None
            or read_limits is not None
        ):
            msg = (
                "export requires a complete bound collection query "
                "without additional query options"
            )
            raise ValueError(msg)
        return artifact
    return inspect(artifact, view=view, select=selected, all=True, read_limits=limits)


def _source(query: SuppliedExport | CollectionView) -> ArrayCollection:
    if isinstance(query, SuppliedExport):
        return query.source
    with reader(query.path) as connection:
        parts = read_parts(connection, query.summary, ReadBudget(query.read_limits))

    def arrays() -> Iterator[RealizedArray]:
        with query.records() as rows:
            for row in rows:
                yield row.realized

    return ArrayCollection(
        tuple(parts.values()), arrays(), query.summary.manifest["provenance"]
    )


def publish_export(
    query: SuppliedExport | CollectionView | CollectionSummary,
    *,
    format_name: str,
    out: str | Path | TextIO,
) -> ExportReceipt:
    """Use complete-file publication; interrupted streams never return a receipt."""
    if format_name == "bundle":
        if not isinstance(out, (str, Path)):
            msg = "collection bundle export requires a directory path"
            raise TypeError(msg)
        return publish_collection(
            _source(query), Path(out).absolute(), query.read_limits
        )
    path = Path(out).absolute() if isinstance(out, (str, Path)) else None
    if path is None and not callable(getattr(out, "write", None)):
        msg = "out must be a path or writable text stream"
        raise TypeError(msg)
    with new_text_file(path) if path is not None else nullcontext(out) as stream:
        if isinstance(query, CollectionSummary):
            stream.write(canonical_json(query.to_dict()) + "\n")
            count, view = 1, "summary"
        elif format_name == "jsonl":
            count, view = (
                write_input(_source(query), stream, query.read_limits),
                "arrays",
            )
        else:
            count, view = _write_records(query, stream, format_name), query.view
    sources = query.sources if isinstance(query, CollectionView) else ()
    files = (
        ()
        if path is None
        else (
            {
                "name": path.name,
                "bytes": path.stat().st_size,
                "sha256": file_digest(path),
            },
        )
    )
    return ExportReceipt(
        str(path) if path is not None else "<stream>",
        format_name,
        view,
        count,
        sources,
        (),
        files,
    )


def _write_records(query: CollectionView, stream: TextIO, format_name: str) -> int:
    count = 0
    writer = None
    if format_name in {"csv", "tsv"}:
        cls = ArraySequence if query.view == "sequences" else ArrayPlacement
        writer = csv.DictWriter(
            stream,
            fieldnames=["schema", *(f.name for f in fields(cls))],
            delimiter="," if format_name == "csv" else "\t",
            lineterminator="\n",
        )
        writer.writeheader()
    elif format_name == "json":
        stream.write(
            canonical_json(
                {
                    "schema": "dense_arrays.record_export.v1",
                    "view": query.view,
                    "sources": list(query.sources),
                }
            )[:-1]
            + ',"records":['
        )
    with query.records() as rows:
        for row in rows:
            value = row.to_dict()
            if writer is not None:
                writer.writerow(value)
            elif format_name == "fasta":
                stream.write(f">{quote(row.array_id, safe='')}\n{row.sequence}\n")
            else:
                stream.write(("," if count else "") + canonical_json(value))
            count += 1
    if format_name == "json":
        stream.write("]}\n")
    return count


def _export_limits(
    read_limits: ReadLimits | None, compare: object, limit: int | None
) -> ReadLimits:
    """Separate admission from source resolution so neither frontend skips it."""
    if read_limits is not None and not isinstance(read_limits, ReadLimits):
        msg = "read_limits must be ReadLimits"
        raise TypeError(msg)
    if compare is not None or limit is not None:
        msg = "array exports do not accept compare or display limits"
        raise ValueError(msg)
    return read_limits or ReadLimits()
