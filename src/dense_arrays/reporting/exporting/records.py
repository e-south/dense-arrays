"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/exporting/records.py

Streaming data formats over the canonical record readers.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import csv
import hashlib
from contextlib import nullcontext
from dataclasses import fields
from pathlib import Path
from typing import TYPE_CHECKING
from urllib.parse import quote

from dense_arrays._record_validation import canonical_json
from dense_arrays.artifacts.publication import new_text_file
from dense_arrays.artifacts.receipts import ExportReceipt
from dense_arrays.reporting.projections import PlacementRecord, SequenceRecord

if TYPE_CHECKING:
    from typing import TextIO

    from .source_bindings import ExportQuery


from .formats import validate_format
from .source_bindings import bind_sources, check_sources


def export_records(
    query: ExportQuery, *, format_name: str, out: str | Path | TextIO
) -> ExportReceipt:
    """Publish a complete query to a new file, or write to a caller-owned stream."""
    validate_format(query.view, format_name)
    if query.limit is not None or query.after is not None:
        msg = "record export requires all=True and a complete query, not a page"
        raise ValueError(msg)
    path = Path(out).absolute() if isinstance(out, (str, Path)) else None
    if path is None and not callable(getattr(out, "write", None)):
        msg = "out must be a file path or writable text stream"
        raise TypeError(msg)
    sources = bind_sources(query)
    with new_text_file(path) if path is not None else nullcontext(out) as stream:
        count, references = _write_records(query, format_name, stream, sources)
        check_sources(query, sources)
    files = ()
    if path is not None:
        with path.open("rb") as stream:
            checksum = hashlib.file_digest(stream, "sha256").hexdigest()
        files = ({"name": path.name, "bytes": path.stat().st_size, "sha256": checksum},)
    return ExportReceipt(
        str(path) if path is not None else "<stream>",
        format_name,
        query.view,
        count,
        sources,
        tuple(references),
        files,
    )


def _write_records(
    query: ExportQuery,
    format_name: str,
    stream: TextIO,
    sources: tuple[dict[str, object], ...],
) -> tuple[int, dict[str, None]]:
    """Stream one complete selection, retaining only receipt identities."""
    references: dict[str, None] = {}
    count = 0
    writer = None
    if format_name in {"csv", "tsv"}:
        cls = SequenceRecord if query.view == "sequences" else PlacementRecord
        writer = csv.DictWriter(
            stream,
            fieldnames=["schema", *(f.name for f in fields(cls))],
            delimiter="," if format_name == "csv" else "\t",
            lineterminator="\n",
        )
        writer.writeheader()
    elif format_name == "json":
        header = canonical_json(
            {
                "schema": "dense_arrays.record_export.v1",
                "view": query.view,
                "sources": list(sources),
            }
        )
        stream.write(header[:-1] + ',"records":[')
    with query.records() as records:
        for record in records:
            value = record.to_dict()
            ref = getattr(record, "design_ref", None) or getattr(
                record, "reference", None
            )
            if ref is not None and ref not in references:
                records.retain_identities()
                references[ref] = None
            if writer is not None:
                writer.writerow(value)
            elif format_name == "fasta":
                stream.write(
                    f">{quote(record.design_ref, safe='/:-._~')} "
                    f"sequence_id={record.sequence_id}\n{record.sequence}\n"
                )
            else:
                stream.write(("," if count else "") + canonical_json(value))
            count += 1
    if format_name == "json":
        stream.write("]}\n")
    return count, references
