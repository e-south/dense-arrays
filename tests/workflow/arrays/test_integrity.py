"""Reject malformed evidence and preserve exact source metadata across exports."""

import json
from dataclasses import replace
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.arrays import ArrayCollection, ArrayFilter
from dense_arrays.artifacts.errors import ArtifactIntegrityError
from dense_arrays.cli import app
from dense_arrays.parts import PartFilter
from dense_arrays.parts.filters import Range


def test_provenance_survives_transport_and_collection_reexport(
    tmp_path: Path, source: ArrayCollection
):
    provenance = {
        "generation": {
            "software": [{"name": "dense-arrays", "version": "0.1.0"}],
            "solver": {"backend": "GUROBI", "version": None},
            "limits": {"attempt_seconds": 5.0, "threads": 12},
            "seeds": {"root": 42, "solver": "14616319817112705343"},
            "evidence": {"configuration_sha256": "a" * 64},
        }
    }
    supplied = replace(source, provenance=provenance)
    transport = tmp_path / "input.jsonl"
    da.export(supplied, view="arrays", all=True, format="jsonl", out=transport)
    da.export(transport, all=True, format="bundle", out=tmp_path / "first")
    da.export(tmp_path / "first", all=True, format="bundle", out=tmp_path / "second")
    summary = da.inspect(tmp_path / "second", verify=True).to_dict()
    assert summary["provenance"] == provenance
    assert summary["exporter"]["version"]
    assert "solver" not in summary["exporter"]


@pytest.mark.parametrize(
    "exporter", [None, {}, {"package": "dense-arrays", "version": ""}]
)
def test_manifest_rejects_invalid_exporter_even_with_recomputed_digest(
    tmp_path: Path, source: ArrayCollection, exporter: object
):
    da.export(source, all=True, format="bundle", out=tmp_path / "data")
    path = tmp_path / "data" / "collection.json"
    value = json.loads(path.read_text())
    value.pop("collection_id")
    value["exporter"] = exporter
    value["collection_id"] = semantic_digest(value)
    path.write_text(canonical_json(value))
    with pytest.raises(ArtifactIntegrityError, match=r"exporter|package|version"):
        da.inspect(path.parent)


@pytest.mark.parametrize(
    "defect", ["missing-footer", "boolean-count", "trailing-record"]
)
def test_incomplete_transport_cannot_publish(
    tmp_path: Path, source: ArrayCollection, defect: str
):
    transport = tmp_path / "input.jsonl"
    da.export(
        replace(source, arrays=tuple(source.arrays)[:1]),
        all=True,
        format="jsonl",
        out=transport,
    )
    lines = transport.read_text().splitlines(keepends=True)
    if defect == "missing-footer":
        lines.pop()
    elif defect == "boolean-count":
        footer = json.loads(lines[-1])
        footer["arrays"] = True
        lines[-1] = canonical_json(footer) + "\n"
    else:
        lines.append(lines[-1])
    transport.write_text("".join(lines))
    with pytest.raises(ValueError, match=r"truncated|completion"):
        da.export(transport, all=True, format="bundle", out=tmp_path / "bad")
    assert not (tmp_path / "bad").exists()


def test_supplied_array_can_render_without_a_generation_plan(
    tmp_path: Path, source: ArrayCollection
):
    da.export(source, all=True, format="bundle", out=tmp_path / "data")
    receipt = da.render(
        tmp_path / "data",
        view="array",
        select=ArrayFilter(array_ids=("array-2",)),
        out=tmp_path / "array.png",
    )
    assert receipt.records == 1
    assert receipt.view == "array"
    assert (tmp_path / "array.png").read_bytes().startswith(b"\x89PNG")
    with pytest.raises(ValueError, match="exactly one"):
        da.render(tmp_path / "data", view="array", out=tmp_path / "ambiguous.png")
    assert not (tmp_path / "ambiguous.png").exists()


def test_cli_render_honors_array_selection(tmp_path: Path, source: ArrayCollection):
    da.export(source, all=True, format="bundle", out=tmp_path / "data")
    result = CliRunner().invoke(
        app,
        [
            "render",
            str(tmp_path / "data"),
            "--view",
            "array",
            "--array-id",
            "array-2",
            "--out",
            str(tmp_path / "cli.png"),
            "--json",
        ],
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout)["view"] == "array"


def test_metric_filter_is_not_silently_discarded(
    tmp_path: Path, source: ArrayCollection
):
    da.export(source, all=True, format="bundle", out=tmp_path / "data")
    with pytest.raises(ValueError, match="metric"):
        da.inspect(
            tmp_path / "data",
            view="parts",
            select=PartFilter(metrics={"length": Range(30, 40)}),
        )
