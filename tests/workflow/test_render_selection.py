"""Render one explicitly selected design from reusable native evidence.

Author: Eric J. South.
"""

import json
import shutil
from dataclasses import replace
from pathlib import Path

import pytest
from PIL import Image
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays.artifacts.records import Design
from dense_arrays.artifacts.store import create_run
from dense_arrays.cli import app


def _request() -> planning.DesignSpec:
    return planning.DesignSpec(
        parts=[parts.Part("a", "AAA"), parts.Part("b", "CCC")],
        length=planning.Length(maximum=6),
        strands="single",
        target=planning.Target(2),
    )


def test_selected_design_python_cli_match_without_generation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    run = da.run(_request(), out=tmp_path / "run")
    with da.inspect(run, view="designs", all=True).records() as rows:
        chosen = list(rows)[1]

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("render must only read persisted design evidence")

    monkeypatch.setattr(da.Optimizer, "solve_report", forbidden)
    monkeypatch.setattr("dense_arrays.generation.acceptance.evaluate", forbidden)
    selected = reporting.DesignFilter(design_ids=(chosen.reference,))
    receipt = da.render(run, select=selected, out=tmp_path / "python.png")
    assert receipt.design_refs == (chosen.reference,)
    with Image.open(tmp_path / "python.png") as image:
        assert chosen.reference in image.info["Description"]
    result = CliRunner().invoke(
        app,
        [
            "render",
            str(run.path),
            "--design-id",
            chosen.reference,
            "--out",
            str(tmp_path / "cli.png"),
            "--json",
        ],
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout)["design_refs"] == [chosen.reference]
    assert json.loads(result.stdout)["sources"] == receipt.to_dict()["sources"]


def test_render_requires_one_match_and_enforces_read_limits(tmp_path: Path):
    run = da.run(_request(), out=tmp_path / "run")
    cases = (
        ({}, "exactly one"),
        (
            {
                "select": reporting.DesignFilter(
                    metrics={"length": reporting.Range(min=20)}
                )
            },
            "exactly one",
        ),
    )
    for index, (options, message) in enumerate(cases):
        out = tmp_path / f"absent-{index}" / "design.png"
        with pytest.raises(ValueError, match=message):
            da.render(run, out=out, **options)
        assert not out.parent.exists()
    with da.inspect(run, view="designs", limit=1).records() as rows:
        chosen = next(rows)
    out = tmp_path / "limited" / "design.png"
    with pytest.raises(reporting.ReadLimitError):
        da.render(
            run,
            select=reporting.DesignFilter(design_ids=(chosen.reference,)),
            read_limits=reporting.ReadLimits(records=1),
            out=out,
        )
    assert not out.parent.exists()


def test_render_selected_matrix_design_from_portable_bundle(tmp_path: Path):
    matrix = planning.MatrixSpec(
        base=_request().with_changes(target=planning.Target(1)),
        axes={
            "variant": {
                "original": planning.Variant(),
                "alternate": planning.Variant(parts=[parts.Part("a", "GGG")]),
            }
        },
        allocation=planning.Allocation(per_cell=1),
        max_cells=2,
    )
    run = da.run(matrix, out=tmp_path / "run")
    with da.inspect(run, view="designs", all=True).records() as rows:
        chosen = list(rows)[1]
    selection = reporting.DesignFilter(design_ids=(chosen.reference,))
    native = da.render(run, select=selection, out=tmp_path / "native.png")
    bundle = tmp_path / "bundle"
    da.export(run, all=True, format="bundle", out=bundle)
    shutil.rmtree(run.path)
    portable = da.render(bundle, select=selection, out=tmp_path / "portable.png")
    assert portable.design_refs == native.design_refs == (chosen.reference,)
    assert portable.sources[0]["bundle_id"] == da.inspect(bundle).bundle_id


def test_saved_selection_render_keeps_revision_after_source_advances(tmp_path: Path):
    reference = da.run(_request(), out=tmp_path / "reference")
    with da.inspect(reference, view="designs", all=True).records() as rows:
        candidates = list(rows)
    evidence = {"solver_status": "optimal", "proof_scope": "offered_packing_model"}
    with create_run(da.inspect(reference, view="plan"), tmp_path / "active") as writer:

        def publish(candidate: Design) -> None:
            ordinal = writer.reserve(active_seconds=0)
            writer.publish(
                ordinal,
                "accepted",
                evidence,
                active_seconds=0,
                design=replace(candidate, run_id=writer.handle.run_id),
            )

        publish(candidates[0])
        snapshot = da.inspect(
            writer.handle,
            view="selection",
            select=reporting.LibrarySelection(take=reporting.Take(count=1)),
        )
        publish(candidates[1])
        rendered = da.render(
            writer.handle, select=snapshot, out=tmp_path / "pinned.png"
        )
        assert rendered.design_refs == tuple(snapshot.references())
        assert rendered.sources[0]["revision"] == snapshot.sources[0].revision
        assert rendered.to_dict()["selection"] == snapshot.summary()
        saved = tmp_path / "selection.json"
        da.export(snapshot, format="selection", out=saved)
        result = CliRunner().invoke(
            app,
            [
                "render",
                str(writer.handle.path),
                "--selection",
                str(saved),
                "--out",
                str(tmp_path / "pinned-cli.png"),
                "--json",
            ],
        )
        assert result.exit_code == 0, result.output
        assert json.loads(result.stdout)["sources"] == rendered.to_dict()["sources"]
        assert json.loads(result.stdout)["design_refs"] == list(snapshot.references())


def test_render_combined_sources_requires_an_unambiguous_selection(tmp_path: Path):
    left = da.run(_request(), out=tmp_path / "left")
    right = da.run(_request(), out=tmp_path / "right")
    with da.inspect(right, view="designs", limit=1).records() as rows:
        chosen = next(rows)
    with pytest.raises(ValueError, match=r"exactly one|ambiguous"):
        da.render(
            [left, right],
            select=reporting.DesignFilter(design_ids=(chosen.design_id,)),
            out=tmp_path / "absent.png",
        )
    rendered = da.render(
        [left, right],
        select=reporting.DesignFilter(design_ids=(chosen.reference,)),
        out=tmp_path / "selected.png",
    )
    assert rendered.design_refs == (chosen.reference,)
    assert rendered.sources[1]["plan_id"] == da.inspect(right).plan_id
