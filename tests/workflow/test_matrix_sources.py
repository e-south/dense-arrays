"""Cells select explicit sources and retain their own collection provenance.

Author: Eric J. South.
"""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.artifacts.run_plans import decode_plan
from dense_arrays.artifacts.store import RunWriter
from dense_arrays.cli import app
from dense_arrays.reporting import ReadLimitError
from dense_arrays.reporting.plans.reading import plan_identities
from dense_arrays.workflow.inputs import read_source


def source_matrix(sources: dict):
    return planning.MatrixSpec(
        planning.DesignSpec(
            [parts.Part("default", "TTT")], planning.Length(maximum=3), strands="single"
        ),
        axes={"pool": {"default": planning.Variant(), "selected": planning.Variant()}},
        allocation=planning.Allocation(per_cell=1),
        max_cells=2,
        sources=sources,
    )


def test_per_cell_pool_selection_preserves_parts_provenance_and_cli(tmp_path: Path):
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence,group\na,AAA,A\nb,CCC,B\n")
    pool = da.prepare(
        parts.PreparationSpec(parts.PartTable(table, "csv")), out=tmp_path / "pool"
    )
    request = source_matrix(
        {"pool=selected": parts.PoolSource(pool, parts.PartFilter(groups=("B",)))}
    )
    plan = da.plan(request)
    first, selected = plan.cells
    assert [p.part_id for p in first.plan.request.parts] == ["default"]
    assert [p.part_id for p in selected.plan.request.parts] == ["b"]
    assert selected.plan.collection_id == pool.pool_id
    assert selected.plan.import_report.selection == parts.PartFilter(groups=("B",))
    assert first.plan.collection_id != selected.plan.collection_id
    saved = tmp_path / "matrix.plan.json"
    plan.write(saved)
    assert read_source(saved).to_dict() == plan.to_dict()
    cli = CliRunner().invoke(
        app, ["run", str(saved), "--out", str(tmp_path / "cli"), "--json"]
    )
    assert cli.exit_code == 0, cli.output
    run = da.run(plan, out=tmp_path / "python")
    assert da.inspect(run, verify=True).accepted == 2
    with da.inspect(run, view="designs", all=True).records() as rows:
        assert {r.cell_id: r.realized.sequence for r in rows} == {
            "pool=default": "TTT",
            "pool=selected": "CCC",
        }
    assert da.inspect(tmp_path / "cli", verify=True).accepted == 2
    assert json.loads(cli.stdout)["plan_id"] == plan.plan_id


def test_sources_validate_cell_references_before_reading_inputs(tmp_path: Path):
    missing = parts.PartTable(tmp_path / "missing.csv", "csv")
    request = source_matrix({"pool=typo": missing})
    request = request.with_changes(base=request.base.with_changes(parts=missing))
    with pytest.raises(ValueError, match=r"sources.*unknown cells"):
        da.plan(request)


def test_source_and_variants_resolve_together(tmp_path: Path):
    request = source_matrix({"pool=selected": [parts.Part("alternative", "AAA")]})
    request = request.with_changes(
        base=request.base.with_changes(
            requirements=[planning.Fixed("anchor", "default", "forward")]
        ),
        axes={
            "pool": {
                "default": planning.Variant(),
                "selected": planning.Variant(
                    parts=[parts.Part("alternative", "CCC")],
                    requirements=[planning.Fixed("anchor", "alternative", "forward")],
                ),
            }
        },
    )
    plan = da.plan(request)
    run = da.run(plan, out=tmp_path / "run")
    assert da.inspect(run, verify=True).accepted == 2
    with da.inspect(run, view="designs", all=True).records() as records:
        assert {d.cell_id: d.realized.sequence for d in records} == {
            "pool=default": "TTT",
            "pool=selected": "CCC",
        }
    invalid = request.with_changes(
        sources={"pool=selected": [parts.Part("unknown", "AAA")]}
    )
    with pytest.raises(ValueError, match="known part"):
        da.plan(invalid)


def test_source_locators_are_relative_and_not_cell_identity(tmp_path: Path):
    original = tmp_path / "original"
    original.mkdir()
    table = original / "parts.csv"
    table.write_text("part_id,sequence\nx,AAA\n")
    request = source_matrix({"pool=selected": parts.PartTable(table, "csv")})
    plan = da.plan(request)
    plan.write(original / "plan.json")
    (original / "request.json").write_text(json.dumps(request.to_dict(base=original)))
    moved = tmp_path / "moved"
    original.rename(moved)
    loaded = read_source(moved / "plan.json")
    assert loaded.plan_id == plan.plan_id
    loaded.verify_inputs()
    assert da.plan(read_source(moved / "request.json")).plan_id == plan.plan_id
    (moved / "parts.csv").write_text("part_id,sequence\nx,CCC\n")
    with pytest.raises(ValueError, match="input changed"):
        da.run(loaded, out=tmp_path / "stale")
    assert not (tmp_path / "stale").exists()
    changed = da.plan(read_source(moved / "request.json"))
    assert changed.cells[1].cell_id == loaded.cells[1].cell_id
    assert changed.cells[1].plan.plan_id != loaded.cells[1].plan.plan_id
    (moved / "parts.csv").unlink()
    assert read_source(moved / "plan.json").plan_id == plan.plan_id


@pytest.mark.parametrize("operation", ["batches", "extension"])
def test_portable_plans_retain_each_source_without_reopening_it(
    tmp_path: Path, operation: str
):
    table = tmp_path / "source.csv"
    table.write_text("part_id,sequence\nx,AAA\ny,CCC\n")
    plan = da.plan(source_matrix({"pool=selected": parts.PartTable(table, "csv")}))
    if operation == "batches":
        saved = tmp_path / "batch.json"
        portable = da.prepare(
            plan, sampling=planning.BatchSampling(size=1, seed=8), out=saved
        )
        expected = 2
    else:
        parent = da.run(plan, out=tmp_path / "parent")
        portable = da.plan(
            planning.ExtensionSpec(
                planning.ParentRun(parent.path),
                {"pool=default": 0, "pool=selected": 1},
                planning.Limits(),
                8,
            )
        )
        saved = tmp_path / "extension.json"
        portable.write(saved)
        parent.path.rename(tmp_path / "moved-parent")
        expected = 1
    table.unlink()
    replay = read_source(saved)
    replay.verify_inputs()
    assert replay.plan_id == portable.plan_id
    assert replay.cells[1].plan.collection_id == plan.cells[1].plan.collection_id
    run = da.run(replay, out=tmp_path / "run")
    assert da.inspect(run, verify=True).accepted == expected
    editable = da.inspect(run, view="request").request
    assert da.plan(editable).plan_id == portable.plan_id
    da.export(run, all=True, format="bundle", out=tmp_path / "bundle")
    assert da.inspect(tmp_path / "bundle", verify=True).designs == expected


def test_stored_sources_are_charged_before_materializing(
    monkeypatch: pytest.MonkeyPatch,
):
    plan = da.plan(source_matrix({"pool=selected": [parts.Part("x", "AAA")]}))
    # Two cell identities, base + two cell parts, and one source snapshot part.
    assert plan_identities(plan) == 6
    assert da.inspect(plan, view="request").identities == 4
    monkeypatch.setattr(
        planning.MatrixPlan,
        "from_dict",
        lambda *_a, **_k: pytest.fail("read cap admitted plan"),
    )
    with pytest.raises(ReadLimitError):
        decode_plan(plan.to_dict(), max_identities=5)


def test_human_preview_explains_combinations_and_selected_parts(tmp_path: Path):
    plan = da.plan(
        source_matrix(
            {"pool=selected": [parts.Part("x", "AAA"), parts.Part("y", "CCC")]}
        )
    )
    saved = tmp_path / "plan.json"
    plan.write(saved)
    result = CliRunner().invoke(app, ["plan", str(saved)])
    assert result.exit_code == 0, result.output
    assert "design combinations" in result.output
    assert "pool=selected: target 1 (active); 2 eligible parts" in result.output
    assert plan.cells[1].plan.collection_id[:12] in result.output


def test_source_specific_resampling_survives_request_export_and_resume(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    request = source_matrix(
        {"pool=selected": [parts.Part("x", "AAA"), parts.Part("y", "CCC")]}
    )
    request = request.with_changes(
        batches={
            "pool=selected": planning.Resampling(
                planning.BatchSampling(size=1, seed=10),
                max_batches=4,
                attempts_per_batch=1,
            )
        }
    )
    plan = da.plan(request)
    assert da.inspect(plan, view="request").identities == 5
    saved_request = tmp_path / "request.json"
    da.export(plan, view="request", out=saved_request)
    assert da.plan(read_source(saved_request)).plan_id == plan.plan_id
    publish = RunWriter.publish

    def interrupt(self: RunWriter, *args: object, **kwargs: object) -> str:
        publish(self, *args, **kwargs)
        raise KeyboardInterrupt

    output = tmp_path / "run"
    with monkeypatch.context() as patch:
        patch.setattr(RunWriter, "publish", interrupt)
        with pytest.raises(KeyboardInterrupt):
            da.run(plan, out=output)
    before = list(da.inspect(output, view="designs", all=True).records())
    run = da.run(resume=output)
    assert da.inspect(run, verify=True).accepted == 2
    with da.inspect(run, view="designs", all=True).records() as rows:
        assert list(rows)[: len(before)] == before


def test_bound_source_tampering_and_unknown_parts_fail_without_output():
    plan = da.plan(source_matrix({"pool=selected": [parts.Part("x", "AAA")]}))
    document = plan.to_dict()
    document["request"]["sources"]["pool=selected"]["parts"][0]["sequence"] = "CCC"
    with pytest.raises(ValueError, match="bound parts identity"):
        planning.MatrixPlan.from_dict(document)
    raw = source_matrix({"pool=selected": [parts.Part("x", "AAA")]})
    with pytest.raises(TypeError, match="frozen BoundParts"):
        planning.MatrixPlan(raw, da.plan(raw.base))
    with pytest.raises(ValueError, match=r"unique|duplicate"):
        source_matrix(
            {"pool=selected": [parts.Part("x", "AAA"), parts.Part("x", "CCC")]}
        )
