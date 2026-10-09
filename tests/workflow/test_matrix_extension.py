"""Additional targets preserve matrix cells and exclude their own accepted lineage.

Author: Eric J. South.
"""

import json
import sqlite3
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.artifacts.store import RunWriter, create_run
from dense_arrays.cli import app


def matrix(*, counts: dict[str, int] | None = None) -> planning.MatrixSpec:
    return planning.MatrixSpec(
        base=planning.DesignSpec(
            parts=[
                parts.Part("a", "AAA"),
                parts.Part("b", "CCC"),
                parts.Part("c", "GGG"),
            ],
            length=planning.Length(maximum=3),
            strands="single",
        ),
        axes={"x": {"a": planning.Variant(), "b": planning.Variant()}},
        allocation=planning.Allocation(counts=counts or {"x=a": 1, "x=b": 1}),
        max_cells=2,
    )


def extension(parent: da.artifacts.RunHandle, counts: dict[str, int]):
    return planning.ExtensionSpec(
        parent=planning.ParentRun(parent.path),
        additional=counts,
        limits=planning.Limits(attempts=30),
        seed=31,
    )


def test_matrix_extensions_keep_cell_targets_and_all_ancestor_exclusions(
    tmp_path: Path,
):
    parent = da.run(matrix(), out=tmp_path / "parent")
    before = (parent.path / "run.sqlite3").read_bytes()
    request = extension(parent, {"x=a": 1, "x=b": 1})
    plan = da.plan(request)
    assert plan.preview["excluded_sequences"] == 2
    child = da.run(plan, out=tmp_path / "child")
    grand = da.run(extension(child, {"x=a": 1, "x=b": 1}), out=tmp_path / "grand")
    for run in (child, grand):
        summary = da.inspect(run, verify=True)
        assert summary.state == "completed"
        assert summary.accepted == summary.target == 2
        assert all(c.accepted == c.target == 1 for c in summary.cells.values())
    designs = list(
        da.inspect([parent, child, grand], view="designs", all=True).records()
    )
    for cell in ("x=a", "x=b"):
        assert len({d.sequence_id for d in designs if d.cell_id == cell}) == 3
    assert (parent.path / "run.sqlite3").read_bytes() == before
    assert planning.MatrixPlan.from_dict(plan.to_dict()).plan_id == plan.plan_id
    recovered = da.inspect(child, view="request").request
    assert da.plan(recovered).plan_id == plan.plan_id


def test_activating_a_cell_does_not_inherit_other_cells_sequences(tmp_path: Path):
    spec = matrix(counts={"x=a": 1, "x=b": 0})
    spec = spec.with_changes(
        base=spec.base.with_changes(parts=[parts.Part("a", "AAA")])
    )
    parent = da.run(spec, out=tmp_path / "parent")
    child = da.run(extension(parent, {"x=a": 0, "x=b": 1}), out=tmp_path / "child")
    summary = da.inspect(child, verify=True)
    assert summary.accepted == 1
    assert summary.cells["x=a"].state == "inactive"
    records = list(da.inspect([parent, child], view="designs", all=True).records())
    assert records[0].sequence_id == records[1].sequence_id
    assert records[0].cell_id != records[1].cell_id
    grand = da.run(extension(child, {"x=a": 1, "x=b": 0}), out=tmp_path / "grand")
    summary = da.inspect(grand, verify=True)
    assert summary.accepted == 0
    assert summary.cells["x=a"].counts["duplicate"] == 1
    assert summary.cells["x=a"].termination_reason == "batch_exhausted"


def test_extension_cli_and_exported_request_retain_frozen_cell_evidence(tmp_path: Path):
    source = tmp_path / "parts.csv"
    source.write_text("part_id,sequence\na,AAA\nb,CCC\nc,GGG\n")
    spec = matrix()
    parent = da.run(
        spec.with_changes(
            base=spec.base.with_changes(parts=parts.PartTable(source, "csv"))
        ),
        out=tmp_path / "parent",
    )
    request = extension(parent, {"x=a": 1, "x=b": 0})
    resolved = da.plan(request)
    input_path = tmp_path / "extend.json"
    input_path.write_text(json.dumps(request.to_dict(base=tmp_path)))
    plan_path = tmp_path / "extension-plan.json"
    runner = CliRunner()
    response = runner.invoke(
        app, ["plan", str(input_path), "--out", str(plan_path), "--json"]
    )
    assert response.exit_code == 0, response.output
    assert json.loads(response.stdout) == resolved.to_dict()
    parent.path.rename(tmp_path / "moved-parent")
    source.unlink()
    output = tmp_path / "child"
    response = runner.invoke(
        app, ["run", str(plan_path), "--out", str(output), "--json"]
    )
    assert response.exit_code == 0, response.output
    assert da.inspect(output, verify=True).accepted == 1
    editable = tmp_path / "revised.json"
    da.export(output, view="request", format="json", out=editable)
    response = runner.invoke(app, ["plan", str(editable), "--json"])
    assert response.exit_code == 0, response.output
    assert json.loads(response.stdout)["plan_id"] == resolved.plan_id
    da.export(output, all=True, format="bundle", out=tmp_path / "bundle")
    assert da.inspect(tmp_path / "bundle", verify=True).designs == 1
    quality = da.inspect(tmp_path / "bundle", view="quality").to_dict()
    assert quality["selection"]["designs"] == 1


@pytest.mark.parametrize("counts", [{"x=a": 1}, {"x=a": 1, "x=b": 0, "x=c": 0}])
def test_additional_targets_require_exactly_the_parent_cell_inventory(
    tmp_path: Path, counts: dict[str, int]
):
    parent = da.run(matrix(), out=tmp_path / "parent")
    before = (parent.path / "run.sqlite3").read_bytes()
    with pytest.raises(ValueError, match="cell"):
        da.plan(extension(parent, counts))
    assert (parent.path / "run.sqlite3").read_bytes() == before


@pytest.mark.parametrize(
    "counts", [{}, {"x=a": 0}, {"x=a": -1}, {"x=a": True}, {"x=a": 1.5}]
)
def test_additional_cell_counts_are_strict_and_request_work(counts: dict[str, object]):
    with pytest.raises((ValueError, TypeError), match="additional"):
        planning.ExtensionSpec(
            planning.ParentRun("unused"), counts, planning.Limits(), 1
        )


def test_extension_limits_and_parent_lock_precede_output(tmp_path: Path):
    resolved = da.plan(matrix())
    live = tmp_path / "live"
    with create_run(resolved, live), pytest.raises(ValueError, match="terminal"):
        da.plan(
            planning.ExtensionSpec(
                planning.ParentRun(live), {"x=a": 1, "x=b": 1}, planning.Limits(), 2
            )
        )
    parent = da.run(resolved, out=tmp_path / "parent")
    request = extension(parent, {"x=a": 1, "x=b": 1})
    for limits in (reporting.ReadLimits(records=1), reporting.ReadLimits(identities=1)):
        with pytest.raises(reporting.ReadLimitError):
            da.plan(request, read_limits=limits)
    with pytest.raises(TypeError, match="per-cell"):
        da.plan(
            planning.ExtensionSpec(
                planning.ParentRun(parent.path), 1, planning.Limits(), 2
            )
        )


def test_duplicate_verification_rejects_a_matching_sequence_from_another_cell(
    tmp_path: Path,
):
    spec = matrix()
    spec = spec.with_changes(
        base=spec.base.with_changes(parts=[parts.Part("a", "AAA")])
    )
    parent = da.run(spec, out=tmp_path / "parent")
    child = da.run(extension(parent, {"x=a": 1, "x=b": 1}), out=tmp_path / "child")
    other_ref = next(
        d.reference
        for d in da.inspect(parent, view="designs", all=True).records()
        if d.cell_id == "x=b"
    )
    with sqlite3.connect(child.path / "run.sqlite3") as connection:
        revision, payload = connection.execute(
            "SELECT revision,payload FROM attempts WHERE "
            "json_extract(payload,'$.outcome')='duplicate' AND "
            "json_extract(payload,'$.cell_id')='x=a'"
        ).fetchone()
        value = json.loads(payload)
        value["evidence"]["matched_design_ref"] = other_ref
        connection.execute(
            "UPDATE attempts SET payload=?,digest=? WHERE revision=?",
            (canonical_json(value), semantic_digest(value), revision),
        )
    with pytest.raises(ValueError, match="parent duplicate"):
        da.inspect(child, verify=True)


def test_resuming_an_extension_preserves_cell_exclusions(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    parent = da.run(matrix(), out=tmp_path / "parent")
    child_path = tmp_path / "child"
    publish = RunWriter.publish

    def interrupted(self: RunWriter, *args: object, **kwargs: object) -> str:
        outcome = publish(self, *args, **kwargs)
        if outcome == "duplicate":
            raise KeyboardInterrupt
        return outcome

    with monkeypatch.context() as patch:
        patch.setattr(RunWriter, "publish", interrupted)
        with pytest.raises(KeyboardInterrupt):
            da.run(extension(parent, {"x=a": 1, "x=b": 1}), out=child_path)
    before = da.inspect(child_path, verify=True)
    assert before.resumable
    child = da.run(resume=child_path)
    assert da.inspect(child, verify=True).accepted == 2
    records = list(da.inspect([parent, child], view="designs", all=True).records())
    assert all(
        len({d.sequence_id for d in records if d.cell_id == c}) == 2
        for c in ("x=a", "x=b")
    )


def test_matrix_base_cannot_silently_broadcast_a_single_cell_parent(tmp_path: Path):
    initial = da.run(matrix().base, out=tmp_path / "single")
    inherited = da.plan(
        planning.ExtensionSpec(
            planning.ParentRun(initial.path), 1, planning.Limits(), 3
        )
    )
    request = matrix().with_changes(base=inherited.request)
    with pytest.raises(ValueError, match=r"base.*exclusion|per-cell"):
        planning.MatrixPlan(request, inherited)


def test_unknown_exclusion_cells_fail_before_source_reads(tmp_path: Path):
    parent = da.run(matrix(), out=tmp_path / "parent")
    request = da.inspect(
        da.plan(extension(parent, {"x=a": 1, "x=b": 1})), view="request"
    ).request
    request = request.with_changes(
        axes={"different": {"a": planning.Variant(), "b": planning.Variant()}},
        base=request.base.with_changes(
            parts=parts.PartTable(tmp_path / "missing.csv", "csv")
        ),
        allocation=planning.Allocation(per_cell=1),
    )
    with pytest.raises(ValueError, match="unknown cells"):
        da.plan(request)
