"""Matrix execution owns one run and records each cell's actual attainment.

Author: Eric J. South.
"""

import json
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.artifacts.errors import ArtifactIntegrityError
from dense_arrays.artifacts.recovery import own_run
from dense_arrays.artifacts.store import RunWriter
from dense_arrays.cli import app
from dense_arrays.workflow import execution, matrices

from .test_recovery_processes import SCRIPT


def matrix(*, attempts: int = 10, per_cell: int = 1):
    return planning.MatrixSpec(
        base=planning.DesignSpec(
            [parts.Part("a", "AAA")],
            planning.Length(maximum=3),
            strands="single",
            limits=planning.Limits(attempts=attempts),
        ),
        axes={
            "condition": {
                "first": planning.Variant(),
                "second": planning.Variant(),
                "inactive": planning.Variant(),
            }
        },
        allocation=planning.Allocation(
            counts={
                "condition=first": per_cell,
                "condition=second": per_cell,
                "condition=inactive": 0,
            }
        ),
        max_cells=3,
    )


def test_matrix_run_reconciles_cells_and_allows_equal_sequences_across_cells(
    tmp_path: Path,
):
    resolved = da.plan(matrix())
    run = da.run(resolved, out=tmp_path / "matrix")
    summary = da.inspect(run, verify=True)
    assert summary.state == "completed"
    assert summary.accepted == summary.target == summary.counts["started"] == 2
    assert summary.plan_id == resolved.plan_id
    assert {name: cell.accepted for name, cell in summary.cells.items()} == {
        "condition=first": 1,
        "condition=second": 1,
        "condition=inactive": 0,
    }
    assert summary.cells["condition=inactive"].counts["started"] == 0
    designs = list(da.inspect(run, view="designs", all=True).records())
    assert [d.cell_id for d in designs] == ["condition=first", "condition=second"]
    assert designs[0].sequence_id == designs[1].sequence_id
    assert designs[0].reference != designs[1].reference
    assert {d.plan_id for d in designs} == {
        c.plan.plan_id for c in resolved.cells if c.active
    }
    selected = da.inspect(
        run, view="designs", select=reporting.DesignFilter(cells=("condition=second",))
    )
    assert [d.reference for d in selected.records()] == [designs[1].reference]


def test_matrix_attempt_limit_is_shared_and_does_not_reallocate_shortfalls(
    tmp_path: Path,
):
    run = da.run(matrix(attempts=1), out=tmp_path / "limited")
    summary = da.inspect(run, verify=True)
    assert summary.state == "stopped"
    assert summary.termination_reason == "attempt_limit"
    assert summary.counts["started"] == summary.accepted == 1
    assert summary.target == 2
    assert summary.cells["condition=first"].accepted == 1
    assert summary.cells["condition=second"].accepted == 0
    assert summary.cells["condition=second"].target == 1


def test_matrix_reports_selection_and_portable_exports_keep_cell_annotations(
    tmp_path: Path,
):
    spec = matrix()
    spec = spec.with_changes(
        base=spec.base.with_changes(parts=[parts.Part("a", "AAA", group="A")]),
        axes={
            "condition": {
                "first": planning.Variant(),
                "second": planning.Variant(parts=[parts.Part("a", "CCC", group="B")]),
                "inactive": planning.Variant(),
            }
        },
    )
    run = da.run(spec, out=tmp_path / "run")
    placements = list(da.inspect(run, view="placements", all=True).records())
    assert [p.group for p in placements] == ["A", "B"]
    selected = reporting.DesignFilter(groups=("B",))
    rows = list(da.inspect(run, view="designs", select=selected).records())
    assert len(rows) == 1
    assert rows[0].cell_id == "condition=second"
    quality = da.inspect(run, view="quality").to_dict()
    assert quality["selection"]["designs"] == 2
    assert quality["supply"]["eligible_parts"] == 2
    assert len(quality["cells"]) == 3
    assert quality["source_runs"][0]["selected_designs"] == 2
    panel = da.inspect(
        run,
        view="selection",
        select=reporting.LibrarySelection(
            take=reporting.Take(
                per_cell={
                    f"{run.run_id}/condition=first": 0,
                    f"{run.run_id}/condition=second": 1,
                    f"{run.run_id}/condition=inactive": 0,
                }
            ),
        ),
    )
    assert list(panel.references()) == [rows[0].reference]
    bundle = tmp_path / "bundle"
    da.export(run, select=panel, format="bundle", out=bundle)
    shutil.rmtree(run.path)
    assert da.inspect(bundle, verify=True).designs == 1
    shared_quality = da.inspect(bundle, view="quality").to_dict()
    assert shared_quality["selection"]["designs"] == 1
    assert shared_quality["search"]["availability"] == "not_included"
    assert len(shared_quality["cells"]) == 3


def test_matrix_cli_and_saved_documents_use_the_same_native_run(tmp_path: Path):
    request = matrix()
    source = tmp_path / "matrix.json"
    source.write_text(json.dumps(request.to_dict()))
    destination = tmp_path / "run"
    runner = CliRunner()
    result = runner.invoke(
        app, ["run", str(source), "--out", str(destination), "--json"]
    )
    assert result.exit_code == 0, result.output
    summary = da.inspect(destination, verify=True)
    assert summary.accepted == 2
    assert (
        json.loads(result.stdout)["cells"]["condition=second"]["counts"]["accepted"]
        == 1
    )
    for view in ("plan", "request", "diagnostics", "quality"):
        response = runner.invoke(
            app, ["inspect", str(destination), "--view", view, "--json"]
        )
        assert response.exit_code == 0, response.output
        assert (
            json.loads(response.stdout) == da.inspect(destination, view=view).to_dict()
        )
    plan_path = tmp_path / "plan.json"
    da.export(destination, view="plan", format="json", out=plan_path)
    assert (
        planning.MatrixPlan.from_dict(json.loads(plan_path.read_text())).plan_id
        == da.plan(request).plan_id
    )
    recovered = da.inspect(destination, view="request").request
    assert da.plan(recovered).plan_id == da.plan(request).plan_id


def test_matrix_exhaustion_is_local_and_round_robin_keeps_other_cells_active(
    tmp_path: Path,
):
    run = da.run(matrix(per_cell=2), out=tmp_path / "exhausted")
    summary = da.inspect(run, verify=True)
    assert summary.state == "stopped"
    assert summary.termination_reason == "cells_exhausted"
    assert summary.accepted == 2
    assert summary.target == 4
    assert summary.counts["started"] == 4
    attempts = list(da.inspect(run, view="attempts", all=True).records())
    assert [a.cell_id for a in attempts] == ["condition=first", "condition=second"] * 2
    assert [a.evidence["cell_attempt"] for a in attempts] == [1, 1, 2, 2]
    assert summary.cells["condition=first"].termination_reason == "batch_exhausted"
    assert summary.cells["condition=second"].termination_reason == "batch_exhausted"


def test_all_inactive_matrix_never_builds_a_model(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("inactive matrix built a solver")

    monkeypatch.setattr(da.Optimizer, "build_model", forbidden)
    run = da.run(matrix(per_cell=0), out=tmp_path / "inactive")
    summary = da.inspect(run, verify=True)
    assert summary.state == "completed"
    assert summary.accepted == summary.target == 0
    assert len(summary.cells) == 3


def test_cell_padding_stream_is_independent_of_matrix_execution_order(tmp_path: Path):
    spec = matrix()
    spec = spec.with_changes(
        base=spec.base.with_changes(
            length=planning.Length(exact=8),
            assembly=planning.Assembly(
                padding=planning.Padding(side="right", max_trials=1)
            ),
        )
    )
    other = spec.with_changes(
        axes={"condition": dict(reversed(tuple(spec.axes["condition"].items())))}
    )
    left = da.run(spec, out=tmp_path / "left")
    right = da.run(other, out=tmp_path / "right")

    def sequences(run: object) -> dict[str, str]:
        return {
            d.cell_id: d.realized.sequence
            for d in da.inspect(run, view="designs", all=True).records()
        }

    assert sequences(left) == sequences(right)


def test_model_build_time_uses_the_shared_budget(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):

    now = [0.0]
    clock = SimpleNamespace(monotonic=lambda: now[0])
    original = execution.build_optimizer

    def slow_build(*args: object, **kwargs: object) -> da.Optimizer:
        result = original(*args, **kwargs)
        now[0] = 2.0
        return result

    monkeypatch.setattr(execution, "time", clock)
    monkeypatch.setattr(matrices, "time", clock)
    monkeypatch.setattr(execution, "build_optimizer", slow_build)
    spec = matrix()
    run = da.run(
        spec.with_changes(
            base=spec.base.with_changes(
                limits=planning.Limits(active_seconds=1),
            )
        ),
        out=tmp_path / "limited",
    )
    summary = da.inspect(run, verify=True)
    assert summary.counts["started"] == 0
    assert summary.active_seconds == 2
    assert summary.termination_reason == "active_time_limit"


def test_interrupted_commit_keeps_completed_cell_when_other_cells_resume(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):

    original = RunWriter._commit  # noqa: SLF001 - inject interruption at a native commit

    def interrupted(writer: RunWriter, state: dict, **kwargs: object) -> None:
        original(writer, state, **kwargs)
        if kwargs.get("design") is not None:
            raise KeyboardInterrupt

    monkeypatch.setattr(RunWriter, "_commit", interrupted)
    path = tmp_path / "interrupted"
    with pytest.raises(KeyboardInterrupt):
        da.run(matrix(), out=path)
    summary = da.inspect(path, verify=True)
    assert summary.accepted == 1
    assert summary.resumable
    assert summary.cells["condition=first"].state == "completed"
    assert summary.cells["condition=second"].state == "stopped"
    monkeypatch.setattr(RunWriter, "_commit", original)
    after = da.inspect(da.run(resume=path), verify=True)
    assert after.state == "completed"
    assert after.accepted == 2
    assert after.cells["condition=first"] == summary.cells["condition=first"]


def test_inactive_cell_provenance_is_checked_before_portable_publication(
    tmp_path: Path,
):
    run = da.run(matrix(), out=tmp_path / "run")
    with sqlite3.connect(run.path / "run.sqlite3") as connection:
        revision, payload = connection.execute(
            "SELECT revision,payload FROM commits ORDER BY revision DESC LIMIT 1"
        ).fetchone()
        value = json.loads(payload)
        value["cells"]["condition=inactive"]["plan_id"] = value["cells"][
            "condition=first"
        ]["plan_id"]
        connection.execute(
            "UPDATE commits SET payload=?,digest=? WHERE revision=?",
            (canonical_json(value), semantic_digest(value), revision),
        )
    with pytest.raises(ArtifactIntegrityError, match="cell plan"):
        da.inspect(run, verify=True)
    with pytest.raises((ArtifactIntegrityError, ValueError), match="cell plan"):
        da.export(run, all=True, format="bundle", out=tmp_path / "bundle")


@pytest.mark.parametrize(
    "mode,accepted,in_progress,code", [("before", 0, 1, 23), ("after", 1, 0, 24)]
)
def test_matrix_process_crash_preserves_the_committed_cell_prefix(
    tmp_path: Path,
    mode: str,
    accepted: int,
    in_progress: int,
    code: int,
):
    script = SCRIPT.replace(
        "try:\n    da.run",
        """request=planning.MatrixSpec(
    base=request.with_changes(target=planning.Target(count=1)),
    axes={'x': {'left': planning.Variant(), 'right': planning.Variant()}},
    allocation=planning.Allocation(per_cell=1), max_cells=2)
try:
    da.run""",
    )
    path = tmp_path / "crashed"
    child = subprocess.run(  # noqa: S603 - bounded local fault-injection process
        [sys.executable, "-c", script, str(path), mode],
        capture_output=True,
        text=True,
        timeout=20,
        check=False,
    )
    assert child.returncode == code, child.stderr
    with own_run(path) as connection:
        connection.execute("SELECT count(*) FROM commits").fetchone()
    summary = da.inspect(path, verify=True)
    assert summary.accepted == accepted
    assert summary.counts["in_progress"] == in_progress
    assert summary.cells["x=left"].accepted == accepted
    assert summary.cells["x=right"].counts["started"] == 0
    assert not summary.resumable


def test_portable_verification_checks_inactive_cell_bindings(tmp_path: Path):
    run = da.run(matrix(), out=tmp_path / "run")
    bundle = tmp_path / "bundle"
    da.export(run, all=True, format="bundle", out=bundle)
    path = bundle / "bundle.json"
    value = json.loads(path.read_text())
    cells = value["source_runs"][0]["cells"]
    cells["condition=inactive"]["plan_id"] = cells["condition=first"]["plan_id"]
    value.pop("bundle_id")
    value["bundle_id"] = semantic_digest(value)
    path.write_text(canonical_json(value) + "\n")
    with pytest.raises(ArtifactIntegrityError, match="cell plan"):
        da.inspect(bundle, verify=True)
