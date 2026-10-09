"""Matrix comparisons retain named choices, exact cell meaning and source locations.

Author: Eric J. South.
"""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays.cli import app


def matrix():
    return planning.MatrixSpec(
        planning.DesignSpec(
            [parts.Part("a", "AAA", group="sites")],
            planning.Length(maximum=6),
            requirements=[
                planning.Occurrences(
                    "copies", parts.PartSelector(groups=("sites",)), min=1
                )
            ],
            strands="single",
        ),
        axes={"pool": {"A": planning.Variant(), "B": planning.Variant()}},
        sources={"pool=B": [parts.Part("b", "CCC", group="sites")]},
        allocation=planning.Allocation(per_cell=1),
        max_cells=4,
    )


def test_matrix_comparison_identifies_resolved_parts_rules_targets_and_declarations():
    request = matrix()
    before = da.plan(request)
    after = da.plan(
        request.with_changes(
            sources={"pool=B": [parts.Part("b", "GGG", group="sites")]},
            allocation=planning.Allocation(counts={"pool=A": 2, "pool=B": 0}),
            axes={
                "pool": {
                    "A": planning.Variant(
                        requirements=[
                            planning.Occurrences(
                                "copies",
                                parts.PartSelector(groups=("sites",)),
                                min=1,
                                max=1,
                            )
                        ]
                    ),
                    "B": planning.Variant(),
                }
            },
        )
    )
    result = da.inspect(before, view="plan", compare=after)
    assert isinstance(result, reporting.PlanComparison)
    changes = {c.path: c for c in result.changes}
    assert changes[("cells", "pool=B", "parts", "b", "sequence")].before == "CCC"
    assert changes[("cells", "pool=B", "parts", "b", "sequence")].after == "GGG"
    assert changes[("cells", "pool=A", "target", "count")].after == 2
    assert changes[("cells", "pool=B", "target", "count")].after == 0
    assert changes[("cells", "pool=A", "requirements", "copies", "max")].after == 1
    assert "matrix" in result.changed_fields
    assert "base" in result.unchanged_fields
    assert result.to_dict()["schema"] == "dense_arrays.matrix_plan_comparison.v1"


def test_reordering_named_choices_is_distinct_from_editing_cell_content():
    request = matrix()
    before = da.plan(request)
    after = da.plan(
        request.with_changes(
            axes={"pool": dict(reversed(request.axes["pool"].items()))}
        )
    )
    difference = da.inspect(before, view="plan", compare=after)
    assert "cells" in difference.unchanged_fields
    assert {c.path for c in difference.changes} == {
        ("cell_order",),
        ("matrix", "choice_order", "pool"),
    }
    assert {c.cell_id: c.plan.plan_id for c in before.cells} == {
        c.cell_id: c.plan.plan_id for c in after.cells
    }


def test_added_and_removed_combinations_are_named_whole_records():
    request = matrix()
    before = da.plan(request)
    after = da.plan(
        request.with_changes(
            axes={"pool": {"A": planning.Variant(), "C": planning.Variant()}},
            sources={"pool=C": [parts.Part("c", "GGG", group="sites")]},
        )
    )
    changes = {
        c.path: c for c in da.inspect(before, view="plan", compare=after).changes
    }
    assert changes[("cells", "pool=B")].kind == "removed"
    assert changes[("cells", "pool=C")].kind == "added"
    assert changes[("cells", "pool=C")].after["parts"]["c"]["sequence"] == "GGG"
    assert not any(c.path[:2] == ("cells", "pool=A") for c in changes.values())
    with pytest.raises(TypeError):
        changes[("cells", "pool=C")].after["target"]["count"] = 9


def test_saved_matrix_comparison_is_source_free_and_matches_cli_and_export(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence,group\nb,CCC,sites\n")
    request = matrix().with_changes(sources={"pool=B": parts.PartTable(table, "csv")})
    before = da.plan(request)
    table.write_text("part_id,sequence,group\nb,GGG,sites\n")
    after = da.plan(request)
    before.write(tmp_path / "before.json")
    after.write(tmp_path / "after.json")
    table.unlink()

    def forbidden(*_a: object, **_kw: object) -> None:
        pytest.fail("comparison invoked a solver")

    monkeypatch.setattr(da.Optimizer, "solve_report", forbidden)
    expected = da.inspect(before, view="plan", compare=after).to_dict()
    runner = CliRunner()
    arguments = [
        "inspect",
        str(tmp_path / "before.json"),
        "--view",
        "plan",
        "--compare",
        str(tmp_path / "after.json"),
    ]
    response = runner.invoke(app, [*arguments, "--json"])
    assert response.exit_code == 0, response.output
    assert json.loads(response.stdout) == expected
    human = runner.invoke(app, arguments)
    assert human.exit_code == 0, human.output
    assert '/cells/pool=B/parts/b/sequence: "CCC" -> "GGG"' in human.stdout
    output = tmp_path / "comparison.json"
    exported = runner.invoke(
        app,
        [
            "export",
            str(tmp_path / "before.json"),
            "--view",
            "plan",
            "--compare",
            str(tmp_path / "after.json"),
            "--out",
            str(output),
        ],
    )
    assert exported.exit_code == 0, exported.output
    assert json.loads(output.read_text()) == expected
    small = tmp_path / "too-small.json"
    with pytest.raises(reporting.ReadLimitError):
        da.export(
            da.inspect(before, view="plan", compare=after),
            out=small,
            read_limits=reporting.ReadLimits(identities=1),
        )
    assert not small.exists()


def test_matrix_locations_do_not_change_design_meaning(tmp_path: Path):
    first, second = tmp_path / "first.csv", tmp_path / "second.csv"
    first.write_text("part_id,sequence,group\nb,CCC,sites\n")
    second.write_bytes(first.read_bytes())
    request = matrix()
    before = da.plan(
        request.with_changes(sources={"pool=B": parts.PartTable(first, "csv")})
    )
    after = da.plan(
        request.with_changes(sources={"pool=B": parts.PartTable(second, "csv")})
    )
    first.unlink()
    second.unlink()
    result = da.inspect(before, view="plan", compare=after)
    assert result.changes == ()
    assert before.plan_id == after.plan_id
    assert result.to_dict()["locations"] == {
        "before": {"base": [], "cells": {"pool=A": [], "pool=B": [str(first)]}},
        "after": {"base": [], "cells": {"pool=A": [], "pool=B": [str(second)]}},
    }


def test_sampling_and_extension_changes_keep_cell_and_ancestor_scope(tmp_path: Path):
    before = da.plan(matrix())
    sampled = da.plan(
        matrix().with_changes(
            batches={
                "pool=B": planning.Resampling(
                    planning.BatchSampling(size=1, seed=2),
                    max_batches=3,
                    attempts_per_batch=2,
                )
            }
        )
    )
    changes = {
        c.path: c for c in da.inspect(before, view="plan", compare=sampled).changes
    }
    assert changes[("cells", "pool=B", "resampling")].kind == "added"
    assert changes[("matrix", "batches", "pool=B")].after["max_batches"] == 3
    parent = da.run(before, out=tmp_path / "parent")
    extension = da.plan(
        planning.ExtensionSpec(
            planning.ParentRun(parent.path),
            {"pool=A": 2, "pool=B": 0},
            planning.Limits(attempts=20),
            6,
        )
    )
    comparison = da.inspect(parent, view="plan", compare=extension)
    exclusions = [
        c
        for c in comparison.changes
        if c.path[0] == "cells" and c.path[2] == "exclusions"
    ]
    assert len(exclusions) == 2
    assert {c.path[1] for c in exclusions} == {"pool=A", "pool=B"}
    assert all(
        c.kind == "added" and c.after["cell_id"] == c.path[1] for c in exclusions
    )


def test_comparison_requires_matching_plan_kinds_and_bounds_both_sources(
    tmp_path: Path,
):
    before = da.plan(matrix())
    with pytest.raises(TypeError, match="matching plan kinds"):
        da.inspect(before, view="plan", compare=before.base)
    with pytest.raises(TypeError, match="matching plan kinds"):
        reporting.PlanComparison(before.base, before)
    with pytest.raises(reporting.ReadLimitError, match="two plan records"):
        da.inspect(
            before,
            view="plan",
            compare=before,
            read_limits=reporting.ReadLimits(records=1),
        )
    with pytest.raises(reporting.ReadLimitError, match="identities"):
        da.inspect(
            before,
            view="plan",
            compare=tmp_path / "unopened.json",
            read_limits=reporting.ReadLimits(identities=6),
        )


def test_new_optional_policies_are_in_changed_field_summary():
    before = da.plan(matrix().base)
    after = da.plan(
        before.request.with_changes(
            resampling=planning.Resampling(
                planning.BatchSampling(size=1, seed=2),
                max_batches=3,
                attempts_per_batch=2,
            )
        )
    )
    result = da.inspect(before, view="plan", compare=after)
    assert result.changed_fields == ("resampling",)
    assert result.changes[0].path == ("resampling",)


@pytest.mark.parametrize("kind", ["matrix", "generation", "evidence", "preparation"])
def test_saved_plan_caps_precede_typed_materialization(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, kind: str
):
    plan = da.plan(matrix())
    if kind == "generation":
        plan = plan.base
    elif kind == "evidence":
        plan = plan.base.evidence
    elif kind == "preparation":
        table = tmp_path / "parts.csv"
        table.write_text("part_id,sequence\na,AAA\nb,CCC\n")
        plan = da.plan(parts.PreparationSpec(parts.PartTable(table, "csv")))
    saved = tmp_path / "plan.json"
    saved.write_text(json.dumps(plan.to_dict()))

    def materialized(*_a: object, **_kw: object) -> None:
        pytest.fail("typed plan was constructed before checking the read cap")

    monkeypatch.setattr(type(plan), "from_dict", materialized)
    with pytest.raises(reporting.ReadLimitError, match="identities"):
        da.inspect(saved, view="plan", read_limits=reporting.ReadLimits(identities=1))


@pytest.mark.parametrize("content", [None, [], False])
def test_plan_read_admission_rejects_nonobject_evidence(
    tmp_path: Path, content: object
):
    path = tmp_path / "invalid.json"
    path.write_text(
        json.dumps({"schema": "dense_arrays.plan_evidence.v1", "content": content})
    )
    with pytest.raises(TypeError, match="plan must be an object"):
        da.inspect(path, view="plan")
