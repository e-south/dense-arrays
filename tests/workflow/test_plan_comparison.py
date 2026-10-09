"""Plan differences identify actionable semantic changes without generation.

Author: Eric J. South.
"""

import hashlib
import json
from dataclasses import replace
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays.cli import app


def plans():
    request = planning.DesignSpec(
        parts=[parts.Part("a", "AAA", group="A"), parts.Part("b", "CCC", group="A")],
        length=planning.Length(maximum=6),
        target=planning.Target(count=2),
        limits=planning.Limits(attempts=10),
        requirements=[
            planning.Occurrences("copies", parts.PartSelector(groups=("A",)), min=1),
            planning.Avoid("no-T", patterns=("TTT",)),
        ],
    )
    changed = replace(
        request,
        parts=[parts.Part("a", "AAT", group="A"), parts.Part("c", "GGG", group="A")],
        length=planning.Length(maximum=8),
        target=planning.Target(count=3),
        limits=planning.Limits(attempts=8),
        seed=12,
        requirements=[
            planning.Occurrences("copies", parts.PartSelector(groups=("A",)), min=2),
            planning.Avoid("no-C", patterns=("CCC",)),
        ],
    )
    return da.plan(request), da.plan(changed)


def test_semantic_changes_are_keyed_by_identity_and_preserve_nulls():
    before, after = plans()
    difference = da.inspect(before, view="plan", compare=after)
    changes = {change.path: change for change in difference.changes}
    assert changes[("parts", "a", "sequence")].before == "AAA"
    assert changes[("parts", "a", "sequence")].after == "AAT"
    assert changes[("parts", "b")].kind == "removed"
    assert changes[("parts", "c")].kind == "added"
    assert changes[("requirements", "no-T")].kind == "removed"
    assert changes[("requirements", "no-C")].kind == "added"
    changed = changes[("requirements", "copies", "min")]
    assert (changed.before, changed.after) == (1, 2)
    assert changes[("target", "count")].after == 3
    assert changes[("limits", "attempts")].after == 8
    assert changes[("seed",)].after == 12
    assert difference.to_dict()["schema"] == "dense_arrays.plan_comparison.v2"
    assert all(
        isinstance(change, reporting.PlanChange) for change in difference.changes
    )
    with pytest.raises(TypeError):
        changes[("parts", "c")].after["sequence"] = "TTT"
    assert "AAT" not in repr(difference)


def test_reordering_does_not_appear_as_changed_part_or_rule_content():
    before, _ = plans()
    after = da.plan(
        replace(
            before.request,
            parts=tuple(reversed(before.request.parts)),
            requirements=tuple(reversed(before.request.requirements)),
        )
    )
    difference = da.inspect(before, view="plan", compare=after)
    assert difference.changed_fields == ("part_order", "requirement_order")
    assert "parts" in difference.unchanged_fields
    assert "requirements" in difference.unchanged_fields
    assert {c.path for c in difference.changes} == {
        ("part_order",),
        ("requirement_order",),
    }


def test_saved_evidence_inspection_export_and_comparison_match_cli(tmp_path: Path):
    before, after = plans()
    evidence = tmp_path / "before.evidence.json"
    evidence.write_text(json.dumps(before.evidence.to_dict()))
    revised = tmp_path / "after.plan.json"
    da.export(after, out=revised)
    result = da.inspect(before.evidence, view="plan", compare=after)
    assert result.to_dict()["locations"] == {"before": None, "after": []}
    same = da.inspect(before, view="plan", compare=before.evidence)
    assert same.changes == ()
    assert same.before.plan_id == same.after.plan_id
    assert da.inspect(evidence, view="plan") == before.evidence
    exported = tmp_path / "exported.evidence.json"
    da.export(before.evidence, out=exported)
    assert json.loads(exported.read_text()) == before.evidence.to_dict()
    response = CliRunner().invoke(
        app,
        [
            "inspect",
            str(evidence),
            "--view",
            "plan",
            "--compare",
            str(revised),
            "--json",
        ],
    )
    assert response.exit_code == 0, response.output
    assert json.loads(response.stdout) == result.to_dict()
    human = CliRunner().invoke(
        app, ["inspect", str(evidence), "--view", "plan", "--compare", str(revised)]
    )
    assert human.exit_code == 0, human.output
    assert "/requirements/copies/min" in human.stdout
    assert "1 -> 2" in human.stdout
    with pytest.raises(reporting.ReadLimitError):
        da.inspect(
            before.evidence,
            view="plan",
            compare=after,
            read_limits=reporting.ReadLimits(records=1),
        )


def test_input_locations_are_separate_and_comparison_does_not_solve(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    first = tmp_path / "first.csv"
    second = tmp_path / "second.csv"
    first.write_text("part_id,sequence\na,AAA\nb,CCC\n")
    second.write_bytes(first.read_bytes())
    request = planning.DesignSpec(
        parts=parts.PartTable(first, "csv"), length=planning.Length(maximum=6)
    )
    before = da.plan(request)
    after = da.plan(replace(request, parts=parts.PartTable(second, "csv")))

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("comparison attempted generation")

    monkeypatch.setattr(da.Optimizer, "solve", forbidden)
    monkeypatch.setattr(da.Optimizer, "solve_report", forbidden)
    result = da.inspect(before, view="plan", compare=after)
    assert result.changes == ()
    assert result.to_dict()["locations"] == {
        "before": [str(first)],
        "after": [str(second)],
    }
    second.write_text("part_id,sequence\na,AAT\nb,CCC\n")
    changed = da.plan(replace(request, parts=parts.PartTable(second, "csv")))
    difference = da.inspect(before, view="plan", compare=changed)
    digests = {c.path[1]: c.kind for c in difference.changes if c.path[0] == "inputs"}
    assert digests == {
        before.inputs[0].sha256: "removed",
        changed.inputs[0].sha256: "added",
    }


def test_extension_exclusions_are_named_and_cannot_be_lost_through_evidence(
    tmp_path: Path,
):
    before, _ = plans()
    run = da.run(before, out=tmp_path / "run")
    database = run.path / "run.sqlite3"
    digest = hashlib.sha256(database.read_bytes()).hexdigest()
    child = da.plan(
        planning.ExtensionSpec(
            parent=planning.ParentRun(run.path),
            additional=1,
            limits=planning.Limits(attempts=10),
            seed=23,
        )
    )
    comparison = da.inspect(before.evidence, view="plan", compare=child.evidence)
    exclusions = [c for c in comparison.changes if c.path[0] == "exclusions"]
    assert len(exclusions) == da.inspect(run).accepted
    assert all(
        c.kind == "added" and c.path[1] == c.after["sequence_id"] for c in exclusions
    )
    assert hashlib.sha256(database.read_bytes()).hexdigest() == digest
    editable = da.inspect(child.evidence, view="request").request
    assert da.plan(editable).exclusions == child.exclusions
    da.export(child.evidence, view="request", out=tmp_path / "request.json")
    with pytest.raises(TypeError):
        da.run(child.evidence, out=tmp_path / "not-executable")
    assert not (tmp_path / "not-executable").exists()


def test_change_values_distinguish_added_null_from_changed_null():
    added = reporting.PlanChange(("metadata", "a~/b"), "added", None, None)
    changed = reporting.PlanChange(("length", "exact"), "changed", None, 12)
    assert added.pointer == "/metadata/a~0~1b"
    assert added.to_dict()["kind"] == "added"
    assert changed.to_dict()["kind"] == "changed"
    assert changed.to_dict()["before"] is None
    with pytest.raises(ValueError, match="added"):
        reporting.PlanChange(("parts",), "added", {"a": "AAA"}, {"b": "CCC"})
    with pytest.raises(ValueError, match="removed"):
        reporting.PlanChange(("parts",), "removed", {}, {})
    with pytest.raises(ValueError, match="changed"):
        reporting.PlanChange(("parts",), "changed", {}, {})


def test_comparison_reserves_retained_plan_state_before_reading_second_source(
    tmp_path: Path,
):
    before, _ = plans()
    with pytest.raises(reporting.ReadLimitError, match="identities"):
        da.inspect(
            before,
            view="plan",
            compare=tmp_path / "unopened.json",
            read_limits=reporting.ReadLimits(identities=4),
        )


def test_metadata_comparison_preserves_json_types_and_signed_zero():
    request = planning.DesignSpec(
        parts=[
            parts.Part(
                "a", "AAA", metadata={"flag": True, "weight": 1, "zero": -0.0, "": None}
            )
        ],
        length=planning.Length(maximum=3),
    )
    after = replace(
        request,
        parts=[
            parts.Part(
                "a", "AAA", metadata={"flag": 1, "weight": 1.0, "zero": 0.0, "": 0}
            )
        ],
    )
    comparison = da.inspect(da.plan(request), view="plan", compare=da.plan(after))
    assert {change.path[-1] for change in comparison.changes} == {
        "flag",
        "weight",
        "zero",
        "",
    }
