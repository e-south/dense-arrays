"""Pure previews bind curated inputs without allocating a packing model.

Author: Eric J. South.
"""

from pathlib import Path

import pytest
from ortools.linear_solver import pywraplp

import dense_arrays as da
from dense_arrays import parts, planning


def test_planning_is_side_effect_free_and_preserves_occurrence_identity(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    def forbidden(*_args: object) -> None:
        pytest.fail("planning allocated a solver")

    monkeypatch.setattr(pywraplp.Solver, "CreateSolver", forbidden)
    monkeypatch.chdir(tmp_path)
    request = planning.DesignSpec(
        parts=[parts.Part("a", "AAA", "A"), parts.Part("b", "AAA", "A")],
        length=planning.Length(maximum=6),
        requirements=[
            planning.Occurrences(
                id="both",
                select=parts.PartSelector(groups=("A",)),
                min=2,
            )
        ],
        target=planning.Target(count=3),
    )
    plan = da.plan(request)
    assert len(plan.request.parts) == 2
    assert plan.request.target.count == 3
    assert plan.preview["oriented_nodes"] == 4
    assert plan.preview["path_variables"] == 20
    assert list(tmp_path.iterdir()) == []


def test_equivalent_requests_have_path_independent_plan_ids(tmp_path: Path):
    first = tmp_path / "first.csv"
    second = tmp_path / "second.csv"
    first.write_text("part_id,sequence,group\na,AAA,A\nb,CCC,B\n")
    second.write_bytes(first.read_bytes())
    plans = [
        da.plan(
            planning.DesignSpec(
                parts=parts.PartTable(path, "csv"),
                length=planning.Length(maximum=6),
            )
        )
        for path in (first, second)
    ]
    assert plans[0].plan_id == plans[1].plan_id
    assert plans[0].inputs[0].path != plans[1].inputs[0].path
    first.write_text("changed\n")
    with pytest.raises(ValueError, match="changed"):
        plans[0].verify_inputs()


def test_native_plan_roundtrip_binds_defaults_and_rejects_tampering():
    plan = da.plan(
        planning.DesignSpec(
            parts=[parts.Part("a", "AAA")],
            length=planning.Length(maximum=3),
        )
    )
    encoded = plan.to_dict()
    assert planning.GenerationPlan.from_dict(encoded) == plan
    encoded["request"]["target"]["count"] = 9
    with pytest.raises(ValueError, match="digest"):
        planning.GenerationPlan.from_dict(encoded)


@pytest.mark.parametrize("minimum", [True, 1.5, -1, 3])
def test_static_count_contradictions_fail_before_execution(minimum: object):
    with pytest.raises((ValueError, TypeError), match=r"minimum|integer|available"):
        da.plan(
            planning.DesignSpec(
                parts=[parts.Part("a", "AAA", "A")],
                length=planning.Length(maximum=3),
                requirements=[
                    planning.Occurrences(
                        id="count",
                        select=parts.PartSelector(groups=("A",)),
                        min=minimum,
                    )
                ],
            )
        )


def test_unavailable_assembly_fails_planning():
    with pytest.raises(ValueError, match=r"exact.*assembly"):
        da.plan(
            planning.DesignSpec(
                parts=[parts.Part("a", "AAA")],
                length=planning.Length(exact=6),
            )
        )


@pytest.mark.parametrize(
    "rule",
    [
        planning.Occurrences("count-A", parts.PartSelector(groups=("A",)), min=2),
        planning.Occurrences("count-missing", parts.PartSelector(groups=("B",)), min=1),
        planning.GroupCoverage("coverage", ("A", "B"), min=1),
    ],
)
def test_resolved_count_errors_identify_the_requirement(rule: object):
    with pytest.raises(ValueError, match=rf"^{rule.id}:"):
        da.plan(
            planning.DesignSpec(
                parts=[parts.Part("a", "AAA", "A")],
                length=planning.Length(maximum=3),
                requirements=[rule],
            )
        )
