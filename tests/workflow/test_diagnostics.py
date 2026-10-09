"""Persist enough final-screen evidence to explain join-created rejection.

Author: Eric J. South.
"""

import json
from pathlib import Path

import pytest
from ortools.linear_solver import pywraplp
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays.artifacts import Attempt
from dense_arrays.cli import app
from dense_arrays.generation.screening import evaluate_screen
from dense_arrays.planning import Avoid
from dense_arrays.realized import Orientation, Placement, PlacementKind, RealizedArray


def test_literal_violation_names_intersecting_part_and_padding_intervals():
    realized = RealizedArray(
        "candidate",
        "ACGAA",
        (Placement("p1", "a", PlacementKind.OTHER, "ACG", 0, Orientation.FORWARD),),
        provenance={"assembly": {"packed_start": 0, "packed_length": 3}},
    )
    evidence = evaluate_screen(
        Avoid("avoid-GAA", patterns=("GAA",), strands="forward"), realized
    )
    match = evidence["observed"][0]
    assert (match["start"], match["end"], match["strand"]) == (2, 5, "forward")
    assert match["intersections"] == [
        {"kind": "part", "part_id": "a", "placement_id": "p1", "start": 0, "end": 3},
        {"kind": "padding", "side": "right", "start": 3, "end": 5},
    ]


def test_persisted_diagnostics_explain_a_real_cbc_rejection_without_solving(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    request = planning.DesignSpec(
        parts=[parts.Part("a", "AA"), parts.Part("b", "CC")],
        length=planning.Length(maximum=4),
        strands="single",
        target=planning.Target(count=2),
        limits=planning.Limits(attempts=1),
        requirements=[
            planning.Fixed("first", "a", "forward", planning.StartWindow(max=0)),
            planning.Fixed("second", "b", "forward"),
            planning.Avoid("no-join", patterns=("AC",), strands="forward"),
        ],
    )
    run = da.run(request, out=tmp_path / "run")

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("diagnostics attempted to create a solver")

    monkeypatch.setattr(pywraplp.Solver, "CreateSolver", forbidden)
    report = da.inspect(run, view="diagnostics", limit=1)
    assert isinstance(report, reporting.DiagnosticReport)
    assert report.cost.records_estimate == 2
    assert report.shortfall == 2
    assert report.diagnostics[0].code == "limit_reached"
    assert report.omitted == 1
    complete = da.inspect(run, view="diagnostics", limit=10)
    assert complete.attempt_counts["rejected"] == 1
    violation = complete.diagnostics[1]
    assert violation.code == "requirement_failed"
    assert violation.requirement_id == "no-join"
    assert violation.expected["patterns"] == ("AC",)
    assert violation.observed["matches"][0]["start"] == 1
    assert [
        i["part_id"] for i in violation.observed["matches"][0]["intersections"]
    ] == ["a", "b"]
    assert complete.proof_scopes == {"offered_packing_model": 1}
    result = CliRunner().invoke(
        app,
        ["inspect", str(run.path), "--view", "diagnostics", "--limit", "10", "--json"],
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout) == complete.to_dict()
    assert "Read cost:" in result.stderr


def test_diagnostic_display_limit_does_not_truncate_reason_totals(tmp_path: Path):
    request = planning.DesignSpec(
        parts=[parts.Part("a", "AAA")],
        length=planning.Length(maximum=3),
        strands="single",
        limits=planning.Limits(attempts=1),
        requirements=[
            planning.Avoid("no-A", patterns=("A",), strands="forward"),
            planning.GC("all-gc", scope="sequence", min=1, max=1),
        ],
    )
    run = da.run(request, out=tmp_path / "run")
    report = da.inspect(run, view="diagnostics", limit=1)
    assert report.reason_counts == {"screening_rejection": 1, "requirement_failed": 2}
    assert report.omitted == 2
    bounded = da.inspect(
        run, view="diagnostics", read_limits=reporting.ReadLimits(records=1)
    )
    with pytest.raises(reporting.ReadLimitError, match="records"):
        bounded.to_dict()


def test_diagnostic_wire_records_validate_identity_and_unknown_fields():
    diagnostic = reporting.Diagnostic(
        "limit_reached",
        "warning",
        "generation",
        {"attempts": 10},
        {"target": 12},
        ("run/revision/21",),
        "Review the persisted attempts.",
        proof_scope="observed_search",
    )
    assert reporting.Diagnostic.from_dict(diagnostic.to_dict()) == diagnostic
    invalid = diagnostic.to_dict()
    invalid["schema"] = "dense_arrays.diagnostic.v99"
    with pytest.raises(ValueError, match="schema"):
        reporting.Diagnostic.from_dict(invalid)
    invalid = diagnostic.to_dict()
    invalid["invented"] = 1
    with pytest.raises(ValueError, match="unknown"):
        reporting.Diagnostic.from_dict(invalid)
    invalid = diagnostic.to_dict()
    invalid["evidence_refs"] = "run/revision/21"
    with pytest.raises(TypeError, match="evidence_refs"):
        reporting.Diagnostic.from_dict(invalid)


@pytest.mark.parametrize("passed", ["false", 1, None])
def test_attempt_evidence_rejects_nonboolean_requirement_results(passed: object):
    with pytest.raises(TypeError, match="passed"):
        Attempt(
            1,
            "default",
            "rejected",
            {
                "requirements": [
                    {"id": "no-join", "observed": [], "passed": passed},
                ]
            },
        )
