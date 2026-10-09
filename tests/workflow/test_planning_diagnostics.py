"""Impossible count requests explain their evidence before allocating a run.

Module Author(s): Eric J. South
"""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays.cli import app


@pytest.mark.parametrize("table_input", [False, True])
def test_insufficient_parts_share_requested_counts_and_evidence(
    tmp_path: Path, table_input: bool
):
    records = (
        parts.Part("a", "ACGTTGCAAGTCCTGA", group="A"),
        parts.Part("other", "ATGCTTAGGACGTTCA", group="B"),
        parts.Part("b", "AGTCCTGATCGTACCG", group="A"),
    )
    source = records
    table = tmp_path / "parts.csv"
    if table_input:
        table.write_text(
            "part_id,sequence,group\na,ACGTTGCAAGTCCTGA,A\n"
            "other,ATGCTTAGGACGTTCA,B\nb,AGTCCTGATCGTACCG,A\n"
        )
        source = parts.PartTable(table, "csv")
    request = planning.DesignSpec(
        parts=source,
        length=planning.Length(maximum=40),
        requirements=(
            planning.Occurrences("three-A", parts.PartSelector(groups=("A",)), min=3),
        ),
    )
    with pytest.raises(planning.PlanningError, match="three-A") as caught:
        da.run(request, out=tmp_path / "python-run")
    error = caught.value
    diagnostic = error.diagnostic
    assert isinstance(error, planning.PlanningError)
    assert isinstance(diagnostic, reporting.Diagnostic)
    assert diagnostic.code == "insufficient_parts"
    assert diagnostic.stage == "planning"
    assert diagnostic.severity == "error"
    assert diagnostic.requirement_id == "three-A"
    assert diagnostic.observed == {"available": 2}
    assert diagnostic.expected == {"minimum": 3}
    assert diagnostic.evidence_refs[:2] == ("parts/0", "parts/2")
    if table_input:
        assert diagnostic.evidence_refs[2:] == (
            f"{table.as_uri()}#row=1",
            f"{table.as_uri()}#row=3",
        )
    else:
        assert len(diagnostic.evidence_refs) == 2
    assert "3" in str(error)
    assert "2" in str(error)
    assert "minimum" in diagnostic.next_action
    assert not (tmp_path / "python-run").exists()
    recipe = tmp_path / "request.json"
    da.export(request, view="request", out=recipe)
    for command in ("plan", "run"):
        out = tmp_path / command
        result = CliRunner().invoke(
            app, [command, str(recipe), "--out", str(out), "--json"]
        )
        assert result.exit_code == 2, result.output
        payload = json.loads(result.stdout)
        assert payload["code"] == "insufficient_parts"
        assert payload["diagnostic"] == diagnostic.to_dict()
        assert not out.exists()
