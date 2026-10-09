"""Input evidence survives normalization and persists independently of source files.

Author: Eric J. South.
"""

from pathlib import Path

import pytest

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.artifacts.store import reader, stored_plan


def test_import_report_survives_planning_execution_and_source_removal(tmp_path: Path):
    source = tmp_path / "curated.csv"
    source.write_text("site,bases,ignored\na, acgt ,note\nb,TTAA,note\n")
    plan = da.plan(
        planning.DesignSpec(
            parts=parts.PartTable(
                source,
                "csv",
                columns={"part_id": "site", "sequence": "bases"},
                normalization=parts.Normalization(
                    uppercase=True, trim_outer_whitespace=True
                ),
            ),
            length=planning.Length(maximum=8),
        )
    )
    report = plan.import_report.to_dict()
    assert report["rows"] == 2
    assert report["ignored_columns"] == ["ignored"]
    assert report["transformations"] == [
        {
            "row": 1,
            "field": "sequence",
            "column": "bases",
            "before": " acgt ",
            "after": "ACGT",
        }
    ]
    restored = planning.GenerationPlan.from_dict(plan.to_dict())
    assert restored.import_report.to_dict() == report
    result = da.run(restored, out=tmp_path / "run")
    source.unlink()
    with reader(result.path) as connection:
        assert stored_plan(connection).import_report.to_dict() == report
    assert da.inspect(result, verify=True).accepted == 1


def test_display_and_repeated_identity_access_do_not_serialize_parts(
    monkeypatch: pytest.MonkeyPatch,
):
    plan = da.plan(
        planning.DesignSpec(
            parts=[parts.Part(f"p{i}", "AAA") for i in range(100)],
            length=planning.Length(maximum=6),
        )
    )
    identity = plan.plan_id

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("display or cached identity rescanned the entire request")

    monkeypatch.setattr("dense_arrays.planning.resolution.request_to_dict", forbidden)
    assert plan.plan_id == identity
    shown = repr(plan)
    assert len(shown) < 300
    assert "100 parts" in shown
    assert "AAA" not in shown


def test_import_report_rejects_unknown_fields():
    plan = da.plan(
        planning.DesignSpec(
            parts=[parts.Part("a", "AAA")], length=planning.Length(maximum=3)
        )
    )
    value = plan.to_dict()
    value["import_report"]["future_policy"] = True
    with pytest.raises(ValueError, match="unknown"):
        planning.GenerationPlan.from_dict(value)
