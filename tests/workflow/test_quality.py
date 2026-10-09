"""Quality metrics describe complete persisted populations with explicit denominators.

Author: Eric J. South.
"""

import importlib
import json
from pathlib import Path

import pytest
from ortools.linear_solver import pywraplp
from PIL import Image
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays.artifacts import Design
from dense_arrays.artifacts.store import create_run
from dense_arrays.cli import app
from dense_arrays.generation.acceptance import evaluate
from dense_arrays.playback.quality import quality_figure
from dense_arrays.realized import Orientation, Placement, PlacementKind, RealizedArray


def shortfall_run(tmp_path: Path):
    sequences = ("ACGTT", "CGTAC", "GTACG", "TACGT", "ACCGT", "CCGTA", "CGTTA", "GTTAC")
    supplied = [
        parts.Part(f"part-{i}", dna, group="A" if i < 4 else "B")
        for i, dna in enumerate(sequences)
    ]
    supplied += [
        parts.Part("unused-a", "AAA", group="C"),
        parts.Part("unused-b", "AA", group="C"),
    ]
    plan = da.plan(
        planning.DesignSpec(
            parts=supplied,
            length=planning.Length(maximum=5),
            strands="single",
            target=planning.Target(count=12),
            limits=planning.Limits(attempts=10),
            requirements=[
                planning.Avoid("no-four-A", patterns=("AAAA",), strands="forward")
            ],
        )
    )
    with create_run(plan, tmp_path / "run") as writer:
        for i, sequence in enumerate((*sequences, sequences[0])):
            attempt = writer.reserve(active_seconds=0)
            design_id = f"d{i + 1}"
            reference = f"{writer.handle.run_id}/default/{design_id}"
            realized = RealizedArray(
                reference,
                sequence,
                (
                    Placement(
                        "p1",
                        f"part-{i % 8}",
                        PlacementKind.OTHER,
                        sequence,
                        0,
                        Orientation.FORWARD,
                    ),
                ),
            )
            writer.publish(
                attempt,
                "accepted",
                {"solver_status": "optimal", "proof_scope": "offered_packing_model"},
                active_seconds=0,
                design=Design(
                    writer.handle.run_id,
                    "default",
                    design_id,
                    plan.plan_id,
                    attempt,
                    realized,
                    evaluate(realized, plan),
                ),
            )
        attempt = writer.reserve(active_seconds=0)
        rejected = RealizedArray(
            "rejected",
            "AAAAA",
            (
                Placement(
                    "p1", "unused-a", PlacementKind.OTHER, "AAA", 0, Orientation.FORWARD
                ),
                Placement(
                    "p2", "unused-b", PlacementKind.OTHER, "AA", 3, Orientation.FORWARD
                ),
            ),
        )
        writer.publish(
            attempt,
            "rejected",
            {
                "solver_status": "optimal",
                "proof_scope": "offered_packing_model",
                "code": "screening_rejection",
                "requirements": evaluate(rejected, plan),
            },
            active_seconds=0,
        )
        writer.finish("stopped", "attempt_limit", active_seconds=1)
        return writer.handle


def test_quality_reconciles_shortfall_and_complete_population_without_solving(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    run = shortfall_run(tmp_path)
    assert da.inspect(run, verify=True).accepted == 8

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("quality inspection invoked a solver")

    monkeypatch.setattr(pywraplp.Solver, "CreateSolver", forbidden)
    report = da.inspect(run, view="quality", limit=2)
    assert isinstance(report, reporting.QualityReport)
    assert report.cost.records_estimate == 19
    value = report.to_dict()
    assert value["status"] == "exact"
    assert value["population"] == "all_accepted_designs_at_revision"
    assert value["attainment"] == {
        "target": 12,
        "accepted": 8,
        "shortfall": 4,
        "distinct_sequences": 8,
    }
    assert value["supply"] == {
        "eligible_parts": 10,
        "eligible_groups": 3,
        "unused_parts": 2,
        "unused_groups": 1,
    }
    assert value["concentration"]["highest_part_occurrence_share"] == 0.125
    assert value["concentration"]["occurrence_denominator"] == 8
    assert value["composition"]["gc_fraction"]["mean"] == 0.5
    assert value["composition"]["gc_fraction"]["histogram"] == [
        {"value": 0.4, "count": 4},
        {"value": 0.6, "count": 4},
    ]
    assert value["composition"]["density"]["mean"] == 1
    assert value["composition"]["compression"]["mean"] == 1
    assert value["search"]["attempt_counts"]["duplicate"] == 1
    assert value["search"]["attempt_counts"]["rejected"] == 1
    assert value["requirements"] == [
        {"id": "no-four-A", "passed": 8, "failed": 0, "missing": 0, "not_applicable": 0}
    ]
    assert len(value["part_usage"]) == 2
    assert value["part_usage"][0]["occurrences"] == 1
    assert value["part_usage"][0]["design_fraction"] == 0.125
    assert value["part_usage"][0]["group"] == "A"
    assert value["examined"] == 19
    cursor = value["next_cursor"]
    remainder = da.inspect(run, view="quality", limit=20, after=cursor).to_dict()
    assert [p["part_id"] for p in remainder["part_usage"]][-2:] == [
        "unused-a",
        "unused-b",
    ]
    assert len(remainder["part_usage"]) == 8
    assert remainder["attainment"] == value["attainment"]
    assert remainder["next_cursor"] is None
    result = CliRunner().invoke(
        app, ["inspect", str(run.path), "--view", "quality", "--limit", "2", "--json"]
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout) == value
    assert "Read cost:" in result.stderr


def test_quality_density_uses_interval_union_and_pre_padding_compression(
    tmp_path: Path,
):
    run = da.run(
        planning.DesignSpec(
            parts=[parts.Part("a", "AAC"), parts.Part("b", "ACG")],
            length=planning.Length(exact=6),
            strands="single",
            assembly=planning.Assembly(
                padding=planning.Padding(side="right", max_trials=1)
            ),
        ),
        out=tmp_path / "run",
    )
    value = da.inspect(run, view="quality").to_dict()
    composition = value["composition"]
    assert composition["density"]["mean"] == pytest.approx(4 / 6)
    assert composition["compression"]["mean"] == 1.5
    assert composition["packed_span"]["mean"] == 4
    assert composition["padding_length"]["mean"] == 2
    assert composition["placement_count"]["mean"] == 2
    assert value["padding_sides"] == {"right": 1}
    assert value["occupancy"] == [
        {"start": 0, "end": 4, "designs": 1, "denominator": 1, "fraction": 1},
        {"start": 4, "end": 6, "designs": 0, "denominator": 1, "fraction": 0},
    ]


def test_empty_quality_uses_null_denominators_and_caps_fail_explicitly(tmp_path: Path):
    run = da.run(
        planning.DesignSpec(
            parts=[parts.Part("a", "AAA")],
            length=planning.Length(maximum=3),
            limits=planning.Limits(attempts=1),
            requirements=[planning.Avoid("no-A", patterns=("A",), strands="both")],
        ),
        out=tmp_path / "empty",
    )
    value = da.inspect(run, view="quality").to_dict()
    assert value["supply"]["unused_parts"] == 1
    assert value["composition"]["gc_fraction"]["mean"] is None
    assert value["composition"]["gc_fraction"]["reason"] == "empty_population"
    assert value["part_usage"][0]["occurrence_share"] is None
    report = da.inspect(
        run, view="quality", read_limits=reporting.ReadLimits(records=1)
    )
    assert "QualityReport" in repr(report)
    assert report.cost.records_estimate == 2
    with pytest.raises(reporting.ReadLimitError, match="records"):
        report.to_dict()
    with pytest.raises(reporting.ReadLimitError, match="identities"):
        da.inspect(
            run, view="quality", read_limits=reporting.ReadLimits(identities=1)
        ).to_dict()


def test_quality_plot_publishes_report_evidence_without_regeneration(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    run = shortfall_run(tmp_path)
    report = da.inspect(run, view="quality").to_dict()
    before = (run.path / "run.sqlite3").read_bytes()

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("quality rendering invoked generation or acceptance")

    monkeypatch.setattr(da.Optimizer, "solve_report", forbidden)
    monkeypatch.setattr("dense_arrays.generation.acceptance.evaluate", forbidden)
    out = tmp_path / "quality.png"
    receipt = da.render(run, view="library-quality", out=out)
    assert receipt.records == 8
    assert receipt.design_refs == ()
    with Image.open(out) as image:
        assert image.width >= 1000
        assert json.loads(image.info["DenseArraysReport"]) == report
    assert (run.path / "run.sqlite3").read_bytes() == before
    response = CliRunner().invoke(
        app,
        [
            "render",
            str(run.path),
            "--view",
            "library-quality",
            "--out",
            str(tmp_path / "cli.png"),
            "--json",
        ],
    )
    assert response.exit_code == 0, response.output
    assert json.loads(response.stdout)["sources"] == receipt.to_dict()["sources"]
    assert "Read cost:" in response.stderr


def test_quality_plot_axes_use_exact_report_counts(tmp_path: Path):
    report = da.inspect(shortfall_run(tmp_path), view="quality").to_dict()
    figure = quality_figure(report)
    axes = {ax.get_label(): ax for ax in figure.axes}
    assert [bar.get_width() for bar in axes["part_usage"].patches] == [1] * 8 + [0] * 2
    for metric, expected in (("gc_fraction", [4, 4]), ("density", [8])):
        axis = axes[metric]
        # One drawing object represents every exact bin as cardinality grows.
        assert len(axis.patches) + len(axis.collections) == 1
        paths = axis.collections[0].get_paths()
        assert [max(path.vertices[:, 1]) for path in paths] == expected
        centers = [
            (min(path.vertices[:, 0]) + max(path.vertices[:, 0])) / 2 for path in paths
        ]
        assert centers == pytest.approx(
            [row["value"] for row in report["composition"][metric]["histogram"]]
        )
    assert [bar.get_width() for bar in axes["outcomes"].patches] == [8, 1, 1]
    figure.clear()


def test_quality_render_work_limit_leaves_destination_uncreated(tmp_path: Path):
    run = shortfall_run(tmp_path)
    out = tmp_path / "uncreated" / "quality.png"
    with pytest.raises(reporting.ReadLimitError, match="records"):
        da.render(
            run,
            view="library-quality",
            out=out,
            read_limits=reporting.ReadLimits(records=1),
        )
    assert not out.parent.exists()
    response = CliRunner().invoke(
        app,
        [
            "render",
            str(run.path),
            "--view",
            "library-quality",
            "--out",
            str(out),
            "--max-read-records",
            "1",
            "--json",
        ],
    )
    assert response.exit_code == 4, response.output
    assert json.loads(response.stdout)["code"] == "read_limit"
    assert not out.parent.exists()


def test_missing_quality_renderer_fails_before_scanning_or_creating_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    run = shortfall_run(tmp_path)
    real_import = importlib.import_module

    def unavailable(name: str) -> object:
        if name == "matplotlib.figure":
            raise ModuleNotFoundError(name)
        return real_import(name)

    monkeypatch.setattr("dense_arrays.reporting.rendering.import_module", unavailable)
    out = tmp_path / "uncreated" / "quality.png"
    with pytest.raises(ValueError, match="playback dependencies"):
        da.render(
            run,
            view="library-quality",
            out=out,
            read_limits=reporting.ReadLimits(records=1),
        )
    assert not out.parent.exists()


def test_background_placements_are_counted_without_implying_motif_hits(tmp_path: Path):
    run = da.run(
        planning.DesignSpec(
            parts=(
                parts.Part("a", "AAA", group="background"),
                parts.Part("b", "CCC", group="background"),
            ),
            length=planning.Length(maximum=6),
            strands="single",
        ),
        out=tmp_path / "background",
    )
    report = da.inspect(run, view="quality").to_dict()
    assert report["composition"]["placement_count"]["mean"] == 2
    assert "motif_count" not in report["composition"]
    selected = reporting.DesignFilter(
        metrics={"placement_count": reporting.Range(min=2, max=2)}
    )
    assert reporting.DesignFilter.from_dict(selected.to_dict()) == selected
    with da.inspect(run, view="designs", select=selected).records() as rows:
        assert len(list(rows)) == 1
    with pytest.raises(ValueError, match="unavailable design metrics"):
        reporting.DesignFilter(metrics={"motif_count": reporting.Range(min=1)})
