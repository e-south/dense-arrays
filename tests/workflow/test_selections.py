"""Saved selections preserve quotas, full identities and source revisions.

Author: Eric J. South.
"""

import io
import json
import shutil
import sqlite3
from dataclasses import replace
from pathlib import Path

import pytest
from PIL import Image
from typer.testing import CliRunner

import dense_arrays as da
import dense_arrays.reporting.selections.materialization as selection_module
from dense_arrays import parts, planning, reporting
from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.artifacts.store import create_run
from dense_arrays.cli import app


def library(path: Path):
    return da.run(
        planning.DesignSpec(
            parts=[parts.Part("a", "AAA", group="A"), parts.Part("b", "CCC")],
            length=planning.Length(maximum=6),
            strands="single",
            target=planning.Target(count=2),
        ),
        out=path,
    )


def test_first_total_and_per_cell_have_distinct_scope_and_no_redistribution(
    tmp_path: Path,
):
    left, right = library(tmp_path / "left"), library(tmp_path / "right")
    designs = list(da.inspect([left, right], view="designs", all=True).records())
    request = reporting.LibrarySelection(take=reporting.Take(count=2))
    total = da.inspect([left, right, left], view="selection", select=request)
    assert isinstance(total, reporting.SelectionSnapshot)
    assert list(total.references()) == [d.reference for d in designs[:2]]
    assert (total.requested, total.available, total.selected, total.status) == (
        2,
        4,
        2,
        "complete",
    )
    quotas = {f"{left.run_id}/default": 1, f"{right.run_id}/default": 0}
    per_cell = da.inspect(
        [left, right],
        view="selection",
        select=reporting.LibrarySelection(take=reporting.Take(per_cell=quotas)),
    )
    assert list(per_cell.references()) == [designs[0].reference]
    assert per_cell.counts[f"{right.run_id}/default"] == {
        "requested": 0,
        "available": 2,
        "selected": 0,
        "shortfall": 0,
    }
    assert len(repr(total)) < 200
    assert designs[0].reference not in repr(total)
    assert reporting.SelectionSnapshot.from_dict(total.to_dict()) == total


def test_random_selection_is_repeatable_without_replacement_and_shortfalls_are_explicit(
    tmp_path: Path,
):
    left, right = library(tmp_path / "left"), library(tmp_path / "right")
    request = reporting.LibrarySelection(
        take=reporting.Take(count=3, policy="random", seed=23)
    )
    first = da.inspect([left, right], view="selection", select=request)
    second = da.inspect([left, right], view="selection", select=request)
    assert first.to_dict() == second.to_dict()
    assert len(set(first.references())) == 3
    assert first.algorithm == "sha256_priority.v1"
    available = list(da.inspect([left, right], view="designs", all=True).records())
    assert list(first.references()) == [
        d.reference for d in available if d.reference in set(first.references())
    ]
    quota = {f"{left.run_id}/default": 3, f"{right.run_id}/default": 1}
    with pytest.raises(
        reporting.SelectionShortfall, match=r"requested 4.*available 4.*selected 3"
    ) as error:
        da.inspect(
            [left, right],
            view="selection",
            select=reporting.LibrarySelection(take=reporting.Take(per_cell=quota)),
        )
    assert error.value.counts[f"{left.run_id}/default"]["shortfall"] == 1
    partial = da.inspect(
        [left, right],
        view="selection",
        select=reporting.LibrarySelection(
            take=reporting.Take(per_cell=quota, shortfall="allow_partial")
        ),
    )
    assert (partial.status, partial.selected, partial.shortfall) == ("partial", 3, 1)
    for bad in ({"default": 0}, {"missing": 0}):
        with pytest.raises(ValueError, match=r"ambiguous|unknown"):
            da.inspect(
                [left, right],
                view="selection",
                select=reporting.LibrarySelection(take=reporting.Take(per_cell=bad)),
            )
    with pytest.raises(reporting.ReadLimitError):
        da.inspect(
            [left, right],
            view="selection",
            select=request,
            read_limits=reporting.ReadLimits(records=1),
        )
    with pytest.raises(ValueError, match="limit"):
        da.inspect(left, view="selection", select=request, limit=1)


def test_saved_membership_reuses_revision_and_exports_without_resampling(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    reference = library(tmp_path / "reference")
    candidates = list(da.inspect(reference, view="designs", all=True).records())
    evidence = {"solver_status": "optimal", "proof_scope": "offered_packing_model"}
    with create_run(da.inspect(reference, view="plan"), tmp_path / "active") as writer:
        for candidate in candidates[:1]:
            ordinal = writer.reserve(active_seconds=0)
            writer.publish(
                ordinal,
                "accepted",
                evidence,
                active_seconds=0,
                design=replace(candidate, run_id=writer.handle.run_id),
            )
        snapshot = da.inspect(
            writer.handle,
            view="selection",
            select=reporting.LibrarySelection(
                take=reporting.Take(count=1, policy="random", seed=23)
            ),
        )
        saved = tmp_path / "panel.selection.json"
        receipt = da.export(snapshot, format="selection", out=saved)
        assert receipt.design_refs == tuple(snapshot.references())
        ordinal = writer.reserve(active_seconds=0)
        writer.publish(
            ordinal,
            "accepted",
            evidence,
            active_seconds=0,
            design=replace(candidates[1], run_id=writer.handle.run_id),
        )

        def forbidden(*_args: object, **_kwargs: object) -> None:
            pytest.fail("saved selections must not sample again")

        monkeypatch.setattr(selection_module, "materialize", forbidden)
        rows = list(
            da.inspect(
                writer.handle, view="designs", select=snapshot, all=True
            ).records()
        )
        assert [r.reference for r in rows] == list(snapshot.references())
        output = io.StringIO()
        exported = da.export(
            writer.handle, view="sequences", select=snapshot, format="fasta", out=output
        )
        assert exported.design_refs == tuple(snapshot.references())
        assert exported.sources[0]["revision"] == snapshot.sources[0].revision
        cli = CliRunner().invoke(
            app,
            [
                "export",
                str(writer.handle.path),
                "--view",
                "sequences",
                "--selection",
                str(saved),
                "--format",
                "fasta",
                "--out",
                "-",
            ],
        )
        assert cli.exit_code == 0, cli.output
        assert cli.stdout == output.getvalue()
        with pytest.raises(ValueError, match="all"):
            da.export(writer.handle, select=snapshot, all=True, out=io.StringIO())


def test_partial_exports_are_qualified_and_default_shortfalls_publish_nothing(
    tmp_path: Path,
):
    run = library(tmp_path / "run")
    request = reporting.LibrarySelection(take=reporting.Take(count=3))
    out = tmp_path / "panel.json"
    with pytest.raises(reporting.SelectionShortfall):
        da.export(run, select=request, format="selection", out=out)
    assert not out.exists()
    selection = tmp_path / "selection.json"
    selection.write_text(
        json.dumps(
            reporting.LibrarySelection(
                take=reporting.Take(count=3, shortfall="allow_partial")
            ).to_dict()
        )
    )
    cli = CliRunner().invoke(
        app,
        [
            "export",
            str(run.path),
            "--selection",
            str(selection),
            "--format",
            "selection",
            "--out",
            str(out),
            "--json",
        ],
    )
    assert cli.exit_code == 3, cli.output
    receipt = json.loads(cli.stdout)
    assert receipt["selection"]["status"] == "partial"
    assert receipt["selection"]["shortfall"] == 1
    assert len(json.loads(out.read_text())["members"]) == 2
    inspected = CliRunner().invoke(
        app,
        [
            "inspect",
            str(run.path),
            "--view",
            "selection",
            "--selection",
            str(selection),
            "--json",
        ],
    )
    assert inspected.exit_code == 0, inspected.output
    assert json.loads(inspected.stdout)["counts"]["total"]["shortfall"] == 1


def test_snapshot_quality_and_render_use_saved_members_with_source_denominators(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):

    left, right = library(tmp_path / "left"), library(tmp_path / "right")
    snapshot = da.inspect(
        [left, right],
        view="selection",
        select=reporting.LibrarySelection(take=reporting.Take(count=1)),
    )

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("rendering must not select again")

    monkeypatch.setattr("dense_arrays.workflow.selections.materialize", forbidden)
    quality = da.inspect([left, right], view="quality", select=snapshot)
    report = quality.to_dict()
    assert report["selection"]["designs"] == 1
    assert report["selection"]["snapshot"]["snapshot_id"] == snapshot.snapshot_id
    assert [s["attainment"]["accepted"] for s in report["source_runs"]] == [2, 2]
    image = tmp_path / "panel.png"
    rendered = da.render(
        [left, right], view="library-quality", select=snapshot, out=image
    )
    assert rendered.selection["snapshot_id"] == snapshot.snapshot_id
    with Image.open(image) as png:
        assert any(snapshot.snapshot_id in str(v) for v in png.info.values())
    exported = da.export(
        [left, right], view="quality", select=snapshot, out=tmp_path / "quality.json"
    )
    assert exported.selection["snapshot_id"] == snapshot.snapshot_id


def test_selected_bundle_moves_and_pinned_content_changes_fail_before_publication(
    tmp_path: Path,
):
    left, right = library(tmp_path / "left"), library(tmp_path / "right")
    snapshot = da.inspect(
        [left, right],
        view="selection",
        select=reporting.LibrarySelection(
            take=reporting.Take(count=3, policy="random", seed=23)
        ),
    )
    bundle = tmp_path / "bundle"
    receipt = da.export([left, right], select=snapshot, format="bundle", out=bundle)
    assert receipt.design_refs == tuple(snapshot.references())
    relocated = tmp_path / "relocated"
    shutil.move(str(bundle), relocated)
    shutil.rmtree(left.path)
    shutil.rmtree(right.path)
    assert da.inspect(relocated, verify=True).verified
    assert [
        d.reference for d in da.inspect(relocated, view="designs", all=True).records()
    ] == list(snapshot.references())
    run = library(tmp_path / "run")
    pinned = da.inspect(
        run,
        view="selection",
        select=reporting.LibrarySelection(take=reporting.Take(count=1)),
    )
    with sqlite3.connect(run.path / "run.sqlite3") as connection:
        payload = connection.execute(
            "SELECT payload FROM designs WHERE ordinal=1"
        ).fetchone()[0]
        value = json.loads(payload)
        value["realized"]["provenance"]["note"] = "changed evidence"
        connection.execute(
            "UPDATE designs SET payload=?,digest=? WHERE ordinal=1",
            (canonical_json(value), semantic_digest(value)),
        )
    out = tmp_path / "changed.fasta"
    with pytest.raises(ValueError, match="content"):
        da.export(run, select=pinned, view="sequences", format="fasta", out=out)
    assert not out.exists()


def test_four_cell_panel_cli_parity_exact_quotas_and_qualified_23_design_shortfall(
    tmp_path: Path,
):
    runs = [
        da.run(
            planning.DesignSpec(
                parts=[
                    parts.Part(label, sequence)
                    for label, sequence in zip(
                        "abcd", ("AAA", "CCC", "GGG", "TTT"), strict=True
                    )
                ],
                length=planning.Length(maximum=12),
                strands="single",
                target=planning.Target(count=5 if name == "short" else 10),
            ),
            out=tmp_path / name,
        )
        for name in ("a", "b", "c", "d", "short")
    ]
    sources = runs[:4]
    request = reporting.LibrarySelection(
        take=reporting.Take(
            per_cell={f"{r.run_id}/default": 6 for r in sources},
            policy="random",
            seed=23,
        )
    )
    snapshot = da.inspect(sources, view="selection", select=request)
    assert snapshot.selected == 24
    assert all(
        c["selected"] == 6 and c["available"] == 10 for c in snapshot.counts.values()
    )
    file = tmp_path / "request.json"
    file.write_text(json.dumps(request.to_dict()))
    response = CliRunner().invoke(
        app,
        [
            "inspect",
            *(str(r.path) for r in sources),
            "--view",
            "selection",
            "--selection",
            str(file),
            "--json",
        ],
    )
    assert response.exit_code == 0, response.output
    assert [m["reference"] for m in json.loads(response.stdout)["members"]] == list(
        snapshot.references()
    )
    first = da.inspect(
        sources,
        view="selection",
        select=reporting.LibrarySelection(take=reporting.Take(count=24)),
    )
    all_designs = list(da.inspect(sources, view="designs", all=True).records())
    assert list(first.references()) == [d.reference for d in all_designs[:24]]
    partial_sources = [*sources[:3], runs[-1]]
    request = reporting.LibrarySelection(
        take=reporting.Take(
            per_cell={f"{r.run_id}/default": 6 for r in partial_sources},
            policy="random",
            seed=23,
            shortfall="allow_partial",
        )
    )
    partial = da.inspect(partial_sources, view="selection", select=request)
    assert (partial.selected, partial.shortfall) == (23, 1)
    output = tmp_path / "panel.tsv"
    receipt = da.export(
        partial_sources, view="placements", select=partial, format="tsv", out=output
    )
    assert receipt.selection["shortfall"] == 1
    assert len(receipt.design_refs) == 23
    with da.inspect(
        sources, view="placements", select=snapshot, limit=3
    ).records() as page:
        prefix = list(page)
        cursor = page.next_cursor
    suffix = list(
        da.inspect(
            sources, view="placements", select=snapshot, after=cursor, all=True
        ).records()
    )
    assert len(prefix + suffix) == 96


def test_selection_scope_limits_and_request_evidence_are_strict(tmp_path: Path):
    run = library(tmp_path / "run")
    snapshot = da.inspect(
        run,
        view="selection",
        select=reporting.LibrarySelection(take=reporting.Take(count=2)),
    )
    for all_rows in ("yes", 1):
        with pytest.raises(TypeError, match="boolean"):
            da.export(run, select=snapshot, all=all_rows, out=io.StringIO())
    with pytest.raises(reporting.ReadLimitError):
        da.export(
            snapshot,
            format="selection",
            read_limits=reporting.ReadLimits(identities=1),
            out=io.StringIO(),
        )
    value = snapshot.to_dict()
    value["request"]["take"]["count"] = 999
    # Even a recomputed outer digest cannot make false allocation evidence valid.
    content = {k: v for k, v in value.items() if k != "snapshot_id"}
    content["sources"] = [
        {k: v for k, v in s.items() if k != "path"} for s in content["sources"]
    ]
    value["snapshot_id"] = semantic_digest(content)
    with pytest.raises(ValueError, match="allocation"):
        reporting.SelectionSnapshot.from_dict(value)


def test_bundle_selection_empty_membership_and_cli_failure_diagnostics(tmp_path: Path):
    run = library(tmp_path / "run")
    bundle = tmp_path / "bundle"
    da.export(run, all=True, format="bundle", out=bundle)
    panel = da.inspect(
        bundle,
        view="selection",
        select=reporting.LibrarySelection(take=reporting.Take(count=1)),
    )
    view = da.inspect(bundle, view="placements", all=True, select=panel)
    with view.records() as rows:
        assert len(list(rows)) == 2
        assert rows.examined <= view.cost.records_estimate
    zero = da.inspect(
        bundle,
        view="selection",
        select=reporting.LibrarySelection(take=reporting.Take(per_cell={})),
    )
    assert zero.selected == zero.requested == 0
    out = tmp_path / "empty"
    da.export(bundle, select=zero, format="bundle", out=out)
    assert da.inspect(out, verify=True).designs == 0
    request = tmp_path / "shortfall.json"
    request.write_text(
        json.dumps(reporting.LibrarySelection(take=reporting.Take(count=3)).to_dict())
    )
    result = CliRunner().invoke(
        app,
        [
            "export",
            str(run.path),
            "--selection",
            str(request),
            "--format",
            "selection",
            "--out",
            str(tmp_path / "absent.json"),
            "--json",
        ],
    )
    assert result.exit_code == 2
    diagnostic = json.loads(result.stdout)
    assert diagnostic["code"] == "selection_shortfall"
    assert diagnostic["counts"]["total"]["shortfall"] == 1
    assert not (tmp_path / "absent.json").exists()
    native_zero = da.inspect(
        run,
        view="selection",
        select=reporting.LibrarySelection(take=reporting.Take(count=0)),
    )
    with sqlite3.connect(run.path / "run.sqlite3") as connection:
        connection.execute(
            "DELETE FROM commits WHERE revision=?", (native_zero.sources[0].revision,)
        )
    with pytest.raises(ValueError, match="no committed record"):
        da.export(run, select=native_zero, out=tmp_path / "lost.json")
    assert not (tmp_path / "lost.json").exists()
    malformed = panel.to_dict()
    malformed.pop("members")
    with pytest.raises(ValueError, match="missing required"):
        reporting.SelectionSnapshot.from_dict(malformed)


def test_bundle_rejects_selection_summary_that_disagrees_with_its_policy(
    tmp_path: Path,
):
    run = library(tmp_path / "run")
    bundle = tmp_path / "bundle"
    da.export(
        run,
        select=reporting.LibrarySelection(take=reporting.Take(count=1)),
        format="bundle",
        out=bundle,
    )
    manifest = bundle / "bundle.json"
    value = json.loads(manifest.read_text())
    value["selection"]["request"]["take"]["count"] = 99
    value.pop("bundle_id")
    value["bundle_id"] = semantic_digest(value)
    manifest.write_text(canonical_json(value) + "\n")
    with pytest.raises(ValueError, match="allocation"):
        da.inspect(bundle)
