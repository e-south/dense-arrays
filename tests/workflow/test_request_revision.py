"""Editable requests preserve declared intent without executing a revision.

Author: Eric J. South.
"""

import json
import shutil
from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays.cli import app
from dense_arrays.workflow.inputs import read_source


def test_revision_retains_parent_without_implying_sequence_exclusion(tmp_path: Path):
    original = planning.DesignSpec(
        parts=(parts.Part("a", "AAA"),),
        length=planning.Length(maximum=3),
        strands="single",
    )
    run = da.run(original, out=tmp_path / "parent")
    report = da.inspect(run, view="request")
    assert isinstance(report, reporting.RequestReport)
    request = report.request
    summary = da.inspect(run)
    assert request.lineage.parent.run_id == run.run_id
    assert request.lineage.parent.plan_id == summary.plan_id
    assert request.lineage.parent.revision == summary.revision
    revised = request.with_changes(length=planning.Length(maximum=4))
    assert request.length.maximum == 3
    assert revised.length.maximum == 4
    assert revised.lineage == request.lineage
    with pytest.raises(FrozenInstanceError):
        request.seed = 99
    with pytest.raises(TypeError, match="seed"):
        request.with_changes(seed=True)
    with pytest.raises(TypeError, match="unexpected"):
        request.with_changes(typo=1)
    output = tmp_path / "handoff" / "revision.json"
    da.export(report, out=output)
    assert read_source(output) == request
    assert json.loads(output.read_text())["schema"] == "dense_arrays.design.v1"
    cli = CliRunner().invoke(
        app, ["inspect", str(run.path), "--view", "request", "--json"]
    )
    assert cli.exit_code == 0, cli.output
    assert json.loads(cli.stdout) == json.loads(output.read_text())
    # Parent lineage is descriptive: repeating the parent's only sequence is valid.
    child = da.run(revised, out=tmp_path / "revision")
    assert da.inspect(child, verify=True).accepted == 1
    with da.inspect(child, view="designs").records() as records:
        assert next(records).realized.sequence == "AAA"
    assert da.inspect(run) == summary


@pytest.mark.parametrize("kind", ["design", "preparation"])
def test_resolved_plan_write_is_portable_and_create_only(tmp_path: Path, kind: str):
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence\na,AAA\nb,CCC\n")
    source = parts.PartTable(table, "csv")
    request = (
        planning.DesignSpec(source, planning.Length(maximum=6))
        if kind == "design"
        else parts.PreparationSpec(source)
    )
    plan = da.plan(request)
    path = tmp_path / "plans" / f"{kind}.json"
    plan.write(path)
    value = json.loads(path.read_text())
    bindings = value["inputs"] if kind == "design" else [value["input"]]
    assert bindings[0]["path"] == "../parts.csv"
    restored = read_source(path)
    assert restored.plan_id == plan.plan_id
    restored.verify_inputs()
    export = path.with_name("export.json")
    da.export(plan, out=export)
    assert path.read_bytes() == export.read_bytes()
    before = path.read_bytes()
    with pytest.raises(FileExistsError):
        plan.write(path)
    assert path.read_bytes() == before


def test_preparation_request_revision_does_not_read_or_retain_mutable_inputs(
    tmp_path: Path,
):
    request = parts.PreparationSpec(parts.PartTable(tmp_path / "missing.csv", "csv"))
    revised = request.with_changes(
        retain=parts.Retention(select=parts.PartFilter(groups=("A",)))
    )
    assert request.retain.select is None
    assert revised.retain.select.groups == ("A",)
    with pytest.raises(TypeError, match="retain"):
        revised.with_changes(retain={})


@pytest.mark.parametrize("source_kind", ["table", "pool"])
def test_revision_preserves_bound_parts_after_request_relocation(
    tmp_path: Path, source_kind: str
):
    folder = tmp_path / "inputs"
    folder.mkdir()
    table = folder / "parts.csv"
    table.write_text("part_id,sequence,group\na,aaa,A\nb,ccc,B\n")
    source = parts.PartTable(
        table, "csv", normalization=parts.Normalization(uppercase=True)
    )
    if source_kind == "pool":
        pool = da.prepare(parts.PreparationSpec(source), out=folder / "pool")
        source = parts.PoolSource(pool)
    original = da.plan(
        planning.DesignSpec(source, planning.Length(maximum=6), strands="single")
    )
    report = da.inspect(original, view="request")
    assert isinstance(report.request.parts, parts.BoundParts)
    revised = report.request.with_changes(target=planning.Target(count=2), seed=23)
    resolved = da.plan(revised)
    assert resolved.collection_id == original.collection_id
    assert resolved.import_report == original.import_report
    assert resolved.inputs == original.inputs
    assert da.inspect(original, view="plan", compare=resolved).changed_fields == (
        "target",
        "seed",
    )
    path = folder / "requests" / "revised.json"
    da.export(revised, out=path)
    value = json.loads(path.read_text())
    assert value["parts"]["inputs"][0]["path"].startswith("../")
    moved = tmp_path / "moved"
    folder.rename(moved)
    restored = read_source(moved / "requests" / "revised.json")
    relocated = da.plan(restored)
    assert relocated.plan_id == resolved.plan_id
    assert relocated.collection_id == original.collection_id
    cli_path = tmp_path / "cli.plan.json"
    cli = CliRunner().invoke(
        app, ["plan", str(moved / "requests" / "revised.json"), "--out", str(cli_path)]
    )
    assert cli.exit_code == 0, cli.output
    cli_plan = read_source(cli_path)
    assert cli_plan.plan_id == relocated.plan_id
    assert [(i.path.resolve(), i.sha256) for i in cli_plan.inputs] == [
        (i.path.resolve(), i.sha256) for i in relocated.inputs
    ]
    human = CliRunner().invoke(app, ["inspect", str(cli_path), "--view", "request"])
    assert human.exit_code == 0, human.output
    assert "2 bound parts" in human.stdout
    assert "source files checked" in human.stdout
    assert (
        da.inspect(da.run(relocated, out=tmp_path / "child"), verify=True).accepted == 2
    )


def test_bound_revision_rejects_changed_sources_and_damaged_evidence(tmp_path: Path):
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence\na,AAA\nb,CCC\n")
    original = da.plan(
        planning.DesignSpec(parts.PartTable(table, "csv"), planning.Length(maximum=6))
    )
    report = da.inspect(original, view="request")
    assert isinstance(report.request.parts, parts.BoundParts)
    output = tmp_path / "uncreated" / "request.json"
    with pytest.raises(reporting.ReadLimitError):
        da.export(report, out=output, read_limits=reporting.ReadLimits(identities=1))
    assert not output.parent.exists()
    malformed = report.to_dict()
    malformed["parts"]["parts"][0]["sequence"] = "TTT"
    request_file = tmp_path / "malformed.json"
    request_file.write_text(json.dumps(malformed))
    with pytest.raises(ValueError, match="identity"):
        read_source(request_file)
    frozen = da.plan(report.request)
    table.write_text("part_id,sequence\na,GGG\nb,CCC\n")
    with pytest.raises(ValueError, match="input changed"):
        da.plan(report.request)
    with pytest.raises(ValueError, match="input changed"):
        da.run(frozen, out=tmp_path / "uncreated-run")
    assert not (tmp_path / "uncreated-run").exists()
    # Replacing the complete source starts a new collection deliberately.
    replacement = da.plan(report.request.with_changes(parts=original.request.parts))
    assert replacement.collection_id != original.collection_id
    assert replacement.import_report.kind == "inline"


def test_portable_request_preserves_recorded_origins_without_file_bindings(
    tmp_path: Path,
):
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence\na,AAA\nb,CCC\n")
    original = da.plan(
        planning.DesignSpec(parts.PartTable(table, "csv"), planning.Length(maximum=6))
    )
    run = da.run(original, out=tmp_path / "run")
    bundle = tmp_path / "bundle"
    da.export(run, format="bundle", all=True, out=bundle)
    table.unlink()
    shutil.rmtree(run.path)
    report = da.inspect(bundle, view="request")
    assert isinstance(report.request.parts, parts.BoundParts)
    assert report.request.parts.locations is None
    revised = da.plan(report.request.with_changes(seed=29))
    assert revised.collection_id == original.collection_id
    assert revised.import_report == original.import_report
    assert revised.evidence.input_digests == original.evidence.input_digests
    assert revised.inputs == ()
    request_path = tmp_path / "portable.request.json"
    da.export(report, out=request_path)
    assert json.loads(request_path.read_text())["parts"]["inputs"][0]["path"] is None
    saved_plan = tmp_path / "portable.plan.json"
    revised.write(saved_plan)
    restored = read_source(saved_plan)
    assert restored == revised
    human = CliRunner().invoke(app, ["inspect", str(saved_plan), "--view", "request"])
    assert human.exit_code == 0, human.output
    assert "embedded parts" in human.stdout
    child = da.run(restored, out=tmp_path / "child")
    assert da.inspect(child, verify=True).accepted == 1
    assert da.inspect(child, view="plan").collection_id == original.collection_id


def test_explicit_exclusion_survives_revision_export_and_parent_removal(tmp_path: Path):
    request = planning.DesignSpec(
        (parts.Part("a", "AAA"),),
        planning.Length(maximum=3),
        strands="single",
    )
    parent = da.run(request, out=tmp_path / "parent")
    exclusion = planning.LibraryExclusion(
        source=planning.ParentRun(parent.path),
        cell_mapping={"default": "default"},
        uniqueness="exact_sequence_per_cell.v1",
    )
    revised = da.inspect(parent, view="request").request.with_changes(
        length=planning.Length(maximum=4), exclude=exclusion
    )
    resolved = da.plan(revised)
    assert resolved.preview["excluded_sequences"] == 1
    assert resolved.preview["exclusion_scope"] == {
        "cell_mapping": {"default": "default"},
        "uniqueness": "exact_sequence_per_cell.v1",
    }
    saved = tmp_path / "handoff" / "request.json"
    da.export(resolved, view="request", out=saved)
    shutil.rmtree(parent.path)
    restored = da.plan(read_source(saved))
    assert restored.exclusions == resolved.exclusions
    rejected = da.run(restored, out=tmp_path / "revised")
    assert da.inspect(rejected, verify=True).accepted == 0
    with da.inspect(rejected, view="attempts", all=True).records() as attempts:
        duplicates = [a for a in attempts if a.outcome == "duplicate"]
    assert len(duplicates) == 1
    assert duplicates[0].evidence["matched_design_ref"].startswith(parent.run_id + "/")
    # An extension's inherited exclusions become explicit when editing its rules.
    other = da.run(request, out=tmp_path / "other")
    extension = da.plan(
        planning.ExtensionSpec(planning.ParentRun(other.path), 1, planning.Limits(), 2)
    )
    edited = da.inspect(extension, view="request").request.with_changes(
        length=planning.Length(maximum=4)
    )
    assert edited.exclude is not None
    assert da.plan(edited).exclusions == extension.exclusions
    assert edited.lineage.parent.run_id == other.run_id


def test_request_caps_and_invalid_exclusion_scope_fail_before_output(tmp_path: Path):
    parent = da.run(
        planning.DesignSpec((parts.Part("a", "AAA"),), planning.Length(maximum=3)),
        out=tmp_path / "parent",
    )
    source = planning.ParentRun(parent.path)
    for mapping, uniqueness in (
        ({}, "exact_sequence_per_cell.v1"),
        ({"other": "default"}, "exact_sequence_per_cell.v1"),
        ({"default": "default"}, "sequence"),
    ):
        with pytest.raises(ValueError, match=r"cell_mapping|uniqueness"):
            planning.LibraryExclusion(source, mapping, uniqueness)
    request = da.inspect(parent, view="request").request.with_changes(
        exclude=planning.LibraryExclusion(
            source, {"default": "default"}, "exact_sequence_per_cell.v1"
        )
    )
    with pytest.raises(reporting.ReadLimitError):
        da.plan(request, read_limits=reporting.ReadLimits(records=1))
    resolved = da.plan(request)
    report = da.inspect(resolved, view="request")
    output = tmp_path / "uncreated" / "request.json"
    with pytest.raises(reporting.ReadLimitError, match="identities"):
        da.export(report, out=output, read_limits=reporting.ReadLimits(identities=1))
    assert not output.parent.exists()
    with pytest.raises(reporting.ReadLimitError, match="identities"):
        da.inspect(
            resolved, view="request", read_limits=reporting.ReadLimits(identities=1)
        )


def test_cli_plan_shares_writer_and_discloses_explicit_exclusion_scope(tmp_path: Path):
    parent = da.run(
        planning.DesignSpec((parts.Part("a", "AAA"),), planning.Length(maximum=3)),
        out=tmp_path / "parent",
    )
    request = da.inspect(parent, view="request").request.with_changes(
        exclude=planning.LibraryExclusion(
            planning.ParentRun(parent.path),
            {"default": "default"},
            "exact_sequence_per_cell.v1",
        )
    )
    path = tmp_path / "requests" / "request.json"
    da.export(request, out=path)
    assert json.loads(path.read_text())["exclude"]["source"]["run"] == "../parent"
    plan = da.plan(read_source(path))
    py_output = tmp_path / "plans" / "python.json"
    plan.write(py_output)
    cli_output = py_output.with_name("cli.json")
    result = CliRunner().invoke(app, ["plan", str(path), "--out", str(cli_output)])
    assert result.exit_code == 0, result.output
    assert cli_output.read_bytes() == py_output.read_bytes()
    assert "1 excluded sequence" in result.stdout
    assert "default -> default" in result.stdout
    assert "exact_sequence_per_cell.v1" in result.stdout


def test_revised_library_exclusions_survive_bundles_and_successive_extension(
    tmp_path: Path,
):
    request = planning.DesignSpec(
        tuple(parts.Part(base, base * 3) for base in "ACG"),
        planning.Length(maximum=9),
        strands="single",
    )
    parent = da.run(request, out=tmp_path / "parent")
    revision = da.inspect(parent, view="request").request.with_changes(
        length=planning.Length(maximum=10),
        exclude=planning.LibraryExclusion(
            planning.ParentRun(parent.path),
            {"default": "default"},
            "exact_sequence_per_cell.v1",
        ),
    )
    child = da.run(revision, out=tmp_path / "child")
    assert da.inspect(child, verify=True).accepted == 1
    child_request = da.inspect(child, view="request").request
    assert child_request.lineage.parent.run_id == child.run_id
    assert child_request.exclude.source.run_id == parent.run_id
    inherited = da.plan(
        planning.ExtensionSpec(
            planning.ParentRun(child.path), 1, planning.Limits(attempts=10), 31
        )
    )
    assert len(inherited.exclusions) == 2
    extension = da.run(inherited, out=tmp_path / "extension")
    assert da.inspect(extension, verify=True).accepted == 1
    refs = [
        d.sequence_id
        for d in da.inspect(
            [parent, child, extension], view="designs", all=True
        ).records()
    ]
    assert len(set(refs)) == 3
    bundle = tmp_path / "bundle"
    da.export(child, all=True, format="bundle", out=bundle)
    child_plan = da.inspect(child, view="plan")
    multi = tmp_path / "multi-bundle"
    da.export([parent, child], all=True, format="bundle", out=multi)
    with pytest.raises(ValueError, match="multiple plans"):
        da.export(multi, view="request", out=tmp_path / "ambiguous.json")
    assert not (tmp_path / "ambiguous.json").exists()
    selected = CliRunner().invoke(
        app,
        [
            "inspect",
            str(multi),
            "--view",
            "request",
            "--plan-id",
            child_plan.plan_id,
            "--json",
        ],
    )
    assert selected.exit_code == 0, selected.output
    assert (
        json.loads(selected.stdout) == da.inspect(child_plan, view="request").to_dict()
    )
    shutil.rmtree(parent.path)
    shutil.rmtree(child.path)
    assert da.inspect(bundle, verify=True).designs == 1
    assert da.inspect(bundle, view="quality").to_dict()["selection"]["designs"] == 1
    evidence = da.inspect(bundle, view="plan")
    assert evidence.exclusions == child_plan.exclusions
    assert da.inspect(evidence, view="request").request == child_plan.request
    assert da.inspect(bundle, view="request").request == child_plan.request
    exported_request = tmp_path / "bundle-request.json"
    da.export(bundle, view="request", out=exported_request)
    assert read_source(exported_request) == child_plan.request
    with pytest.raises(reporting.ReadLimitError, match="identities"):
        da.inspect(bundle, view="plan", read_limits=reporting.ReadLimits(identities=3))


def test_lineage_edits_are_visible_and_strict_without_changing_source_plan(
    tmp_path: Path,
):
    plan = da.plan(
        planning.DesignSpec((parts.Part("a", "AAA"),), planning.Length(maximum=3))
    )
    revised = plan.request.with_changes(
        lineage=planning.Lineage(planning.RunReference("origin", plan.plan_id, 2))
    )
    changes = da.inspect(plan, view="plan", compare=da.plan(revised))
    assert changes.changed_fields == ("lineage",)
    assert plan.request.lineage is None
    path = tmp_path / "request.json"
    da.export(revised, out=path)
    value = json.loads(path.read_text())
    value["lineage"]["parent"]["typo"] = 1
    path.write_text(json.dumps(value))
    with pytest.raises(ValueError, match="unknown"):
        read_source(path)
    with pytest.raises(TypeError, match="revision"):
        planning.RunReference("origin", plan.plan_id, revision=True)
