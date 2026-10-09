"""Native exports bind the committed manifests read during publication.

Author: Eric J. South.
"""

import hashlib
import json
import sqlite3
from pathlib import Path

import pytest

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays.artifacts import RunHandle
from dense_arrays.artifacts.store import create_run
from dense_arrays.reporting.exporting import records as exporting


def digest(value: dict) -> str:
    """Calculate the specified wire fingerprint independently of package helpers."""
    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=False
        ).encode()
    ).hexdigest()


def run(path: Path):
    return da.run(
        planning.DesignSpec(
            [parts.Part("a", "AAA")], planning.Length(maximum=3), strands="single"
        ),
        out=path,
    )


def manifest(path: Path, kind: str, revision: int = 0) -> dict:
    if kind == "bundle":
        return json.loads((path / "bundle.json").read_text())
    with sqlite3.connect(path / f"{kind}.sqlite3") as connection:
        sql, params = (
            ("SELECT payload FROM commits WHERE revision=?", (revision,))
            if kind == "run"
            else ("SELECT payload FROM manifest WHERE id=1", ())
        )
        return json.loads(connection.execute(sql, params).fetchone()[0])


@pytest.mark.parametrize("kind", ["run", "pool", "bundle", "union", "selection"])
def test_record_exports_bind_each_canonical_source_manifest(tmp_path: Path, kind: str):
    first = run(tmp_path / "run")
    if kind == "pool":
        table = tmp_path / "parts.csv"
        table.write_text("part_id,sequence\na,AAA\n")
        source = da.prepare(
            parts.PreparationSpec(parts.PartTable(table, "csv")), out=tmp_path / "pool"
        )
        view, sources = "parts", [(source.path, "pool", 0)]
    elif kind == "bundle":
        source = tmp_path / "bundle"
        da.export(first, all=True, format="bundle", out=source)
        view, sources = "designs", [(source, "bundle", 0)]
    else:
        source = [first, run(tmp_path / "other"), first] if kind == "union" else first
        handles = source if isinstance(source, list) else [source]
        sources = [(h.path, "run", da.inspect(h).revision) for h in handles]
        view = "designs"
        if kind == "selection":
            source = da.inspect(
                source,
                view="selection",
                select=reporting.LibrarySelection(take=reporting.Take(count=1)),
            )
    expected = [digest(manifest(p, k, r)) for p, k, r in sources]
    before = [manifest(p, k, r) for p, k, r in sources]
    output = tmp_path / "export.json"
    receipt = da.export(
        source, view=view, all=kind != "selection", format="json", out=output
    )
    bindings = receipt.to_dict()["sources"]
    assert [s["manifest_digest"] for s in bindings] == expected
    assert json.loads(output.read_text())["sources"] == bindings
    assert [manifest(p, k, r) for p, k, r in sources] == before
    if kind == "selection":
        assert [s.manifest_digest for s in source.sources] == expected
    assert (
        receipt.to_dict()["files"][0]["sha256"]
        == hashlib.sha256(output.read_bytes()).hexdigest()
    )


def test_export_rechecks_manifest_before_publishing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    source = run(tmp_path / "run")
    original = exporting._write_records  # noqa: SLF001 - inject after completed stream

    def changed(*args: object, **kwargs: object) -> tuple[int, dict[str, None]]:
        result = original(*args, **kwargs)
        with sqlite3.connect(source.path / "run.sqlite3") as connection:
            revision = da.inspect(source).revision
            value = manifest(source.path, "run", revision)
            value["active_seconds"] += 1
            connection.execute(
                "UPDATE commits SET payload=?,digest=? WHERE revision=?",
                (
                    json.dumps(value, sort_keys=True, separators=(",", ":")),
                    digest(value),
                    revision,
                ),
            )
        return result

    monkeypatch.setattr(exporting, "_write_records", changed)
    output = tmp_path / "must-not-publish.json"
    with pytest.raises(ValueError, match=r"manifest.*changed|source.*changed"):
        da.export(source, all=True, format="json", out=output)
    assert not output.exists()


def test_later_commit_keeps_pinned_export_digest(tmp_path: Path):
    plan = da.plan(
        planning.DesignSpec([parts.Part("a", "AAA")], planning.Length(maximum=3))
    )
    with create_run(plan, tmp_path / "active") as writer:
        attempt = writer.reserve(active_seconds=0)
        writer.publish(
            attempt, "rejected", {"code": "screening_rejection"}, active_seconds=0
        )
        revision = da.inspect(writer.handle).revision
        pinned = RunHandle(writer.handle.path, writer.handle.run_id, revision)
        expected = digest(manifest(writer.handle.path, "run", revision))
        attempt = writer.reserve(active_seconds=0)
        writer.publish(
            attempt, "rejected", {"code": "screening_rejection"}, active_seconds=0
        )
        output = tmp_path / "attempts.json"
        receipt = da.export(pinned, view="attempts", all=True, out=output)
        assert receipt.sources[0]["manifest_digest"] == expected
        assert receipt.sources[0]["revision"] == revision
        assert receipt.records == 1


def test_manifest_binding_does_not_consume_projected_record_allowance(tmp_path: Path):
    source = run(tmp_path / "run")
    query = da.inspect(
        source,
        view="sequences",
        all=True,
        read_limits=reporting.ReadLimits(records=1, identities=1),
    )
    assert query.cost.records_estimate == 1
    receipt = da.export(query, all=True, format="json", out=tmp_path / "one.json")
    assert receipt.records == 1
    assert len(receipt.sources) == 1
    assert len(receipt.sources[0]["manifest_digest"]) == 64
    with pytest.raises(reporting.ReadLimitError, match="identities"):
        da.export(
            [source, source],
            view="sequences",
            all=True,
            read_limits=reporting.ReadLimits(identities=1),
            out=tmp_path / "bounded.json",
        )
    assert not (tmp_path / "bounded.json").exists()


def test_saved_selection_digest_is_not_rebound_after_source_change(tmp_path: Path):
    source = run(tmp_path / "run")
    selection = da.inspect(
        source,
        view="selection",
        select=reporting.LibrarySelection(take=reporting.Take(count=1)),
    )
    revision = selection.sources[0].revision
    value = manifest(source.path, "run", revision)
    value["active_seconds"] += 1
    with sqlite3.connect(source.path / "run.sqlite3") as connection:
        connection.execute(
            "UPDATE commits SET payload=?,digest=? WHERE revision=?",
            (
                json.dumps(value, sort_keys=True, separators=(",", ":")),
                digest(value),
                revision,
            ),
        )
    with pytest.raises(ValueError, match=r"source.*changed"):
        da.export(selection, out=tmp_path / "selection.json")
    assert not (tmp_path / "selection.json").exists()


@pytest.mark.parametrize("kind", ["pool", "bundle"])
def test_immutable_source_change_during_export_does_not_publish(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, kind: str
):
    if kind == "bundle":
        source = tmp_path / "bundle"
        da.export(run(tmp_path / "run"), all=True, format="bundle", out=source)
        view = "designs"
    else:
        table = tmp_path / "parts.csv"
        table.write_text("part_id,sequence\na,AAA\n")
        source = da.prepare(
            parts.PreparationSpec(parts.PartTable(table, "csv")), out=tmp_path / "pool"
        ).path
        view = "parts"
    original = exporting._write_records  # noqa: SLF001 - inject after completed stream

    def changed(*args: object, **kwargs: object) -> tuple[int, dict[str, None]]:
        result = original(*args, **kwargs)
        value = manifest(source, kind)
        if kind == "bundle":
            value["source_runs"][0]["active_seconds"] += 1
            value.pop("bundle_id")
            value["bundle_id"] = digest(value)
            (source / "bundle.json").write_text(
                json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n"
            )
        else:
            value["producer"]["package_version"] = "changed"
            with sqlite3.connect(source / "pool.sqlite3") as connection:
                connection.execute(
                    "UPDATE manifest SET payload=?,digest=? WHERE id=1",
                    (
                        json.dumps(value, sort_keys=True, separators=(",", ":")),
                        digest(value),
                    ),
                )
        return result

    monkeypatch.setattr(exporting, "_write_records", changed)
    with pytest.raises(ValueError, match=r"source.*changed"):
        da.export(source, view=view, all=True, out=tmp_path / "output.json")
    assert not (tmp_path / "output.json").exists()
