"""Faults cannot publish incoherent native evidence or silently weaken reading.

Author: Eric J. South.
"""

import json
import sqlite3
from pathlib import Path

import pytest

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.artifacts import Design
from dense_arrays.artifacts.store import create_run
from dense_arrays.realized import Orientation, Placement, PlacementKind, RealizedArray


def plan() -> planning.GenerationPlan:
    """Use a one-part, fully bounded fixture."""
    return da.plan(
        planning.DesignSpec(
            parts=[parts.Part("a", "AAA")], length=planning.Length(maximum=3)
        )
    )


def test_transaction_failure_preserves_committed_prefix(tmp_path: Path):
    resolved = plan()
    out = tmp_path / "run"
    with create_run(resolved, out) as writer:
        attempt = writer.reserve(active_seconds=0.01)
        design_id = "d00000001"
        reference = f"{writer.handle.run_id}/default/{design_id}"
        realized = RealizedArray(
            reference,
            "AAA",
            (Placement("p1", "a", PlacementKind.OTHER, "AAA", 0, Orientation.FORWARD),),
        )
        design = Design(
            writer.handle.run_id,
            "default",
            design_id,
            resolved.plan_id,
            attempt,
            realized,
            (),
        )
        writer.connection.executescript("""
            CREATE TRIGGER reject_commit BEFORE INSERT ON commits
            WHEN NEW.revision=2 BEGIN SELECT RAISE(ABORT, 'injected commit fault'); END;
        """)
        with pytest.raises(sqlite3.IntegrityError, match="injected commit fault"):
            writer.publish(
                attempt,
                "accepted",
                {"solver_status": "optimal", "proof_scope": "offered_packing_model"},
                active_seconds=0.1,
                design=design,
            )
        report = da.inspect(writer.handle, verify=True)
        assert report.accepted == 0
        assert report.counts["in_progress"] == 1
        assert report.counts["started"] == 1
        assert writer.state["counts"]["accepted"] == 0
        assert writer.connection.execute("SELECT count(*) FROM designs").fetchone() == (
            0,
        )


def test_semantic_corruption_fails_even_when_checksum_is_recomputed(tmp_path: Path):
    result = da.run(plan(), out=tmp_path / "run")
    connection = sqlite3.connect(result.path / "run.sqlite3")
    try:
        record = json.loads(
            connection.execute("SELECT payload FROM designs").fetchone()[0]
        )
        record["realized"]["placements"][0]["feature_id"] = "unknown"
        with connection:
            connection.execute(
                "UPDATE designs SET payload=?,digest=?",
                (canonical_json(record), semantic_digest(record)),
            )
    finally:
        connection.close()
    with pytest.raises(ValueError, match="unknown part"):
        da.inspect(result, verify=True)


def test_unknown_attempt_fields_are_not_silently_dropped(tmp_path: Path):
    result = da.run(plan(), out=tmp_path / "run")
    connection = sqlite3.connect(result.path / "run.sqlite3")
    try:
        revision, payload = connection.execute(
            "SELECT revision,payload FROM attempts ORDER BY revision DESC LIMIT 1"
        ).fetchone()
        record = json.loads(payload)
        record["future_field"] = True
        with connection:
            connection.execute(
                "UPDATE attempts SET payload=?,digest=? WHERE revision=?",
                (canonical_json(record), semantic_digest(record), revision),
            )
    finally:
        connection.close()
    with (
        pytest.raises(ValueError, match="unknown"),
        da.inspect(result, view="attempts").records() as records,
    ):
        next(records)


def test_changed_planned_input_does_not_create_output(tmp_path: Path):
    source = tmp_path / "parts.csv"
    source.write_text("part_id,sequence\na,AAA\n")
    resolved = da.plan(
        planning.DesignSpec(
            parts=parts.PartTable(source, "csv"), length=planning.Length(maximum=3)
        )
    )
    source.write_text("part_id,sequence\na,CCC\n")
    out = tmp_path / "run"
    with pytest.raises(ValueError, match="changed"):
        da.run(resolved, out=out)
    assert not out.exists()


def test_summary_does_not_scan_records(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    result = da.run(plan(), out=tmp_path / "run")

    def forbidden(*_args: object) -> None:
        pytest.fail("summary scanned record tables")

    monkeypatch.setattr("dense_arrays.reporting.readers.reader", forbidden)
    assert da.inspect(result).accepted == 1
