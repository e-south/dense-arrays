"""Runtime feedback and recovery preserve the meaning of committed outcomes.

Author: Eric J. South.
"""

import json
import sqlite3
from pathlib import Path

import pytest

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.artifacts.errors import ArtifactIntegrityError
from dense_arrays.artifacts.store import RunWriter
from dense_arrays.generation.batches import sampling


def recipe(*, attempts: int = 1, size: int = 1) -> planning.DesignSpec:
    return planning.DesignSpec(
        (parts.Part("a", "AAA", group="A"), parts.Part("b", "CCC", group="B")),
        planning.Length(maximum=3),
        target=planning.Target(5),
        strands="single",
        resampling=planning.Resampling(
            planning.BatchSampling(size, seed=3),
            3,
            attempts,
            feedback=planning.FeedbackPolicy(failure_alpha=2),
        ),
    )


def decisions(path: Path) -> list:
    return list(da.inspect(path, view="batches", all=True).records())


def test_rejections_feed_back_but_duplicates_and_exhaustion_do_not(tmp_path: Path):
    rejected = recipe(size=2).with_changes(
        requirements=(planning.Avoid("all", ("AAA", "CCC")),)
    )
    path = tmp_path / "rejected"
    da.run(rejected, out=path)
    assert da.inspect(path, verify=True).counts["rejected"] == 3
    assert [dict(d.batch.feedback.failed) for d in decisions(path)] == [
        {},
        {"a": 1, "b": 1},
        {"a": 2, "b": 2},
    ]

    path = tmp_path / "exhausted"
    da.run(recipe(attempts=4, size=2), out=path)
    summary = da.inspect(path, verify=True)
    assert summary.accepted == 2
    assert summary.counts["duplicate"] > 0
    assert summary.counts["no_candidate"] == 3
    assert all(not d.batch.feedback.failed for d in decisions(path))
    assert dict(decisions(path)[1].batch.feedback.used) == {"a": 1, "b": 1}


def test_initial_infeasibility_marks_the_offered_parts(tmp_path: Path):
    request = recipe().with_changes(
        requirements=(planning.GroupCoverage("both", ("A", "B"), min=2),)
    )
    path = tmp_path / "run"
    da.run(request, out=path)
    assert da.inspect(path, verify=True).accepted == 0
    saved = decisions(path)
    assert not saved[0].batch.feedback.failed
    assert sum(saved[1].batch.feedback.failed.values()) == 1
    assert sum(saved[2].batch.feedback.failed.values()) == 2


def test_resume_inside_a_batch_never_samples_that_batch_again(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    request = recipe(attempts=4, size=2).with_changes(target=planning.Target(2))
    original = RunWriter.publish

    def interrupt(self: RunWriter, *args: object, **kwargs: object) -> object:
        result = original(self, *args, **kwargs)
        if self.state["counts"]["accepted"] == 1:
            raise KeyboardInterrupt
        return result

    path = tmp_path / "run"
    with monkeypatch.context() as patch:
        patch.setattr(RunWriter, "publish", interrupt)
        with pytest.raises(KeyboardInterrupt):
            da.run(request, out=path)
    before = decisions(path)

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("recovery sampled the saved batch again")

    monkeypatch.setattr(sampling, "sample_batch", forbidden)
    monkeypatch.setattr("dense_arrays.workflow.batch_cursor.sample_batch", forbidden)
    da.run(resume=path)
    assert da.inspect(path, verify=True).accepted == 2
    assert decisions(path) == before


def test_changed_feedback_cannot_be_hidden_by_rehashing(tmp_path: Path):
    path = tmp_path / "run"
    da.run(recipe(), out=path)
    with sqlite3.connect(path / "run.sqlite3") as connection:
        raw = connection.execute(
            "SELECT payload FROM batches WHERE ordinal=2"
        ).fetchone()[0]
        payload = json.loads(raw)
        prior_id = payload["batch"]["batch_id"]
        payload["batch"]["feedback"]["used"] = {"a": 99}
        content = dict(payload["batch"])
        content.pop("batch_id")
        payload["batch"]["batch_id"] = semantic_digest(content)
        connection.execute(
            "UPDATE batches SET payload=?,digest=? WHERE ordinal=2",
            (canonical_json(payload), semantic_digest(payload)),
        )
        for table in ("attempts", "designs"):
            for rowid, raw in connection.execute(
                f"SELECT rowid,payload FROM {table}"  # noqa: S608 - two literal fixture tables
            ).fetchall():
                data = json.loads(raw)
                evidence = data["evidence"] if table == "attempts" else data
                if evidence.get("batch_id") == prior_id:
                    evidence["batch_id"] = payload["batch"]["batch_id"]
                    connection.execute(
                        f"UPDATE {table} SET payload=?,digest=? WHERE rowid=?",  # noqa: S608
                        (canonical_json(data), semantic_digest(data), rowid),
                    )
    with pytest.raises(ArtifactIntegrityError, match="feedback"):
        da.inspect(path, verify=True)
