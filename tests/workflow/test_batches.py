"""Prepared candidate batches keep selection, packing proof and replay distinct.

Author: Eric J. South.
"""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.cli import app
from dense_arrays.generation.batches import sampling
from dense_arrays.workflow.inputs import read_source


def recipe():
    return planning.DesignSpec(
        parts=[
            parts.Part("a", "AAA", group="x"),
            parts.Part("b", "CCC", group="x"),
            parts.Part("c", "GGG", group="y"),
            parts.Part("d", "TTT", group="y"),
        ],
        length=planning.Length(maximum=3),
        strands="single",
        target=planning.Target(3),
    )


def test_prepared_batch_is_portable_executable_and_bounded(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    source = da.plan(recipe())
    sampled = da.prepare(
        source,
        sampling=planning.BatchSampling(size=2, seed=7),
        out=tmp_path / "batch.json",
    )
    batch = sampled.request.batch
    assert len(batch.part_ids) == 2
    assert batch.collection_id == source.collection_id
    assert sampled.preview["parts"] == 4
    assert sampled.preview["offered_parts"] == 2
    assert sampled.preview["oriented_nodes"] == 2
    loaded = read_source(tmp_path / "batch.json")
    assert loaded.plan_id == sampled.plan_id

    monkeypatch.setattr(
        sampling, "sample_batch", lambda *_a, **_kw: pytest.fail("replay sampled again")
    )
    result = da.run(loaded, out=tmp_path / "run")
    summary = da.inspect(result, verify=True)
    assert summary.accepted == 2
    assert summary.termination_reason == "batch_exhausted"
    assert summary.state == "stopped"
    attempts = list(da.inspect(result, view="attempts", all=True).records())
    assert {a.evidence["batch_id"] for a in attempts} == {batch.batch_id}
    designs = list(da.inspect(result, view="designs", all=True).records())
    assert {p.feature_id for d in designs for p in d.realized.placements} == set(
        batch.part_ids
    )
    replay = da.plan(da.inspect(result, view="request").request)
    assert replay.request.batch == sampled.request.batch
    assert replay.request.lineage.parent.plan_id == sampled.plan_id


def test_explicit_batch_preserves_identity_and_local_infeasibility(tmp_path: Path):
    base = da.plan(recipe())
    batch = planning.CandidateBatch(
        part_ids=("b", "a"), collection_id=base.collection_id
    )
    request = recipe().with_changes(
        batch=batch,
        requirements=(
            planning.Occurrences("need_y", parts.PartSelector(groups=("y",)), min=1),
        ),
    )
    result = da.run(request, out=tmp_path / "impossible")
    summary = da.inspect(result, verify=True)
    assert summary.accepted == 0
    assert summary.termination_reason == "batch_infeasible"
    assert summary.counts["no_candidate"] == 1
    assert (
        da.inspect(result, view="attempts", all=True)
        .records()
        .__next__()
        .evidence["solver_status"]
        == "infeasible"
    )
    # A different offered batch can satisfy the unchanged packing requirements.
    other = planning.CandidateBatch(part_ids=("c",), collection_id=base.collection_id)
    good = da.run(
        request.with_changes(batch=other, target=planning.Target(1)),
        out=tmp_path / "good",
    )
    assert da.inspect(good, verify=True).accepted == 1


def test_group_balancing_fixed_parts_and_matrix_streams(tmp_path: Path):
    spec = recipe().with_changes(
        requirements=(planning.Fixed("anchor", "a", "forward"),)
    )
    matrix = da.plan(
        planning.MatrixSpec(
            base=spec.with_changes(target=planning.Target()),
            axes={"cell": {"one": planning.Variant(), "two": planning.Variant()}},
            allocation=planning.Allocation(per_cell=1),
            max_cells=2,
        )
    )
    policy = planning.BatchSampling(size=2, strategy="group_balanced", seed=10)
    prepared = da.prepare(matrix, sampling=policy, out=tmp_path / "matrix.json")
    for cell in prepared.cells:
        batch = cell.plan.request.batch
        assert "a" in batch.part_ids
        offered = [p for p in cell.plan.request.parts if p.part_id in batch.part_ids]
        assert {p.group for p in offered} == {"x", "y"}
    assert (
        prepared.cells[0].plan.request.batch.batch_id
        != prepared.cells[1].plan.request.batch.batch_id
    )
    result = da.run(read_source(tmp_path / "matrix.json"), out=tmp_path / "run")
    assert da.inspect(result, verify=True).accepted == 2
    assert (
        da.plan(da.inspect(result, view="request").request).plan_id == prepared.plan_id
    )


def test_prepare_cli_matches_python_and_never_replaces_output(tmp_path: Path):
    source = da.plan(recipe())
    source.write(tmp_path / "plan.json")
    expected = da.prepare(
        source,
        sampling=planning.BatchSampling(size=2, seed=7),
        out=tmp_path / "api.json",
    )
    runner = CliRunner()
    args = [
        "prepare",
        str(tmp_path / "plan.json"),
        "--out",
        str(tmp_path / "cli.json"),
        "--batch-size",
        "2",
        "--batch-seed",
        "7",
        "--json",
    ]
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout)["plan_id"] == expected.plan_id
    before = (tmp_path / "cli.json").read_bytes()
    assert runner.invoke(app, args).exit_code != 0
    assert (tmp_path / "cli.json").read_bytes() == before


def test_batch_membership_binding_and_policy_errors_fail_before_output(tmp_path: Path):
    base = da.plan(recipe())
    for ids, collection in [(("missing",), base.collection_id), (("a",), "0" * 64)]:
        batch = planning.CandidateBatch(part_ids=ids, collection_id=collection)
        with pytest.raises(ValueError, match="batch"):
            da.run(recipe().with_changes(batch=batch), out=tmp_path / "bad")
        assert not (tmp_path / "bad").exists()
    with pytest.raises(ValueError, match="size"):
        da.prepare(
            base, sampling=planning.BatchSampling(size=5), out=tmp_path / "too_big"
        )
    assert not (tmp_path / "too_big").exists()
    with pytest.raises(ValueError, match="unique"):
        planning.CandidateBatch(part_ids=("a", "a"), collection_id=base.collection_id)


@pytest.mark.parametrize("matrix", [False, True])
@pytest.mark.parametrize("boundary", ["reserved", "accepted"])
def test_batch_recovery_retains_committed_selection_and_prefix(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, matrix: bool, boundary: str
):
    from dense_arrays.artifacts.store import RunWriter  # noqa: PLC0415

    spec = recipe().with_changes(target=planning.Target(2))
    source = (
        da.plan(
            planning.MatrixSpec(
                base=spec.with_changes(target=planning.Target()),
                axes={"x": {"one": planning.Variant(), "two": planning.Variant()}},
                allocation=planning.Allocation(per_cell=2),
                max_cells=2,
            )
        )
        if matrix
        else da.plan(spec)
    )
    prepared = da.prepare(
        source, sampling=planning.BatchSampling(3, seed=5), out=tmp_path / "plan.json"
    )
    method = "reserve" if boundary == "reserved" else "publish"
    original = getattr(RunWriter, method)

    def interrupt(self: RunWriter, *args: object, **kwargs: object) -> object:
        result = original(self, *args, **kwargs)
        if boundary == "reserved" or result == "accepted":
            raise KeyboardInterrupt
        return result

    with monkeypatch.context() as patch:
        patch.setattr(RunWriter, method, interrupt)
        with pytest.raises(KeyboardInterrupt):
            da.run(prepared, out=tmp_path / "run")
    before = da.inspect(tmp_path / "run", verify=True)
    assert before.resumable
    prefix = [
        r.to_dict()
        for r in da.inspect(tmp_path / "run", view="designs", all=True).records()
    ]
    monkeypatch.setattr(
        sampling, "sample_batch", lambda *_a, **_kw: pytest.fail("resampled")
    )
    result = da.run(resume=tmp_path / "run")
    summary = da.inspect(result, verify=True)
    assert summary.accepted == (4 if matrix else 2)
    assert summary.state == "completed"
    assert [
        r.to_dict() for r in da.inspect(result, view="designs", all=True).records()
    ][: len(prefix)] == prefix


def test_compact_model_keeps_identical_sequences_distinct(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    from dense_arrays.optimizer import Optimizer  # noqa: PLC0415

    spec = recipe().with_changes(
        parts=[parts.Part("first", "AAA"), parts.Part("second", "AAA")],
        target=planning.Target(1),
    )
    source = da.plan(spec)
    bound = da.plan(
        spec.with_changes(
            batch=planning.CandidateBatch(("second",), source.collection_id)
        )
    )
    counts = []
    original = Optimizer.__init__

    def capture(
        self: Optimizer, library: list[str], *args: object, **kwargs: object
    ) -> None:
        counts.append(len(library))
        original(self, library, *args, **kwargs)

    monkeypatch.setattr(Optimizer, "__init__", capture)
    result = da.run(bound, out=tmp_path / "run")
    design = next(da.inspect(result, view="designs", all=True).records())
    assert design.realized.placements[0].feature_id == "second"
    assert design.realized.placements[0].placement_id == "p2"
    assert counts == [1]
    assert da.inspect(result, verify=True).accepted == 1


def test_sampling_order_is_stable_and_fixed_occurrences_count_toward_balance(
    tmp_path: Path,
):
    request = recipe().with_changes(
        requirements=(planning.Fixed("anchor", "a", "forward"),)
    )
    policy = planning.BatchSampling(2, "group_balanced", 2)
    a = da.prepare(da.plan(request), sampling=policy, out=tmp_path / "one.json")
    b = da.prepare(
        da.plan(request.with_changes(parts=tuple(reversed(request.parts)))),
        sampling=policy,
        out=tmp_path / "two.json",
    )
    assert a.request.batch.part_ids == b.request.batch.part_ids
    assert a.request.batch.part_ids[0] == "a"
    assert set(a.request.batch.part_ids) & {"c", "d"}
    with pytest.raises(ValueError, match="already"):
        da.prepare(a, sampling=policy, out=tmp_path / "again.json")
    assert not (tmp_path / "again.json").exists()


def test_prepared_file_survives_source_removal_and_checks_batch_digest(tmp_path: Path):
    source = tmp_path / "parts.csv"
    source.write_text("part_id,sequence\na,AAA\nb,CCC\nc,GGG\n")
    plan = da.plan(recipe().with_changes(parts=parts.PartTable(source, "csv")))
    prepared = da.prepare(
        plan, sampling=planning.BatchSampling(2), out=tmp_path / "prepared.json"
    )
    source.unlink()
    result = da.run(read_source(tmp_path / "prepared.json"), out=tmp_path / "run")
    assert da.inspect(result, verify=True).accepted == 2
    wire = prepared.to_dict()
    wire["request"]["batch"]["part_ids"].reverse()
    with pytest.raises(ValueError, match="batch digest"):
        planning.GenerationPlan.from_dict(wire)


def test_plan_read_limits_count_prepared_membership(tmp_path: Path):
    prepared = da.prepare(
        da.plan(recipe()),
        sampling=planning.BatchSampling(2),
        out=tmp_path / "plan.json",
    )
    with pytest.raises(da.reporting.ReadLimitError, match="identities"):
        da.inspect(
            prepared, view="plan", read_limits=da.reporting.ReadLimits(identities=4)
        )


def test_group_fill_fixture_records_behavior_preserved_and_corrected(tmp_path: Path):
    from collections import Counter  # noqa: PLC0415

    fixture = json.loads(
        (
            Path(__file__).parents[1] / "fixtures/workflow/densegen-batch-fill-v1.json"
        ).read_text()
    )
    for i, case in enumerate(fixture["cases"]):
        fixed = tuple(
            planning.Fixed(f"fixed_{key}", key, "forward") for key in case["fixed"]
        )
        source = da.plan(recipe().with_changes(requirements=fixed))
        prepared = da.prepare(
            source,
            sampling=planning.BatchSampling(
                case["size"], "group_balanced", case["seed"]
            ),
            out=tmp_path / f"{i}.json",
        )
        offered = set(prepared.request.batch.part_ids)
        counts = dict(
            Counter(p.group for p in prepared.request.parts if p.part_id in offered)
        )
        assert counts == {"x": 1, "y": 1}
        if not case["fixed"]:
            assert counts == case["legacy_group_counts"]
        else:
            assert case["legacy_group_counts"] == {"x": 2}
        assert (
            prepared.request.batch.sampling.to_dict()["policy"]
            == fixture["successor_policy"]
        )


def test_interrupted_reservation_does_not_imply_batch_enumeration(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    from dense_arrays.artifacts.store import RunWriter  # noqa: PLC0415

    request = recipe().with_changes(
        requirements=(planning.GroupCoverage("both", ("x", "y"), min=2),)
    )
    source = da.plan(request)
    request = request.with_changes(
        batch=planning.CandidateBatch(("a",), source.collection_id)
    )
    original = RunWriter.reserve

    def interrupt(self: RunWriter, *args: object, **kwargs: object) -> int:
        original(self, *args, **kwargs)
        raise KeyboardInterrupt

    with monkeypatch.context() as patch:
        patch.setattr(RunWriter, "reserve", interrupt)
        with pytest.raises(KeyboardInterrupt):
            da.run(request, out=tmp_path / "run")
    resumed = da.run(resume=tmp_path / "run")
    summary = da.inspect(resumed, verify=True)
    assert summary.termination_reason == "batch_infeasible"
    assert summary.counts["interrupted_unresolved"] == 1
    assert summary.counts["no_candidate"] == 1
