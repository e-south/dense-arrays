"""Model admission bounds offered work without invalidating saved evidence.

Author: Eric J. South.
"""

import json
from dataclasses import replace
from itertools import islice, product
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.cli import app
from dense_arrays.generation import packing
from dense_arrays.planning.serialization import request_from_dict, request_to_dict


def recipe(count: int = 4, **changes: object):
    return planning.DesignSpec(
        parts=tuple(
            parts.Part(f"p{i}", "".join(sequence))
            for i, sequence in enumerate(islice(product("ACGT", repeat=5), count))
        ),
        length=planning.Length(maximum=5),
        strands="single",
        **changes,
    )


def test_default_model_admission_fails_before_destination_or_optimizer(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("oversized work reached optimizer allocation")

    monkeypatch.setattr(packing, "Optimizer", forbidden)
    request = recipe(501)
    saved = da.plan(request)
    assert saved.preview["oriented_nodes"] == 501
    restored = planning.GenerationPlan.from_dict(saved.to_dict())
    assert restored.plan_id == saved.plan_id  # Reading does not allocate a model.
    with pytest.raises(ValueError, match=r"limits\.model_pairs"):
        da.run(restored, out=tmp_path / "run")
    assert not (tmp_path / "run").exists()


@pytest.mark.parametrize("search", ["exact", "greedy"])
def test_model_admission_counts_orientations_and_preserves_boundary(search: str):
    bounded = planning.Limits(model_pairs=16)
    request = recipe(limits=bounded, search=search)
    assert da.plan(request).preview["oriented_nodes"] == 4
    with pytest.raises(ValueError, match=r"64.*limits.model_pairs.*16"):
        da.plan(request.with_changes(strands="double")).admit_work()
    assert da.plan(
        request.with_changes(limits=replace(bounded, model_pairs=64), strands="double")
    )


def test_offered_batches_and_inactive_matrix_cells_control_admission():
    source = da.plan(recipe())
    small = planning.CandidateBatch(("p0", "p1"), source.collection_id)
    request = recipe(limits=planning.Limits(model_pairs=4))
    assert da.plan(request.with_changes(batch=small)).preview["oriented_nodes"] == 2
    matrix = planning.MatrixSpec(
        base=request,
        axes={
            "condition": {"active": planning.Variant(), "inactive": planning.Variant()}
        },
        allocation=planning.Allocation(
            counts={"condition=active": 1, "condition=inactive": 0}
        ),
        max_cells=2,
        batches={"condition=active": small},
    )
    assert len(da.plan(matrix).cells) == 2
    with pytest.raises(ValueError, match=r"limits\.model_pairs"):
        da.plan(
            matrix.with_changes(
                allocation=planning.Allocation(total=2, policy="balanced")
            )
        ).admit_work()


def test_declared_model_cap_roundtrips_with_matching_cli_errors(tmp_path: Path):
    request = recipe(limits=planning.Limits(model_pairs=15))
    document = request_to_dict(request)
    assert document["limits"]["model_pairs"] == 15
    assert request_from_dict(document) == request
    path = tmp_path / "request.json"
    path.write_text(json.dumps(document))
    output = tmp_path / "run"
    result = CliRunner().invoke(app, ["run", str(path), "--out", str(output), "--json"])
    assert result.exit_code == 2, result.output
    assert "limits.model_pairs" in json.loads(result.stdout)["message"]
    assert not output.exists()
    # New defaults preserve the existing canonical request encoding.
    assert request_to_dict(recipe())["limits"] == {
        "attempts": 1000,
        "active_seconds": 300,
        "solver_seconds": 30,
    }
    extension = planning.ExtensionSpec(
        planning.ParentRun("unused"), 1, request.limits, 0
    )
    assert (
        planning.ExtensionSpec.from_dict(extension.to_dict()).limits == request.limits
    )


@pytest.mark.parametrize("cap", [True, 0, -1, 1.5, None])
def test_model_cap_requires_a_positive_integer(cap: object):
    with pytest.raises((TypeError, ValueError), match="model_pairs"):
        planning.Limits(model_pairs=cap)


def test_model_builder_rechecks_admission_without_allocating(
    monkeypatch: pytest.MonkeyPatch,
):
    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("oversized work reached optimizer allocation")

    monkeypatch.setattr(packing, "Optimizer", forbidden)
    saved = planning.GenerationPlan(recipe(limits=planning.Limits(model_pairs=15)))
    with pytest.raises(ValueError, match=r"limits\.model_pairs"):
        packing.build_optimizer(saved, seconds=0.001)


def test_scheduled_and_resampled_work_use_largest_offered_batch():
    source = da.plan(recipe())
    small = planning.CandidateBatch(("p0", "p1"), source.collection_id)
    large = planning.CandidateBatch(("p0", "p1", "p2"), source.collection_id)
    request = recipe(limits=planning.Limits(model_pairs=4))
    schedule = planning.BatchSchedule((small,), attempts_per_batch=1)
    assert da.plan(request.with_changes(schedule=schedule))
    with pytest.raises(ValueError, match=r"9.*limits.model_pairs.*4"):
        da.plan(
            request.with_changes(schedule=replace(schedule, batches=(small, large)))
        ).admit_work()
    resampling = planning.Resampling(
        planning.BatchSampling(size=2, seed=7), max_batches=2, attempts_per_batch=1
    )
    assert da.plan(request.with_changes(resampling=resampling))
    with pytest.raises(ValueError, match=r"9.*limits.model_pairs.*4"):
        da.plan(
            request.with_changes(
                resampling=replace(
                    resampling, sampling=replace(resampling.sampling, size=3)
                )
            )
        ).admit_work()


def test_large_collection_can_be_planned_then_prepared_as_a_bounded_batch(
    tmp_path: Path,
):
    source = da.plan(recipe(501))
    selected = da.prepare(
        source,
        sampling=planning.BatchSampling(size=4, seed=7),
        out=tmp_path / "batch.json",
    )
    assert selected.preview["oriented_nodes"] == 4
    result = da.run(selected, out=tmp_path / "run")
    assert da.inspect(result, verify=True).accepted == 1
