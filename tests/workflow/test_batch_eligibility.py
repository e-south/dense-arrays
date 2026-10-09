"""Batch identity policies compose without sacrificing feasible selections.

Author: Eric J. South.
"""

import json
from collections import Counter
from itertools import combinations
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.cli import app
from dense_arrays.planning.batches.validation import validate_selected
from dense_arrays.workflow.inputs import read_source


def source() -> planning.GenerationPlan:
    return da.plan(
        planning.DesignSpec(
            parts=(
                parts.Part(
                    "a",
                    "AAA",
                    group="G",
                    core_start=0,
                    core_end=1,
                    core_orientation="forward",
                ),
                parts.Part(
                    "b",
                    "AAA",
                    group="G",
                    core_start=0,
                    core_end=2,
                    core_orientation="forward",
                ),
                parts.Part(
                    "c",
                    "AAC",
                    group="G",
                    core_start=0,
                    core_end=1,
                    core_orientation="forward",
                ),
            ),
            length=planning.Length(maximum=3),
            strands="single",
            target=planning.Target(2),
        )
    )


@pytest.mark.parametrize("strategy", ["uniform", "group_balanced"])
def test_joint_uniqueness_finds_full_selection_instead_of_greedy_dead_end(
    tmp_path: Path, strategy: str
):
    plan = source()
    for seed in range(6):
        policy = planning.BatchSampling(
            2, strategy=strategy, seed=seed, unique_sequences=True, unique_cores=True
        )
        prepared = da.prepare(plan, sampling=policy, out=tmp_path / f"{seed}.json")
        assert set(prepared.request.batch.part_ids) == {"b", "c"}
        replay = read_source(tmp_path / f"{seed}.json")
        assert replay.plan_id == prepared.plan_id
    result = da.run(replay, out=tmp_path / "run")
    assert da.inspect(result, verify=True).accepted == 2


def test_group_cap_does_not_drop_sequence_alternatives_prematurely(tmp_path: Path):
    plan = source()
    candidates = tuple(
        parts.Part(p.part_id, p.sequence, group="B" if p.part_id == "b" else "A")
        for p in plan.request.parts
    )
    plan = da.plan(plan.request.with_changes(parts=candidates))
    for seed in range(6):
        prepared = da.prepare(
            plan,
            sampling=planning.BatchSampling(
                2, seed=seed, unique_sequences=True, max_per_group=1
            ),
            out=tmp_path / f"{seed}.json",
        )
        assert set(prepared.request.batch.part_ids) == {"b", "c"}


def test_core_identity_uses_declared_orientation_and_group(tmp_path: Path):
    candidates = (
        parts.Part(
            "forward",
            "CAA",
            group="A",
            core_start=1,
            core_end=3,
            core_orientation="forward",
        ),
        parts.Part(
            "reverse",
            "ATT",
            group="A",
            core_start=1,
            core_end=3,
            core_orientation="reverse",
        ),
        parts.Part(
            "other",
            "GAA",
            group="B",
            core_start=1,
            core_end=3,
            core_orientation="forward",
        ),
    )
    assert candidates[0].core_sequence == candidates[1].core_sequence == "AA"
    plan = da.plan(source().request.with_changes(parts=candidates))
    prepared = da.prepare(
        plan,
        sampling=planning.BatchSampling(2, unique_cores=True),
        out=tmp_path / "plan.json",
    )
    selected = set(prepared.request.batch.part_ids)
    assert "other" in selected
    assert len(selected & {"forward", "reverse"}) == 1


def test_fixed_parts_count_toward_uniqueness_and_caps(tmp_path: Path):
    plan = da.plan(
        source().request.with_changes(
            requirements=(planning.Fixed("anchor", "a", "forward"),)
        )
    )
    with pytest.raises(ValueError, match=r"batch size.*uniqueness.*cap"):
        da.prepare(
            plan,
            sampling=planning.BatchSampling(
                2, unique_sequences=True, unique_cores=True
            ),
            out=tmp_path / "impossible.json",
        )
    assert not (tmp_path / "impossible.json").exists()
    with pytest.raises(ValueError, match=r"batch size.*cap"):
        da.prepare(
            plan,
            sampling=planning.BatchSampling(2, max_per_group=1),
            out=tmp_path / "cap.json",
        )
    assert not (tmp_path / "cap.json").exists()


def test_missing_core_annotations_fail_before_sampling(tmp_path: Path):
    plan = da.plan(
        source().request.with_changes(parts=(parts.Part("a", "AAA", group="A"),))
    )
    with pytest.raises(ValueError, match="core annotations"):
        da.prepare(
            plan,
            sampling=planning.BatchSampling(1, unique_cores=True),
            out=tmp_path / "plan.json",
        )
    assert not (tmp_path / "plan.json").exists()


def test_prepared_membership_must_obey_declared_sampling_policy():
    plan = source()
    forged = planning.CandidateBatch(
        ("a", "b"),
        plan.collection_id,
        sampling=planning.BatchSampling(2, unique_sequences=True),
    )
    with pytest.raises(ValueError, match="sequence uniqueness"):
        da.plan(plan.request.with_changes(batch=forged))


def test_cli_policy_equals_python_and_preserves_unconstrained_wire(tmp_path: Path):
    plan = source()
    plan.write(tmp_path / "source.json")
    policy = planning.BatchSampling(
        2, unique_sequences=True, unique_cores=True, max_per_group=2
    )
    prepared = da.prepare(plan, sampling=policy, out=tmp_path / "api.json")
    response = CliRunner().invoke(
        app,
        [
            "prepare",
            str(tmp_path / "source.json"),
            "--batch-size",
            "2",
            "--unique-sequences",
            "--unique-cores",
            "--max-per-group",
            "2",
            "--out",
            str(tmp_path / "cli.json"),
            "--json",
        ],
    )
    assert response.exit_code == 0, response.output
    assert read_source(tmp_path / "cli.json").plan_id == prepared.plan_id
    assert planning.BatchSampling(2).to_dict() == {
        "policy": "part_priority_sha256.v1",
        "size": 2,
        "strategy": "uniform",
        "seed": 0,
    }
    assert policy.to_dict()["policy"] != planning.BatchSampling(2).to_dict()["policy"]


@pytest.mark.parametrize(
    "kwargs",
    [
        {"unique_sequences": 1},
        {"unique_cores": "true"},
        {"max_per_group": 0},
        {"max_per_group": True},
    ],
)
def test_invalid_eligibility_policies_fail(kwargs: dict[str, object]):
    with pytest.raises((ValueError, TypeError)):
        planning.BatchSampling(2, **kwargs)


def test_group_balance_matches_enumerated_feasible_subsets_and_counts_fixed(
    tmp_path: Path,
):
    candidates = tuple(
        parts.Part(f"p{i}", seq, group=group)
        for i, (seq, group) in enumerate(
            zip(
                ("AAA", "CCC", "AAA", "GGG", "TTT", "GGG", "AAC"),
                ("A", "A", "B", "B", "B", "C", "C"),
                strict=True,
            )
        )
    )
    plan = da.plan(
        source().request.with_changes(
            parts=candidates,
            requirements=(planning.Fixed("anchor", "p0", "forward"),),
        )
    )
    feasible = [
        combo
        for combo in combinations(candidates, 4)
        if "p0" in {p.part_id for p in combo}
        and len({p.sequence for p in combo}) == 4
        and max(Counter(p.group for p in combo).values()) <= 2
    ]
    best_balance = min(
        sum(n * n for n in Counter(p.group for p in c).values()) for c in feasible
    )
    for seed in range(6):
        prepared = da.prepare(
            plan,
            sampling=planning.BatchSampling(
                4, "group_balanced", seed, unique_sequences=True, max_per_group=2
            ),
            out=tmp_path / f"{seed}.json",
        )
        selected = [
            p for p in candidates if p.part_id in prepared.request.batch.part_ids
        ]
        assert "p0" in prepared.request.batch.part_ids
        assert len(selected) == len({p.sequence for p in selected}) == 4
        counts = Counter(p.group for p in selected)
        assert sum(n * n for n in counts.values()) == best_balance
        assert max(counts.values()) <= 2


def test_constrained_selection_is_invariant_to_input_order(tmp_path: Path):
    plan = source()
    reversed_plan = da.plan(
        plan.request.with_changes(parts=tuple(reversed(plan.request.parts)))
    )
    policy = planning.BatchSampling(2, unique_cores=True, max_per_group=2)
    left = da.prepare(plan, sampling=policy, out=tmp_path / "left.json")
    right = da.prepare(reversed_plan, sampling=policy, out=tmp_path / "right.json")
    assert left.request.batch.part_ids == right.request.batch.part_ids


def test_fixed_membership_cannot_violate_policy(tmp_path: Path):
    request = source().request.with_changes(
        requirements=(
            planning.Fixed("first", "a", "forward"),
            planning.Fixed("second", "b", "forward"),
        ),
        length=planning.Length(maximum=6),
    )
    with pytest.raises(ValueError, match="sequence uniqueness"):
        da.prepare(
            da.plan(request),
            sampling=planning.BatchSampling(2, unique_sequences=True),
            out=tmp_path / "duplicate.json",
        )
    with pytest.raises(ValueError, match="max_per_group"):
        da.prepare(
            da.plan(request),
            sampling=planning.BatchSampling(2, max_per_group=1),
            out=tmp_path / "cap.json",
        )


def test_backend_failure_is_not_reported_as_an_impossible_batch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    from dense_arrays.generation.batches import constrained  # noqa: PLC0415

    monkeypatch.setattr(
        constrained.SimpleMinCostFlow, "solve", lambda self: self.BAD_COST_RANGE
    )
    with pytest.raises(RuntimeError, match="backend failed"):
        da.prepare(
            source(),
            sampling=planning.BatchSampling(2, unique_sequences=True),
            out=tmp_path / "failed.json",
        )
    assert not (tmp_path / "failed.json").exists()


def test_uniqueness_membership_matches_pinned_append_helper():
    fixture = json.loads(
        (
            Path(__file__).parents[1]
            / "fixtures/workflow/densegen-batch-uniqueness-v1.json"
        ).read_text()
    )
    candidates = [
        parts.Part(
            str(i),
            r["tfbs"],
            group=r["tf"],
            core_start=r["tfbs"].index(r["tfbs_core"]),
            core_end=r["tfbs"].index(r["tfbs_core"]) + len(r["tfbs_core"]),
            core_orientation="forward",
        )
        for i, r in enumerate(fixture["rows"])
    ]
    for case in fixture["cases"]:
        policy = planning.BatchSampling(
            5,
            unique_sequences=case["unique_sequences"],
            unique_cores=case["unique_cores"],
        )
        selected, accepted = [], []
        for part in candidates:
            try:
                validate_selected([*selected, part], policy)
            except ValueError:
                accepted.append(False)
            else:
                selected.append(part)
                accepted.append(True)
        assert accepted == case["accepted"]


def test_saved_constrained_batch_and_bundle_do_not_sample_again(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    from dense_arrays.generation.batches import constrained  # noqa: PLC0415

    plan = da.prepare(
        source(),
        sampling=planning.BatchSampling(2, unique_sequences=True, unique_cores=True),
        out=tmp_path / "plan.json",
    )
    monkeypatch.setattr(
        constrained.SimpleMinCostFlow,
        "solve",
        lambda _self: pytest.fail("replay must not select again"),
    )
    result = da.run(read_source(tmp_path / "plan.json"), out=tmp_path / "run")
    assert da.inspect(result, verify=True).accepted == 2
    da.export(result, all=True, format="bundle", out=tmp_path / "bundle")
    assert da.inspect(tmp_path / "bundle", verify=True).designs == 2
    assert da.inspect(result, view="plan").plan_id == plan.plan_id
