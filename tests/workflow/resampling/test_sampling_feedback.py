"""Feedback weights preserve declared observations and version random selection.

Author: Eric J. South.
"""

import json
import math
from collections import Counter
from pathlib import Path

import pytest

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.artifacts.reading import ReadLimitError
from dense_arrays.artifacts.run_plans import decode_plan
from dense_arrays.generation.batches.sampling import sample_batch
from dense_arrays.planning.batches.bindings import membership_size


def source() -> planning.GenerationPlan:
    return da.plan(
        planning.DesignSpec(
            parts=(
                parts.Part("a", "AAA", group="A"),
                parts.Part("alias", "AAA", group="A"),
                parts.Part("b", "CCC", group="B"),
            ),
            length=planning.Length(maximum=3),
            strands="single",
        )
    )


def test_feedback_weights_share_group_sequence_observations():
    feedback = planning.FeedbackSnapshot(
        planning.FeedbackPolicy(coverage_alpha=3, failure_alpha=2),
        used={"a": 3},
        failed={"b": 2},
    )
    assert feedback.weights(source().request.parts) == pytest.approx(
        {"a": 1.75, "alias": 1.75, "b": 0.8}
    )
    assert planning.FeedbackSnapshot.from_dict(feedback.to_dict()) == feedback


def test_weighted_batch_retains_feedback_and_favors_unused_parts():
    feedback = planning.FeedbackSnapshot(
        planning.FeedbackPolicy(coverage_alpha=1e9), used={"a": 1000000}
    )
    batch = sample_batch(
        source(),
        planning.BatchSampling(1, seed=7),
        stream="default/batch/2",
        feedback=feedback,
    )
    assert batch.part_ids == ("b",)
    assert batch.feedback == feedback
    assert planning.CandidateBatch.from_dict(batch.to_dict()) == batch
    assert da.plan(source().request.with_changes(batch=batch)).request.batch == batch


def test_failure_penalty_is_applied_without_removing_candidates():
    feedback = planning.FeedbackSnapshot(
        planning.FeedbackPolicy(coverage_alpha=0, failure_alpha=1e9),
        failed={"a": 1000000},
    )
    assert all(value > 0 for value in feedback.weights(source().request.parts).values())
    batch = sample_batch(
        source(),
        planning.BatchSampling(1, seed=7),
        stream="default/batch/2",
        feedback=feedback,
    )
    assert batch.part_ids == ("b",)


def test_empty_feedback_keeps_part_priority_order():
    plan = source()
    policy = planning.BatchSampling(2, seed=41, unique_sequences=True)
    ordinary = sample_batch(plan, policy, stream="same")
    weighted = sample_batch(
        plan,
        policy,
        stream="same",
        feedback=planning.FeedbackSnapshot(
            planning.FeedbackPolicy(coverage_alpha=0, failure_alpha=0)
        ),
    )
    assert weighted.part_ids == ordinary.part_ids
    assert weighted.batch_id != ordinary.batch_id


def test_unknown_feedback_parts_fail_before_selection():
    feedback = planning.FeedbackSnapshot(planning.FeedbackPolicy(), used={"typo": 1})
    with pytest.raises(ValueError, match=r"unknown.*feedback"):
        sample_batch(source(), planning.BatchSampling(1), stream="x", feedback=feedback)


@pytest.mark.parametrize("value", [-1, True, float("inf"), float("nan")])
def test_feedback_policy_rejects_invalid_numeric_values(value: object):
    with pytest.raises((ValueError, TypeError)):
        planning.FeedbackPolicy(coverage_alpha=value)


def test_feedback_snapshot_detaches_caller_counts_and_rejects_non_counts():
    values = {"a": 2}
    feedback = planning.FeedbackSnapshot(planning.FeedbackPolicy(), used=values)
    values["a"] = 99
    assert feedback.used["a"] == 2
    for invalid in (-1, 0.5, True):
        with pytest.raises((ValueError, TypeError)):
            planning.FeedbackSnapshot(planning.FeedbackPolicy(), failed={"a": invalid})


def test_numeric_underflow_is_explicit():
    feedback = planning.FeedbackSnapshot(
        planning.FeedbackPolicy(failure_alpha=1e200, failure_power=2),
        failed={"a": 10000},
    )
    with pytest.raises(ValueError, match="numeric range"):
        feedback.weights(source().request.parts)
    assert math.isfinite(planning.FeedbackPolicy().weight(1, 1))


def test_disabled_weight_terms_do_not_overflow():
    policy = planning.FeedbackPolicy(
        coverage_alpha=0, coverage_power=1e308, failure_alpha=0, failure_power=1e308
    )
    assert policy.weight(10, 10) == 1


def test_prepared_feedback_is_charged_to_plan_read_limits():

    plan = source()
    batch = sample_batch(
        plan,
        planning.BatchSampling(1),
        stream="one",
        feedback=planning.FeedbackSnapshot(
            planning.FeedbackPolicy(), used={"a": 1, "b": 2}, failed={"alias": 1}
        ),
    )
    bound = da.plan(plan.request.with_changes(batch=batch))
    assert membership_size(bound.request) == 4
    with pytest.raises(ReadLimitError):
        decode_plan(bound.to_dict(), 6)


def test_weight_formula_matches_pinned_source_fixture():
    fixture = json.loads(
        (
            Path(__file__).parents[2] / "fixtures/workflow/densegen-feedback-v1.json"
        ).read_text()
    )
    rows = fixture["parts"]
    items = tuple(
        parts.Part(r["part_id"], r["sequence"], group=r["group"]) for r in rows
    )
    snapshot = planning.FeedbackSnapshot(
        planning.FeedbackPolicy(**fixture["parameters"]),
        used={r["part_id"]: r["used"] for r in rows},
        failed={r["part_id"]: r["failed"] for r in rows},
    )
    weights = snapshot.weights(items)
    grouped = Counter()
    for part in items:
        grouped[part.group] += weights[part.part_id]
    assert dict(grouped) == pytest.approx(fixture["group_weights"])
