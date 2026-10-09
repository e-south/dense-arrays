"""Exact assembly translates geometry and records finite deterministic effort.

Author: Eric J. South.
"""

from dataclasses import replace
from pathlib import Path

import pytest

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.generation.acceptance import evaluate


@pytest.mark.parametrize("side,start", [("left", 4), ("right", 0)])
def test_padding_preserves_final_windows_and_replays(
    side: str, start: int, tmp_path: Path
):
    request = planning.DesignSpec(
        parts=[parts.Part("a", "AACC"), parts.Part("b", "CCGT")],
        length=planning.Length(exact=10),
        assembly=planning.Assembly(padding=planning.Padding(side=side, max_trials=5)),
        strands="single",
        seed=7,
        requirements=[
            planning.Fixed(
                "anchor", "a", "forward", planning.StartWindow(min=start, max=start)
            ),
            planning.Fixed("second", "b", "forward"),
            planning.Spacing("overlap", "a", "b", min=-2, max=-2),
        ],
    )
    plan = da.plan(request)
    assert planning.GenerationPlan.from_dict(plan.to_dict()) == plan
    sequences = []
    for name in ("first", "second"):
        run = da.run(plan, out=tmp_path / name)
        assert da.inspect(run, verify=True).accepted == 1
        with da.inspect(run, view="designs").records() as records:
            design = next(records)
        assert len(design.realized.sequence) == 10
        assert design.realized.placements[0].start == start
        assembly = design.realized.provenance["assembly"]
        assert assembly["padding_length"] == 4
        assert assembly["packed_start"] == start
        assert assembly["trial"] == 1
        sequences.append(design.realized.sequence)
    assert sequences[0] == sequences[1]


def test_exact_without_padding_never_fills_a_short_packing(tmp_path: Path):
    request = planning.DesignSpec(
        parts=[parts.Part("a", "AAA")],
        length=planning.Length(exact=4),
        assembly=planning.Assembly(),
        strands="single",
    )
    report = da.inspect(da.run(request, out=tmp_path / "run"), verify=True)
    assert report.accepted == 0
    assert report.termination_reason == "batch_infeasible"


@pytest.mark.parametrize("trials", [True, 0, -1, 1.5])
def test_padding_requires_a_finite_integral_trial_bound(trials: object):
    with pytest.raises((TypeError, ValueError), match="integer"):
        planning.Padding(side="right", max_trials=trials)


def test_changed_assembly_geometry_fails_independent_validation(tmp_path: Path):
    plan = da.plan(
        planning.DesignSpec(
            parts=[parts.Part("a", "AAA")],
            length=planning.Length(exact=5),
            assembly=planning.Assembly(
                padding=planning.Padding(side="left", max_trials=3)
            ),
            strands="single",
        )
    )
    run = da.run(plan, out=tmp_path / "run")
    with da.inspect(run, view="designs").records() as records:
        realized = next(records).realized
    altered = dict(realized.provenance["assembly"], packed_start=0)
    corrupted = replace(realized, provenance={"assembly": altered})
    with pytest.raises(ValueError, match="assembly"):
        evaluate(corrupted, plan)
