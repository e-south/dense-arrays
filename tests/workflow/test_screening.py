"""Final sequence screening includes junctions, padding and interval exceptions.

Author: Eric J. South.
"""

from pathlib import Path

import pytest

import dense_arrays as da
from dense_arrays import parts, planning


def test_join_created_pattern_rejects_an_optimal_packing(tmp_path: Path):
    request = planning.DesignSpec(
        parts=[parts.Part("a", "AA"), parts.Part("b", "CC")],
        length=planning.Length(maximum=4),
        strands="single",
        limits=planning.Limits(attempts=1),
        requirements=[
            planning.Fixed("a", "a", "forward", planning.StartWindow(max=0)),
            planning.Fixed("b", "b", "forward"),
            planning.Avoid(
                "no-join",
                patterns=("AC",),
                strands="both",
                except_placements=("a", "b"),
            ),
        ],
    )
    result = da.run(request, out=tmp_path / "run")
    summary = da.inspect(result, verify=True)
    assert summary.accepted == 0
    assert summary.counts["rejected"] == 1
    with da.inspect(result, view="attempts").records() as records:
        attempt = next(records)
    assert attempt.evidence["solver_status"] == "optimal"
    assert attempt.evidence["code"] == "screening_rejection"
    rule = attempt.evidence["requirements"][-1]
    assert rule["observed"][0]["start"] == 1
    assert rule["observed"][0]["end"] == 3


@pytest.mark.parametrize("exempt", [False, True])
def test_reverse_pattern_exception_is_tied_to_a_fixed_interval(
    exempt: bool, tmp_path: Path
):
    request = planning.DesignSpec(
        parts=[parts.Part("a", "AAA")],
        length=planning.Length(maximum=3),
        strands="single",
        limits=planning.Limits(attempts=1),
        requirements=[
            planning.Fixed("anchor", "a", "forward"),
            planning.Avoid(
                "no-TT",
                patterns=("TT",),
                strands="both",
                except_placements=("a",) if exempt else (),
            ),
        ],
    )
    report = da.inspect(da.run(request, out=tmp_path / "run"), verify=True)
    assert report.accepted == int(exempt)


def test_padding_exhaustion_is_bounded_and_is_not_infeasibility(tmp_path: Path):
    request = planning.DesignSpec(
        parts=[parts.Part("a", "AAA")],
        length=planning.Length(exact=4),
        assembly=planning.Assembly(
            padding=planning.Padding(side="right", max_trials=3)
        ),
        strands="single",
        limits=planning.Limits(attempts=1),
        requirements=[
            planning.Fixed("anchor", "a", "forward"),
            planning.Avoid(
                "no-extra",
                patterns=("A", "C", "G", "T"),
                strands="forward",
                except_placements=("a",),
            ),
        ],
    )
    result = da.run(request, out=tmp_path / "run")
    summary = da.inspect(result, verify=True)
    assert summary.counts["started"] == summary.counts["rejected"] == 1
    assert summary.termination_reason == "attempt_limit"
    with da.inspect(result, view="attempts").records() as records:
        attempt = next(records)
    assert attempt.evidence["assembly_trials"] == 3
    assert attempt.evidence["code"] == "padding_trials_exhausted"
    assert attempt.evidence["requirements"][-1]["observed"][0]["start"] == 3


def test_zero_padding_gc_is_not_applicable(tmp_path: Path):
    request = planning.DesignSpec(
        parts=[parts.Part("a", "ACG")],
        length=planning.Length(exact=3),
        assembly=planning.Assembly(
            padding=planning.Padding(side="right", max_trials=3)
        ),
        strands="single",
        requirements=[planning.GC("padding-gc", scope="padding", min=1, max=1)],
    )
    result = da.run(request, out=tmp_path / "run")
    assert da.inspect(result, verify=True).accepted == 1
    with da.inspect(result, view="designs").records() as records:
        design = next(records)
    assert design.requirements[0]["observed"] is None
    assert design.requirements[0]["status"] == "not_applicable"


def test_gc_uses_discrete_counts_without_rounding_into_acceptance(tmp_path: Path):
    request = planning.DesignSpec(
        parts=[parts.Part("a", "ACG")],
        length=planning.Length(maximum=3),
        strands="single",
        limits=planning.Limits(attempts=1),
        requirements=[planning.GC("gc", scope="sequence", min=0.66, max=0.66)],
    )
    summary = da.inspect(da.run(request, out=tmp_path / "run"), verify=True)
    assert summary.accepted == 0
    assert summary.counts["rejected"] == 1
