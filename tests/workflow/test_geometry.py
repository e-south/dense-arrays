"""Final-coordinate requirements preserve occurrence identity through the workflow.

Author: Eric J. South.
"""

from pathlib import Path

import pytest

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.generation.acceptance import evaluate
from dense_arrays.generation.randomness import PADDING_POLICY
from dense_arrays.realized import Orientation, Placement, PlacementKind, RealizedArray


def test_native_fixed_pair_roundtrips_and_recounts_signed_spacing(tmp_path: Path):
    request = planning.DesignSpec(
        parts=[parts.Part("a", "AACC"), parts.Part("b", "CCGT")],
        length=planning.Length(maximum=6),
        strands="single",
        requirements=[
            planning.Fixed(
                "a-fixed", "a", "forward", planning.StartWindow(min=0, max=0)
            ),
            planning.Fixed("b-fixed", "b", "forward"),
            planning.Spacing("overlap", "a", "b", min=-2, max=-2),
        ],
    )
    plan = da.plan(request)
    assert planning.GenerationPlan.from_dict(plan.to_dict()) == plan
    run = da.run(plan, out=tmp_path / "run")
    assert da.inspect(run, verify=True).accepted == 1
    with da.inspect(run, view="designs").records() as records:
        design = next(records)
    assert design.realized.sequence == "AACCGT"
    assert design.requirements[-1]["observed"] == -2
    assert all(row["passed"] for row in design.requirements)


def test_native_reverse_fixed_identity_excludes_equal_string(tmp_path: Path):
    request = planning.DesignSpec(
        parts=[parts.Part("omit", "AAA"), parts.Part("fixed", "AAA")],
        length=planning.Length(maximum=3),
        requirements=[
            planning.Fixed("anchor", "fixed", "reverse"),
            planning.Occurrences("omit", parts.PartSelector(part_ids=("omit",)), max=0),
        ],
    )
    run = da.run(request, out=tmp_path / "run")
    with da.inspect(run, view="designs").records() as records:
        design = next(records)
    assert design.realized.sequence == "TTT"
    assert design.realized.placements[0].feature_id == "fixed"
    assert da.inspect(run, verify=True).accepted == 1


def test_fixed_requirement_static_conflict_fails_before_output(tmp_path: Path):
    request = planning.DesignSpec(
        parts=[parts.Part("a", "AAA")],
        length=planning.Length(maximum=3),
        requirements=[
            planning.Fixed("anchor", "a", "forward"),
            planning.Occurrences("omit", parts.PartSelector(part_ids=("a",)), max=0),
        ],
    )
    with pytest.raises(ValueError, match=r"fixed.*maximum"):
        da.run(request, out=tmp_path / "run")
    assert not (tmp_path / "run").exists()


@pytest.mark.parametrize(
    "sequence,starts,assembly",
    [
        ("GACGTTGCAAGTCTGCAGTACCGAT", (1, 13), False),
        ("ACGTTGCAAGTCTGCAGTACCGATG", (0, 12), False),
        ("ACGTTGCAAGTCGTGCAGTACCGAT", (0, 13), False),
        ("ACGTTGCAAGTCGTGCAGTACCGAT", (0, 13), True),
    ],
)
def test_final_verification_rejects_uncovered_packing_bases(
    sequence: str, starts: tuple[int, int], assembly: bool
):
    """The packed interval must be entirely explained by selected occurrences."""
    model = planning.PlanEvidence(
        planning.DesignSpec(
            parts=(parts.Part("a", "ACGTTGCAAGTC"), parts.Part("b", "TGCAGTACCGAT")),
            length=planning.Length(exact=25)
            if assembly
            else planning.Length(maximum=25),
            assembly=planning.Assembly() if assembly else None,
            strands="single",
        )
    )
    provenance = (
        {
            "assembly": {
                "policy": PADDING_POLICY,
                "side": None,
                "padding_length": 0,
                "packed_start": 0,
                "packed_length": 25,
                "trial": 1,
                "stream_id": "0" * 64,
            }
        }
        if assembly
        else {}
    )
    realized = RealizedArray(
        "probe/default/d1",
        sequence,
        tuple(
            Placement(
                f"p{i}", label, PlacementKind.OTHER, dna, start, Orientation.FORWARD
            )
            for i, (label, dna, start) in enumerate(
                zip(("a", "b"), ("ACGTTGCAAGTC", "TGCAGTACCGAT"), starts, strict=True),
                1,
            )
        ),
        provenance=provenance,
    )
    with pytest.raises(ValueError, match=r"cover.*packed interval"):
        evaluate(realized, model)
