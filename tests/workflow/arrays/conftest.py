"""Part identities and overlapping, reverse-strand supplied geometries."""

import pytest

from dense_arrays.arrays import ArrayCollection
from dense_arrays.parts import Part
from dense_arrays.realized import Orientation, Placement, PlacementKind, RealizedArray
from dense_arrays.sequence import reverse_complement


@pytest.fixture
def source() -> ArrayCollection:
    """Include equal sequences, repeated parts and reverse-strand core composition."""
    sequence = "ACGTTGCAAGTCCTGA"
    parts = (
        Part(
            "part-a",
            sequence,
            group="regulator",
            core_start=2,
            core_end=13,
            core_orientation="forward",
        ),
        Part("part-b", sequence, group="background"),
        Part(
            "part-r",
            reverse_complement(sequence),
            group="regulator",
            core_start=2,
            core_end=13,
            core_orientation="reverse",
        ),
    )
    first = RealizedArray(
        "array-1",
        "TT" + sequence + "AA",
        (
            Placement(
                "one", "part-a", PlacementKind.TFBS, sequence, 2, Orientation.FORWARD
            ),
            Placement(
                "two", "part-b", PlacementKind.OTHER, sequence, 2, Orientation.FORWARD
            ),
            Placement(
                "repeat", "part-a", PlacementKind.TFBS, sequence, 2, Orientation.FORWARD
            ),
        ),
        provenance={"family": "example"},
    )
    second = RealizedArray(
        "array-2",
        first.sequence,
        (
            Placement(
                "reverse",
                "part-r",
                PlacementKind.TFBS,
                sequence,
                2,
                Orientation.REVERSE,
            ),
        ),
    )
    return ArrayCollection(
        parts=parts, arrays=(first, second), provenance={"dataset": "supplied-arrays"}
    )
