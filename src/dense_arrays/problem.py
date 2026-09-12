"""Validated, immutable motif-packing problem inputs.

Module Author(s): Eric J. South
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import TYPE_CHECKING

from .sequence import VALID_BASES, adjacency_matrix, reverse_complement

if TYPE_CHECKING:
    from collections.abc import Sequence


def discrete_integer(value: object, name: str, *, minimum: int | None = None) -> int:
    """Normalize an integer without truncation or boolean coercion.

    Returns
    -------
    int
        The validated integer.

    Raises
    ------
    ValueError
        If the value is not an integer or falls below the minimum.
    """
    if isinstance(value, bool) or not isinstance(value, Integral):
        msg = f"{name} must be an integer (not bool)"
        raise ValueError(msg)  # noqa: TRY004 - invalid configuration is a ValueError
    normalized = int(value)
    if minimum is not None and normalized < minimum:
        msg = f"{name} must be >= {minimum}"
        raise ValueError(msg)
    return normalized


def motif_library(library: Sequence[str]) -> tuple[str, ...]:
    """Take a validated snapshot of uppercase DNA motif entries.

    Returns
    -------
    tuple[str, ...]
        The immutable library, retaining duplicate entry identities.

    Raises
    ------
    ValueError
        If the library or a motif is empty or has an invalid alphabet.
    """
    if isinstance(library, (str, bytes)) or not library:
        msg = "library must contain at least one motif"
        raise ValueError(msg)
    for i, motif in enumerate(library):
        if not isinstance(motif, str) or not motif:
            msg = f"motif at index {i} must be a non-empty string"
            raise ValueError(msg)
        invalid = set(motif) - VALID_BASES
        if invalid:
            msg = (
                f"motif at index {i} contains invalid bases: {sorted(invalid)}. "
                "Use uppercase A/C/G/T only."
            )
            raise ValueError(msg)
    return tuple(library)


@dataclass(frozen=True)
class PackingProblem:
    """The immutable sequence library, length bound, and strand policy."""

    library: tuple[str, ...]
    sequence_length: int
    strands: str

    @classmethod
    def create(
        cls, library: Sequence[str], sequence_length: int, strands: str
    ) -> PackingProblem:
        """Validate caller inputs before model allocation.

        Returns
        -------
        PackingProblem
            A snapshot independent of caller-owned collections.

        Raises
        ------
        ValueError
            If a length, motif, or strand policy is invalid.
        """
        length = discrete_integer(sequence_length, "sequence_length", minimum=1)
        if not isinstance(strands, str) or strands not in {"single", "double"}:
            msg = "strands must be single or double"
            raise ValueError(msg)
        return cls(motif_library(library), length, strands)

    @property
    def oriented_library(self) -> tuple[str, ...]:
        """Motifs in node order: forward entries, then reverse complements."""
        if self.strands == "single":
            return self.library
        return self.library + tuple(reverse_complement(m) for m in self.library)

    @property
    def adjacency(self) -> tuple[tuple[int, ...], ...]:
        """Pairwise path-entry shifts, retaining repeated-entry semantics."""
        return tuple(tuple(row) for row in adjacency_matrix(self.oriented_library))
