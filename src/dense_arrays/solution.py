"""Immutable, contiguous motif-placement results for dense arrays.

Module Author(s): Virgile Andreani, Eric J. South
Dunlop Lab
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Self

from .problem import discrete_integer, motif_library
from .sequence import COMPLEMENT, dispatch_labels, reverse_complement


@dataclass(frozen=True)
class DenseArray:
    """An immutable, contiguous collection of compatible motif placements.

    Parameters
    ----------
    library
        Nonempty uppercase A/C/G/T motif entries, copied on construction.
    sequence_length
        Positive integer maximum sequence length, excluding booleans.
    offsets_fwd, offsets_rev
        One nonnegative integer offset or None per library entry. An entry may
        select at most one orientation. Public views return defensive copies.

    Notes
    -----
    Compatible contained and non-maximal overlaps are allowed. This representation
    describes placements; it does not itself assert an exact solver path.

    Raises
    ------
    ValueError
        If motifs, lengths, offsets, or orientation choices are invalid, or if
        placements conflict, leave gaps, exceed the bound, or select no entries.
    """

    _library: tuple[str, ...]
    sequence_length: int
    sequence: str
    _offsets_fwd: tuple[int | None, ...]
    _offsets_rev: tuple[int | None, ...]

    def __init__(  # noqa: C901, PLR0912
        self: Self,
        library: list[str],
        sequence_length: int,
        offsets_fwd: list[int | None],
        offsets_rev: list[int | None],
    ) -> None:
        sequence_length = discrete_integer(
            sequence_length, "sequence_length", minimum=1
        )
        library = list(motif_library(library))
        if len(offsets_fwd) != len(library) or len(offsets_rev) != len(library):
            msg = "offsets_fwd and offsets_rev must match library length"
            raise ValueError(msg)

        normalized_fwd = tuple(
            None if o is None else discrete_integer(o, "offsets", minimum=0)
            for o in offsets_fwd
        )
        normalized_rev = tuple(
            None if o is None else discrete_integer(o, "offsets", minimum=0)
            for o in offsets_rev
        )
        object.__setattr__(self, "_library", tuple(library))
        object.__setattr__(self, "sequence_length", sequence_length)
        object.__setattr__(self, "_offsets_fwd", normalized_fwd)
        object.__setattr__(self, "_offsets_rev", normalized_rev)

        for fwd, rev in zip(offsets_fwd, offsets_rev, strict=True):
            if fwd is not None and rev is not None:
                msg = "A library entry can have only one selected orientation"
                raise ValueError(msg)

        placements: list[tuple[int, str]] = []
        for i, offset in enumerate(self.offsets_fwd):
            if offset is None:
                continue
            placements.append((offset, self.library[i]))
        for i, offset in enumerate(self.offsets_rev):
            if offset is None:
                continue
            placements.append((offset, reverse_complement(self.library[i])))

        if not placements:
            msg = "solution must contain at least one motif"
            raise ValueError(msg)

        max_end = 0
        for offset, motif in placements:
            end = offset + len(motif)
            if end > sequence_length:
                msg = "motif extends beyond sequence_length"
                raise ValueError(msg)
            max_end = max(max_end, end)

        sequence_chars: list[str | None] = [None] * max_end
        for offset, motif in placements:
            for i, base in enumerate(motif):
                pos = offset + i
                existing = sequence_chars[pos]
                if existing is None:
                    sequence_chars[pos] = base
                elif existing != base:
                    msg = f"Offsets conflict at position {pos}"
                    raise ValueError(msg)

        if any(base is None for base in sequence_chars):
            msg = "Offsets leave gaps; sequence must be contiguous"
            raise ValueError(msg)

        object.__setattr__(self, "sequence", "".join(sequence_chars))

    @property
    def library(self) -> list[str]:
        """A defensive copy of the original motif library."""
        return list(self._library)

    @property
    def offsets_fwd(self) -> list[int | None]:
        """A defensive copy of selected forward offsets."""
        return list(self._offsets_fwd)

    @property
    def offsets_rev(self) -> list[int | None]:
        """A defensive copy of selected reverse-complement offsets."""
        return list(self._offsets_rev)

    def offset_indices_in_order(self: Self) -> list[tuple[int, int]]:
        """
        List the motifs in the solution by ascending offset.

        Returns
        -------
        offset_indices :
            Each element represents `(offset, index)` where `offset` is the
            offset where the motif starts and `index` is its index in the motif library.
        """
        order_fwd = [
            (offset, i)
            for i, offset in enumerate(self.offsets_fwd)
            if offset is not None
        ]
        order_rev = [
            (offset, i + len(self.library))
            for i, offset in enumerate(self.offsets_rev)
            if offset is not None
        ]
        # We sort by offset first and then by motif length:
        # If two motifs have the same index, we want the shortest first
        return sorted(
            order_fwd + order_rev,
            key=lambda o_i: (o_i[0], len(self.library[o_i[1] % len(self.library)])),
        )

    @property
    def nb_motifs(self: Self) -> int:
        """Number of motifs that fit in this solution."""
        nb_fwd = sum(offset is not None for offset in self.offsets_fwd)
        nb_rev = sum(offset is not None for offset in self.offsets_rev)
        return nb_fwd + nb_rev

    @property
    def compression_ratio(self: Self) -> float:
        """Compression ratio, i.e. `length of motifs in solution / sequence length`."""
        total_length = sum(
            len(motif)
            for motif, fwd, rev in zip(
                self.library,
                self.offsets_fwd,
                self.offsets_rev,
                strict=True,
            )
            if fwd is not None or rev is not None
        )
        return total_length / self.sequence_length

    def __str__(self: Self) -> str:
        """
        Build a string that visually represents the solution.

        Returns
        -------
        s : str
            The solution as a string.
        """
        sequence = self.sequence + "-" * (self.sequence_length - len(self.sequence))
        seq_rev = "".join(COMPLEMENT[c] for c in sequence)
        lines_fwd = dispatch_labels(self.library, self.offsets_fwd, rev=False)
        lines_rev = dispatch_labels(self.library, self.offsets_rev, rev=True)

        s_fwd = "--> " + "\n--> ".join([*lines_fwd[::-1], sequence])
        s_rev = "<-- " + "\n<-- ".join([seq_rev, *lines_rev])

        return s_fwd + "\n" + s_rev
