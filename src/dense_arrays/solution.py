"""Solution representation for dense-arrays.

--------------------------------------------------------------------------------
<dense-array project>

Solution representation for dense-arrays.

Module Author(s): Virgile Andreani, Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Self

from .sequence import COMPLEMENT, dispatch_labels, reverse_complement


@dataclass
class DenseArray:
    """Representation of a solution."""

    library: list[str]
    sequence_length: int
    sequence: str
    offsets_fwd: list[int | None]
    offsets_rev: list[int | None]

    def __init__(  # noqa: C901, PLR0912
        self: Self,
        library: list[str],
        sequence_length: int,
        offsets_fwd: list[int | None],
        offsets_rev: list[int | None],
    ) -> None:
        if sequence_length <= 0:
            msg = "sequence_length must be > 0"
            raise ValueError(msg)
        if len(offsets_fwd) != len(library) or len(offsets_rev) != len(library):
            msg = "offsets_fwd and offsets_rev must match library length"
            raise ValueError(msg)

        self.library = list(library)
        self.sequence_length = sequence_length
        self.offsets_fwd = list(offsets_fwd)
        self.offsets_rev = list(offsets_rev)

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
            if not isinstance(offset, int):
                msg = "offsets must be integers or None"
                raise TypeError(msg)
            if offset < 0:
                msg = "offsets must be >= 0"
                raise ValueError(msg)
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

        self.sequence = "".join(sequence_chars)

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
