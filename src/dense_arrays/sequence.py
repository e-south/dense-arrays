"""Sequence utilities for dense-arrays.

--------------------------------------------------------------------------------
<dense-array project>

Sequence utilities for dense-arrays.

Module Author(s): Virgile Andreani, Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

from __future__ import annotations

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "-": "-"}
VALID_BASES = {"A", "C", "G", "T"}


def shift_metric(motifa: str, motifb: str) -> int:
    """Compute how much we have to shift `motifb` to match the end of `motifa`.

    Example
    -------
    shift_metric("ATGCATTA", "CATTATG") == 3 because

        motifa: ATGCATTA
        motifb:    CATTATG
        shift : 0123

    and shift_metric("ATGCATTA", "TATGA") == 6 because

        motifa: ATGCATTA
        motifb:       TATGA
        shift : 0123456

    Note: we only consider shifts such that the shifted `motifb` overhangs from
    `motifa`.  If `motifb` is contained inside `motifa`, it will not be counted:

    shift_metric("ATGTTAACT", "TTAA") == 8 because

        motifa: ATGTTAACT
        motifb:    TTAA         # not allowed
        motifb:         TTAA    # okay
        shift : 012345678

    Returns
    -------
    shift : int
        The result.
    """
    if motifa == motifb:
        return len(motifa)
    for shift in range(max(len(motifa) - len(motifb), 0), len(motifa)):
        if motifa[shift:] == motifb[: len(motifa) - shift]:
            return shift
    return len(motifa)


def adjacency_matrix(motifs: list[str]) -> list[list[int]]:
    """Return the matrix A_ij such that A_ij = shift_metric(motifs[i], motifs[j]).

    Parameters
    ----------
    motifs
        List of motifs.

    Returns
    -------
    adj :
        Adjacency matrix.
    """
    return [[shift_metric(motifa, motifb) for motifb in motifs] for motifa in motifs]


def reverse_complement(sequence: str) -> str:
    """
    Return the reverse complement of the sequence.

    Parameters
    ----------
    sequence
        String composed of ATGC characters.

    Returns
    -------
    rev_comp :
        Reverse complement of the sequence.
    """
    return "".join(COMPLEMENT[c] for c in sequence[::-1])


def dispatch_labels(
    library: list[str],
    offsets: list[int | None],
    *,
    rev: bool,
) -> list[str]:
    """Render motif labels into aligned rows for printing.

    Returns
    -------
    lines
        One or more lines representing the motif placements.
    """
    lines: list[str] = []
    order = sorted((o, i) for i, o in enumerate(offsets) if o is not None)
    for offset, i in order:
        motif = library[i][::-1] if rev else library[i]
        for iline, line in enumerate(lines):
            if not line or len(line) < offset:
                lines[iline] += " " * (offset - len(line)) + motif
                break
        else:
            line = " " * offset + motif
            lines.append(line)
    return lines


def take_best_run(
    sequence: str,
    sequence_length: int,
    library: list[str],
    strands: str,
) -> str:
    """Return the subsequence maximizing motif coverage.

    Raises
    ------
    AssertionError
        If no subsequence can be selected.

    Returns
    -------
    subsequence
        The best-scoring subsequence.
    """
    max_nb_motifs = -1
    offset_max_nb_motifs = None
    for offset in range(max(1, len(sequence) - sequence_length + 1)):
        subseq = sequence[offset : offset + sequence_length]
        if strands == "single":
            nb_motifs = sum(motif in subseq for motif in library)
        else:
            nb_motifs = sum(
                (motif in subseq) or (reverse_complement(motif) in subseq)
                for motif in library
            )
        if nb_motifs > max_nb_motifs:
            max_nb_motifs = nb_motifs
            offset_max_nb_motifs = offset
    if offset_max_nb_motifs is None:
        msg = "This should not have happened."
        raise AssertionError(msg)
    subseq = sequence[offset_max_nb_motifs : offset_max_nb_motifs + sequence_length]
    return subseq
