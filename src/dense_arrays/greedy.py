"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/greedy.py

Greedy packing with explicit entry, orientation, and occurrence identity.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import time
from typing import TYPE_CHECKING

from .errors import InfeasibleError
from .solution import DenseArray

if TYPE_CHECKING:
    from .problem import PackingProblem


def realize_greedy(
    problem: PackingProblem, *, deadline: float | None = None
) -> DenseArray:
    """Choose the best of deterministic greedy paths, starting at each fitting node.

    Each extension uses the exact path-entry shift. A library entry is selected
    at most once, on one strand. Incidental substrings never count as entries.
    The heuristic maximizes selected entries across its starts and breaks ties
    by shorter realized length; it does not prove an optimal packing.

    An optional monotonic-clock deadline is checked between starts and
    extensions. Adjacency construction and one extension scan are cooperative
    units, not hard wall-clock bounds.

    Returns
    -------
    DenseArray
        A contiguous realization of the chosen path.

    Raises
    ------
    InfeasibleError
        If no library entry fits the length bound on either allowed strand.
    TimeoutError
        If the cooperative deadline is reached before all starts complete.
    """
    _check_deadline(deadline)
    nb_motifs = len(problem.library)
    adjacency = problem.adjacency
    library = problem.oriented_library
    candidates = sorted(
        (i for i, motif in enumerate(library) if len(motif) <= problem.sequence_length),
        key=lambda i: (len(library[i]), i),
    )
    if not candidates:
        msg = "No feasible solution was found."
        raise InfeasibleError(msg)
    best_path: list[tuple[int, int]] = []
    best_length = problem.sequence_length + 1
    for start in candidates:
        _check_deadline(deadline)
        path = [(0, start)]
        used = {start % nb_motifs}
        length = len(library[start])
        while True:
            _check_deadline(deadline)
            last_offset, last = path[-1]
            extensions = [
                (last_offset + adjacency[last][node] + len(library[node]), node)
                for node in candidates
                if node % nb_motifs not in used
                and last_offset + adjacency[last][node] + len(library[node])
                <= problem.sequence_length
            ]
            if not extensions:
                break
            length, node = min(extensions)
            path.append((length - len(library[node]), node))
            used.add(node % nb_motifs)
        if (len(path), -length) > (len(best_path), -best_length):
            best_path, best_length = path, length
    offsets_fwd = [None] * nb_motifs
    offsets_rev = [None] * nb_motifs
    for offset, node in best_path:
        offsets = offsets_fwd if node < nb_motifs else offsets_rev
        offsets[node % nb_motifs] = offset
    return DenseArray(
        problem.library, problem.sequence_length, offsets_fwd, offsets_rev
    )


def _check_deadline(deadline: float | None) -> None:
    if deadline is not None and time.monotonic() >= deadline:
        msg = "greedy search reached its cooperative deadline"
        raise TimeoutError(msg)
