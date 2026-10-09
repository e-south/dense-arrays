"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/background/automaton.py

Deterministic forbidden-word recognition using bounded prefix/failure states.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from collections import deque

from .contracts import ConstructionWork


def transitions(
    patterns: tuple[str, ...], work: ConstructionWork
) -> tuple[tuple[int, ...], ...]:
    """Compile overlaps and suffix matches; -1 marks a forbidden transition."""
    children, failure, forbidden = [{}], [0], [False]
    work.admit("automaton_states")
    for pattern in patterns:
        node = 0
        for base in pattern:
            work.check_time()
            if base not in children[node]:
                work.admit("automaton_states")
                children[node][base] = len(children)
                children.append({})
                failure.append(0)
                forbidden.append(False)
            node = children[node][base]
        forbidden[node] = True
    rows = [[0] * 4 for _ in children]
    pending = deque([0])
    while pending:
        work.check_time()
        node = pending.popleft()
        for index, base in enumerate("ACGT"):
            if base in children[node]:
                child = children[node][base]
                failure[child] = rows[failure[node]][index] if node else 0
                forbidden[child] |= forbidden[failure[child]]
                rows[node][index] = child
                pending.append(child)
            else:
                rows[node][index] = rows[failure[node]][index] if node else 0
    return tuple(
        tuple(-1 if forbidden[child] else child for child in row) for row in rows
    )
