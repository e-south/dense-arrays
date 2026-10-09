"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/generation/batches/constrained.py

Resolve simultaneous sequence/core uniqueness and group capacity by flow.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from collections import Counter
from dataclasses import dataclass, field

from ortools.graph.python.min_cost_flow import SimpleMinCostFlow

from dense_arrays.parts import Part
from dense_arrays.planning import BatchSampling
from dense_arrays.planning.batches.validation import validate_selected


@dataclass
class _Network:
    flow: SimpleMinCostFlow = field(default_factory=SimpleMinCostFlow)
    nodes: dict[tuple[str, object], int] = field(default_factory=dict)

    def node(self, kind: str, key: object) -> int:
        return self.nodes.setdefault((kind, key), len(self.nodes) + 2)


def select_constrained(
    ordered: list[Part], fixed: list[Part], policy: BatchSampling
) -> list[Part]:
    """Select a feasible full-size subset, then prefer balance and seeded ranks."""
    validate_selected(fixed, policy)
    need = policy.size - len(fixed)
    if need == 0:
        return fixed
    ids = {p.part_id for p in fixed}
    sequences = {p.sequence for p in fixed} if policy.unique_sequences else set()
    cores = (
        {(p.group, p.core_sequence) for p in fixed} if policy.unique_cores else set()
    )
    candidates = [
        p
        for p in ordered
        if p.part_id not in ids
        and p.sequence not in sequences
        and (p.group, p.core_sequence) not in cores
    ]
    network = _Network()
    arcs = []
    left_seen, right_seen = set(), set()
    groups = Counter()
    for rank, part in enumerate(candidates):
        left = network.node(
            "sequence", part.sequence if policy.unique_sequences else part.part_id
        )
        right = network.node(
            "core",
            (part.group, part.core_sequence) if policy.unique_cores else part.part_id,
        )
        group = network.node("group", part.group)
        if left not in left_seen:
            network.flow.add_arc_with_capacity_and_unit_cost(0, left, 1, 0)
            left_seen.add(left)
        arcs.append(
            network.flow.add_arc_with_capacity_and_unit_cost(left, right, 1, rank)
        )
        if right not in right_seen:
            network.flow.add_arc_with_capacity_and_unit_cost(right, group, 1, 0)
            right_seen.add(right)
            groups[part.group] += 1
    _group_capacities(
        network, groups, Counter(p.group for p in fixed), policy, len(candidates)
    )
    network.flow.set_node_supply(0, need)
    network.flow.set_node_supply(1, -need)
    status = network.flow.solve()
    if status == network.flow.INFEASIBLE:
        msg = "batch size cannot be satisfied under declared uniqueness and group caps"
        raise ValueError(msg)
    if status != network.flow.OPTIMAL:
        msg = f"batch selection backend failed: {status}"
        raise RuntimeError(msg)
    selected = fixed + [
        p for p, arc in zip(candidates, arcs, strict=True) if network.flow.flow(arc)
    ]
    validate_selected(selected, policy)
    return selected


def _group_capacities(
    network: _Network,
    groups: Counter,
    fixed: Counter,
    policy: BatchSampling,
    candidates: int,
) -> None:
    """Apply group caps and convex balance costs; rank costs break objective ties."""
    scale = candidates * policy.size + 1
    for group, available in groups.items():
        capacity = (
            available
            if policy.max_per_group is None
            else min(available, policy.max_per_group - fixed[group])
        )
        node = network.node("group", group)
        if policy.strategy == "uniform":
            network.flow.add_arc_with_capacity_and_unit_cost(node, 1, capacity, 0)
        else:
            for i in range(capacity):
                network.flow.add_arc_with_capacity_and_unit_cost(
                    node, 1, 1, (2 * (fixed[group] + i) + 1) * scale
                )
