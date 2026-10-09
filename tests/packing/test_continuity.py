"""Independent rank-feasibility checks for the subtour-elimination subsystem.

Author: Eric J. South.
"""

from itertools import permutations, product

from ortools.linear_solver import linear_solver_pb2

from dense_arrays import Optimizer


def test_continuity_rows_admit_exactly_acyclic_directed_edge_sets():
    """A selected edge must increase rank; a directed cycle has no valid ranking."""
    optimizer = Optimizer(["AA", "CC", "GG"], 6, strands="single")
    optimizer.build_model("CBC")
    model = linear_solver_pb2.MPModelProto()
    optimizer.model.ExportModelToProto(model)
    names = {index: variable.name for index, variable in enumerate(model.variable)}
    ranks = {
        index: int(variable.name[2:-1])
        for index, variable in enumerate(model.variable)
        if variable.name.startswith("u[")
    }
    rows = [
        row
        for row in model.constraint
        if any(index in ranks for index in row.var_index)
    ]
    assert len(ranks) == 3
    assert len(rows) == 6
    assert all(
        model.variable[index].is_integer
        and model.variable[index].lower_bound == 1
        and model.variable[index].upper_bound == 3
        for index in ranks
    )
    edges = tuple(permutations(range(3), 2))
    for enabled in product((0, 1), repeat=len(edges)):
        selected = {
            edge for edge, present in zip(edges, enabled, strict=True) if present
        }
        # A finite graph is acyclic exactly when some topological order admits it.
        expected = any(
            all(order.index(i) < order.index(j) for i, j in selected)
            for order in permutations(range(3))
        )
        edge_values = {
            f"X[{i},{j}]": present
            for (i, j), present in zip(edges, enabled, strict=True)
        }
        feasible = False
        for assignment in product(range(1, 4), repeat=3):
            values = edge_values | {
                f"u[{node}]": rank for node, rank in enumerate(assignment)
            }
            if all(
                row.lower_bound
                <= sum(
                    coefficient * values[names[index]]
                    for index, coefficient in zip(
                        row.var_index, row.coefficient, strict=True
                    )
                )
                <= row.upper_bound
                for row in rows
            ):
                feasible = True
                break
        assert feasible == expected, selected
