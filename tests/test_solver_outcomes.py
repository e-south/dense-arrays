"""Solver outcome contracts.

Author: Eric J. South.
"""

from types import SimpleNamespace

import pytest
from ortools.linear_solver import pywraplp

from dense_arrays import DenseArray, Optimizer


@pytest.mark.parametrize(
    "status, error_name",
    [
        (pywraplp.Solver.INFEASIBLE, "InfeasibleError"),
        (pywraplp.Solver.FEASIBLE, "UnprovenSolutionError"),
        (pywraplp.Solver.ABNORMAL, "SolverBackendError"),
        (pywraplp.Solver.UNBOUNDED, "SolverBackendError"),
        (pywraplp.Solver.NOT_SOLVED, "SolverBackendError"),
        (99, "SolverBackendError"),
    ],
)
def test_solver_outcome_types(status: int, error_name: str):
    optimizer = Optimizer(["AAA"], 3, "single")
    optimizer.model = SimpleNamespace(Solve=lambda: status)
    with pytest.raises((RuntimeError, ValueError)) as caught:
        optimizer.solve()
    assert type(caught.value).__name__ == error_name


@pytest.mark.parametrize("method", ["solutions", "solutions_diverse", "optimal"])
def test_abnormal_status_propagates(monkeypatch: pytest.MonkeyPatch, method: str):
    optimizer = Optimizer(["AAA"], 3, "single")
    monkeypatch.setattr(optimizer, "build_model", lambda *_a, **_kw: None)
    optimizer.model = SimpleNamespace(Solve=lambda: pywraplp.Solver.ABNORMAL)
    invoke = (
        optimizer.optimal
        if method == "optimal"
        else lambda: next(getattr(optimizer, method)())
    )
    with pytest.raises(RuntimeError, match="abnormal"):
        invoke()


@pytest.mark.parametrize("method", ["solutions", "solutions_diverse"])
def test_failure_after_one_result_propagates(
    monkeypatch: pytest.MonkeyPatch, method: str
):
    optimizer = Optimizer(["AAA"], 3, "single")
    monkeypatch.setattr(optimizer, "build_model", lambda *_a, **_kw: None)
    monkeypatch.setattr(optimizer, "forbid", lambda _solution: None)
    solution = DenseArray(["AAA"], 3, [0], [None])
    results = iter([solution, ValueError("invalid result construction")])

    def solve() -> DenseArray:
        result = next(results)
        if isinstance(result, Exception):
            raise result
        return result

    monkeypatch.setattr(optimizer, "solve", solve)
    iterator = getattr(optimizer, method)()
    assert next(iterator) is solution
    with pytest.raises(ValueError, match="invalid result construction"):
        next(iterator)


def test_rejected_solver_options_leave_model_unchanged():
    optimizer = Optimizer(["AAA"], 3, "single")
    optimizer.build_model("CBC")
    original = optimizer.model
    with pytest.raises(ValueError, match="option"):
        optimizer.build_model("CBC", solver_options=["unsupported=1"])
    assert optimizer.model is original


def test_cbc_success_and_exhaustion():
    optimizer = Optimizer(["AAA"], 3, "single")
    assert [solution.sequence for solution in optimizer.solutions("CBC")] == ["AAA"]
    assert list(Optimizer(["AAA"], 2, "single").solutions("CBC")) == []


def test_reconstruction_errors_are_typed(monkeypatch: pytest.MonkeyPatch):
    optimizer = Optimizer(["AAA"], 3, "single")
    optimizer.build_model("CBC")

    def invalid(*_args: object) -> None:
        msg = "broken reconstructed offsets"
        raise ValueError(msg)

    monkeypatch.setattr("dense_arrays.optimizer.DenseArray", invalid)
    with pytest.raises(RuntimeError, match="broken reconstructed offsets") as caught:
        optimizer.solve()
    assert type(caught.value).__name__ == "InvalidSolverResultError"


def test_backend_execution_error_is_typed():
    optimizer = Optimizer(["AAA"], 3, "single")

    def fail() -> int:
        msg = "backend execution failed"
        raise ValueError(msg)

    optimizer.model = SimpleNamespace(Solve=fail)
    with pytest.raises(RuntimeError, match="backend execution failed") as caught:
        optimizer.solve()
    assert type(caught.value).__name__ == "SolverBackendError"


@pytest.mark.parametrize(
    "selected",
    [
        {(-1, 0)},
        {(-1, 0), (0, -1), (-1, 1), (1, -1)},
        {(-1, 0), (0, -1), (1, 2), (2, 1)},
    ],
)
def test_invalid_solver_path_is_typed(selected: set[tuple[int, int]]):
    optimizer = Optimizer(["AAA", "CCC", "GGG"], 9, "single")
    arcs = {
        (i, j): SimpleNamespace(
            solution_value=lambda edge=(i, j): float(edge in selected)
        )
        for i in range(-1, 3)
        for j in range(-1, 3)
        if i != j
    }
    optimizer.model = SimpleNamespace(Solve=lambda: pywraplp.Solver.OPTIMAL, X=arcs)
    with pytest.raises(RuntimeError) as caught:
        optimizer.solve()
    assert type(caught.value).__name__ == "InvalidSolverResultError"


def test_solver_option_collection_is_explicit(monkeypatch: pytest.MonkeyPatch):
    optimizer = Optimizer(["AAA"], 3, "single")
    allocations = []

    def create(_solver: str) -> None:
        allocations.append(True)

    monkeypatch.setattr(pywraplp.Solver, "CreateSolver", create)
    with pytest.raises(ValueError, match="solver_options"):
        optimizer.build_model("CBC", solver_options="threads=1")
    assert not allocations


@pytest.mark.parametrize("method", ["solutions", "solutions_diverse"])
@pytest.mark.parametrize(
    "status",
    [pywraplp.Solver.FEASIBLE, pywraplp.Solver.ABNORMAL, pywraplp.Solver.NOT_SOLVED],
)
def test_backend_status_after_a_real_result_remains_visible(
    monkeypatch: pytest.MonkeyPatch, method: str, status: int
):
    optimizer = Optimizer(["AAA"], 3, "single")
    iterator = getattr(optimizer, method)("CBC")
    assert next(iterator).sequence == "AAA"
    monkeypatch.setattr(optimizer.model, "Solve", lambda: status)
    with pytest.raises(RuntimeError):
        next(iterator)


def test_backend_creation_error_is_typed(monkeypatch: pytest.MonkeyPatch):
    optimizer = Optimizer(["AAA"], 3, "single")

    def fail(_solver: str) -> None:
        msg = "backend initialization failed"
        raise ValueError(msg)

    monkeypatch.setattr(pywraplp.Solver, "CreateSolver", fail)
    with pytest.raises(
        RuntimeError, match=r"CBC.*backend initialization failed"
    ) as caught:
        optimizer.build_model("CBC")
    assert type(caught.value).__name__ == "SolverBackendError"
    assert optimizer.model is None
