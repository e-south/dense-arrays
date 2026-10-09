"""Supported controls and typed evidence at the packing boundary.

Author: Eric J. South.
"""

from types import SimpleNamespace

import pytest
from ortools.linear_solver import pywraplp

from dense_arrays import Optimizer
from dense_arrays.solver import SolverControls, SolveStatus


def test_real_cbc_reports_proof_and_exhaustion():
    optimizer = Optimizer(["AAA"], 3, "single")
    optimizer.build_model(controls=SolverControls(time_limit_seconds=1))
    report = optimizer.solve_report()
    assert report.status is SolveStatus.OPTIMAL
    assert report.solution.sequence == "AAA"
    assert report.proof_scope == "offered_packing_model"
    optimizer.forbid(report.solution)
    exhausted = optimizer.solve_report()
    assert exhausted.status is SolveStatus.INFEASIBLE
    assert exhausted.solution is None


@pytest.mark.parametrize(
    "raw,expected",
    [
        (pywraplp.Solver.FEASIBLE, SolveStatus.UNPROVEN),
        (pywraplp.Solver.NOT_SOLVED, SolveStatus.UNKNOWN),
        (pywraplp.Solver.ABNORMAL, SolveStatus.BACKEND_ERROR),
        (99, SolveStatus.UNKNOWN),
    ],
)
def test_reports_do_not_invent_a_timeout_cause(raw: int, expected: SolveStatus):
    optimizer = Optimizer(["AAA"], 3, "single")
    optimizer.model = SimpleNamespace(Solve=lambda: raw)
    report = optimizer.solve_report()
    assert report.status is expected
    assert report.backend_status == raw
    assert report.termination_reason == "unknown"
    assert report.solution is None
    assert report.proof_scope is None


def test_backend_exception_is_a_report():
    def fail() -> int:
        msg = "external failure"
        raise OSError(msg)

    optimizer = Optimizer(["AAA"], 3, "single")
    optimizer.model = SimpleNamespace(Solve=fail)
    report = optimizer.solve_report()
    assert report.status is SolveStatus.BACKEND_ERROR
    assert report.backend_status is None
    assert report.termination_reason == "backend_exception"


def test_controls_rejected_without_replacing_model():
    optimizer = Optimizer(["AAA"], 3, "single")
    optimizer.build_model()
    original = optimizer.model
    with pytest.raises(ValueError, match=r"threads.*CBC"):
        optimizer.build_model(controls=SolverControls(threads=2))
    assert optimizer.model is original
    assert optimizer.solve().sequence == "AAA"


@pytest.mark.parametrize("value", [True, 0, -1, float("inf"), float("nan"), "1"])
def test_invalid_time_limits_fail_before_solver_creation(value: object):
    with pytest.raises(ValueError, match="time_limit_seconds"):
        SolverControls(time_limit_seconds=value)


@pytest.mark.parametrize("value", [True, 0, -1, 1.5, "2"])
def test_invalid_thread_counts(value: object):
    with pytest.raises(ValueError, match="threads"):
        SolverControls(threads=value)


def test_time_control_reaches_backend(monkeypatch: pytest.MonkeyPatch):
    model = pywraplp.Solver.CreateSolver("CBC")
    observed = []
    monkeypatch.setattr(model, "SetTimeLimit", observed.append)
    monkeypatch.setattr(pywraplp.Solver, "CreateSolver", lambda _name: model)
    Optimizer(["AAA"], 3).build_model(controls=SolverControls(time_limit_seconds=0.125))
    assert observed == [125]


@pytest.mark.parametrize("method", ["optimal", "solutions", "solutions_diverse"])
def test_convenience_methods_apply_time_controls(
    method: str, monkeypatch: pytest.MonkeyPatch
):
    """A bounded request must reach the actual backend through each public method."""
    model = pywraplp.Solver.CreateSolver("CBC")
    observed = []
    monkeypatch.setattr(model, "SetTimeLimit", observed.append)
    monkeypatch.setattr(pywraplp.Solver, "CreateSolver", lambda _name: model)
    optimizer = Optimizer(["ACGTTGCAAGTCCTGA"], 16, "single")
    result = getattr(optimizer, method)(controls=SolverControls(time_limit_seconds=2))
    solution = result if method == "optimal" else next(result)
    assert solution.sequence == "ACGTTGCAAGTCCTGA"
    assert observed == [2000]


@pytest.mark.parametrize("method", ["optimal", "solutions", "solutions_diverse"])
def test_convenience_methods_reject_unsupported_controls(method: str):
    """CBC thread controls must fail instead of being ignored by a wrapper."""
    optimizer = Optimizer(["ACGTTGCAAGTCCTGA"], 16, "single")
    controls = SolverControls(threads=2)
    if method == "optimal":
        with pytest.raises(ValueError, match=r"threads.*CBC"):
            optimizer.optimal(controls=controls)
    else:
        iterator = getattr(optimizer, method)(controls=controls)
        with pytest.raises(ValueError, match=r"threads.*CBC"):
            next(iterator)
    assert optimizer.model is None
