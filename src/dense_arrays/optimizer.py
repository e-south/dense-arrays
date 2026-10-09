"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/optimizer.py

Optimization and enumeration of dense motif-packing paths.

Module Author(s): Virgile Andreani, Eric J. South
Maintainer(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import itertools as it
import math
from numbers import Real
from typing import TYPE_CHECKING, Self

if TYPE_CHECKING:
    from collections.abc import Iterator

from ortools.linear_solver import pywraplp

from .constraints import (
    CountConstraint,
    CoverageConstraint,
    FixedOccurrence,
    PromoterConstraint,
    RegulatorRequirements,
    SpacingConstraint,
    _interval,
    _normalize_min_counts,
    _normalize_min_required,
    _normalize_regulator_mapping,
    count_bounds,
    occurrence_indices,
)
from .errors import (
    InfeasibleError,
    InvalidSolverResultError,
    SolverBackendError,
    UnprovenSolutionError,
)
from .greedy import realize_greedy
from .model import ModelConfiguration, build_solver_model, part_usage_weights
from .problem import PackingProblem, discrete_integer
from .solution import DenseArray
from .solver import SolverControls, SolveReport, SolverIdentity, SolveStatus

_INTEGER_TOLERANCE = 1e-6


class Optimizer:
    """Configure and solve a motif-packing problem.

    Parameters
    ----------
    library
        Nonempty uppercase A/C/G/T motifs. Repeated strings remain distinct
        library entries. Caller collections and public views are copied.
    sequence_length
        Positive integer length bound; booleans and fractional values are invalid.
    strands
        ``single`` selects forward entries; ``double`` permits either orientation
        of each entry, at most once. This problem configuration is immutable.
    length_mode
        ``maximum`` retains the packing bound. ``exact`` requires the packed
        sequence itself to attain it, without generating padding.

    Raises
    ------
    ValueError
        If a motif, length, or strand policy is invalid.
    """

    def __init__(
        self: Self,
        library: list[str],
        sequence_length: int,
        strands: str = "double",
        *,
        length_mode: str = "maximum",
    ) -> None:
        self._problem = PackingProblem.create(
            library, sequence_length, strands, length_mode=length_mode
        )
        self._adjacency = self._problem.adjacency
        self._promoters: tuple[PromoterConstraint, ...] = ()
        self._regulator_constraints: RegulatorRequirements | None = None
        self._count_constraints: tuple[CountConstraint, ...] = ()
        self._coverage_constraints: tuple[CoverageConstraint, ...] = ()
        self._fixed_occurrences: tuple[FixedOccurrence, ...] = ()
        self._spacing_constraints: tuple[SpacingConstraint, ...] = ()
        self.model = None
        self._model_modified = False
        self._ilefts: tuple[int, ...] = ()
        self._irights: tuple[int, ...] = ()

    @property
    def library(self) -> list[str]:
        """A defensive copy of the immutable motif library."""
        return list(self._problem.library)

    @property
    def sequence_length(self) -> int:
        """The immutable maximum array length."""
        return self._problem.sequence_length

    @property
    def strands(self) -> str:
        """The immutable strand policy."""
        return self._problem.strands

    @property
    def length_mode(self) -> str:
        """The immutable maximum/exact interpretation of the packing length."""
        return self._problem.length_mode

    @property
    def adjacency_matrix(self) -> list[list[int]]:
        """A defensive view of cached path-entry shifts."""
        return [list(row) for row in self._adjacency]

    @property
    def promoters(self) -> list[PromoterConstraint]:
        """A defensive view of immutable promoter requirements."""
        return list(self._promoters)

    @property
    def ilefts(self) -> list[int]:
        """A defensive view of entries with left-side preference."""
        return list(self._ilefts)

    @property
    def irights(self) -> list[int]:
        """A defensive view of entries with right-side preference."""
        return list(self._irights)

    def _ensure_model_not_built(self: Self, action: str) -> None:
        if self.model is not None:
            msg = (
                f"{action} must be added before build_model(); "
                "create a new Optimizer to change them."
            )
            raise RuntimeError(msg)

    def add_promoter_constraints(
        self: Self,
        *,
        upstream: str,
        downstream: str,
        upstream_pos: int | tuple[int | None, int | None] | None = None,
        downstream_pos: int | tuple[int | None, int | None] | None = None,
        spacer_length: int | tuple[int | None, int | None] | None = None,
    ) -> None:
        """
        Add a promoter constraint to the optimization problem.

        Parameters
        ----------
        upstream
            The upstream element (typically -35). Must appear in the library.
        downstream
            The downstream element (typically -10). Must appear in the library.
        upstream_pos
            Position for the upstream element, or tuple (min, max).
        downstream_pos
            Position for the downstream element, or tuple (min, max).
        spacer_length
            Integer spacer or ordered two-bound tuple. Negative spacers permit
            intentional overlap. Position bounds must be nonnegative integers;
            a None bound leaves that side unbounded. Booleans are invalid.

        Raises
        ------
        ValueError
            If motifs or position/spacer values are invalid.
        RuntimeError
            If the model has already been built.
        """
        self._ensure_model_not_built("Promoter constraints")
        upstream_index = self._find_motif_index(upstream)
        downstream_index = self._find_motif_index(downstream, avoid=upstream_index)

        constraint = PromoterConstraint(
            upstream_index=upstream_index,
            downstream_index=downstream_index,
            upstream_pos=upstream_pos,
            downstream_pos=downstream_pos,
            spacer_length=spacer_length,
        )
        self._promoters += (constraint,)

    def _find_motif_index(self: Self, motif: str, avoid: int | None = None) -> int:
        upstream_indices = {p.upstream_index for p in self.promoters}
        downstream_indices = {p.downstream_index for p in self.promoters}
        all_indices = upstream_indices | downstream_indices | {avoid}
        try:
            start = 0
            while True:
                index = self.library.index(motif, start)
                if index not in all_indices:
                    return index
                start = index + 1
        except ValueError as err:
            if start == 0:
                msg = "Promoter elements must be present in the library."
            else:
                msg = (
                    "If a promoter element is reused, "
                    "it must appear several times in the library."
                )
            raise ValueError(msg) from err

    def add_side_biases(
        self: Self, *, left: list[str] | None = None, right: list[str] | None = None
    ) -> None:
        """
        Add side biases for motifs.

        Everything else being equal, the motifs specified in `left` will
        prefer being as much on the left as possible, and likewise for
        those specified in `right`.

        Parameters
        ----------
        left
            List of motifs that should preferentially appear on the left.
        right
            List of motifs that should preferentially appear on the right.

        Raises
        ------
        ValueError
            If the left or right motifs don't belong to the initial library.
        """
        self._ensure_model_not_built("Side biases")
        if any(
            value is not None and not isinstance(value, list) for value in (left, right)
        ):
            msg = "Side biases must be lists of library motifs"
            raise ValueError(msg)
        try:
            ilefts = tuple(self.library.index(motif) for motif in left) if left else ()
            irights = (
                tuple(self.library.index(motif) for motif in right) if right else ()
            )
        except ValueError as err:
            msg = "All motifs must belong to the initial library."
            raise ValueError(msg) from err
        self._ilefts, self._irights = ilefts, irights

    def add_regulator_constraints(
        self: Self,
        regulator_by_index: list[str] | dict[int, str],
        *,
        required: set[str] | None = None,
        min_count_by_regulator: dict[str, int] | None = None,
        min_required_regulators: int | None = None,
    ) -> None:
        """
        Add regulator coverage constraints to the optimization problem.

        Parameters
        ----------
        regulator_by_index
            Mapping from motif index to regulator label (len == nb_motifs).
        required
            Regulators that must appear at least once (>=1).
        min_count_by_regulator
            Per-regulator minimum counts (>=1).
        min_required_regulators
            Require at least this many unique regulators to appear (k-of-n).

        Raises
        ------
        ValueError
            If constraints are invalid or infeasible given the motif library.
        RuntimeError
            If regulator constraints are already set.
        """
        self._ensure_model_not_built("Regulator constraints")
        if self._regulator_constraints is not None:
            msg = (
                "Regulator constraints already set; create a new Optimizer "
                "to change them."
            )
            raise RuntimeError(msg)

        if (
            required is None
            and not min_count_by_regulator
            and min_required_regulators is None
        ):
            msg = "At least one regulator constraint must be provided."
            raise ValueError(msg)

        if required is not None and not isinstance(required, (set, frozenset)):
            msg = "required regulators must be a set of labels"
            raise ValueError(msg)
        if min_count_by_regulator is not None and not isinstance(
            min_count_by_regulator, dict
        ):
            msg = "min_count_by_regulator must be a dict of labels to counts"
            raise ValueError(msg)
        mapping = _normalize_regulator_mapping(self.nb_motifs, regulator_by_index)
        available = set(mapping.values())
        required_set = set(required or [])
        if not required_set.issubset(available):
            missing = sorted(required_set - available)
            msg = f"Required regulators missing from mapping: {missing}"
            raise ValueError(msg)

        min_counts = _normalize_min_counts(
            mapping, required_set, min_count_by_regulator
        )
        min_required_regulators = _normalize_min_required(
            min_required_regulators,
            available,
        )

        self._regulator_constraints = RegulatorRequirements(
            tuple(mapping.items()), tuple(min_counts.items()), min_required_regulators
        )

    @property
    def nb_motifs(self: Self) -> int:
        """The number of motifs in the library (not counting reverse duplicates)."""
        return len(self.library)

    @property
    def nb_nodes(self: Self) -> int:
        """
        The number of nodes in the library (ignoring the starting and end nodes).

        It is equal to `nb_motifs` for single-stranded optimization
        and `2 * nb_motifs` for double-stranded optimization.
        """
        return self.nb_motifs * {"single": 1, "double": 2}[self.strands]

    def add_count_constraint(
        self,
        indices: list[int],
        *,
        minimum: int | None = None,
        maximum: int | None = None,
    ) -> None:
        """Bound selected supplied occurrences independently of string equality.

        Bounds are inclusive, nonnegative integers. At least one is required.
        Overlapping index sets impose separate requirements. Configure before
        building the model; greedy solving cannot enforce these constraints.
        """
        self._ensure_model_not_built("Count constraints")
        selected = occurrence_indices(indices, self.nb_motifs)
        low, high = count_bounds(minimum, maximum, len(selected))
        self._count_constraints += (CountConstraint(selected, low, high),)

    def add_group_coverage(self, groups: list[list[int]], *, minimum: int) -> None:
        """Require a minimum number of represented, disjoint occurrence groups."""
        self._ensure_model_not_built("Group coverage constraints")
        selected = tuple(occurrence_indices(group, self.nb_motifs) for group in groups)
        count = discrete_integer(minimum, "minimum groups", minimum=1)
        flattened = [index for group in selected for index in group]
        if count > len(selected) or len(flattened) != len(set(flattened)):
            msg = "coverage requires distinct groups and enough available groups"
            raise ValueError(msg)
        self._coverage_constraints += (CoverageConstraint(selected, count),)

    def add_fixed_occurrence(
        self,
        index: int,
        *,
        orientation: str,
        start: int | tuple[int | None, int | None] | None = None,
        origin: str = "start",
    ) -> None:
        """Require an exact occurrence/strand without matching by sequence string.

        Start windows are inclusive. ``origin="start"`` uses zero-based packing
        coordinates; ``origin="end"`` uses signed positions relative to the
        packed sequence's exclusive end. The generation layer translates final
        assembly coordinates explicitly.
        """
        self._ensure_model_not_built("Fixed occurrences")
        selected = occurrence_indices([index], self.nb_motifs)[0]
        if orientation not in {"forward", "reverse"}:
            msg = "fixed orientation must be forward or reverse"
            raise ValueError(msg)
        if orientation == "reverse" and self.strands == "single":
            msg = "reverse fixed occurrence requires double-strand eligibility"
            raise ValueError(msg)
        if any(item.index == selected for item in self._fixed_occurrences):
            msg = "a supplied occurrence can be fixed only once"
            raise ValueError(msg)
        if origin not in {"start", "end"}:
            msg = "fixed coordinate origin must be start or end"
            raise ValueError(msg)
        interval = _interval(
            start, "fixed start", minimum=0 if origin == "start" else None
        )
        if (
            origin == "start"
            and interval[0] is not None
            and interval[0] + len(self.library[selected]) > self.sequence_length
        ):
            msg = "fixed start exceeds the sequence length bound"
            raise ValueError(msg)
        self._fixed_occurrences += (
            FixedOccurrence(selected, orientation, interval, origin),
        )

    def add_spacing_constraint(
        self,
        upstream: int,
        downstream: int,
        *,
        minimum: int,
        maximum: int,
    ) -> None:
        """Constrain one declared fixed pair using signed end-to-start spacing."""
        self._ensure_model_not_built("Spacing constraints")
        pair = occurrence_indices([upstream, downstream], self.nb_motifs)
        if not set(pair) <= {item.index for item in self._fixed_occurrences}:
            msg = "spacing requires two declared fixed occurrences"
            raise ValueError(msg)
        if self._spacing_constraints:
            msg = "only one fixed spacing pair is currently supported"
            raise ValueError(msg)
        low = discrete_integer(minimum, "spacing minimum")
        high = discrete_integer(maximum, "spacing maximum")
        interval = _interval((low, high), "spacing")
        self._spacing_constraints = (SpacingConstraint(pair[0], pair[1], interval),)

    def build_model(
        self: Self,
        solver: str = "CBC",
        solver_options: list[str] | None = None,
        *,
        controls: SolverControls | None = None,
    ) -> None:
        """
        Create the solver instance and build the linear model.

        This method belongs to the advanced API: most users should not build the model
        themselves, but rather use functions which build it automatically, such as
        `optimal`, `solutions` or `solutions_diverse`.

        Raises
        ------
        SolverBackendError
            If the requested backend could not create the model.
        ValueError
            If a solver option is invalid or rejected by the backend. The
            existing model is retained when a replacement build is rejected.
        """
        configuration = ModelConfiguration(
            self._problem,
            self._promoters,
            self._regulator_constraints,
            self._ilefts,
            self._irights,
            self._count_constraints,
            self._coverage_constraints,
            self._fixed_occurrences,
            self._spacing_constraints,
        )
        self.model = build_solver_model(configuration, solver, solver_options, controls)
        self._solver_name = solver
        self._model_modified = False

    @property
    def solver_identity(self) -> SolverIdentity | None:
        """Describe the built backend without constructing or solving a model.

        Returns
        -------
        SolverIdentity or None
            Requested name and the backend's reported version; None before build.
        """
        if self.model is None:
            return None
        return SolverIdentity(self._solver_name, self.model.SolverVersion())

    def solve(self: Self) -> DenseArray:
        """
        Solve the currently built model and return its optimal solution.

        This belongs to the advanced API: most users should rather use `optimal()`
        instead, which builds the model automatically.

        Raises
        ------
        RuntimeError
            If the model has not been built yet (`build_model` should be called).
        InfeasibleError
            If the solver proves that no feasible path exists.
        UnprovenSolutionError
            If a feasible incumbent exists but optimality is not proved.
        SolverBackendError
            If the backend fails or reports an unsupported result status.
        InvalidSolverResultError
            If the selected arcs or reconstructed array violate their contracts.

        Returns
        -------
        solution : DenseArray
            The optimal solution.
        """
        report = self.solve_report()
        if report.solution is not None:
            return report.solution
        error_type = {
            SolveStatus.INFEASIBLE: InfeasibleError,
            SolveStatus.UNPROVEN: UnprovenSolutionError,
            SolveStatus.INVALID_RESULT: InvalidSolverResultError,
        }.get(report.status, SolverBackendError)
        raise error_type(report.detail)

    def solve_report(
        self: Self, *, time_limit_seconds: float | None = None
    ) -> SolveReport:
        """Solve once and return typed evidence, preserving unknown causes.

        Returns
        -------
        SolveReport
            Backend status and, only after validation, an optimal solution.

        Raises
        ------
        RuntimeError
            If no model has been built.
        """
        if self.model is None:
            msg = "Model not built: call `build_model(solver)` first"
            raise RuntimeError(msg)

        if time_limit_seconds is not None:
            controls = SolverControls(time_limit_seconds=time_limit_seconds)
            self.model.SetTimeLimit(controls.time_limit_ms)

        # Solve the problem
        try:
            status = self.model.Solve()
        except Exception as err:  # noqa: BLE001 - preserve backend failures as evidence
            msg = f"Solver backend execution failed: {err}"
            return SolveReport(
                SolveStatus.BACKEND_ERROR,
                None,
                termination_reason="backend_exception",
                detail=msg,
            )

        if status != pywraplp.Solver.OPTIMAL:
            status_messages = {
                pywraplp.Solver.FEASIBLE: "A feasible solution was found, but not necessarily optimal.",  # noqa: E501
                pywraplp.Solver.INFEASIBLE: "No feasible solution was found.",
                pywraplp.Solver.UNBOUNDED: "The model is unbounded.",
                pywraplp.Solver.ABNORMAL: "The model is abnormal.",
                pywraplp.Solver.NOT_SOLVED: "The model has not been solved.",
            }
            msg = status_messages.get(
                status, f"Solver ended with unknown status: {status}."
            )
            outcome = {
                pywraplp.Solver.INFEASIBLE: SolveStatus.INFEASIBLE,
                pywraplp.Solver.FEASIBLE: SolveStatus.UNPROVEN,
                pywraplp.Solver.ABNORMAL: SolveStatus.BACKEND_ERROR,
                pywraplp.Solver.UNBOUNDED: SolveStatus.BACKEND_ERROR,
            }.get(status, SolveStatus.UNKNOWN)
            return SolveReport(
                outcome,
                status,
                proof_scope=(
                    "offered_packing_model"
                    if outcome is SolveStatus.INFEASIBLE
                    else None
                ),
                termination_reason=(
                    "proven_infeasible"
                    if outcome is SolveStatus.INFEASIBLE
                    else "unknown"
                ),
                detail=msg,
            )

        try:
            path = self._selected_path()
            offsets_fwd = [None] * self.nb_motifs
            offsets_rev = [None] * self.nb_motifs
            offset = 0
            previous = None
            for node in path:
                if previous is not None:
                    offset += self._adjacency[previous][node]
                offsets = offsets_fwd if node < self.nb_motifs else offsets_rev
                offsets[node % self.nb_motifs] = offset
                previous = node
            solution = DenseArray(
                self.library, self.sequence_length, offsets_fwd, offsets_rev
            )
        except (
            ValueError,
            TypeError,
            AttributeError,
            KeyError,
            IndexError,
            InvalidSolverResultError,
        ) as err:
            msg = f"Invalid solver result: {err}"
            return SolveReport(SolveStatus.INVALID_RESULT, status, detail=msg)
        return SolveReport(
            SolveStatus.OPTIMAL,
            status,
            solution,
            proof_scope="offered_packing_model",
            termination_reason="proven_optimal",
        )

    def _selected_path(self) -> list[int]:
        selected = set()
        for edge, variable in self.model.X.items():
            value = variable.solution_value()
            if (
                not math.isfinite(value)
                or min(abs(value), abs(value - 1)) > _INTEGER_TOLERANCE
            ):
                msg = "Solver returned a non-binary path variable"
                raise InvalidSolverResultError(msg)
            if round(value) == 1:
                selected.add(edge)
        remaining = set(selected)
        path = []
        used_entries = set()
        current = -1
        while True:
            successors = [j for i, j in remaining if i == current]
            if len(successors) != 1:
                msg = "Solver returned an invalid path: expected one successor"
                raise InvalidSolverResultError(msg)
            node = successors[0]
            remaining.remove((current, node))
            if node == -1:
                break
            if not 0 <= node < self.nb_nodes or node % self.nb_motifs in used_entries:
                msg = "Solver returned an invalid path: repeated or unknown entry"
                raise InvalidSolverResultError(msg)
            used_entries.add(node % self.nb_motifs)
            path.append(node)
            current = node
        if remaining or not path:
            msg = "Solver returned an invalid path: disconnected or empty result"
            raise InvalidSolverResultError(msg)
        return path

    def forbid(self: Self, solution: DenseArray) -> None:
        """
        Add a constraint to the model to forbid a given solution.

        Parameters
        ----------
        solution
            A result with this library, length bound, and permitted orientations,
            whose ordered offsets describe an exact packing path.

        Raises
        ------
        RuntimeError
            If the model has not been built yet (`build_model` should be called).
        ValueError
            If the result is foreign or its placements do not form a valid path.
            Rejected results leave the model unchanged.
        """
        if self.model is None:
            msg = "Model not built: call `build_model(solver)` first"
            raise RuntimeError(msg)

        if not isinstance(solution, DenseArray) or (
            solution.library != self.library
            or solution.sequence_length != self.sequence_length
        ):
            msg = "Solution does not belong to this packing problem"
            raise ValueError(msg)
        if (
            self.length_mode == "exact"
            and len(solution.sequence) != self.sequence_length
        ):
            msg = "Solution violates this problem's exact packing length"
            raise ValueError(msg)
        if self.strands == "single" and any(
            o is not None for o in solution.offsets_rev
        ):
            msg = "Solution violates this problem's single-strand policy"
            raise ValueError(msg)
        ordered = solution.offset_indices_in_order()
        expected_offset = 0
        previous = None
        for offset, index in ordered:
            if previous is not None:
                expected_offset += self._adjacency[previous][index]
            if offset != expected_offset:
                msg = "Solution placements do not describe an exact packing path"
                raise ValueError(msg)
            previous = index
        sol = [-1, *(i for _, i in ordered), -1]
        sum_on_path = sum(self.model.X[i, j] for i, j in it.pairwise(sol))
        self.model.Add(sum_on_path <= solution.nb_motifs)
        self._model_modified = True

    def set_motif_weight(self: Self, imotif: int, weight: float) -> None:
        """
        Set the weight of a particular motif in the score.

        Parameters
        ----------
        imotif
            Index of the motif.
        weight
            Finite real weight, representable as a solver coefficient. Negative
            weights are permitted; booleans and nonfinite values are invalid.

        Raises
        ------
        RuntimeError
            If the model has not been built yet (`build_model` should be called).
        ValueError
            If the index is not an original library entry or the weight is
            invalid. Rejected updates leave all coefficients unchanged.
        """
        if self.model is None:
            msg = "Model not built: call `build_model(solver)` first"
            raise RuntimeError(msg)

        imotif = discrete_integer(imotif, "imotif", minimum=0)
        if imotif >= self.nb_motifs:
            msg = "imotif must identify an entry in the motif library"
            raise ValueError(msg)
        if isinstance(weight, bool) or not isinstance(weight, Real):
            msg = "weight must be a finite real number (not bool)"
            raise ValueError(msg)  # noqa: TRY004 - invalid configuration is a ValueError
        try:
            normalized_weight = float(weight)
        except (OverflowError, ValueError) as err:
            msg = "weight must be representable as a finite real number"
            raise ValueError(msg) from err
        if not math.isfinite(normalized_weight):
            msg = "weight must be a finite real number"
            raise ValueError(msg)
        objective = self.model.Objective()

        for i in range(-1, self.nb_nodes):
            if i != imotif:
                objective.SetCoefficient(self.model.X[i, imotif], normalized_weight)
            imotif2 = imotif + self.nb_motifs
            if self.strands == "double" and i != imotif2:
                objective.SetCoefficient(self.model.X[i, imotif2], normalized_weight)
        self._model_modified = True

    def solutions(
        self: Self,
        solver: str = "CBC",
        solver_options: list[str] | None = None,
        *,
        controls: SolverControls | None = None,
    ) -> Iterator[DenseArray]:
        """
        Iterate over solutions in decreasing order of score.

        Only proven infeasibility ends iteration normally. Unproven, backend,
        and invalid-result failures propagate, including after a yielded result.

        Note that this function (re)builds the model automatically.

        Parameters
        ----------
        solver
            Solver name given to OrTools.
        solver_options
            List of strings passed to the solver
            with `SetSolverSpecificParametersAsString`.
        controls
            Cooperative limits applied to each backend solve. They do not bound
            total enumeration time. Unsupported controls raise `ValueError`.

        Yields
        ------
        solution :
            Solutions in decreasing order of score.
        """
        self.build_model(solver, solver_options=solver_options, controls=controls)

        while True:
            try:
                sol = self.solve()
            except InfeasibleError:
                break
            yield sol
            self.forbid(sol)

    def solutions_diverse(
        self: Self,
        solver: str = "CBC",
        solver_options: list[str] | None = None,
        *,
        controls: SolverControls | None = None,
    ) -> Iterator[DenseArray]:
        """
        Yield optimal packings under weights that favor underused motif entries.

        Entry usage updates the objective after each yield. Each result is an
        optimum for that iteration's weighted model; scores across iterations
        need not decrease, and equal cumulative usage is not guaranteed.

        Only proven infeasibility ends iteration normally. Other solve failures
        propagate, including after a yielded result.

        Note that this function (re)builds the model automatically.

        Parameters
        ----------
        solver
            Solver name given to OrTools.
        solver_options
            List of strings passed to the solver
            with `SetSolverSpecificParametersAsString`.
        controls
            Cooperative limits applied to each backend solve. They do not bound
            total enumeration time. Unsupported controls raise `ValueError`.

        Yields
        ------
        solution :
            A packing proven optimal under the current motif weights.
        """
        self.build_model(solver, solver_options=solver_options, controls=controls)

        motifs = [0] * self.nb_motifs
        weights = (1.0,) * self.nb_motifs
        while True:
            try:
                sol = self.solve()
            except InfeasibleError:
                break
            yield sol
            # Forbid the solution
            self.forbid(sol)
            # Tally up the motifs
            for i, (fwd, rev) in enumerate(
                zip(sol.offsets_fwd, sol.offsets_rev, strict=True),
            ):
                if fwd is not None or rev is not None:
                    motifs[i] += 1
            # Update motif weights
            updated = part_usage_weights(motifs)
            for index, (before, after) in enumerate(zip(weights, updated, strict=True)):
                if before != after:
                    self.set_motif_weight(index, after)
            weights = updated

    def optimal(
        self: Self,
        solver: str = "CBC",
        solver_options: list[str] | None = None,
        *,
        controls: SolverControls | None = None,
    ) -> DenseArray:
        """
        Return the optimal solution.

        Note that this function (re)builds the model automatically.

        Parameters
        ----------
        solver
            Solver name given to OrTools.
        solver_options
            List of strings passed to the solver
            with `SetSolverSpecificParametersAsString`.
        controls
            Cooperative limits applied to each backend solve. They do not bound
            total enumeration time. Unsupported controls raise `ValueError`.

        Returns
        -------
        solution :
            Optimal solution.

        Raises
        ------
        InfeasibleError
            If no feasible solution exists.
        """
        try:
            return next(
                self.solutions(solver, solver_options=solver_options, controls=controls)
            )
        except StopIteration as err:
            msg = "No feasible solution was found."
            raise InfeasibleError(msg) from err

    def approximate(self: Self) -> DenseArray:
        """
        Return the best path found by deterministic greedy starts.

        Every selected entry owns one occurrence on one strand. Duplicate entries
        need separate occurrences; incidental contained substrings are not added.
        The heuristic does not prove optimality and does not call a solver.

        Raises
        ------
        ValueError
            If exact length, promoter/regulator requirements, count/coverage
            bounds, fixed occurrences, spacing, side biases, or model changes
            from ``forbid`` or ``set_motif_weight`` are configured.
        InfeasibleError
            If no motif fits the sequence length bound.

        Returns
        -------
        solution :
            Approximate solution.
        """
        if (
            self._promoters
            or self._regulator_constraints
            or self._count_constraints
            or self._coverage_constraints
            or self._fixed_occurrences
            or self._spacing_constraints
            or self.length_mode == "exact"
            or self._ilefts
            or self._irights
            or self._model_modified
        ):
            msg = (
                "approximate() does not support configured constraints, "
                "side biases, or model changes"
            )
            raise ValueError(msg)
        return realize_greedy(self._problem)
