"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/model.py

Construct the integer path model from a validated problem snapshot.

Module Author(s): Virgile Andreani, Eric J. South
Maintainer(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from ortools.linear_solver import pywraplp

from .errors import SolverBackendError
from .problem import discrete_integer
from .solver import SolverControls

if TYPE_CHECKING:
    from collections.abc import Sequence

    from .constraints import (
        CountConstraint,
        CoverageConstraint,
        FixedOccurrence,
        PromoterConstraint,
        RegulatorRequirements,
        SpacingConstraint,
    )
    from .problem import PackingProblem


def part_usage_weights(counts: Sequence[int]) -> tuple[float, ...]:
    """Favor below-mean proposed usage without trading away one occurrence.

    At most N selected parts each receive 0.5/N bonus, so the total preference
    contribution is at most 0.5. Integer arithmetic identifies below-mean usage.
    """
    if not counts:
        msg = "part usage requires a nonempty collection"
        raise ValueError(msg)
    values = tuple(discrete_integer(c, "part usage", minimum=0) for c in counts)
    total, size = sum(values), len(values)
    return tuple(1 + 0.5 / size if c * size < total else 1.0 for c in values)


@dataclass(frozen=True)
class ModelConfiguration:
    """Immutable problem and configured constraints for one model build."""

    problem: PackingProblem
    promoters: tuple[PromoterConstraint, ...]
    regulator_constraints: RegulatorRequirements | None
    left_bias: tuple[int, ...]
    right_bias: tuple[int, ...]
    counts: tuple[CountConstraint, ...] = ()
    coverage: tuple[CoverageConstraint, ...] = ()
    fixed: tuple[FixedOccurrence, ...] = ()
    spacing: tuple[SpacingConstraint, ...] = ()


def build_solver_model(
    configuration: ModelConfiguration,
    solver: str,
    solver_options: list[str] | None,
    controls: SolverControls | None = None,
) -> pywraplp.Solver:
    """Build a fresh model, exposing it only after all options are accepted.

    Returns
    -------
    pywraplp.Solver
        The initialized backend model with path variables and constraints.

    Raises
    ------
    SolverBackendError
        If the requested backend cannot be created.
    ValueError
        If a solver option is invalid or rejected by the backend.
    """
    if solver_options is not None and not isinstance(solver_options, list):
        msg = "solver_options must be a list of non-empty strings"
        raise ValueError(msg)
    if solver_options is not None and any(
        not isinstance(o, str) or not o.strip() for o in solver_options
    ):
        msg = "solver_options must be a list of non-empty strings"
        raise ValueError(msg)
    controls = _validated_controls(controls, solver)
    try:
        model = pywraplp.Solver.CreateSolver(solver)
    except Exception as err:
        msg = f"Could not create model for solver {solver!r}: {err}"
        raise SolverBackendError(msg) from err
    if model is None:
        msg = f"Could not create model for solver {solver!r}; backend unavailable."
        raise SolverBackendError(msg)
    if controls.time_limit_ms is not None:
        model.SetTimeLimit(controls.time_limit_ms)
    if controls.threads is not None and not model.SetNumThreads(controls.threads):
        msg = f"Solver {solver!r} rejected threads={controls.threads}"
        raise ValueError(msg)
    if solver_options:
        for option in solver_options:
            if not model.SetSolverSpecificParametersAsString(option):
                msg = f"Solver {solver!r} rejected option: {option!r}"
                raise ValueError(msg)
    _ModelBuilder(configuration, model).build()
    return model


def _validated_controls(controls: SolverControls | None, solver: str) -> SolverControls:
    if controls is None:
        return SolverControls()
    if not isinstance(controls, SolverControls):
        msg = "controls must be SolverControls"
        raise TypeError(msg)
    if controls.threads is not None and solver.upper() not in {
        "SCIP",
        "SCIP_MIXED_INTEGER_PROGRAMMING",
    }:
        msg = f"Explicit threads are unsupported for {solver}; use SCIP or omit them"
        raise ValueError(msg)
    return controls


class _ModelBuilder:
    def __init__(
        self, configuration: ModelConfiguration, model: pywraplp.Solver
    ) -> None:
        self.model = model
        self.library = configuration.problem.library
        self.sequence_length = configuration.problem.sequence_length
        self.length_mode = configuration.problem.length_mode
        self.strands = configuration.problem.strands
        self.nb_motifs = len(self.library)
        self.nb_nodes = len(configuration.problem.oriented_library)
        self._adjacency = configuration.problem.adjacency
        self.promoters = configuration.promoters
        self._regulator_constraints = configuration.regulator_constraints
        self.counts = configuration.counts
        self.coverage = configuration.coverage
        self.fixed = configuration.fixed
        self.spacing = configuration.spacing
        self.ilefts = configuration.left_bias
        self.irights = configuration.right_bias

    def build(self) -> None:
        """Add the objective and constraints for a single path."""
        # X_ij are binary variables. X_ij == 1 means that motif #j directly follows
        # (and possibly overlaps) motif #i in the sequence.
        start = {
            (-1, j): self.model.BoolVar(f"X[-1,{j}]") for j in range(self.nb_nodes)
        }
        end = {(i, -1): self.model.BoolVar(f"X[{i},-1]") for i in range(self.nb_nodes)}
        middle = {
            (i, j): self.model.BoolVar(f"X[{i},{j}]")
            for i in range(self.nb_nodes)
            for j in range(self.nb_nodes)
            if i != j
        }
        X = start | end | middle  # noqa: N806
        self.model.X = X

        # Path starts at the start
        self.model.Add(sum(X[-1, j] for j in range(self.nb_nodes)) == 1)

        # Path ends at the end
        self.model.Add(sum(X[i, -1] for i in range(self.nb_nodes)) == 1)

        # Conservation of flow
        for k in range(self.nb_nodes):
            enter_direct = sum(X[i, k] for i in range(-1, self.nb_nodes) if i != k)
            exit_direct = sum(X[k, j] for j in range(-1, self.nb_nodes) if k != j)
            self.model.Add(enter_direct == exit_direct)

        # Don't include any motif more than once
        for k in range(self.nb_motifs):
            enter_direct = sum(X[i, k] for i in range(-1, self.nb_nodes) if i != k)
            exit_direct = sum(X[k, j] for j in range(-1, self.nb_nodes) if k != j)
            if self.strands == "single":
                self.model.Add(enter_direct <= 1)
                self.model.Add(exit_direct <= 1)
                continue
            # krev is the index of the reverse complement of motif k
            krev = k + self.nb_motifs
            enter_rev = sum(X[i, krev] for i in range(-1, self.nb_nodes) if i != krev)
            exit_rev = sum(X[krev, j] for j in range(-1, self.nb_nodes) if krev != j)
            self.model.Add(enter_direct + enter_rev <= 1)
            self.model.Add(exit_direct + exit_rev <= 1)

        # Global length constraint
        size_inside = sum(
            self._adjacency[i][j] * X[i, j]
            for i in range(self.nb_nodes)
            for j in range(self.nb_nodes)
            if i != j
        )
        size_terminal = sum(
            len(self.library[i % self.nb_motifs]) * X[i, -1]
            for i in range(self.nb_nodes)
        )
        self.packed_span = size_inside + size_terminal
        if self.length_mode == "exact":
            self.model.Add(self.packed_span == self.sequence_length)
        else:
            self.model.Add(self.packed_span <= self.sequence_length)

        # Subtour elimination variables
        self._add_continuity_variables()

        # Apply user-defined distance constraints
        self._add_promoter_constraints()

        # Apply regulator coverage constraints
        self._add_regulator_constraints()
        self._add_occurrence_constraints()
        self._add_fixed_geometry()

        # Objective
        self.model.Maximize(
            sum(
                X[i, j]
                for i in range(-1, self.nb_nodes)
                for j in range(self.nb_nodes)
                if i != j
            ),
        )

        # Apply user-defined side biases
        # (needs to be after the objective definition because it modifies it)
        self._add_side_biases()

    def _occurrence_total(self, indices: tuple[int, ...]) -> pywraplp.LinearExpr:
        nodes = list(indices)
        if self.strands == "double":
            nodes.extend(i + self.nb_motifs for i in indices)
        return sum(
            self.model.X[i, node]
            for node in nodes
            for i in range(-1, self.nb_nodes)
            if i != node
        )

    def _add_occurrence_constraints(self) -> None:
        for constraint in self.counts:
            total = self._occurrence_total(constraint.indices)
            if constraint.minimum is not None:
                self.model.Add(total >= constraint.minimum)
            if constraint.maximum is not None:
                self.model.Add(total <= constraint.maximum)
        for number, constraint in enumerate(self.coverage):
            represented = []
            for ordinal, indices in enumerate(constraint.groups):
                used = self.model.BoolVar(f"coverage[{number},{ordinal}]")
                total = self._occurrence_total(indices)
                self.model.Add(used <= total)
                self.model.Add(total <= len(indices) * used)
                represented.append(used)
            self.model.Add(sum(represented) >= constraint.minimum)

    def _add_fixed_geometry(self) -> None:
        if not self.fixed:
            return
        self._add_position_variables()
        nodes = {}
        for fixed in self.fixed:
            node = fixed.index + (
                self.nb_motifs if fixed.orientation == "reverse" else 0
            )
            nodes[fixed.index] = node
            incoming = sum(
                self.model.X[i, node] for i in range(-1, self.nb_nodes) if i != node
            )
            self.model.Add(incoming == 1)
            minimum, maximum = fixed.start
            position = self.model.position[node]
            if fixed.origin == "end":
                position = position - self.packed_span
            if minimum is not None:
                self.model.Add(position >= minimum)
            if maximum is not None:
                self.model.Add(position <= maximum)
        for spacing in self.spacing:
            distance = (
                self.model.position[nodes[spacing.downstream]]
                - self.model.position[nodes[spacing.upstream]]
                - len(self.library[spacing.upstream])
            )
            self.model.Add(distance >= spacing.interval[0])
            self.model.Add(distance <= spacing.interval[1])

    def _add_continuity_variables(self) -> None:
        """Implement subtour elimination variables and constraints into the model."""
        try:
            self.model.cont  # noqa: B018
        except AttributeError:
            pass
        else:
            # Continuity variables already exist
            return

        self.model.cont = [
            self.model.IntVar(1, self.nb_nodes, f"u[{i}]") for i in range(self.nb_nodes)
        ]

        for i in range(self.nb_nodes):
            for j in range(self.nb_nodes):
                if i == j:
                    continue
                # A selected edge increases rank: u_j >= u_i + 1. Direct
                # coefficients avoid temporary expression trees for every edge.
                row = self.model.Constraint(-self.model.infinity(), self.nb_nodes - 1)
                row.SetCoefficient(self.model.cont[i], 1)
                row.SetCoefficient(self.model.cont[j], -1)
                row.SetCoefficient(self.model.X[i, j], self.nb_nodes)

    def _add_position_variables(self) -> None:
        """Implement position variables and constraints into the model."""
        try:
            self.model.position  # noqa: B018
        except AttributeError:
            pass
        else:
            # Position variables already exist
            return

        # Initialize position variables
        self.model.position = [
            self.model.IntVar(0, self.sequence_length - 1, f"position[{i}]")
            for i in range(self.nb_nodes)
        ]

        # Define position for each node
        for i in range(-1, self.nb_nodes):
            for j in range(self.nb_nodes):
                if i == j:
                    continue
                shift = 0 if i == -1 else self._adjacency[i][j]
                base_pos = 0 if i == -1 else self.model.position[i]
                distance_i_j = self.model.position[j] - base_pos
                slack = (self.sequence_length - 1) * (1 - self.model.X[i, j])
                self.model.Add(shift * self.model.X[i, j] - slack <= distance_i_j)
                self.model.Add(distance_i_j <= shift * self.model.X[i, j] + slack)

    def _add_promoter_constraints(self) -> None:
        """Implement promoter constraints into the model."""
        if not self.promoters:
            return

        self._add_position_variables()

        for constraint in self.promoters:
            # Both upstream and downstream elements must appear in the sequence
            for k in [constraint.upstream_index, constraint.downstream_index]:
                self.model.Add(
                    sum(self.model.X[i, k] for i in range(-1, self.nb_nodes) if i != k)
                    >= 1
                )

            # Position both upstream and downstream elements
            spacer_length = (
                self.model.position[constraint.downstream_index]
                - self.model.position[constraint.upstream_index]
                - len(self.library[constraint.upstream_index])
            )
            for pos_or_len, (min_val, max_val) in [
                (
                    self.model.position[constraint.upstream_index],
                    constraint.upstream_pos,
                ),
                (
                    self.model.position[constraint.downstream_index],
                    constraint.downstream_pos,
                ),
                (spacer_length, constraint.spacer_length),
            ]:
                if min_val is not None:
                    self.model.Add(min_val <= pos_or_len)
                if max_val is not None:
                    self.model.Add(pos_or_len <= max_val)

    def _add_side_biases(self) -> None:
        """Implement the side biases into the model."""
        if not self.ilefts and not self.irights:
            return

        self._add_position_variables()

        objective = self.model.Objective()

        weight = 0.5 / (self.nb_motifs * self.sequence_length)

        for i in self.ilefts:
            objective.SetCoefficient(self.model.position[i], -weight)
            if self.strands == "double":
                irev = i + self.nb_motifs
                objective.SetCoefficient(self.model.position[irev], -weight)
        for i in self.irights:
            objective.SetCoefficient(self.model.position[i], weight)
            if self.strands == "double":
                irev = i + self.nb_motifs
                objective.SetCoefficient(self.model.position[irev], weight)

    def _add_regulator_constraints(self) -> None:
        """Implement regulator coverage constraints into the model."""
        if not self._regulator_constraints:
            return

        mapping = dict(self._regulator_constraints.mapping)
        min_counts = dict(self._regulator_constraints.min_counts)
        min_required = self._regulator_constraints.min_required

        self.model.selected = [
            self.model.BoolVar(f"selected[{i}]") for i in range(self.nb_motifs)
        ]

        def _incoming(node: int) -> pywraplp.LinearExpr:
            return sum(
                self.model.X[i, node] for i in range(-1, self.nb_nodes) if i != node
            )

        for i in range(self.nb_motifs):
            used_fwd = _incoming(i)
            if self.strands == "double":
                used_rev = _incoming(i + self.nb_motifs)
                used_total = used_fwd + used_rev
            else:
                used_total = used_fwd
            self.model.Add(used_total <= self.model.selected[i])
            self.model.Add(self.model.selected[i] <= used_total)

        groups: dict[str, list[int]] = {}
        for idx, label in mapping.items():
            groups.setdefault(label, []).append(idx)

        for label, indices in groups.items():
            total = sum(self.model.selected[i] for i in indices)
            min_count = min_counts.get(label)
            if min_count is not None:
                self.model.Add(total >= min_count)

        if min_required is not None:
            covered_flags = []
            for ridx, indices in enumerate(groups.values()):
                covered = self.model.BoolVar(f"covered[{ridx}]")
                total = sum(self.model.selected[i] for i in indices)
                self.model.Add(total >= covered)
                self.model.Add(total <= len(indices) * covered)
                covered_flags.append(covered)
            self.model.Add(sum(covered_flags) >= min_required)
