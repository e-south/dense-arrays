"""Optimization model and solver for dense-arrays.

--------------------------------------------------------------------------------
<dense-array project>

Optimization model and solver for dense-arrays.

Module Author(s): Virgile Andreani, Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

from __future__ import annotations

import itertools as it
from typing import TYPE_CHECKING, Self

if TYPE_CHECKING:
    from collections.abc import Iterator

from ortools.linear_solver import pywraplp

from .constraints import (
    PromoterConstraint,
    _normalize_min_counts,
    _normalize_min_required,
    _normalize_regulator_mapping,
)
from .sequence import VALID_BASES, adjacency_matrix, reverse_complement, take_best_run
from .solution import DenseArray


class Optimizer:
    """Optimizer."""

    def __init__(
        self: Self,
        library: list[str],
        sequence_length: int,
        strands: str = "double",
    ) -> None:
        if sequence_length <= 0:
            msg = "sequence_length must be > 0"
            raise ValueError(msg)
        if strands not in {"single", "double"}:
            msg = "strands must be single or double"
            raise ValueError(msg)
        if not library:
            msg = "library must contain at least one motif"
            raise ValueError(msg)
        for i, motif in enumerate(library):
            if not isinstance(motif, str) or not motif:
                msg = f"motif at index {i} must be a non-empty string"
                raise ValueError(msg)
            invalid = set(motif) - VALID_BASES
            if invalid:
                msg = (
                    f"motif at index {i} contains invalid bases: {sorted(invalid)}. "
                    "Use uppercase A/C/G/T only."
                )
                raise ValueError(msg)

        self.library = list(library)
        self.sequence_length = sequence_length
        self.strands = strands
        self.promoters: list[PromoterConstraint] = []
        self._regulator_constraints: dict | None = None
        if strands == "double":
            library = library + [reverse_complement(motif) for motif in library]  # noqa: PLR6104
        self.adjacency_matrix = adjacency_matrix(library)
        self.model = None
        self.ilefts: list[int] = []
        self.irights: list[int] = []

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
            Length of the spacer between both elements, or tuple (min, max).
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
        self.promoters.append(constraint)

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
        try:
            self.ilefts = [self.library.index(motif) for motif in left] if left else []
            self.irights = (
                [self.library.index(motif) for motif in right] if right else []
            )
        except ValueError as err:
            msg = "All motifs must belong to the initial library."
            raise ValueError(msg) from err

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

        self._regulator_constraints = {
            "mapping": mapping,
            "required": required_set,
            "min_counts": min_counts,
            "min_required": min_required_regulators,
        }

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

    def build_model(
        self: Self, solver: str = "CBC", solver_options: list[str] | None = None
    ) -> None:
        """
        Create the solver instance and build the linear model.

        This method belongs to the advanced API: most users should not build the model
        themselves, but rather use functions which build it automatically, such as
        `optimal`, `solutions` or `solutions_diverse`.

        Raises
        ------
        RuntimeError
            If the backend could not create the model.
        """
        self.model = pywraplp.Solver.CreateSolver(solver)

        if self.model is None:
            msg = "Could not create model. There is a problem with the backend."
            raise RuntimeError(msg)

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
            self.adjacency_matrix[i][j] * X[i, j]
            for i in range(self.nb_nodes)
            for j in range(self.nb_nodes)
            if i != j
        )
        size_terminal = sum(
            len(self.library[i % self.nb_motifs]) * X[i, -1]
            for i in range(self.nb_nodes)
        )
        self.model.Add(size_inside + size_terminal <= self.sequence_length)

        # Subtour elimination variables
        self._add_continuity_variables()

        # Apply user-defined distance constraints
        self._add_promoter_constraints()

        # Apply regulator coverage constraints
        self._add_regulator_constraints()

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

        if solver_options:
            for option in solver_options:
                self.model.SetSolverSpecificParametersAsString(option)

    def _add_continuity_variables(self: Self) -> None:
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
                distance_i_j = self.model.cont[j] - self.model.cont[i]
                slack = self.nb_nodes * (1 - self.model.X[i, j])
                self.model.Add(-distance_i_j + 1 <= slack)

    def _add_position_variables(self: Self) -> None:
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
                shift = 0 if i == -1 else self.adjacency_matrix[i][j]
                base_pos = 0 if i == -1 else self.model.position[i]
                distance_i_j = self.model.position[j] - base_pos
                slack = (self.sequence_length - 1) * (1 - self.model.X[i, j])
                self.model.Add(shift * self.model.X[i, j] - slack <= distance_i_j)
                self.model.Add(distance_i_j <= shift * self.model.X[i, j] + slack)

    def _add_promoter_constraints(self: Self) -> None:
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

    def _add_side_biases(self: Self) -> None:
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

    def _add_regulator_constraints(self: Self) -> None:
        """Implement regulator coverage constraints into the model."""
        if not self._regulator_constraints:
            return

        mapping: dict[int, str] = self._regulator_constraints["mapping"]
        min_counts: dict[str, int] = self._regulator_constraints["min_counts"]
        min_required = self._regulator_constraints["min_required"]

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

    def solve(self: Self) -> DenseArray:  # noqa: C901, PLR0912
        """
        Solve the currently built model and return its optimal solution.

        This belongs to the advanced API: most users should rather use `optimal()`
        instead, which builds the model automatically.

        Raises
        ------
        RuntimeError
            If the model has not been built yet (`build_model` should be called).
        ValueError
            If the model could not be solved optimally.

        Returns
        -------
        solution : DenseArray
            The optimal solution.
        """
        if self.model is None:
            msg = "Model not built: call `build_model(solver)` first"
            raise RuntimeError(msg)

        # Solve the problem
        status = self.model.Solve()

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
            raise ValueError(msg)

        # Extract the solution
        sol = [-1]
        offset = 0
        offsets_fwd = [None] * self.nb_motifs
        offsets_rev = [None] * self.nb_motifs
        while sol[-1] >= 0 or len(sol) == 1:
            current = sol[-1]
            if current == -1:
                candidates = range(self.nb_nodes)
            else:
                candidates = range(-1, self.nb_nodes)
            for j in candidates:
                if current >= 0 and j == current:
                    continue
                if round(self.model.X[current, j].solution_value()) == 1:
                    if len(sol) > 1:
                        offset += self.adjacency_matrix[current][j]
                    if j >= self.nb_motifs:
                        offsets_rev[j % self.nb_motifs] = offset
                    elif j >= 0:
                        offsets_fwd[j] = offset
                    sol.append(j)
                    break
            else:
                msg = "Solver returned an invalid path."
                raise RuntimeError(msg)
        sol = sol[1:-1]
        if not sol:
            msg = "No feasible solution was found."
            raise ValueError(msg)

        return DenseArray(
            self.library,
            self.sequence_length,
            offsets_fwd,
            offsets_rev,
        )

    def forbid(self: Self, solution: DenseArray) -> None:
        """
        Add a constraint to the model to forbid a given solution.

        Parameters
        ----------
        solution
            The solution to forbid.

        Raises
        ------
        RuntimeError
            If the model has not been built yet (`build_model` should be called).
        """
        if self.model is None:
            msg = "Model not built: call `build_model(solver)` first"
            raise RuntimeError(msg)

        sol = [-1, *(i for _, i in solution.offset_indices_in_order()), -1]
        sum_on_path = sum(self.model.X[i, j] for i, j in it.pairwise(sol))
        self.model.Add(sum_on_path <= solution.nb_motifs)

    def set_motif_weight(self: Self, imotif: int, weight: float) -> None:
        """
        Set the weight of a particular motif in the score.

        Parameters
        ----------
        imotif
            Index of the motif.
        weight
            Weight of the motif.

        Raises
        ------
        RuntimeError
            If the model has not been built yet (`build_model` should be called).
        """
        if self.model is None:
            msg = "Model not built: call `build_model(solver)` first"
            raise RuntimeError(msg)

        objective = self.model.Objective()

        for i in range(-1, self.nb_nodes):
            if i != imotif:
                objective.SetCoefficient(self.model.X[i, imotif], weight)
            imotif2 = imotif + self.nb_motifs
            if self.strands == "double" and i != imotif2:
                objective.SetCoefficient(self.model.X[i, imotif2], weight)

    def solutions(
        self: Self, solver: str = "CBC", solver_options: list[str] | None = None
    ) -> Iterator[DenseArray]:
        """
        Iterate over solutions in decreasing order of score.

        Note that this function (re)builds the model automatically.

        Parameters
        ----------
        solver
            Solver name given to OrTools.
        solver_options
            List of strings passed to the solver
            with `SetSolverSpecificParametersAsString`.

        Yields
        ------
        solution :
            Solutions in decreasing order of score.
        """
        self.build_model(solver, solver_options=solver_options)

        while True:
            try:
                sol = self.solve()
            except ValueError:
                break
            yield sol
            self.forbid(sol)

    def solutions_diverse(
        self: Self, solver: str = "CBC", solver_options: list[str] | None = None
    ) -> Iterator[DenseArray]:
        """
        Return an iterator of optimal solutions trying to minimize the bias in motifs.

        Note that this function (re)builds the model automatically.

        Parameters
        ----------
        solver
            Solver name given to OrTools.
        solver_options
            List of strings passed to the solver
            with `SetSolverSpecificParametersAsString`.

        Yields
        ------
        solution :
            Solutions in decreasing order of score.
        """
        self.build_model(solver, solver_options=solver_options)

        epsilon = 0.5 / self.nb_motifs

        motifs = [0] * self.nb_motifs
        imins: list[int] = []
        while True:
            try:
                sol = self.solve()
            except ValueError:
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
            for imin in imins:
                self.set_motif_weight(imin, 1)
            avg_abundance = sum(motifs) / len(motifs)
            imins = [i for i, qty in enumerate(motifs) if qty < avg_abundance]
            for imin in imins:
                self.set_motif_weight(imin, 1 + epsilon)

    def optimal(
        self: Self, solver: str = "CBC", solver_options: list[str] | None = None
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

        Returns
        -------
        solution :
            Optimal solution.

        Raises
        ------
        ValueError
            If no feasible solution exists.
        """
        try:
            return next(self.solutions(solver, solver_options=solver_options))
        except StopIteration as err:
            msg = "No feasible solution was found."
            raise ValueError(msg) from err

    def approximate(self: Self) -> DenseArray:
        """
        Return a solution approximated with a greedy algorithm.

        Raises
        ------
        ValueError
            If no feasible solution exists.

        Returns
        -------
        solution :
            Approximate solution.
        """
        library = list(self.library)
        if self.strands == "double":
            library += [reverse_complement(motif) for motif in library]

        while len(library) > {"single": 1, "double": 2}[self.strands]:
            adj = adjacency_matrix(library)
            _min_dist, i, j = min(
                (adj[i][j], i, j)
                for i in range(len(library))
                for j in range(len(library))
                if i != j
                and not (self.strands == "double" and abs(i - j) == len(library) // 2)
            )
            library[i] = library[i][: adj[i][j]] + library[j]
            if self.strands == "double":
                library[(i + len(library) // 2) % len(library)] = reverse_complement(
                    library[i],
                )
                del library[max(j, (j + len(library) // 2) % len(library))]
                del library[min(j, (j + len(library) // 2) % len(library))]
            else:
                del library[j]
        sequence = take_best_run(
            library[0],
            self.sequence_length,
            self.library,
            self.strands,
        )
        offsets_fwd = [
            sequence.index(motif) if motif in sequence else None
            for motif in self.library
        ]
        if self.strands == "double":
            offsets_rev = [
                sequence.index(reverse_complement(motif))
                if reverse_complement(motif) in sequence
                else None
                for motif in self.library
            ]
            for i in range(self.nb_motifs):
                if offsets_fwd[i] is not None and offsets_rev[i] is not None:
                    offsets_rev[i] = None
        else:
            offsets_rev = [None] * self.nb_motifs
        if all(offset is None for offset in offsets_fwd) and all(
            offset is None for offset in offsets_rev
        ):
            msg = "No feasible solution was found."
            raise ValueError(msg)
        return DenseArray(self.library, self.sequence_length, offsets_fwd, offsets_rev)
