"""Optimization module."""

import itertools as it
from collections.abc import Iterator
from dataclasses import dataclass
from typing import Self

from ortools.linear_solver import pywraplp

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "-": "-"}

__all__ = ["shift_metric", "DenseArray", "Optimizer"]


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
    """
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
    """
    return [[shift_metric(motifa, motifb) for motifb in motifs] for motifa in motifs]


def reverse_complement(sequence: str) -> str:
    return "".join(COMPLEMENT[c] for c in sequence[::-1])


def dispatch_labels(
    library: list[str],
    offsets: list[int | None],
    *,
    rev: bool,
) -> list[str]:
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


@dataclass
class DenseArray:
    """Representation of a solution."""

    library: list[str]
    sequence_length: int
    sequence: str
    offsets_fwd: list[int | None]
    offsets_rev: list[int | None]

    def __init__(
        self: Self,
        library: list[str],
        sequence_length: int,
        offsets_fwd: list[int | None],
        offsets_rev: list[int | None],
    ) -> None:
        self.library = library
        self.sequence_length = sequence_length
        self.offsets_fwd = offsets_fwd
        self.offsets_rev = offsets_rev
        sequence = ""
        for offset, i in self.offset_indices_in_order():
            motif = library[i % len(library)]
            if i >= len(library):
                motif = reverse_complement(motif)
            sequence = sequence[:offset] + motif
        self.sequence = sequence

    def offset_indices_in_order(self: Self) -> list[tuple[int, int]]:
        """
        List the motifs in the solution by ascending offset.

        Returns
        -------
        offset_indices :
            Each element represents `(offset, index)` where `offset` is the
            offset where the motif starts and `index` is its index in the motif library.
        """
        order_fwd = [
            (offset, i)
            for i, offset in enumerate(self.offsets_fwd)
            if offset is not None
        ]
        order_rev = [
            (offset, i + len(self.library))
            for i, offset in enumerate(self.offsets_rev)
            if offset is not None
        ]
        # We sort by offset first and then by motif length:
        # If two motifs have the same index, we want the shortest first
        return sorted(
            order_fwd + order_rev,
            key=lambda o_i: (o_i[0], len(self.library[o_i[1] % len(self.library)])),
        )

    @property
    def nb_motifs(self: Self) -> int:
        """Number of motifs that fit in this solution."""
        nb_fwd = sum(offset is not None for offset in self.offsets_fwd)
        nb_rev = sum(offset is not None for offset in self.offsets_rev)
        return nb_fwd + nb_rev

    @property
    def compression_ratio(self: Self) -> float:
        """Compression ratio, i.e. `length of motifs in solution / sequence length`."""
        total_length = sum(
            len(motif)
            for motif, fwd, rev in zip(
                self.library,
                self.offsets_fwd,
                self.offsets_rev,
                strict=True,
            )
            if fwd is not None or rev is not None
        )
        return total_length / self.sequence_length

    def __str__(self: Self) -> str:
        """Build a string that visually represents the solution."""
        sequence = self.sequence + "-" * (self.sequence_length - len(self.sequence))
        seq_rev = "".join(COMPLEMENT[c] for c in sequence)
        lines_fwd = dispatch_labels(self.library, self.offsets_fwd, rev=False)
        lines_rev = dispatch_labels(self.library, self.offsets_rev, rev=True)

        s_fwd = "--> " + "\n--> ".join(lines_fwd[::-1] + [sequence])
        s_rev = "<-- " + "\n<-- ".join([seq_rev, *lines_rev])

        return s_fwd + "\n" + s_rev


@dataclass(init=False)
class PromoterConstraint:
    upstream_index: int
    downstream_index: int
    upstream_pos: tuple[int | None, int | None]
    downstream_pos: tuple[int | None, int | None]
    spacer_length: tuple[int | None, int | None]

    def __init__(  # noqa: PLR0913
        self: Self,
        *,
        upstream_index: int,
        downstream_index: int,
        upstream_pos: int | tuple[int | None, int | None] | None = None,
        downstream_pos: int | tuple[int | None, int | None] | None = None,
        spacer_length: int | tuple[int | None, int | None] | None = None,
    ) -> None:
        self.upstream_index = upstream_index
        self.downstream_index = downstream_index
        self.upstream_pos = (
            upstream_pos
            if isinstance(upstream_pos, tuple)
            else (upstream_pos, upstream_pos)
        )
        self.downstream_pos = (
            downstream_pos
            if isinstance(downstream_pos, tuple)
            else (downstream_pos, downstream_pos)
        )
        self.spacer_length = (
            spacer_length
            if isinstance(spacer_length, tuple)
            else (spacer_length, spacer_length)
        )


class Optimizer:
    """Optimizer."""

    def __init__(
        self: Self,
        library: list[str],
        sequence_length: int,
        strands: str = "double",
    ) -> None:
        if strands not in {"single", "double"}:
            msg = "strands must be single or double"
            raise ValueError(msg)

        self.library = list(library)
        self.sequence_length = sequence_length
        self.strands = strands
        self.promoters: list[PromoterConstraint] = []
        if strands == "double":
            library = library + [reverse_complement(motif) for motif in library]
        self.adjacency_matrix = adjacency_matrix(library)
        self.model = None

    def add_promoter_constraints(  # noqa: PLR0913
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
        try:
            upstream_index = self.library.index(upstream)
            downstream_index = self.library.index(downstream)
        except ValueError as err:
            msg = "Both promoter elements must be present in the library."
            raise ValueError(msg) from err

        constraint = PromoterConstraint(
            upstream_index=upstream_index,
            downstream_index=downstream_index,
            upstream_pos=upstream_pos,
            downstream_pos=downstream_pos,
            spacer_length=spacer_length,
        )
        self.promoters.append(constraint)

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

    def _build_model(self: Self, solver: str = "CBC") -> None:
        self.model = pywraplp.Solver.CreateSolver(solver)

        if self.model is None:
            msg = "Could not create model. There is a problem with the backend."
            raise RuntimeError(msg)

        # X_ij are binary variables. X_ij == 1 means that motif #j directly follows
        # (and possibly overlaps) motif #i in the sequence.
        start = {
            (-1, j): self.model.BoolVar(f"X[-1,{j}]") for j in range(self.nb_nodes)
        }
        end = {
            (i, -1): self.model.BoolVar(f"X[{i},-1]") for i in range(-1, self.nb_nodes)
        }
        middle = {
            (i, j): self.model.BoolVar(f"X[{i},{j}]")
            for i in range(self.nb_nodes)
            for j in range(self.nb_nodes)
            if i != j
        }
        X = start | end | middle  # noqa: N806
        self.model.X = X

        # Path starts at the start
        self.model.Add(sum(X[-1, j] for j in range(-1, self.nb_nodes)) == 1)

        # Path ends at the end
        self.model.Add(sum(X[i, -1] for i in range(-1, self.nb_nodes)) == 1)

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

        # Continuity constraints
        cont = [
            self.model.IntVar(1, self.nb_nodes, f"u[{i}]") for i in range(self.nb_nodes)
        ]
        self.model.cont = cont

        for i in range(self.nb_nodes):
            for j in range(self.nb_nodes):
                if i == j:
                    continue
                self.model.Add(cont[i] - cont[j] + 1 <= self.nb_nodes * (1 - X[i, j]))

        # Apply user-defined distance constraints if any
        if self.promoters:
            self._add_promoter_constraints()

        # Objective
        self.model.Maximize(
            sum(
                X[i, j]
                for i in range(-1, self.nb_nodes)
                for j in range(self.nb_nodes)
                if i != j
            ),
        )

    def _add_promoter_constraints(self: Self) -> None:
        """Add promoter constraints to the problem."""
        if self.model is None:
            msg = "The model must exist first (call _build_model)."
            raise AssertionError(msg)

        # Initialize cumulative length variables
        position = [
            self.model.IntVar(0, self.sequence_length - 1, f"position[{i}]")
            for i in range(self.nb_nodes)
        ]
        self.model.position = position

        # Set cumulative length for the start point
        self.model.Add(position[-1] == 0)

        # Define cumulative length for each node
        for i in range(-1, self.nb_nodes):
            for j in range(self.nb_nodes):
                if i == j:
                    continue
                shift = 0 if i == -1 else self.adjacency_matrix[i][j]
                distance_i_j = position[j] - position[i]
                slack = (self.sequence_length - 1) * (1 - self.model.X[i, j])
                self.model.Add(shift * self.model.X[i, j] - slack <= distance_i_j)
                self.model.Add(distance_i_j <= shift * self.model.X[i, j] + slack)

        for constraint in self.promoters:
            # Both upstream and downstream elements must appear in the sequence
            for k in [constraint.upstream_index, constraint.downstream_index]:
                self.model.Add(
                    sum(self.model.X[i, k] for i in range(-1, self.nb_nodes) if i != k)
                    >= 1
                )

            # Position both upstream and downstream elements
            spacer_length = (
                position[constraint.downstream_index]
                - position[constraint.upstream_index]
                - len(self.library[constraint.upstream_index])
            )
            for pos_or_len, (min_val, max_val) in [
                (position[constraint.upstream_index], constraint.upstream_pos),
                (position[constraint.downstream_index], constraint.downstream_pos),
                (spacer_length, constraint.spacer_length),
            ]:
                if min_val is not None:
                    self.model.Add(min_val <= pos_or_len)
                if max_val is not None:
                    self.model.Add(pos_or_len <= max_val)

    def _solve(self: Self) -> DenseArray:  # noqa: C901
        if self.model is None:
            msg = "Model not built: call `_build_model(solver)` first"
            raise RuntimeError(msg)

        # Solve the problem
        status = self.model.Solve()

        if status != pywraplp.Solver.OPTIMAL:
            match status:
                case pywraplp.Solver.FEASIBLE:
                    msg = "A feasible solution was found, but not necessarily optimal."
                case pywraplp.Solver.INFEASIBLE:
                    msg = "No feasible solution was found."
                case pywraplp.Solver.UNBOUNDED:
                    msg = "The model is unbounded."
                case pywraplp.Solver.ABNORMAL:
                    msg = "The model is abnormal."
                case pywraplp.Solver.NOT_SOLVED:
                    msg = "The model has not been solved."
                case _:
                    msg = f"Solver ended with unknown status: {status}."
            raise ValueError(msg)

        # Extract the solution
        sol = [-1]
        offset = 0
        offsets_fwd = [None] * self.nb_motifs
        offsets_rev = [None] * self.nb_motifs
        while sol[-1] >= 0 or len(sol) == 1:
            for j in range(-1, self.nb_nodes):
                if sol[-1] >= 0 and j == sol[-1]:
                    continue
                if round(self.model.X[sol[-1], j].solution_value()) == 1:
                    if len(sol) > 1:
                        offset += self.adjacency_matrix[sol[-1]][j]
                    if j >= self.nb_motifs:
                        offsets_rev[j % self.nb_motifs] = offset
                    elif j >= 0:
                        offsets_fwd[j] = offset
                    sol.append(j)
                    break
        sol = sol[1:-1]

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
        """
        sol = [-1, *(i for _, i in solution.offset_indices_in_order()), -1]
        sum_on_path = sum(self.model.X[i, j] for i, j in it.pairwise(sol))
        self.model.Add(sum_on_path <= solution.nb_motifs)

    def solutions(self: Self, solver: str = "CBC") -> Iterator[DenseArray]:
        """
        Iterate over solutions in decreasing order of score.

        Parameters
        ----------
        solver
            Solver name given to OrTools.

        Yields
        ------
        solution :
            Solutions in decreasing order of score.
        """
        self._build_model(solver)

        while True:
            try:
                sol = self._solve()
            except ValueError:
                break
            yield sol
            self.forbid(sol)

    def set_motif_weight(self: Self, imotif: int, weight: float) -> None:
        """
        Set the weight of a particular motif in the score.

        Parameters
        ----------
        imotif
            Index of the motif.
        weight
            Weight of the motif.
        """
        objective = self.model.Objective()

        for i in range(-1, self.nb_nodes):
            if i != imotif:
                objective.SetCoefficient(self.model.X[i, imotif], weight)
            imotif2 = imotif + self.nb_motifs
            if self.strands == "double" and i != imotif2:
                objective.SetCoefficient(self.model.X[i, imotif2], weight)

    def solutions_diverse(self: Self, solver: str = "CBC") -> Iterator[DenseArray]:
        """
        Return an iterator of optimal solutions trying to minimize the bias in motifs.

        Parameters
        ----------
        solver
            Solver name given to OrTools.

        Yields
        ------
        solution :
            Solutions in decreasing order of score.
        """
        self._build_model(solver)

        epsilon = 0.5 / self.nb_motifs

        motifs = [0] * self.nb_motifs
        imins: list[int] = []
        while True:
            try:
                sol = self._solve()
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

    def optimal(self: Self, solver: str = "CBC") -> DenseArray:
        """
        Return the optimal solution.

        Returns
        -------
        solution :
            Optimal solution.
        """
        return next(self.solutions(solver))

    def approximate(self: Self) -> DenseArray:
        """
        Return a solution approximated with a greedy algorithm.

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
        return DenseArray(self.library, self.sequence_length, offsets_fwd, offsets_rev)


def take_best_run(
    sequence: str,
    sequence_length: int,
    library: list[str],
    strands: str,
) -> str:
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


# def optimize_basepairs(
#    motifs: list[str], sequence_length: int, double: bool = False, solver: str = "CBC"
# ):
#    solver = pywraplp.Solver.CreateSolver(solver)
#
#    ATGC = list("ATGC")
#
#    sequence = [
#        [solver.IntVar(0, 1, f"{c}[{i}]") for c in ATGC]
#        for i in range(sequence_length)
#    ]
#
#    motif_present = [
#        [
#            solver.IntVar(0, 1, f"d[{i},{j}]")
#            for j in range(sequence_length - len(motifs[i]) + 1)
#        ]
#        for i in range(len(motifs))
#    ]
#
#    motif_present_summary = [solver.IntVar(0, 1, f"d[{i}]")
#                             for i in range(len(motifs))]
#
#    # There needs to be one basepair selected on every position of the sequence
#    for basepair in sequence:
#        solver.Add(sum(basepair) == 1)
#
#    for imotif, motif in enumerate(motifs):
#        for offset in range(sequence_length - len(motif) + 1):
#            solver.Add(
#                len(motif) * motif_present[imotif][offset]
#                <= sum(
#                    sequence[offset + ibp][{"A": 0, "T": 1, "G": 2, "C": 3}[bp]]
#                    for ibp, bp in enumerate(motif)
#                )
#            )
#        solver.Add(motif_present_summary[imotif] <= sum(motif_present[imotif]))
#
#    solver.Maximize(sum(motif_present_summary))
#
#    print("Number of variables =", solver.NumVariables())
#    print("Number of constraints =", solver.NumConstraints())
#
#    status = solver.Solve()
#
#    if status == pywraplp.Solver.OPTIMAL:
#        print("Optimal solution found!")
#        print("Objective value =", round(solver.Objective().Value()))
#        solution = []
#        for basepair in sequence:
#            for var, bp in zip(basepair, ATGC):
#                if round(var.solution_value()) == 1:
#                    solution.append(bp)
#                    break
#        return "".join(solution)
#
#    print("The problem does not have an optimal solution.")
