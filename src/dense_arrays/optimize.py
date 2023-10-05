"""Optimization module."""

from collections.abc import Iterator
from dataclasses import dataclass

from ortools.linear_solver import pywraplp

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "-": "-"}

__all__ = ["DenseArray", "Optimizer"]


def shift_metric(motifa: str, motifb: str) -> int:
    """
    Compute how much we have to shift `motifb` to match the end of `motifa`.

    Example
    -------
    shift_metric("ATGCATTA", "CATTATG") == 3 because

        motifa: ATGCATTA
        motifb:    CATTATG
        shift : 0123

    and shift_metric("ATGCATTA", "TATGA") == 8 because

        motifa: ATGCATTA
        motifb:         TATGA
        shift : 012345678

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
    """
    Return the matrix A_ij such that A_ij = shift_metric(motifs[i], motifs[j]).

    Parameters
    ----------
    motifs : list[str]
        List of motifs.
    """
    return [[shift_metric(motifa, motifb) for motifb in motifs] for motifa in motifs]


def reverse_complement(sequence: str) -> str:
    return "".join(COMPLEMENT[c] for c in sequence[::-1])


def dispatch_labels(
    library: list[str], offsets: list[int | None], rev: bool
) -> list[str]:
    lines = []
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
    sequence_size: int
    sequence: str
    offsets_fwd: list[int | None]
    offsets_rev: list[int | None]

    def __init__(
        self,
        library: list[str],
        sequence_size: int,
        offsets_fwd: list[int | None],
        offsets_rev: list[int | None],
    ):
        self.library = library
        self.sequence_size = sequence_size
        self.offsets_fwd = offsets_fwd
        self.offsets_rev = offsets_rev
        sequence = ""
        for offset, i in self._offset_indices_in_order():
            motif = library[i % len(library)]
            if i >= len(library):
                motif = reverse_complement(motif)
            sequence = sequence[:offset] + motif
        self.sequence = sequence

    def _offset_indices_in_order(self):
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
        return sorted(order_fwd + order_rev)

    @property
    def nb_motifs(self):
        """Number of motifs that fit in this solution."""
        nb_fwd = sum(offset is not None for offset in self.offsets_fwd)
        nb_rev = sum(offset is not None for offset in self.offsets_rev)
        return nb_fwd + nb_rev

    @property
    def compression_ratio(self):
        """Compression ratio, i.e. total length of motifs / solution size."""
        total_length = sum(
            len(motif)
            for motif, fwd, rev in zip(self.library, self.offsets_fwd, self.offsets_rev)
            if fwd is not None or rev is not None
        )
        return total_length / self.sequence_size

    def __str__(self) -> str:
        """Str dunder."""
        sequence = self.sequence + "-" * (self.sequence_size - len(self.sequence))
        seq_rev = "".join(COMPLEMENT[c] for c in sequence)
        lines_fwd = dispatch_labels(self.library, self.offsets_fwd, False)
        lines_rev = dispatch_labels(self.library, self.offsets_rev, True)

        s_fwd = "--> " + "\n--> ".join(lines_fwd[::-1] + [sequence])
        s_rev = "<-- " + "\n<-- ".join([seq_rev] + lines_rev)

        return s_fwd + "\n" + s_rev


class Optimizer:
    """Optimizer."""

    def __init__(self, library: list[str], sequence_size: int, strands: str = "double"):
        if strands not in ["single", "double"]:
            return ValueError("strands must be single or double")

        self.library = list(library)
        self.sequence_size = sequence_size
        self.strands = strands
        if strands == "double":
            library = library + [reverse_complement(motif) for motif in library]
        self.adjacency_matrix = adjacency_matrix(library)
        self.model = None

    def _build_model(self, solver: str = "CBC") -> None:
        nb_motifs = len(self.library)
        nb_nodes = nb_motifs if self.strands == "single" else 2 * nb_motifs

        self.model = pywraplp.Solver.CreateSolver(solver)

        # X_ij are binary variables. X_ij == 1 means that motif #j directly follows
        # (and possibly overlaps) motif #i in the sequence.
        start = {(-1, j): self.model.BoolVar(f"X[-1,{j}]") for j in range(nb_nodes)}
        end = {(i, -1): self.model.BoolVar(f"X[{i},-1]") for i in range(-1, nb_nodes)}
        middle = {
            (i, j): self.model.BoolVar(f"X[{i},{j}]")
            for i in range(nb_nodes)
            for j in range(nb_nodes)
            if i != j
        }
        X = start | end | middle
        self.model.X = X

        # Path starts at the start
        self.model.Add(sum(X[-1, j] for j in range(-1, nb_nodes)) == 1)

        # Path ends at the end
        self.model.Add(sum(X[i, -1] for i in range(-1, nb_nodes)) == 1)

        # Conservation of flow
        for k in range(nb_nodes):
            self.model.Add(
                sum(X[i, k] for i in range(-1, nb_nodes) if i != k)
                == sum(X[k, j] for j in range(-1, nb_nodes) if j != k)
            )

        # Don't include any motif more than once
        if self.strands == "single":
            for k in range(nb_nodes):
                self.model.Add(sum(X[i, k] for i in range(-1, nb_nodes) if i != k) <= 1)
                self.model.Add(sum(X[k, j] for j in range(-1, nb_nodes) if j != k) <= 1)
        else:
            for k in range(nb_nodes // 2):
                krev = k + nb_nodes // 2
                enter_direct = sum(X[i, k] for i in range(-1, nb_nodes) if i != k)
                enter_rev = sum(X[i, krev] for i in range(-1, nb_nodes) if i != krev)
                self.model.Add(enter_direct + enter_rev <= 1)
                exit_direct = sum(X[k, j] for j in range(-1, nb_nodes) if k != j)
                exit_rev = sum(X[krev, j] for j in range(-1, nb_nodes) if krev != j)
                self.model.Add(exit_direct + exit_rev <= 1)

        # Global length constraint
        size_inside = sum(
            self.adjacency_matrix[i][j] * X[i, j]
            for i in range(nb_nodes)
            for j in range(nb_nodes)
            if i != j
        )
        size_terminal = sum(
            len(self.library[i % nb_motifs]) * X[i, -1] for i in range(nb_nodes)
        )
        self.model.Add(size_inside + size_terminal <= self.sequence_size)

        # Continuity constraints
        cont = [self.model.IntVar(1, nb_nodes, f"u[{i}]") for i in range(nb_nodes)]
        self.model.cont = cont

        for i in range(nb_nodes):
            for j in range(nb_nodes):
                if i == j:
                    continue
                self.model.Add(cont[i] - cont[j] + 1 <= nb_nodes * (1 - X[i, j]))

        # Objective
        self.model.Maximize(
            sum(X[i, j] for i in range(-1, nb_nodes) for j in range(nb_nodes) if i != j)
        )

    def _solve(self) -> DenseArray:
        if self.model is None:
            raise RuntimeError("Model not built: call _build_model(solver) first")

        nb_motifs = len(self.library)
        nb_nodes = nb_motifs if self.strands == "single" else 2 * nb_motifs

        # Solve the problem
        status = self.model.Solve()

        if status == pywraplp.Solver.OPTIMAL:
            # Extract solution
            sol = [-1]
            offset = 0
            offsets_fwd = [None] * nb_motifs
            offsets_rev = [None] * nb_motifs
            while sol[-1] >= 0 or len(sol) == 1:
                for j in range(-1, nb_nodes):
                    if sol[-1] >= 0 and j == sol[-1]:
                        continue
                    if round(self.model.X[sol[-1], j].solution_value()) == 1:
                        if len(sol) > 1:
                            offset += self.adjacency_matrix[sol[-1]][j]
                        if j >= len(self.library):
                            offsets_rev[j % nb_motifs] = offset
                        elif j >= 0:
                            offsets_fwd[j] = offset
                        sol.append(j)
                        break
            sol = sol[1:-1]
            return DenseArray(
                self.library, self.sequence_size, offsets_fwd, offsets_rev
            )

    def forbid(self, solution: DenseArray) -> None:
        """Add a constraint to the model to forbid a given solution."""
        sol = [-1, *(i for _, i in solution._offset_indices_in_order()), -1]
        self.model.Add(
            sum(self.model.X[i, j] for i, j in zip(sol, sol[1:])) <= solution.nb_motifs
        )

    def solutions(self, solver: str = "CBC") -> Iterator[DenseArray]:
        """Iterate over solutions in decreasing order of score."""
        self._build_model(solver)

        sol = self._solve()
        while sol is not None:
            yield sol
            self.forbid(sol)
            sol = self._solve()

    def solutions_diverse(self, solver: str = "CBC") -> Iterator[DenseArray]:
        """Return an iterator of optimal solutions trying to minimize the bias in motifs."""
        self._build_model(solver)

        nb_nodes = len(self.library)
        if self.strands == "double":
            nb_nodes *= 2

        motifs = [0] * len(self.library)
        # Define a dummy constraint, clear it and add it to the solver (I don't
        # know how to generate an empty constraint otherwise)
        constraint = self.model.Constraint()
        constraint.SetBounds(0, 0)
        imax = 0
        while True:
            print(motifs)
            sol = self._solve()
            if sol is None:
                break
            yield sol
            # Forbid the solution
            self.forbid(sol)
            # Tally up the motifs
            for i, (fwd, rev) in enumerate(zip(sol.offsets_fwd, sol.offsets_rev)):
                if fwd is not None or rev is not None:
                    motifs[i] += 1
            # Update the constraint to forbid the most common motif
            for i in range(-1, nb_nodes):
                if i != imax:
                    constraint.SetCoefficient(self.model.X[i, imax], 0)
                    constraint.SetCoefficient(self.model.X[imax, i], 0)
                imax2 = imax + len(self.library)
                if self.strands == "double" and i != imax2:
                    constraint.SetCoefficient(self.model.X[i, imax2], 0)
                    constraint.SetCoefficient(self.model.X[imax2, i], 0)

            imax = motifs.index(max(motifs))
            for i in range(-1, nb_nodes):
                if i != imax:
                    constraint.SetCoefficient(self.model.X[i, imax], 1)
                    constraint.SetCoefficient(self.model.X[imax, i], 1)
                imax2 = imax + len(self.library)
                if self.strands == "double" and i != imax2:
                    constraint.SetCoefficient(self.model.X[i, imax2], 1)
                    constraint.SetCoefficient(self.model.X[imax2, i], 1)

    def optimal(self, solver: str = "CBC") -> DenseArray:
        """Return the optimal solution."""
        return next(self.solutions(solver))


## TODO: put into Optimizer
# def approximate(motifs: list[str], sequence_size: int, double: bool = False) -> list[int]:
#    assert not double
#
#    old_motifs = list(motifs)
#    motifs = list(motifs)
#
#    while len(motifs) > 1:
#        adj = adjacency_matrix(motifs)
#        min_dist = max(len(motif) for motif in motifs)
#        index_min_dist = None
#        for i in range(len(motifs)):
#            for j in range(len(motifs)):
#                if i == j:
#                    continue
#                if adj[i][j] < min_dist:
#                    min_dist = adj[i][j]
#                    index_min_dist = (i, j)
#        i, j = index_min_dist
#        motifs[i] = motifs[i][: adj[i][j]] + motifs[j]
#        del motifs[j]
#    print(motifs)
#    max_nb_motifs = 0
#    offset_max_nb_motifs = None
#    for offset in range(len(motifs[0]) - sequence_size + 1):
#        subseq = motifs[0][offset : offset + sequence_size]
#        nb_motifs = sum(motif in subseq for motif in old_motifs)
#        if nb_motifs > max_nb_motifs:
#            max_nb_motifs = nb_motifs
#            offset_max_nb_motifs = offset
#    subseq = motifs[0][offset_max_nb_motifs : offset_max_nb_motifs + sequence_size]
#    present = [imotif for imotif, motif in enumerate(old_motifs) if motif in subseq]
#    indices = [(subseq.index(old_motifs[imotif]), imotif) for imotif in present]
#    indices.sort()
#    return [imotif for _, imotif in indices]
#
#
# def optimize_basepairs(
#    motifs: list[str], sequence_size: int, double: bool = False, solver: str = "CBC"
# ):
#    solver = pywraplp.Solver.CreateSolver(solver)
#
#    ATGC = list("ATGC")
#
#    sequence = [
#        [solver.IntVar(0, 1, f"{c}[{i}]") for c in ATGC] for i in range(sequence_size)
#    ]
#
#    motif_present = [
#        [
#            solver.IntVar(0, 1, f"d[{i},{j}]")
#            for j in range(sequence_size - len(motifs[i]) + 1)
#        ]
#        for i in range(len(motifs))
#    ]
#
#    motif_present_summary = [solver.IntVar(0, 1, f"d[{i}]") for i in range(len(motifs))]
#
#    # There needs to be one basepair selected on every position of the sequence
#    for basepair in sequence:
#        solver.Add(sum(basepair) == 1)
#
#    for imotif, motif in enumerate(motifs):
#        for offset in range(sequence_size - len(motif) + 1):
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
