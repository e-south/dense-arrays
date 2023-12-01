"""
Dense Arrays allows to create densely packed string arrays from a library of motifs.

It models the String Packing Problem as an Integer Linear Programming problem.
"""

import dense_arrays as da

motifs = [
    "ATGC",
    "CGT",
    "ATTA",
    "TTATTA",
    "GGGGGGG",
    "CCCCGCC",
    "CATGAGAT",
    "CAGGAAA",
    "ACCGGA",
    "CGGCATTA",
    "TATCCCG",
]


opt = da.Optimizer(
    library=motifs,
    sequence_length=30,
    strands="double",
    special_motif_a="CCCCGCC",
    special_motif_b="GGGGGGG",
    a_start_min=5,
    a_start_max=10,
    a_b_min_distance=10,
    a_b_max_distance=20,
)


best = opt.optimal()
print(f"Optimal solution, score {best.nb_motifs}")
print(best)

# print("List of all solutions")
# for solution in opt.solutions():
#     print(solution)
