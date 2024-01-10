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
)

opt.add_promoter_constraints(
    upstream="CCCCGCC",
    downstream="GGGGGGG",
    upstream_pos=(5, 10),
    spacer_length=(3, 13),
)


best = opt.optimal()
print(f"Optimal solution, score {best.nb_motifs}")
print(best)

# print("List of all solutions")
# for solution in opt.solutions():
#     print(solution)
