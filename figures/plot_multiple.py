"""
Plot "multiple_promoters.svg" from the file "multiple_promoters.csv".
"""
from pathlib import Path

import dense_arrays as da
import matplotlib.pyplot as plt
import numpy as np

motifs = np.zeros(100, dtype=np.float64)
upstream = np.zeros((3, 100), dtype=np.float64)
downstream = np.zeros((3, 100), dtype=np.float64)

n = 0
library_size = 20
motif_min_size = 5
motif_max_size = 15
sequence_length = 100
upstreams = ["TTGACA", "TTGACA", "TGGCAGG"]
downstreams = ["TATAAT", "TATACT", "TTGCA"]
ATGC = list("ATGC")

with Path("../benchmarks/multiple_promoters.csv").open("r") as f:
    lines = iter(f.readlines())
    while True:
        try:
            library = next(lines).rstrip().split(",")[1:]
            positions = next(lines).rstrip().split(",")
        except StopIteration:
            break
        n += 1

        offsets_fwd = [int(x) if x != "None" else None for x in positions[1:27]]
        offsets_rev = [int(x) if x != "None" else None for x in positions[27:]]

        sol = da.DenseArray(library, sequence_length, offsets_fwd, offsets_rev)
        for offset, index in sol.offset_indices_in_order():
            if 20 <= index % 26 < 23:  # noqa: PLR2004
                iupstream = index % 26 - 20
                upstream[iupstream, offset : offset + len(upstreams[iupstream])] += 1
            elif 23 <= index % 26 < 26:  # noqa: PLR2004
                idownstream = index % 26 - 23
                downstream[
                    idownstream, offset : offset + len(upstreams[idownstream])
                ] += 1
            else:
                motifs[offset : offset + len(library[index % 26])] += 1

fig, ax = plt.subplots(figsize=(7, 4))

bps = range(sequence_length)

ax.plot(bps, motifs, label="normal motifs", color="C0")
for ipromoter, (up, down) in enumerate(zip(upstream, downstream, strict=True), start=1):
    ax.plot(bps, up, label=f"#{ipromoter} up", color=f"C{ipromoter}")
    ax.plot(bps, down, "--", label=f"#{ipromoter} down", color=f"C{ipromoter}")

ax.set_xlabel("Sequence position (bp)")
ax.set_ylabel("Number of motifs\ncovering this position")
ax.legend(loc="lower left")
plt.savefig("multiple_promoters.png", bbox_inches="tight")
plt.savefig("multiple_promoters.svg", bbox_inches="tight")
