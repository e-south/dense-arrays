---
title: Packing method
description: Understand compatible overlaps and the integer optimization formulation.
---

# How motifs share sequence space

Dense Arrays formulates motif packing as an integer optimization problem.
Given exact DNA strings and a sequence-length limit, compatible suffixes and
prefixes can share bases. Selecting an order with useful overlaps allows more
motif entries to fit into that limit. Double-strand optimization also admits
reverse-complement orientations.

For example, `CAG` overlaps `AGC` by two bases, and `AGC` overlaps `CGT` by one.
Together they occupy `CAGCGT`: six bases for three three-base motifs. The
banner uses this synthetic arrangement. It illustrates compatible placement,
without asserting a biological effect or a recorded solver trajectory.

![String packing formulation: motif library and length limit, pairwise shifts, an oriented graph, and example packed sequences](assets/SPP_overview.png)

The figure connects the nucleotide String Packing Problem to an Orienteering
Problem: motifs become oriented graph nodes and transitions account for the
sequence span needed to place successive motifs. The optimizer chooses a
feasible arrangement within the requested length. The Python result retains
the realized sequence and motif offsets.

The exact model counts selected entries along this path. It does not count
every motif that happens to occur as a substring: an entry wholly contained
inside another does not get an additional placement through the overlap metric.
Repeated identical strings also require distinct path placements. For example,
`ACGTA` and `CGT` with a five-base limit yield one selected entry in the exact
model. Distinguish selected entries from all sequence matches when reporting a
motif count.

See the [associated paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012276)
for the formulation and scientific context. When citing software results,
also record the package version or commit used.

## Additional requirements

The Python API can enforce positional relationships and regulator coverage,
or express side preferences. These act on the supplied strings, labels, and
coordinates. Their [guide](constraints.md) explains what each requirement means
and supplies runnable examples.

## Playback is a separate view

A producer can supply saved sequence placements to the
[playback interface](playback.md). Reconstruction orders those placements by
coordinates; it does not recover an unrecorded optimizer trace. Producer
translation, study labels, and biological interpretation remain with their
respective owners. See [architecture](architecture/README.md) for the boundaries.
