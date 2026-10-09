---
title: Assemble an exact-length design
description: Place named binding sites, add bounded padding and verify the final sequence.
author: Eric J. South
---

# Assemble an exact-length design

Install the [playback extra](../../installation.md#optional-features) for the figure.
Run this example in a new directory; it creates `runs/exact` and two PNG files.

```python
import dense_arrays as da  # Generate, inspect and render saved records.
from dense_arrays import (
    parts,
    planning,
    reporting,
)  # Declare parts, rules and selection.
```

## Assemble and render an exact-length design

Choose a final length and declare whether padding is allowed. `Assembly()`
without padding requires the packing itself to fill that length. It never
invents missing bases. This example adds bounded right padding:

```python
# Pack three 16-base sites before adding bounded right padding.
exact_request = planning.DesignSpec(
    parts=(
        parts.Part("u", "ACGTTGCAAGTCCTGA"),
        parts.Part("bridge", "AGTCCTGATCGTACCG"),
        parts.Part("d", "TCGTACCGATGCTTAG"),
    ),
    length=planning.Length(exact=40),
    assembly=planning.Assembly(padding=planning.Padding(side="right", max_trials=20)),
    strands="single",
    requirements=(
        planning.Fixed("upstream", "u", "forward", planning.StartWindow(max=0)),
        planning.Fixed("downstream", "d", "forward"),
        planning.Spacing("adjacent-anchors", "u", "d", min=0, max=0),
        planning.Avoid("no-G-run", patterns=("GGGG",), strands="both"),
        planning.GC("final-gc", scope="sequence", min=0.2, max=0.8),
    ),
    seed=7,  # Fix seeded sampling and padding streams.
)
exact_run = da.run(
    exact_request, out="runs/exact"
)  # Persist placements and final screens.
assert da.inspect(exact_run, verify=True).accepted == 1
with da.inspect(exact_run, view="designs").records() as records:
    exact_design = next(records)
assert len(exact_design.realized.sequence) == 40  # Check the final, padded length.
assert exact_design.realized.provenance["assembly"]["padding_length"] == 8
```

Fixed requirements name supplied part IDs and a forward/reverse orientation.
Start windows are inclusive and refer to the final sequence. Left padding
translates those windows into packing coordinates, then rechecks them after
assembly. Spacing is downstream start minus upstream end, regardless of strand;
negative spacing means intentional overlap. One fixed pair is supported.

Final checks include junctions and padding. `Avoid` accepts literal A/C/G/T
patterns and checks both strands unless configured otherwise. Its optional
`except_placements` names fixed part IDs; only matches wholly inside those
realized intervals are exempt. GC limits use fractions, without rounding a
candidate into acceptance. `scope="padding"` requires a padding policy; zero
added bases are reported as `not_applicable`.

For PNG rendering, install the [playback extra](../../installation.md#optional-features). Rendering reads persisted placements and requirement
evidence; it does not solve or repeat final screening:

```python
# Draw the selected saved design without solving again.
receipt = da.render(
    exact_run,
    select=reporting.DesignFilter(design_ids=(exact_design.reference,)),
    out="exact-design.png",
)
assert receipt.records == 1
assert receipt.design_refs == (exact_design.reference,)
```

```bash
# Render saved evidence to a new image.
dense-arrays render runs/exact --out exact-design-cli.png --json
```

The `design` view requires exactly one selected design and a `.png` destination.
For a larger library, add `--design-id REF` using a full reference from
`inspect --view designs --json`. Rendering accepts the same design filters,
bounded selections and saved panels as inspection. Zero or multiple matches
fail before output; the command never chooses an implicit first design.
Native runs, portable bundles and combined sources retain the selected design's
own cell plan. A saved panel also retains its original source revision.

Use `library-quality` below for aggregate reporting. Negative spacing is valid
generation evidence, but playback v1 cannot represent it; rendering fails before
creating output. Existing [playback APIs](../../playback.md) remain available for their
supported saved-record and media formats.
