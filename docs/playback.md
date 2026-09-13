---
title: Saved-placement playback
description: Build a realized-array record and render its placements as a PNG, MP4, or GIF.
---

# Render saved feature placements

Turn saved feature placements into a PNG, MP4, or GIF to inspect their positions
and overlaps. This guide follows the four 16-base motifs from the
[first-array example](quickstart.md) across their 40-base sequence.

## Watch four overlapping motifs

<p id="playback-example-description">The opening frame shows the complete graph
and DNA duplex in gray. Color then follows the four placements in coordinate
order: the first motif covers 16 bases, and each later motif adds eight.
Overlapping portions share the same sequence positions throughout.</p>

<video controls preload="metadata" poster="../assets/playback-opening.png"
       aria-describedby="playback-example-description" style="width: 100%; height: auto;">
  <source src="../assets/playback-example.mp4" type="video/mp4">
  Your browser cannot display this video. Use the download link below.
</video>

[Download the MP4](assets/playback-example.mp4),
[view the opening frame](assets/playback-opening.png), or
[inspect the completed poster](assets/playback-poster.png).

The video explains saved coordinates; it does not rerun optimization or show
the solver's search history. Reproduce the placements and exports below.

## Create a PNG example

From the [installed checkout](quickstart.md#install-from-source), install the
playback dependencies:

```bash
uv sync --frozen --extra playback
```

Then run this Python example with `uv run python`. It describes the same four
synthetic placements as the first-array example, without running a solver or
reading an external data file. Outputs go to a new temporary directory whose
path is printed:

```python
from pathlib import Path
from tempfile import mkdtemp

from dense_arrays.playback import (
    PlaybackDocument,
    dumps_playback_plan,
    dumps_realized_array,
    reconstruct_playback,
)
from dense_arrays.playback.matplotlib_renderer import render_collection_poster_png
from dense_arrays.realized import Orientation, Placement, PlacementKind, RealizedArray

realized = RealizedArray(
    source_id="synthetic:first-array",
    sequence="ACGTTGCAAGTCCTGATCGTACCGATGCTTAGGACGTTCA",
    placements=tuple(
        Placement(
            placement_id=identifier,
            feature_id=identifier,
            kind=PlacementKind.OTHER,
            sequence=motif,
            start=start,
            orientation=Orientation.FORWARD,
        )
        for identifier, motif, start in (
            ("motif-1", "ACGTTGCAAGTCCTGA", 0),
            ("motif-2", "AGTCCTGATCGTACCG", 8),
            ("motif-3", "TCGTACCGATGCTTAG", 16),
            ("motif-4", "ATGCTTAGGACGTTCA", 24),
        )
    ),
)
plan = reconstruct_playback(realized)
output = Path(mkdtemp(prefix="dense-arrays-playback-"))
(output / "realized.json").write_text(dumps_realized_array(realized), encoding="utf-8")
(output / "plan.json").write_text(dumps_playback_plan(plan), encoding="utf-8")
document = PlaybackDocument(plan=plan, title="Four overlapping motifs")
poster = render_collection_poster_png((document,), output / "poster.png")
assert poster.is_file() and poster.stat().st_size > 0
print(output)
print(plan.authority.value)  # placement_reconstructed
assert plan.ordering_status.value == "unique"
assert tuple((span.start, span.end) for span in plan.steps[1].added_spans) == (
    (16, 24),
)
```

Open `poster.png` to inspect the four overlapping motifs. Coordinates are
zero-based and half-open:

| Placement | Sequence | Occupied span |
| --- | --- | --- |
| motif-1 | `ACGTTGCAAGTCCTGA` | `[0, 16)` |
| motif-2 | `AGTCCTGATCGTACCG` | `[8, 24)` |
| motif-3 | `TCGTACCGATGCTTAG` | `[16, 32)` |
| motif-4 | `ATGCTTAGGACGTTCA` | `[24, 40)` |

Each placement sequence is already oriented to the realized sequence. The
poster and video exports share the NetworkX layout and Matplotlib renderer.

## Interpret the result

Every v1 plan uses `placement_reconstructed` authority: its order comes from
coordinates. The order can be:

- `unique`: a strict coordinate order.
- `ambiguous`: equal starts or containment require a deterministic tie-break.
- `layout_only`: internal uncovered spans prevent a complete placement chain;
  the renderer shows the layout without an active traversal chain.

Reconstruction also checks declared distances. A layout that violates a
requirement retains that result as `passed=False`; rendering it does not make
the requirement pass. Media keeps the reconstructed authority, ordering
qualifications, and failed requirements visible. Long failure details are
retained in native media metadata. See
[how to read the full evidence](reference/playback-presentation.md#read-the-evidence).

## Render serialized input

From the repository root, activate the installed environment:

```bash
source .venv/bin/activate
```

Then change to the output directory printed by the Python example. Render
either `RealizedArray` or `PlaybackPlan` JSON:

```bash
dense-arrays-playback realized.json --poster rendered.png
```

For saved files, the Python entrypoints are `loads_realized_array()` and
`loads_playback_plan()`; both accept the file's JSON text. The matching
`dumps_*()` functions return JSON text. Python constructors and JSON loaders
share semantic validation, including placement bounds, sequence agreement,
references, and reveal geometry. Invalid input is rejected before any requested
export is published. See [validation details](reference/playback.md#interpretation-and-validation).

## Export a still or video

The command accepts `--poster`, `--mp4`, and `--gif`; request at least one.
The playback extra supplies the PNG and GIF dependencies. MP4 additionally
requires a local FFmpeg executable on `PATH`:

```bash
dense-arrays-playback realized.json --mp4 playback.mp4
```

Request a poster and an animated GIF together with:

```bash
dense-arrays-playback realized.json --poster poster-copy.png --gif playback.gif
```

Existing outputs require `--replace`. Input/output aliases and colliding
destinations are rejected. The command renders every requested format before
publishing any destination, so a rendering failure leaves prior outputs
untouched. Publication is atomic per file; a filesystem failure during publication
reports which files were already written. See [CLI export behavior](reference/cli.md#render-saved-placements).
Run these render commands from the input directory as above, or pass explicit
input and output paths. `dense-arrays-playback --help` lists all export options.

## Adapt your own records

Supply the sequence and placements as a `RealizedArray`, then call
`reconstruct_playback()` to build the `PlaybackPlan` used by the renderer.
The optimizer's `DenseArray` result is a separate interface; producer adapters
own the translation to saved placements.

If an adapter recovered coordinates, pass evidence of that procedure through
`reconstruct_playback(realized, notices=(notice,))`. Metadata names alone do not
establish how coordinates were obtained. The
[playback reference](reference/playback.md) defines accepted records and notices.

Use [media presentation settings](reference/playback-presentation.md) to choose
labels, colors, graph detail, and optional notice summaries. Follow the
[ownership and evidence rules](architecture/solution-playback.md) when adapting
producer data or adding publication captions.
