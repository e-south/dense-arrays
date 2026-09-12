# Render saved feature placements

Playback explains an existing sequence and its persisted feature placements.
The public seam is `RealizedArray` → `PlaybackPlan` → renderer. It does not
rerun optimization. Producer adapters translate their own records into this
contract; the optimizer's `DenseArray` result is a separate interface.

## Create a self-contained example

Run this Python example in the [installed checkout](quickstart.md#install-from-source).
It describes the same three synthetic placements as the first-array example,
without depending on a solver or an external data file. Outputs go to a new
temporary directory whose path is printed:

```python
from pathlib import Path
from tempfile import mkdtemp

from dense_arrays.playback import (
    dumps_playback_plan,
    dumps_realized_array,
    reconstruct_playback,
    render_playback_html,
)
from dense_arrays.realized import Orientation, Placement, PlacementKind, RealizedArray

realized = RealizedArray(
    source_id="synthetic:first-array",
    sequence="CAGCGT",
    placements=tuple(
        Placement(
            placement_id=identifier,
            feature_id=identifier,
            kind=PlacementKind.OTHER,
            sequence=motif,
            start=start,
            orientation=Orientation.FORWARD,
            label=motif,
        )
        for identifier, motif, start in (
            ("motif-1", "CAG", 0),
            ("motif-2", "AGC", 1),
            ("motif-3", "CGT", 3),
        )
    ),
)
plan = reconstruct_playback(realized)
output = Path(mkdtemp(prefix="dense-arrays-playback-"))
(output / "realized.json").write_text(dumps_realized_array(realized), encoding="utf-8")
(output / "plan.json").write_text(dumps_playback_plan(plan), encoding="utf-8")
(output / "playback.html").write_text(
    render_playback_html(plan, title="Three overlapping motifs"), encoding="utf-8"
)
print(output)
print(plan.authority.value)  # placement_reconstructed
```

Open `playback.html` in a browser. The HTML is self-contained and requires only
the core package. Coordinates are zero-based and half-open: `CAG` occupies
`[0, 3)`, `AGC` occupies `[1, 4)`, and `CGT` occupies `[3, 6)`.
Each placement sequence is already oriented to the realized sequence.

## Render serialized input

From the repository root, activate the installed environment:

```bash
source .venv/bin/activate
```

Then change to the output directory printed by the Python example. Render
either strict `RealizedArray` or `PlaybackPlan` JSON:

```bash
dense-arrays-playback realized.json --html rendered.html
```

For saved files, the Python entrypoints are `loads_realized_array()` and
`loads_playback_plan()`; both accept the file's JSON text. The matching
`dumps_*()` functions return JSON text. Unknown schema fields are rejected.
See the [playback API](api.md#realized-arrays-and-playback).

## Export a still or video

From the repository root, install the optional renderer dependencies:

```bash
uv sync --frozen --extra playback
```

The installed `dense-arrays-playback` command can now write a PNG poster:

```bash
dense-arrays-playback realized.json --html rendered.html --poster poster.png
```

MP4 export additionally requires a local FFmpeg executable on `PATH`:

```bash
dense-arrays-playback realized.json --html rendered.html --mp4 playback.mp4
```

Run these render commands from the input directory as above, or pass explicit
input and output paths. `dense-arrays-playback --help` lists all export options.

## Interpret the result

Reconstruction checks placement bounds, identities, sequence agreement, and
constraint references. Declared distance constraints are evaluated in the plan.
The plan always reports `placement_reconstructed` authority. Its ordering status
distinguishes a unique coordinate order, an ambiguous order requiring a
deterministic tie-break, and a layout with internal uncovered spans.

The displayed order is a coordinate explanation. It is not the optimizer's
recorded search or selected path. `solver_selected` authority is reserved for
future exact traces. Preserve these distinctions when adding captions or
adapting producer data; the [playback contract](architecture/solution-playback.md)
is the authority for integration details.
