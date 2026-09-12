---
title: Media presentation
description: Set labels, colors, and producer sequence frames for the established playback renderer.
author: Eric J. South
---

# Media presentation

PNG, MP4, and GIF share the NetworkX graph layout and Matplotlib scene renderer.
A `PlaybackDocument` combines the validated plan with artifact metadata,
placement labels, colors, and presentation choices. Install the `playback`
extra before rendering.

Use the `plan` and output directory from the [PNG example](../playback.md):

```python
from dense_arrays.playback import PlaybackDocument
from dense_arrays.playback.matplotlib_renderer import render_collection_poster_png
from dense_arrays.playback.theme import LegendEntry, PlaybackPresentation

presentation = PlaybackPresentation(
    graph_detail="reduced",
    graph_fraction=0.3,
    show_authority_notice=True,
    legend_entries=(LegendEntry("selected", "Selected motif", "#365E80"),),
)
document = PlaybackDocument(
    plan=plan,
    title="Three overlapping motifs",
    subtitle="Coordinate reconstruction from persisted placements",
    label_overrides={"motif-1": "First selected motif"},
    color_overrides={"motif-1": "#365E80"},
    presentation=presentation,
)
poster = render_collection_poster_png((document,), output / "labeled-poster.png")
assert poster.is_file()
```

Labels and colors are keyed by **placement ID**. Unknown IDs are rejected;
colors use opaque `#RRGGBB` notation. Mappings are immutable snapshots, so
later changes to caller dictionaries cannot change the document.

## Scene settings

The scene preserves resolved labels and a compact evidence summary. Authority,
ordering, and failed requirements remain visible when optional notices are
disabled. Titles and subtitles provide artifact metadata.

| Setting | Accepted values and effect |
| --- | --- |
| `color_profile` | `categorical` cycles generic colors; `uniform` uses one color; `constraints` distinguishes declared fixed elements |
| `color_overrides` | Explicit placement colors on the document, taking precedence over the profile |
| `legend_entries` | Caller-authored `LegendEntry(key, label, color)` records with unique keys |
| `graph_detail` | `full` includes context and traversal relations; `reduced` includes traversal relations at the chosen graph fraction; `none` hides the graph |
| `graph_fraction` | Finite number from 0 to 0.5; it must be zero exactly when `graph_detail="none"` |
| `show_edge_costs` | Boolean controlling edge-cost labels |
| `show_authority_notice` | Boolean including plan notices in the visible summary and full evidence metadata |
| `show_distance_bracket` | `never`, `when_declared`, or `always`; brackets describe declared constraints, and `always` reports when none exist |

A `layout_only` document never animates a complete placement chain, regardless
of graph detail. No profile infers biological identity from a label or ID.
Study-specific profiles such as `secg` are unsupported; supply caller colors
and legend entries instead.

## Read the evidence

The visible summary uses at most three lines and reserves no more than 25% of
the figure height. It always identifies reconstructed authority and coordinate
order. A short failed requirement shows its actual and required distance;
longer failure text becomes a `FAILED` count directing readers to metadata.
Long optional notices are excerpted. This keeps the figure compact while
native metadata retains every declared distance result, including passed results,
and the full text of enabled notices. The native sequence panel summarizes more
than two distance brackets, or overlong bracket labels, in the same way. Producer
callbacks retain ownership of their own annotations.

| Format | Native metadata |
| --- | --- |
| PNG | `Title` and `Description` hold the first scene's title and subtitle; `PlaybackEvidence` holds its full evidence |
| MP4 | `title` joins scene titles; `comment` contains one title, subtitle, and full evidence block per scene |
| GIF | The native `comment` contains the same scene blocks, encoded as UTF-8 |

Read the PNG evidence from the example above with Pillow:

```python
from PIL import Image

with Image.open(poster) as image:
    evidence = image.info["PlaybackEvidence"]
print(evidence)
assert "Reconstructed from placements" in evidence
```

For GIF, read `image.info["comment"].decode("utf-8")` instead. For MP4, inspect
the container's `comment` tag with a media metadata tool such as FFprobe.
Titles and subtitles do not add visual headings to the scenes in any format.

## Export and timing

Public media functions live in `dense_arrays.playback.matplotlib_renderer`:

```text
render_collection_poster_png(documents, output_path, ...)
render_collection_mp4(documents, output_path, ...)
render_collection_gif(documents, output_path, ...)
```

The PNG poster shows the completed first scene. MP4 and GIF render the
collection in order. All formats require the `playback` extra; MP4 also
requires FFmpeg on `PATH`.

Media timing requires a positive integer `fps` and finite positive
`seconds_per_step`. Lead, hold, and scene-transition durations are finite and
non-negative. Zero lead or hold means zero frames for that phase. MP4 and GIF
use the same frame schedule. Writers stage their output and close figures on
failure. The CLI adds destination preflight and stages every requested format;
see [publication behavior](cli.md#render-saved-placements).

## Producer-owned duplex frames

A caller can supply `duplex_frame_renderer(document, step_index)` to the media
functions. It returns a nonempty NumPy `uint8` RGB or RGBA image for every
step, with constant image shape within each scene. Frames are validated as they
are requested during rendering, with at most two images cached for transitions.
The callback must support every requested step;
a poster requests the completed first scene. A late invalid frame fails the
staged export and leaves any existing destination file untouched.
The graph geometry, routing, and playback timing remain Dense Arrays-owned.

For a bound callback, its owner may declare `native_nucleotide_cap_height_px`
and `preferred_figure_height_inches` as finite positive values. These metrics
align the graph typography and figure dimensions with the supplied duplex.
Without a callback, the same media renderer draws the sequence placements
natively.

If the producer draws its own distance brackets, its callback owner declares
`renders_distance_brackets=True`. Otherwise the scene renderer draws them
according to `show_distance_bracket`. This capability defaults to `False` and
must be boolean; explicit ownership prevents duplicate span annotations.

The existing DenseGen publication examples use this callback with
`BaseRenderDuplexProjection.render_rgba`. Their producer-owned recipes in
`dnadesign` are:

- `src/dnadesign/densegen/workspaces/demo_dense_array_showcase/playback.yaml`
  for generic overlap packing.
- `src/dnadesign/densegen/workspaces/demo_dense_array_showcase/playback-constraints.yaml`
  for fixed anchors and the RNAP illustration.

Those recipes publish MP4/poster bundles under their workspace's
`outputs/publication/playback/`. They anchor the refined graph-and-duplex
presentation. Dense Arrays does not depend on `dnadesign`; producer translation,
BaseRender frames, biological labels, and recipe publication remain with their
respective owners. Review the exported media at its intended size; the
[product brief](../architecture/animation-product-spec.md) records visual goals.

## Signatures

::: dense_arrays.playback.presentation.PlaybackDocument

::: dense_arrays.playback.theme
    options:
      show_root_heading: false
      members:
        - PlaybackPresentation
        - LegendEntry

::: dense_arrays.playback.matplotlib_renderer
    options:
      show_root_heading: false
      members:
        - render_collection_poster_png
        - render_collection_mp4
        - render_collection_gif
