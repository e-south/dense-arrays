---
title: Write documentation
description: Give each page one reader task, verified examples, and a clear route to detail.
---

# Write documentation

Help a reader complete a task or make a decision. For Dense Arrays, that can
mean fitting a motif library into a length limit, interpreting returned
offsets, or checking whether a playback preserves the supplied evidence.
State that outcome before explaining the implementation.

## Choose the page that owns the information

| Reader need | Owner |
| --- | --- |
| Decide whether the tool fits the task | README and documentation index |
| Complete a first successful example | First-array tutorial |
| Add requirements or render saved data | Constraints and playback guides |
| Look up inputs, outputs, defaults, and failures | One reference page per interface |
| Understand the formulation or ownership | Method and architecture pages |
| Change and verify the package | Code map and development guide |

Keep one home for each contract. Link to it instead of copying its full text.
Split a page when its reader task or owner changes, not when it reaches an
arbitrary number of lines. Check the rendered page too: a short mkdocstrings
directive can expand into a large reference.

## Write the page

Use `title` and `description` front matter for docs-site metadata, followed by
one descriptive H1 and a short opening that says what the reader can do.
Keep prerequisites beside the first command. State its working directory,
expected output, and relevant failure behavior. Put optional detail after the
successful path and link to the next task.

Use concrete verbs and the package's actual names. Replace phrases such as
“the public seam” with “the data flow” or name the inputs and outputs.
Explain a limitation where it affects a decision. Avoid repeated claims that
the documentation is canonical, modern, robust, or easy to use; demonstrate
those properties through a clear route and a working example.

Preserve domain terms, uncertainty, citations, and all existing authorship.
Keep Virgile Andreani's credits wherever they occur. Attribute newly authored
material to Eric J. South, without rewriting joint credits or adding
machine-generated authors. Source-derived technical claims must stay within
what the implementation and tests establish.

## Dogfood the route

Give a fresh reader or agent a task, not the answer's file path:

1. “Create the tutorial's 40-base array from four 16-base motifs using CBC,
   explain an overlap and the returned offsets, and find the associated paper.”
2. “Require two entries from one motif group and one from another group.”
3. “Render saved placements as a PNG, find the MP4 export command, and
   determine whether the order was solver-recorded.”
4. From `src/dense_arrays/playback/`: “Find the owner and tests for rejecting a
   malformed serialized placement before rendering.”

Record the actual pages read, commands attempted, missing prerequisites, and
the proposed edit/test scope. A route passes when the reader finds the task,
owner, relevant limitations, and verification command without guessing an
undocumented dependency. Distinguish a delegated reader check from a fresh
installation or an independently launched runtime.

Run changed examples, build docs strictly, inspect internal links and anchors,
and check the built site at narrow and wide widths. Review plain Markdown as
well as the site. Keep the banner's accessible description; use text labels
alongside color. Complete the [development gate](../development.md).

## Maintain the teaching media

The playback guide maintains one MP4, its opening PNG, and its completed
poster in `docs/assets/`. Regenerate them from the guide's Python example
using `render_collection_mp4((document,), output / "playback.mp4")` with the
default timing. Use the guide's `poster.png` for the completed still, and
decode the first video frame for the opening image:

```bash
ffmpeg -i playback.mp4 -frames:v 1 playback-opening.png
```

Review the encoded opening, middle, and final frames before replacing the
three assets. Keep transient producer runs outside tracked source. The
editable process figure is `docs/assets/motif-packing-process.svg`; its
graph uses start-to-start shifts and a final-motif cost, while playback
labels count newly covered bases.

The reader/task separation follows [Diátaxis](https://diataxis.fr/start-here/).
The focus on reader benefit is informed by McEnerney's
[writing workshop on value](https://calendar.fsu.edu/event/invited-speaker-scholarly-writing-meeting-the-demand-for-value-graduate-student-workshop-gsrc).
Here, the practical measure is whether a reader can act correctly from the page.
