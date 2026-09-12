---
title: Command-line reference
description: Select optimization or playback commands and interpret their limits and failures.
---

# Command-line reference

Run these help commands from an [installed source checkout](../quickstart.md):

```bash
uv run dense-arrays optimize --help
uv run dense-arrays solutions --help
uv run dense-arrays-playback --help
```

## Optimize motifs

| Option | Meaning |
| --- | --- |
| `--motif` | One motif; repeat for each library entry |
| `--motifs-file` | One motif per line; replaces `--motif` |
| `--length` | Required positive integer sequence-length limit |
| `--strands` | `single` or `double`; defaults to `double` |
| `--solver` | Backend name passed to OR-Tools; defaults to `CBC` |
| `--max-solutions` | Maximum displayed results for `solutions`; defaults to 10 |
| `--diverse` | Bias `solutions` toward less represented motif entries |

Positional and regulator constraints use the [Python API](../constraints.md).
There is no CLI solve-time limit; `--max-solutions` limits result count only.
Output is a terminal display. It is not the persisted placement JSON expected
by playback.

Bad options and missing motif input produce usage errors. Invalid DNA and a
solve that returns no results exit nonzero. If `solutions` has already printed
a result, a later solve failure can be swallowed and the command can exit zero.
Backend creation failures can currently show a Python traceback;
[solver limitations](optimizer.md#current-solver-limitations) explain why a
solve ending early is not always proof of infeasibility.

## Render saved placements

Supply a realized-array or playback-plan JSON file and the required `--html`
path. `--poster`, `--mp4`, and `--gif` add exports. `--title` labels the HTML
document; `--subtitle` is accepted but is not visibly rendered. Raster exports
do not currently draw either field. See [the export guide](../playback.md#export-a-still-or-video)
for dependencies and working-directory instructions.

Some malformed inputs and missing renderer dependencies currently produce
tracebacks. HTML is written before optional exports, so a failed media export
can leave an HTML file. Output paths are overwritten; use a fresh output
directory and keep inputs separate from outputs.
