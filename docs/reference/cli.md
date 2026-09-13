---
title: Command-line reference
description: Choose optimization or playback commands and interpret failures and export behavior.
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

Bad options, unreadable motif files, malformed inputs, no feasible first
result, and solver failures produce errors on stderr and a nonzero exit.
If a failure follows an already printed result, the command still exits
nonzero; preceding output does not imply that enumeration completed.
See [solver outcomes](optimizer.md#solver-outcomes) for the Python exception types.

## Render saved placements

Supply a realized-array or playback-plan JSON file and at least one output:
`--poster` for PNG, `--mp4`, or `--gif`. Multiple formats can be requested in
one command. `--title` and `--subtitle` supply artifact metadata; the
[evidence reference](playback-presentation.md#read-the-evidence) explains where
each format stores it. All formats require the playback extra, and MP4 also
requires FFmpeg. See the
[export guide](../playback.md#export-a-still-or-video) for working-directory
instructions and the [media presentation reference](playback-presentation.md)
for Python settings and producer frame callbacks.

Inputs must pass schema and semantic validation before export. The command
rejects input/output aliases, colliding destinations, symlink output paths,
and existing files unless `--replace` is given. Invalid input, missing media
dependencies, and export failures produce concise errors on stderr and exit
nonzero.

Every requested format is rendered to temporary files before any destination
is published. A rendering failure leaves existing destination files untouched.
Publication then occurs atomically **per file**, not as one filesystem
transaction across all formats. If publication fails partway through, the
error lists the files already published. Successful commands print each
written path.
