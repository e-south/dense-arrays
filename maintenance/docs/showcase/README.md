# Rebuild the documentation assets

This recipe renders the same selected 100-bp and 60-bp records used in the
published demonstrations. It runs the public BaseRender and Dense Arrays
interfaces in a locked environment. It has no dependency on the dissertation
workspace or a full upstream collection.

From the repository root, with Python 3.12, uv, Arial, Inkscape and FFmpeg
available, choose a fresh output directory.

```bash
uv run --locked --project maintenance/docs/showcase \
  python maintenance/docs/showcase/render.py \
  --output build/showcase-candidate
```

The 100-bp movie uses reconstructed placement order. It does not display solver
search history. Traversal and matching color reveals are 50% faster than the
previous reviewed sequence. Holds and the 60-bp movie keep their timing.
The End edge completes with the final part and has no separate animation.

| Change | File |
| --- | --- |
| Selected realized records and drawing records | `inputs/records.json` |
| Origins and coordinate-recovery qualifications | `inputs/origins.json` |
| Fonts, palette, sequence geometry, graph style and timing | `rendering.json` |
| Workflow artwork, centered captions and arrows | `inputs/workflow-context.svg` |
| Producer commits and dependencies | `pyproject.toml` and `uv.lock` |

Drawing records are checked against the realized sequences and placements before
rendering. New inputs must preserve identities and supply both representations.
The workflow SVG is a retained closed subset with no editor paths or external
resources. The sequence in that workflow is rendered from the selected record
rather than copied from the illustration. Its video uses the same geometry.
Only empty video margins are clipped, with eight pixels of additional clearance.

The README plays the full packing loop inline on GitHub. Its `picture` selects
`packing-loop.webp` when the reader permits motion and uses `packing-card.png`
as a static fallback. PyPI's renderer drops the motion source and retains the
PNG. Both delivery images use immutable Git commit URLs, so they also load from
feature branches. No Play label or external player link is added.

The animation is encoded from the same MP4, at the same pace, within the shared
rounded artwork frame. It has a 1 MiB delivery budget. The documentation keeps
SVGs and native video. Commit new delivery assets first, then select their
commit in the README. Retain that asset commit in the merged history. Verify
actual GitHub playback, the public asset responses and the packaged long
description before merging. Do not select a new asset on `main` prematurely.

Review the candidate in the built documentation before copying its files into
`docs/assets/showcase/`. Keep the player JavaScript there. The generated CSS
owns media placement. `receipt.json` records input and code hashes, producer
commits, record IDs and output hashes. Changed inputs during a render fail
before the receipt is written.

BaseRender owns sequence artwork and native output. Dense Arrays owns the graph
and playback schedule. This recipe owns the selected examples and documentation
composition. It does not add a runtime dependency between the two packages.

The secret-scan baseline records reviewed DNA strings and provenance hashes.
When inputs or receipts change, inspect each new detection before updating that
baseline. Do not exclude the fixture or receipt files from scanning.
