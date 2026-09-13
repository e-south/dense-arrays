"""Preserve native GIF comments through Matplotlib's Pillow writer.

Module Author(s): Eric J. South
"""

from __future__ import annotations

from matplotlib.animation import PillowWriter


class EvidencePillowWriter(PillowWriter):
    """Pass scene evidence through the existing Pillow save operation."""

    def finish(self) -> None:
        """Attach the native comment without re-encoding the completed GIF."""
        if self._frames:
            self._frames[0].info["comment"] = self.metadata["comment"].encode("utf-8")
        super().finish()
