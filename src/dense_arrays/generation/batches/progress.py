"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/generation/batches/progress.py

Logical search identity for one attempt within a prepared batch schedule.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dataclasses import dataclass


@dataclass(frozen=True)
class BatchWork:
    """Attempt coordinates supplied by the schedule, independent of solver calls."""

    index: int
    attempt: int
    enumerated: bool
    on_unproven: str = "stop"

    def to_dict(self) -> dict[str, int]:
        """Reserve both ordinals before starting search."""
        return {"batch_index": self.index, "batch_attempt": self.attempt}
