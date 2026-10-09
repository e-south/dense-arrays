"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/selections/membership.py

Check saved identity, content and ordering without evaluating a selection policy.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

from dense_arrays._record_validation import semantic_digest

if TYPE_CHECKING:
    from dense_arrays.artifacts.reading import ReadBudget
    from dense_arrays.artifacts.records import Design

    from .snapshots import SelectionSnapshot


class Membership:
    """A bounded lookup shared by record projection and quality accumulation."""

    def __init__(self, snapshot: SelectionSnapshot, budget: ReadBudget) -> None:
        """Charge retained membership before indexing full design references."""
        budget.retain(snapshot.selected)
        self.members = {
            m.reference: (index, m.record_digest)
            for index, m in enumerate(snapshot.members)
        }
        self.selected = 0

    def matches(self, design: Design) -> bool:
        """Validate the next saved member before exposing or accumulating it."""
        if design.reference not in self.members:
            return False
        index, fingerprint = self.members[design.reference]
        if index != self.selected or fingerprint != semantic_digest(design.to_dict()):
            msg = "saved selection design content or order changed"
            raise ValueError(msg)
        self.selected += 1
        return True

    def finish(self) -> None:
        """Reject lost membership after a complete source scan."""
        if self.selected != len(self.members):
            msg = "saved selection members are missing from the pinned source revisions"
            raise ValueError(msg)
