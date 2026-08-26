"""Public, renderer-independent dense-array playback contracts."""

from .html import (
    PlaybackDocument,
    render_playback_collection_html,
    render_playback_html,
)
from .models import (
    PLAYBACK_PLAN_SCHEMA_VERSION,
    ConstraintResult,
    CoordinateSpan,
    NoticeLevel,
    OrderingStatus,
    PlaybackAuthority,
    PlaybackNotice,
    PlaybackPlan,
    PlaybackStep,
)
from .reconstruction import reconstruct_playback
from .serialization import (
    dumps_playback_plan,
    dumps_realized_array,
    loads_playback_plan,
    loads_realized_array,
    playback_plan_from_dict,
    playback_plan_to_dict,
    realized_array_from_dict,
    realized_array_to_dict,
)

__all__ = [
    "PLAYBACK_PLAN_SCHEMA_VERSION",
    "ConstraintResult",
    "CoordinateSpan",
    "NoticeLevel",
    "OrderingStatus",
    "PlaybackAuthority",
    "PlaybackDocument",
    "PlaybackNotice",
    "PlaybackPlan",
    "PlaybackStep",
    "dumps_playback_plan",
    "dumps_realized_array",
    "loads_playback_plan",
    "loads_realized_array",
    "playback_plan_from_dict",
    "playback_plan_to_dict",
    "realized_array_from_dict",
    "realized_array_to_dict",
    "reconstruct_playback",
    "render_playback_collection_html",
    "render_playback_html",
]
