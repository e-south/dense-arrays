"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/generation/geometry.py

Validate final placement coordinates against the declared assembly transform.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dense_arrays._record_validation import digest, integer, object_fields
from dense_arrays.generation.randomness import PADDING_POLICY
from dense_arrays.planning import GenerationPlan, PlanEvidence
from dense_arrays.realized import RealizedArray


def validate_assembly(
    realized: RealizedArray, plan: GenerationPlan | PlanEvidence
) -> None:
    """Check geometry and policy evidence without drawing or regenerating padding."""
    assembly = plan.request.assembly
    value = realized.provenance.get("assembly")
    if assembly is None:
        if value is not None:
            msg = "unexpected assembly provenance without an assembly policy"
            raise ValueError(msg)
        _validate_coverage(realized, 0, len(realized.sequence))
        return
    keys = {
        "policy",
        "side",
        "padding_length",
        "packed_start",
        "packed_length",
        "trial",
        "stream_id",
    }
    data = object_fields(value, keys, "assembly provenance")
    if set(data) != keys or data["policy"] != PADDING_POLICY:
        msg = "unsupported or incomplete assembly provenance"
        raise ValueError(msg)
    for name in ("padding_length", "packed_start", "packed_length", "trial"):
        integer(
            data[name],
            field_name=f"assembly.{name}",
            minimum=1 if name in {"trial", "packed_length"} else 0,
        )
    digest(data["stream_id"], field_name="assembly.stream_id")
    padding = assembly.padding
    side = None if padding is None else padding.side
    total = len(realized.sequence)
    shift = data["padding_length"] if side == "left" else 0
    if (
        data["side"] != side
        or (padding is None and data["padding_length"] != 0)
        or data["packed_length"] + data["padding_length"] != total
        or data["packed_start"] != shift
        or min(p.start for p in realized.placements) != shift
        or max(p.end for p in realized.placements) != shift + data["packed_length"]
        or data["trial"] > (1 if padding is None else padding.max_trials)
    ):
        msg = "assembly provenance disagrees with policy or final placement geometry"
        raise ValueError(msg)
    _validate_coverage(realized, shift, shift + data["packed_length"])


def _validate_coverage(realized: RealizedArray, start: int, end: int) -> None:
    """Check the interval union without allocating one object per sequence base."""
    cursor = start
    for placement in sorted(realized.placements, key=lambda item: item.start):
        if placement.start < start or placement.start > cursor or placement.end > end:
            break
        cursor = max(cursor, placement.end)
    else:
        if cursor == end:
            return
    msg = "selected placements must continuously cover the packed interval"
    raise ValueError(msg)
