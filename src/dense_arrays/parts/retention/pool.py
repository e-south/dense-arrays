"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/retention/pool.py

Resolve target-relative retention effort under a mandatory candidate cap.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dataclasses import asdict, dataclass
from fractions import Fraction
from math import ceil

from dense_arrays._record_validation import integer, object_fields
from dense_arrays.parts.scoring.configuration import finite

POOL_SIZE_POLICY = "retained_multiplier.v1"


@dataclass(frozen=True)
class PoolSize:
    """Size an MMR choice pool relative to the requested retained count.

    The multiplier is at least one. The maximum is a hard effort cap; it can
    produce a retained-count shortfall when smaller than the requested count.
    Score admission remains independent and can further reduce the actual pool.
    """

    per_retained: float
    maximum: int

    def __post_init__(self) -> None:
        """Require finite explicit scaling and a positive integer cap."""
        factor = finite(self.per_retained, "pool_size.per_retained")
        if factor < 1:
            msg = "pool_size.per_retained must be at least one"
            raise ValueError(msg)
        integer(self.maximum, field_name="pool_size.maximum", minimum=1)
        object.__setattr__(self, "per_retained", factor)

    def requested(self, retained: int) -> int:
        """Round the decimal multiplier upward before applying the cap."""
        integer(retained, field_name="requested retention", minimum=0)
        return ceil(Fraction(str(self.per_retained)) * retained)

    def resolve(self, retained: int) -> int:
        """Return the maximum admitted pool size for this retained target."""
        return min(self.maximum, self.requested(retained))

    def to_dict(self) -> dict[str, object]:
        """Declare both sizing inputs and their arithmetic policy."""
        return {"policy": POOL_SIZE_POLICY, **asdict(self)}

    @classmethod
    def from_dict(cls, value: object) -> "PoolSize":
        """Reject unknown sizing policies and incomplete declarations."""
        data = object_fields(value, {"policy", "per_retained", "maximum"}, "pool size")
        if data.pop("policy", None) != POOL_SIZE_POLICY:
            msg = "unsupported pool sizing policy"
            raise ValueError(msg)
        return cls(**data)
