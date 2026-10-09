"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/retention/models.py

Explicit MMR policy and per-candidate retention evidence.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dataclasses import asdict, dataclass

from dense_arrays._record_validation import integer, object_fields
from dense_arrays.parts.scoring.configuration import finite

from .pool import PoolSize


@dataclass(frozen=True)
class MMR:
    """Bound greedy core-diversity selection and name its numerical meaning."""

    pool_size: int | PoolSize
    relevance_weight: float
    score_scaling: str
    minimum_fraction_of_max: float | None = None
    distance: str = "pwm_tolerant_hamming"
    tie_break: str = "score_core_sequence"
    algorithm: str = "greedy_mmr.v1"

    def __post_init__(self) -> None:
        """Reject unsupported algorithms, implicit scales and unbounded pools."""
        if not isinstance(self.pool_size, PoolSize):
            integer(self.pool_size, field_name="mmr.pool_size", minimum=1)
        if not 0 < finite(self.relevance_weight, "mmr.relevance_weight") <= 1:
            msg = "mmr.relevance_weight must be in (0, 1]"
            raise ValueError(msg)
        if self.score_scaling not in {"fraction_of_max_clipped", "score_percentile"}:
            msg = (
                "mmr.score_scaling requires fraction_of_max_clipped or score_percentile"
            )
            raise ValueError(msg)
        if (
            self.minimum_fraction_of_max is not None
            and not 0
            < finite(self.minimum_fraction_of_max, "mmr.minimum_fraction_of_max")
            <= 1
        ):
            msg = "mmr.minimum_fraction_of_max must be in (0, 1]"
            raise ValueError(msg)
        if (
            self.distance != "pwm_tolerant_hamming"
            or self.tie_break != "score_core_sequence"
            or self.algorithm != "greedy_mmr.v1"
        ):
            msg = "unsupported MMR distance, tie_break or algorithm"
            raise ValueError(msg)

    def to_dict(self) -> dict[str, object]:
        """Encode every effective policy, including defaults and version."""
        data = asdict(self)
        if isinstance(self.pool_size, PoolSize):
            data["pool_size"] = self.pool_size.to_dict()
        return data

    def pool_limit(self, retained: int) -> int:
        """Resolve the admitted pool bound from one declared size policy."""
        return (
            self.pool_size.resolve(retained)
            if isinstance(self.pool_size, PoolSize)
            else self.pool_size
        )

    @classmethod
    def from_dict(cls, value: object) -> "MMR":
        """Parse only named policy fields; construction validates their meaning."""
        data = object_fields(value, set(cls.__dataclass_fields__), "mmr")
        if isinstance(data.get("pool_size"), dict):
            data["pool_size"] = PoolSize.from_dict(data["pool_size"])
        return cls(**data)


@dataclass(frozen=True)
class MMRDecision:
    """Pool admission and greedy choice evidence; absent utility was not computed."""

    pool_status: str
    relevance: float | None = None
    utility: float | None = None
    nearest_distance: float | None = None
    nearest_similarity: float | None = None

    def __post_init__(self) -> None:
        """Keep bounded pool admission distinct from eligibility and selection."""
        if self.pool_status not in {"included", "below_score", "beyond_limit"}:
            msg = "unsupported MMR pool status"
            raise ValueError(msg)
        for name in ("relevance", "utility", "nearest_distance", "nearest_similarity"):
            value = getattr(self, name)
            if value is not None:
                finite(value, f"mmr decision {name}")
        if self.pool_status == "included":
            if self.relevance is None or not 0 <= self.relevance <= 1:
                msg = "included MMR candidates require relevance in [0, 1]"
                raise ValueError(msg)
        elif any(
            value is not None
            for value in (
                self.relevance,
                self.utility,
                self.nearest_distance,
                self.nearest_similarity,
            )
        ):
            msg = "excluded MMR candidates cannot carry greedy choice evidence"
            raise ValueError(msg)
        if self.utility is None and (
            self.nearest_distance is not None or self.nearest_similarity is not None
        ):
            msg = "MMR distance evidence requires a computed utility"
            raise ValueError(msg)
        if self.nearest_distance is not None and self.nearest_distance < 0:
            msg = "MMR nearest distance must be nonnegative"
            raise ValueError(msg)
        if (
            self.nearest_similarity is not None
            and not 0 <= self.nearest_similarity <= 1
        ):
            msg = "MMR nearest similarity must be in [0, 1]"
            raise ValueError(msg)

    def to_dict(self) -> dict[str, object]:
        """Publish explicit absent values without reconstructing a decision."""
        return asdict(self)

    @classmethod
    def from_dict(cls, value: object) -> "MMRDecision":
        """Require complete evidence fields from a stored decision."""
        keys = set(cls.__dataclass_fields__)
        data = object_fields(value, keys, "MMR decision")
        if set(data) != keys:
            msg = "incomplete MMR decision"
            raise ValueError(msg)
        return cls(**data)
