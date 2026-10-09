"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/batches/models.py

Explicit offered-part identities and versioned finite sampling requests.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dataclasses import dataclass, field

from dense_arrays._record_validation import (
    digest,
    integer,
    object_fields,
    required_text,
    semantic_digest,
)

from .feedback import FeedbackSnapshot

BATCH_SCHEMA = "dense_arrays.candidate_batch.v1"
SAMPLING_POLICY = "part_priority_sha256.v1"
CONSTRAINED_SAMPLING_POLICY = "part_flow_priority_sha256.v1"


def validate_unproven_policy(value: str) -> None:
    """Require an explicit action without accepting unproven candidates."""
    if not isinstance(value, str) or value not in {"stop", "next_batch"}:
        msg = "on_unproven must be stop or next_batch"
        raise ValueError(msg)


@dataclass(frozen=True)
class BatchSampling:
    """Sample distinct part identities; group balancing accounts for fixed parts."""

    size: int
    strategy: str = "uniform"
    seed: int = 0
    unique_sequences: bool = field(default=False, kw_only=True)
    unique_cores: bool = field(default=False, kw_only=True)
    max_per_group: int | None = field(default=None, kw_only=True)

    def __post_init__(self) -> None:
        """Reject unbounded sizes, implicit strategies and invalid random seeds."""
        integer(self.size, field_name="batch.size", minimum=1)
        integer(self.seed, field_name="batch.seed", minimum=0)
        if self.strategy not in {"uniform", "group_balanced"}:
            msg = "batch strategy must be uniform or group_balanced"
            raise ValueError(msg)
        if not isinstance(self.unique_sequences, bool) or not isinstance(
            self.unique_cores, bool
        ):
            msg = "batch uniqueness flags must be booleans"
            raise TypeError(msg)
        if self.max_per_group is not None:
            integer(self.max_per_group, field_name="batch.max_per_group", minimum=1)

    @property
    def policy_id(self) -> str:
        """Version constrained selection separately from unconstrained priorities."""
        return (
            CONSTRAINED_SAMPLING_POLICY
            if self.unique_sequences
            or self.unique_cores
            or self.max_per_group is not None
            else SAMPLING_POLICY
        )

    def to_dict(self) -> dict[str, object]:
        """Persist algorithm ownership with every sampling parameter."""
        return {
            "policy": self.policy_id,
            "size": self.size,
            "strategy": self.strategy,
            "seed": self.seed,
            **(
                {
                    "unique_sequences": self.unique_sequences,
                    "unique_cores": self.unique_cores,
                    "max_per_group": self.max_per_group,
                }
                if self.policy_id == CONSTRAINED_SAMPLING_POLICY
                else {}
            ),
        }

    @classmethod
    def from_dict(cls, value: object) -> "BatchSampling":
        """Require the explicit supported algorithm version on saved selections."""
        keys = {"policy", "size", "strategy", "seed"}
        constraints = {"unique_sequences", "unique_cores", "max_per_group"}
        data = object_fields(value, keys | constraints, "batch sampling")
        policy = data.get("policy")
        expected = keys | constraints if policy == CONSTRAINED_SAMPLING_POLICY else keys
        if set(data) != expected or policy not in {
            SAMPLING_POLICY,
            CONSTRAINED_SAMPLING_POLICY,
        }:
            msg = "unsupported or incomplete batch sampling policy"
            raise ValueError(msg)
        data.pop("policy")
        result = cls(**data)
        if result.policy_id != policy:
            msg = "batch sampling policy does not match its settings"
            raise ValueError(msg)
        return result


@dataclass(frozen=True)
class CandidateBatch:
    """An ordered subset of one eligible collection, sufficient for exact replay."""

    part_ids: tuple[str, ...]
    collection_id: str
    sampling: BatchSampling | None = None
    stream: str = "default"
    feedback: FeedbackSnapshot | None = field(default=None, kw_only=True)
    batch_id: str = field(init=False)

    def __post_init__(self) -> None:
        """Freeze ordered identities without redrawing a saved selection."""
        digest(self.collection_id, field_name="batch.collection_id")
        required_text(self.stream, field_name="batch.stream")
        if not isinstance(self.part_ids, (list, tuple)) or not self.part_ids:
            msg = "batch.part_ids must be a nonempty ordered array"
            raise ValueError(msg)
        for part in self.part_ids:
            required_text(part, field_name="batch.part_ids")
        if len(set(self.part_ids)) != len(self.part_ids):
            msg = "batch.part_ids must be unique"
            raise ValueError(msg)
        object.__setattr__(self, "part_ids", tuple(self.part_ids))
        if self.sampling is not None:
            if not isinstance(self.sampling, BatchSampling):
                msg = "batch.sampling must be BatchSampling"
                raise TypeError(msg)
            if self.sampling.size != len(self.part_ids):
                msg = "batch membership does not match sampling size"
                raise ValueError(msg)
        if self.feedback is not None and (
            not isinstance(self.feedback, FeedbackSnapshot) or self.sampling is None
        ):
            msg = "weighted batches require a FeedbackSnapshot and sampling policy"
            raise TypeError(msg)
        object.__setattr__(self, "batch_id", semantic_digest(self.content()))

    def content(self) -> dict[str, object]:
        """Bind selection meaning without source paths or runtime state."""
        return {
            "schema": BATCH_SCHEMA,
            "part_ids": list(self.part_ids),
            "collection_id": self.collection_id,
            "sampling": None if self.sampling is None else self.sampling.to_dict(),
            "stream": self.stream,
            **(
                {"feedback": self.feedback.to_dict()}
                if self.feedback is not None
                else {}
            ),
        }

    def to_dict(self) -> dict[str, object]:
        """Encode a self-checking batch with explicit sampled or supplied provenance."""
        return {**self.content(), "batch_id": self.batch_id}

    @classmethod
    def from_dict(cls, value: object) -> "CandidateBatch":
        """Validate a prepared batch without invoking its sampling algorithm."""
        keys = {"schema", "part_ids", "collection_id", "sampling", "stream", "batch_id"}
        data = object_fields(value, keys | {"feedback"}, "candidate batch")
        if not keys <= data.keys() or data.pop("schema") != BATCH_SCHEMA:
            msg = "unsupported or incomplete candidate batch"
            raise ValueError(msg)
        identity = data.pop("batch_id")
        if data["sampling"] is not None:
            data["sampling"] = BatchSampling.from_dict(data["sampling"])
        if "feedback" in data:
            data["feedback"] = FeedbackSnapshot.from_dict(data["feedback"])
        result = cls(**data)
        if result.batch_id != identity:
            msg = "candidate batch digest mismatch"
            raise ValueError(msg)
        return result
