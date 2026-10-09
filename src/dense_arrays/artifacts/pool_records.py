"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/pool_records.py

Versioned immutable pool summaries and collection-scoped part records.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dataclasses import dataclass, field

from dense_arrays._record_validation import digest, integer, object_fields
from dense_arrays.artifacts.preparation.records import PoolAccounting
from dense_arrays.artifacts.preparation.sets import SetAccounting, read_accounting
from dense_arrays.artifacts.provenance import Producer
from dense_arrays.artifacts.reading import ReadCost, ReadLimits, Verification
from dense_arrays.parts.models import Part
from dense_arrays.parts.serialization import part_from_dict, part_to_dict

POOL_IDENTITY_SCHEMA = "dense_arrays.pool.v1"
LEGACY_POOL_SCHEMA = "dense_arrays.pool.v1"
POOL_SCHEMA = "dense_arrays.pool.v2"
SAMPLED_POOL_SCHEMA = "dense_arrays.pool.v3"


@dataclass(frozen=True)
class PoolSummary:
    """Small committed preparation accounting, independent of part-table size."""

    pool_id: str
    plan_id: str
    source_parts: int
    retained_parts: int
    verified: bool = False
    revision: int = 0
    state: str = "completed"
    read_limits: ReadLimits = field(
        default_factory=ReadLimits, compare=False, repr=False
    )
    verification: Verification | None = field(default=None, compare=False)
    producer: Producer | None = None
    preparation: PoolAccounting | SetAccounting | None = None

    @property
    def cost(self) -> ReadCost:
        """Describe the single pool manifest read."""
        return ReadCost(self.pool_id, 0, "manifest", "summary", 1, self.read_limits)

    @property
    def verification_cost(self) -> ReadCost:
        """Preview preparation and retained-record verification work."""
        return ReadCost(
            self.pool_id,
            0,
            "scan",
            "verification",
            1 + self.retained_parts + (self.source_parts if self.preparation else 0),
            self.read_limits,
        )

    def __post_init__(self) -> None:
        """Reject impossible counts or unsupported mutable-pool states."""
        if self.producer is not None and not isinstance(self.producer, Producer):
            msg = "producer must be a recorded Producer or null"
            raise TypeError(msg)
        digest(self.pool_id, field_name="pool_id")
        digest(self.plan_id, field_name="plan_id")
        for name in ("source_parts", "retained_parts"):
            integer(getattr(self, name), field_name=name, minimum=0)
        if self.preparation is not None and (
            not isinstance(self.preparation, (PoolAccounting, SetAccounting))
            or self.preparation.counts["processed"] != self.source_parts
            or self.preparation.counts["retained"] != self.retained_parts
            or self.preparation.state != self.state
            or self.producer is None
        ):
            msg = "sampled pool accounting disagrees with its manifest"
            raise ValueError(msg)
        if (
            self.retained_parts > self.source_parts
            or self.revision != 0
            or (self.preparation is None and self.state != "completed")
        ):
            msg = "invalid immutable pool accounting or revision"
            raise ValueError(msg)
        if not isinstance(self.verified, bool):
            msg = "verified must be boolean"
            raise TypeError(msg)

    def to_dict(self) -> dict[str, object]:
        """Publish explicit identity, completion and verification scope."""
        return {
            "schema": SAMPLED_POOL_SCHEMA
            if self.preparation is not None
            else LEGACY_POOL_SCHEMA
            if self.producer is None
            else POOL_SCHEMA,
            **({"preparation": self.preparation.to_dict()} if self.preparation else {}),
            "pool_id": self.pool_id,
            "plan_id": self.plan_id,
            "source_parts": self.source_parts,
            "retained_parts": self.retained_parts,
            "verified": self.verified,
            "revision": self.revision,
            "state": self.state,
            **(
                {"producer": self.producer.to_dict()}
                if self.producer is not None
                else {}
            ),
        }

    @classmethod
    def from_dict(cls, value: object) -> "PoolSummary":
        """Validate the complete immutable manifest without loading source parts."""
        keys = {
            "schema",
            "pool_id",
            "plan_id",
            "source_parts",
            "retained_parts",
            "verified",
            "revision",
            "state",
        }
        data = object_fields(value, keys | {"producer", "preparation"}, "pool summary")
        schema = data.get("schema")
        if schema not in {LEGACY_POOL_SCHEMA, POOL_SCHEMA, SAMPLED_POOL_SCHEMA}:
            msg = (
                f"unsupported pool schema {schema!r}; supported: "
                f"{LEGACY_POOL_SCHEMA}, {POOL_SCHEMA}"
            )
            raise ValueError(msg)
        if schema in {POOL_SCHEMA, SAMPLED_POOL_SCHEMA}:
            keys.add("producer")
        if schema == SAMPLED_POOL_SCHEMA:
            keys.add("preparation")
        if set(data) != keys:
            msg = "unsupported or incomplete pool schema"
            raise ValueError(msg)
        data.pop("schema")
        if "producer" in data:
            data["producer"] = Producer.from_dict(data["producer"])
        if "preparation" in data:
            data["preparation"] = read_accounting(data["preparation"])
        return cls(**data)


@dataclass(frozen=True)
class PoolPart:
    """A retained supplied occurrence qualified by its immutable collection."""

    pool_id: str
    ordinal: int
    part: Part

    def __post_init__(self) -> None:
        """Keep record identity separate from DNA sequence equivalence."""
        digest(self.pool_id, field_name="pool_id")
        integer(self.ordinal, field_name="ordinal", minimum=1)
        if not isinstance(self.part, Part):
            msg = "pool part requires Part"
            raise TypeError(msg)

    @property
    def part_id(self) -> str:
        """Expose the original supplied ID within this pool."""
        return self.part.part_id

    def to_dict(self) -> dict[str, object]:
        """Serialize one normalized occurrence without erasing its collection."""
        return {
            "schema": "dense_arrays.pool_part.v1",
            "pool_id": self.pool_id,
            "ordinal": self.ordinal,
            "part": part_to_dict(self.part),
        }

    @classmethod
    def from_dict(cls, value: object) -> "PoolPart":
        """Reject unknown part-record fields and noncanonical contents."""
        keys = {"schema", "pool_id", "ordinal", "part"}
        data = object_fields(value, keys, "pool part")
        if set(data) != keys or data.pop("schema") != "dense_arrays.pool_part.v1":
            msg = "unsupported or incomplete pool part schema"
            raise ValueError(msg)
        data["part"] = part_from_dict(data["part"])
        result = cls(**data)
        if result.to_dict() != value:
            msg = "pool part requires complete canonical fields"
            raise ValueError(msg)
        return result
