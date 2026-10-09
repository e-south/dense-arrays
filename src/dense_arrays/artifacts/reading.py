"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/reading.py

Work bounds and descriptors for native artifact reads.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dataclasses import dataclass, field

from dense_arrays._record_validation import integer, required_text


class ReadLimitError(RuntimeError):
    """Requested inspection exceeded its explicit work bound."""


@dataclass(frozen=True)
class ReadLimits:
    """Positive work caps, independent of output pagination and generation effort."""

    records: int = 100_000
    pairs: int = 100_000
    identities: int = 100_000

    def __post_init__(self) -> None:
        """Reject booleans, fractions and unlimited sentinel values."""
        for name in ("records", "pairs", "identities"):
            integer(getattr(self, name), field_name=f"read_limits.{name}", minimum=1)

    def to_dict(self) -> dict[str, int]:
        """Expose resolved caps rather than hidden front-end defaults."""
        return {
            name: getattr(self, name) for name in ("records", "pairs", "identities")
        }


@dataclass(frozen=True)
class ReadCost:
    """Describe an immutable record query without opening an iterator."""

    source_id: str
    revision: int
    mode: str
    projection: str
    records_estimate: int | None
    limits: ReadLimits
    bytes_estimate: int | None = None

    def __post_init__(self) -> None:
        """Keep unknown estimates explicit and constrain declared access modes."""
        required_text(self.source_id, field_name="source_id")
        integer(self.revision, field_name="revision", minimum=0)
        required_text(self.projection, field_name="projection")
        if self.mode not in {"manifest", "indexed", "scan"}:
            msg = "unsupported read cost mode"
            raise ValueError(msg)
        if not isinstance(self.limits, ReadLimits):
            msg = "read cost requires ReadLimits"
            raise TypeError(msg)
        for name in ("records_estimate", "bytes_estimate"):
            if getattr(self, name) is not None:
                integer(getattr(self, name), field_name=name, minimum=0)

    def to_dict(self) -> dict[str, object]:
        """Serialize conservative scan estimates and caller-visible caps."""
        return {
            "schema": "dense_arrays.read_cost.v1",
            "source_id": self.source_id,
            "revision": self.revision,
            "mode": self.mode,
            "projection": self.projection,
            "records_estimate": self.records_estimate,
            "estimate_kind": "upper_bound",
            "bytes_estimate": self.bytes_estimate,
            "limits": self.limits.to_dict(),
        }


@dataclass
class ReadBudget:
    """Per-iterator counters; never shared between independent record reads."""

    limits: ReadLimits = field(default_factory=ReadLimits)
    examined: int = 0
    returned: int = 0
    position: int = 0
    offset: int = 0
    identities: int = 0
    bytes_checked: int = 0
    pair_evaluations: int = 0

    def compare(self, count: int) -> None:
        """Admit pair work cumulatively before evaluating comparisons."""
        integer(count, field_name="pair work", minimum=0)
        if self.pair_evaluations + count > self.limits.pairs:
            msg = (
                f"read_limits.pairs={self.limits.pairs} reached; "
                "explicitly increase the comparison allowance"
            )
            raise ReadLimitError(msg)
        self.pair_evaluations += count

    def retain(self, count: int = 1) -> None:
        """Charge persistent lookup entries before allocating additional state."""
        if self.identities + count > self.limits.identities:
            msg = (
                f"read_limits.identities={self.limits.identities} reached; "
                "narrow the population or explicitly increase the state bound"
            )
            raise ReadLimitError(msg)
        self.identities += count

    def examine(self, payload: str | None = None) -> None:
        """Charge one data record before decoding it or evaluating a predicate."""
        if self.examined >= self.limits.records:
            msg = (
                f"read_limits.records={self.limits.records} reached; "
                "narrow the query or explicitly increase the read bound"
            )
            raise ReadLimitError(msg)
        self.examined += 1
        if payload is not None:
            self.bytes_checked += len(payload.encode("utf-8"))


@dataclass(frozen=True)
class Verification:
    """The checked native JSON records, excluding physical database and index bytes."""

    boundary: str
    records_checked: int
    bytes_checked: int

    def to_dict(self) -> dict[str, object]:
        """Expose the evidence boundary without claiming execution replay."""
        return {
            "boundary": self.boundary,
            "records_checked": self.records_checked,
            "bytes_checked": self.bytes_checked,
            "byte_scope": "native_record_json_utf8",
        }
