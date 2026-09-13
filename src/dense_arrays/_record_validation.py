"""Validate scalar fields and immutable JSON snapshots for persisted records.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from enum import StrEnum
from types import MappingProxyType

_IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")
_SHA256_LENGTH = 64


def required_text(value: object, *, field_name: str) -> str:
    """Return nonblank text without coercing identities."""
    if not isinstance(value, str):
        msg = f"{field_name} must be a non-empty string"
        raise TypeError(msg)
    if not value.strip():
        msg = f"{field_name} must be a non-empty string"
        raise ValueError(msg)
    return value


def normalized_dna(value: object, *, field_name: str) -> str:
    """Return uppercase DNA after checking the IUPAC alphabet."""
    sequence = required_text(value, field_name=field_name).strip().upper()
    invalid = sorted(set(sequence) - _IUPAC_DNA)
    if invalid:
        msg = f"{field_name} contains non-IUPAC DNA symbols: {invalid}"
        raise ValueError(msg)
    return sequence


def integer(value: object, *, field_name: str, minimum: int | None = None) -> int:
    """Return an integer within its declared domain, excluding booleans."""
    if isinstance(value, bool) or not isinstance(value, int):
        msg = f"{field_name} must be an integer"
        raise TypeError(msg)
    if minimum is not None and value < minimum:
        msg = f"{field_name} must be an integer >= {minimum}"
        raise ValueError(msg)
    return value


def enum_value[T: StrEnum](value: object, enum: type[T], *, field_name: str) -> T:
    """Return a supported enum value from an enum or exact string."""
    required_text(value, field_name=field_name)
    return enum(value)


def digest(value: object, *, field_name: str) -> str:
    """Return a canonical SHA-256 hexadecimal digest."""
    text = required_text(value, field_name=field_name).strip().lower()
    if len(text) != _SHA256_LENGTH or any(
        char not in "0123456789abcdef" for char in text
    ):
        msg = f"{field_name} must be a SHA-256 hex digest"
        raise ValueError(msg)
    return text


def records[T](
    value: object, record_type: type[T], *, field_name: str
) -> tuple[T, ...]:
    """Freeze a sequence of records after checking every member's type."""
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        msg = f"{field_name} must be a sequence of {record_type.__name__} records"
        raise TypeError(msg)
    result = tuple(value)
    if any(not isinstance(item, record_type) for item in result):
        msg = f"{field_name} must contain only {record_type.__name__} records"
        raise TypeError(msg)
    return result


def _freeze_json(value: object, ancestors: frozenset[int]) -> object:
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if math.isfinite(value):
            return value
        msg = "JSON provenance numbers must be finite"
        raise ValueError(msg)
    if id(value) in ancestors:
        msg = "JSON provenance must not contain cyclic containers"
        raise ValueError(msg)
    parents = ancestors | {id(value)}
    if isinstance(value, Mapping):
        if any(not isinstance(key, str) for key in value):
            msg = "JSON provenance object keys must be strings"
            raise TypeError(msg)
        return MappingProxyType(
            {key: _freeze_json(item, parents) for key, item in value.items()}
        )
    if isinstance(value, (list, tuple)):
        return tuple(_freeze_json(item, parents) for item in value)
    msg = "provenance and metadata must contain only JSON values"
    raise TypeError(msg)


def immutable_json_mapping(value: object) -> Mapping[str, object]:
    """Return a recursively immutable, detached JSON object snapshot."""
    if not isinstance(value, Mapping):
        msg = "provenance and metadata must be JSON objects"
        raise TypeError(msg)
    return _freeze_json(value, frozenset())


def mutable_json(value: object) -> object:
    """Return independent JSON dictionaries and arrays from a frozen snapshot."""
    if isinstance(value, Mapping):
        return {key: mutable_json(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [mutable_json(item) for item in value]
    return value


def validate_placement_sequence(
    *, placement_id: str, start: int, end: int, sequence: str, realized_sequence: str
) -> None:
    """Check placement bounds and exact alignment with the realized sequence."""
    if end > len(realized_sequence):
        msg = (
            f"placement {placement_id!r} ends at {end}, "
            f"beyond sequence length {len(realized_sequence)}"
        )
        raise ValueError(msg)
    observed = realized_sequence[start:end]
    if observed != sequence:
        msg = (
            f"placement {placement_id!r} is sequence-inconsistent: "
            f"expected {sequence!r}, observed {observed!r}"
        )
        raise ValueError(msg)
