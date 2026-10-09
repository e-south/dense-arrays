"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/cursors.py

Opaque, versioned continuation tokens bound to one native query snapshot.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import base64
import binascii
import json
from dataclasses import dataclass, replace

from dense_arrays._record_validation import (
    canonical_json,
    digest,
    integer,
    object_fields,
    required_text,
    semantic_digest,
)

CURSOR_SCHEMA = "dense_arrays.cursor.v1"
MAX_CURSOR_CHARACTERS = 4096
_BINDING_FIELDS = 2


@dataclass(frozen=True)
class Cursor:
    """An ordered continuation point, independent of page size and read limits."""

    source_id: str
    revision: int
    query_id: str
    ordinal: int = 0
    offset: int = 0
    bindings: tuple[tuple[str, int], ...] = ()

    def __post_init__(self) -> None:
        """Reject invalid identities or positions without opening an artifact."""
        required_text(self.source_id, field_name="cursor.source_id")
        digest(self.query_id, field_name="cursor.query_id")
        integer(self.revision, field_name="cursor.revision", minimum=0)
        integer(self.ordinal, field_name="cursor.ordinal", minimum=0)
        integer(self.offset, field_name="cursor.offset", minimum=0)
        if not isinstance(self.bindings, (tuple, list)):
            msg = "cursor bindings must be source/revision pairs"
            raise TypeError(msg)
        for binding in self.bindings:
            if (
                not isinstance(binding, (tuple, list))
                or len(binding) != _BINDING_FIELDS
            ):
                msg = "cursor bindings must be source/revision pairs"
                raise ValueError(msg)
            required_text(binding[0], field_name="cursor binding source")
            integer(binding[1], field_name="cursor binding revision", minimum=0)
        object.__setattr__(self, "bindings", tuple(tuple(b) for b in self.bindings))

    def advance(self, ordinal: int, offset: int = 0) -> Cursor:
        """Keep the query and snapshot while advancing to a returned record."""
        return replace(self, ordinal=ordinal, offset=offset)

    def token(self) -> str:
        """Encode an integrity-checked token; this is not an authorization token."""
        payload = {
            "schema": CURSOR_SCHEMA,
            "source_id": self.source_id,
            "revision": self.revision,
            "query_id": self.query_id,
            "ordinal": self.ordinal,
        }
        if self.offset:
            payload["offset"] = self.offset
        if self.bindings:
            payload["bindings"] = [list(b) for b in self.bindings]
        envelope = {"cursor": payload, "digest": semantic_digest(payload)}
        token = (
            base64.urlsafe_b64encode(canonical_json(envelope).encode())
            .decode()
            .rstrip("=")
        )
        if len(token) > MAX_CURSOR_CHARACTERS:
            msg = "cursor exceeds 4096 characters; use fewer sources or all=True"
            raise ValueError(msg)
        return token

    @classmethod
    def from_token(cls, token: str) -> Cursor:
        """Reject unknown, altered or oversized tokens before native reads."""
        if (
            not isinstance(token, str)
            or not token
            or len(token) > MAX_CURSOR_CHARACTERS
        ):
            msg = "cursor must be a nonempty token of at most 4096 characters"
            raise ValueError(msg)
        try:
            value = json.loads(
                base64.b64decode(
                    token + "=" * (-len(token) % 4), altchars=b"-_", validate=True
                )
            )
        except (binascii.Error, ValueError, UnicodeError) as err:
            msg = "malformed cursor token"
            raise ValueError(msg) from err
        envelope = object_fields(value, {"cursor", "digest"}, "cursor envelope")
        fields = {"schema", "source_id", "revision", "query_id", "ordinal"}
        payload = object_fields(
            envelope.get("cursor"), fields | {"offset", "bindings"}, "cursor"
        )
        if (
            not fields <= set(payload)
            or payload["schema"] != CURSOR_SCHEMA
            or envelope.get("digest") != semantic_digest(payload)
        ):
            msg = "unsupported or corrupted cursor schema"
            raise ValueError(msg)
        payload.pop("schema")
        cursor = cls(**payload)
        if cursor.token() != token:
            msg = "cursor requires its canonical encoding"
            raise ValueError(msg)
        return cursor
