"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/provenance.py

Recorded execution environment; readers never rediscover a producer's versions.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import platform
from dataclasses import asdict, dataclass
from importlib.metadata import PackageNotFoundError, version

from dense_arrays._record_validation import object_fields, required_text
from dense_arrays.solver import SolverIdentity

PRODUCER_SCHEMA = "dense_arrays.producer.v1"


@dataclass(frozen=True)
class Producer:
    """Software/runtime observations, separate from part and plan identities.

    Missing distribution metadata stays null. A solver identity is supplied only
    after a model was built; preparation and pre-build failures never infer one.
    Versions describe the runtime and do not attest an unpublished source tree.
    """

    package_version: str | None
    python_implementation: str
    python_version: str
    platform_system: str
    platform_machine: str
    ortools_version: str | None
    solver: SolverIdentity | None = None

    def __post_init__(self) -> None:
        """Freeze explicit observations with no hostname, source path or timestamp."""
        for name in ("package_version", "ortools_version"):
            if (value := getattr(self, name)) is not None:
                required_text(value, field_name=f"producer.{name}")
        for name in (
            "python_implementation",
            "python_version",
            "platform_system",
            "platform_machine",
        ):
            required_text(getattr(self, name), field_name=f"producer.{name}")
        if self.solver is not None and not isinstance(self.solver, SolverIdentity):
            msg = "producer.solver must be SolverIdentity or null"
            raise TypeError(msg)

    @classmethod
    def capture(cls) -> Producer:
        """Capture local version metadata without loading a backend or invoking Git."""
        return cls(
            _version("dense-arrays"),
            platform.python_implementation(),
            platform.python_version(),
            platform.system(),
            platform.machine(),
            _version("ortools"),
        )

    def to_dict(self) -> dict[str, object]:
        """Serialize a stable observation; absent values remain explicit nulls."""
        return {"schema": PRODUCER_SCHEMA, **asdict(self)}

    @classmethod
    def from_dict(cls, value: object) -> Producer:
        """Validate recorded evidence without querying the inspecting environment."""
        keys = {
            "schema",
            "package_version",
            "python_implementation",
            "python_version",
            "platform_system",
            "platform_machine",
            "ortools_version",
            "solver",
        }
        data = object_fields(value, keys, "producer")
        if set(data) != keys or data.pop("schema") != PRODUCER_SCHEMA:
            msg = (
                "unsupported or incomplete producer schema; "
                f"supported: {PRODUCER_SCHEMA}"
            )
            raise ValueError(msg)
        if data["solver"] is not None:
            data["solver"] = SolverIdentity.from_dict(data["solver"])
        return cls(**data)


def _version(distribution: str) -> str | None:
    try:
        return version(distribution)
    except PackageNotFoundError:
        return None
