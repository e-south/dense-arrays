"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/documents.py

Read strict YAML/JSON objects without interpreting a producer schema.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING, override

import yaml

if TYPE_CHECKING:
    from pathlib import Path


class _StrictLoader(yaml.SafeLoader):
    """Reject duplicate YAML/JSON object fields instead of overwriting them."""

    @override
    def construct_mapping(
        self, node: yaml.MappingNode, deep: bool = False
    ) -> dict[str, object]:
        """Construct one mapping without losing repeated fields."""
        self.flatten_mapping(node)
        result = {}
        for key_node, value_node in node.value:
            key = self.construct_object(key_node, deep=deep)
            if not isinstance(key, str):
                msg = "request object keys must be strings"
                raise TypeError(msg)
            if key in result:
                msg = (
                    f"duplicate request field {key!r} "
                    f"at line {key_node.start_mark.line + 1}"
                )
                raise ValueError(msg)
            result[key] = self.construct_object(value_node, deep=deep)
        return result


def read_document(path: Path) -> object:
    """Read one strict YAML/JSON object without interpreting a filename."""
    try:
        value = yaml.load(path.read_text(encoding="utf-8"), Loader=_StrictLoader)  # noqa: S506 - SafeLoader subclass
    except yaml.YAMLError as err:
        msg = f"{path}: malformed request: {err}"
        raise ValueError(msg) from err
    return value
