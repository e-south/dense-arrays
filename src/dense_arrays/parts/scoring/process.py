"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/scoring/process.py

Bound external scoring time and output without collecting unbounded pipes.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import os
import selectors
import subprocess
import time
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

    from .configuration import ScoringLimits


class ScoringError(RuntimeError):
    """An execution failure, distinct from a candidate failing a score cutoff."""

    def __init__(self, reason: str, detail: str) -> None:
        """Preserve the machine-readable reason and bounded backend context."""
        self.reason = reason
        super().__init__(f"FIMO {reason}: {detail}")


@dataclass(frozen=True)
class ToolOutput:
    """Captured standard output and the total bytes consumed on both pipes."""

    stdout: bytes
    total_bytes: int


def remaining_seconds(deadline: float) -> float:
    """Refuse work when a shared scoring deadline has elapsed."""
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        reason = "timeout"
        raise ScoringError(reason, "scoring time limit reached")
    return remaining


def invoke(command: list[str], *, cwd: Path, limits: ScoringLimits) -> ToolOutput:
    """Reap the child on timeout, excess output, cancellation or backend failure."""
    deadline = time.monotonic() + limits.seconds
    try:
        child = subprocess.Popen(  # noqa: S603 - resolved executable and argument list
            command,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except OSError as err:
        msg = "unavailable"
        raise ScoringError(msg, str(err)) from err
    stdout, stderr = bytearray(), bytearray()
    try:
        with selectors.DefaultSelector() as selector:
            selector.register(child.stdout, selectors.EVENT_READ, stdout)
            selector.register(child.stderr, selectors.EVENT_READ, stderr)
            while selector.get_map():
                remaining = remaining_seconds(deadline)
                for key, _ in selector.select(timeout=remaining):
                    chunk = os.read(key.fileobj.fileno(), 65536)
                    if not chunk:
                        selector.unregister(key.fileobj)
                        continue
                    if len(stdout) + len(stderr) + len(chunk) > limits.output_bytes:
                        msg = "output_limit"
                        raise ScoringError(msg, "tool output byte limit reached")
                    key.data.extend(chunk)
            remaining = remaining_seconds(deadline)
            try:
                code = child.wait(timeout=remaining)
            except subprocess.TimeoutExpired as err:
                msg = "timeout"
                raise ScoringError(msg, "scoring time limit reached") from err
        if code:
            detail = stderr[:4096].decode("utf-8", errors="replace").strip()
            msg = "backend"
            raise ScoringError(msg, f"exit {code}: {detail}")
        return ToolOutput(bytes(stdout), len(stdout) + len(stderr))
    finally:
        if child.poll() is None:
            child.kill()
        child.wait()
        child.stdout.close()
        child.stderr.close()
