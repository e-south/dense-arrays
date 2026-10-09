"""Separate processes cannot share a destination or lose a reserved attempt.

Author: Eric J. South.
"""

import subprocess
import sys
from pathlib import Path

import dense_arrays as da


def test_competing_creators_have_one_owner(tmp_path: Path):
    out = tmp_path / "run"
    command = [
        sys.executable,
        "-m",
        "dense_arrays.cli",
        "run",
        "--motif",
        "AAA",
        "--length",
        "3",
        "--out",
        str(out),
        "--json",
    ]
    children = [
        subprocess.Popen(  # noqa: S603 - fixed local command, no shell
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        for _ in range(2)
    ]
    try:
        results = [child.communicate(timeout=20) for child in children]
        assert sorted(child.returncode for child in children) == [0, 2], results
        assert da.inspect(out, verify=True).accepted == 1
        assert (out / ".writer.lock").is_file()
    finally:
        for child in children:
            if child.poll() is None:
                child.kill()
            child.wait(timeout=10)


def test_process_exit_preserves_inspectable_reserved_attempt(tmp_path: Path):
    script = """
import os, sys
from pathlib import Path
import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.artifacts.store import create_run
plan = da.plan(planning.DesignSpec(
    parts=[parts.Part('a', 'AAA')], length=planning.Length(maximum=3)))
with create_run(plan, Path(sys.argv[1])) as writer:
    writer.reserve(active_seconds=0.1)
    os._exit(23)
"""
    out = tmp_path / "interrupted"
    result = subprocess.run(  # noqa: S603 - fixed local script, no shell
        [sys.executable, "-c", script, str(out)],
        capture_output=True,
        text=True,
        timeout=20,
        check=False,
    )
    assert result.returncode == 23, result.stderr
    report = da.inspect(out, verify=True)
    assert report.counts["started"] == report.counts["in_progress"] == 1
    assert report.accepted == 0
    assert report.resumable is False
