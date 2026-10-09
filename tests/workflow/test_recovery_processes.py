"""Real process interruption and lock handoff preserve native committed evidence.

Author: Eric J. South.
"""

import json
import selectors
import signal
import subprocess
import sys
from pathlib import Path

import pytest

import dense_arrays as da
from dense_arrays.artifacts.recovery import own_run

SCRIPT = """
import os,signal,sys
from pathlib import Path
import dense_arrays as da
from dense_arrays import parts,planning
from dense_arrays.artifacts.store import RunWriter
mode=sys.argv[2]
original=RunWriter._commit
def commit(self,state,**kwargs):
    if mode=='before' and kwargs.get('design') is not None:
        self.connection.create_function('abort_process',0,lambda: os._exit(23))
        self.connection.execute(
            'CREATE TEMP TRIGGER abort_commit BEFORE INSERT ON commits '
            'BEGIN SELECT abort_process(); END')
    original(self,state,**kwargs)
    if kwargs.get('design') is not None:
        if mode=='after': os._exit(24)
        if mode=='interrupt':
            print('ready',flush=True)
            signal.pause()
RunWriter._commit=commit
request=planning.DesignSpec(
    parts=[parts.Part('a','AAA'),parts.Part('b','CCC'),parts.Part('c','GGG')],
    length=planning.Length(maximum=6),strands='single',target=planning.Target(count=3))
try:
    da.run(request,out=Path(sys.argv[1]))
except KeyboardInterrupt:
    sys.exit(130)
"""


@pytest.mark.parametrize("matrix_run", [False, True])
def test_signal_interruption_and_lock_handoff(tmp_path: Path, matrix_run: bool):
    path = tmp_path / "run"
    script = matrix_script() if matrix_run else SCRIPT
    child = subprocess.Popen(  # noqa: S603 - fixed local script and isolated destination
        [sys.executable, "-c", script, str(path), "interrupt"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    try:
        with selectors.DefaultSelector() as selector:
            selector.register(child.stdout, selectors.EVENT_READ)
            assert selector.select(timeout=20), "child never committed a design"
        assert child.stdout.readline().strip() == "ready"
        inode = (path / ".writer.lock").stat().st_ino
        busy = subprocess.run(  # noqa: S603 - fixed local CLI, no shell
            [
                sys.executable,
                "-m",
                "dense_arrays.cli",
                "run",
                "--resume",
                str(path),
                "--json",
            ],
            capture_output=True,
            text=True,
            timeout=20,
            check=False,
        )
        assert busy.returncode == 4, busy.stderr
        assert json.loads(busy.stdout)["code"] == "writer_busy"
        child.send_signal(signal.SIGINT)
        _, stderr = child.communicate(timeout=20)
        assert child.returncode == 130, stderr
        before = da.inspect(path, verify=True)
        assert before.accepted == 1
        assert before.resumable
        with own_run(path), pytest.raises(ValueError, match="active writer"):
            da.run(resume=path)
        result = da.run(resume=path)
        assert da.inspect(result, verify=True).accepted == (4 if matrix_run else 3)
        assert (path / ".writer.lock").stat().st_ino == inode
    finally:
        if child.poll() is None:
            child.kill()
        child.wait(timeout=10)


@pytest.mark.parametrize(
    "mode,expected", [("before", (0, 1, 23)), ("after", (1, 1, 24))]
)
@pytest.mark.parametrize("matrix_run", [False, True])
def test_abrupt_exit_preserves_prefix_and_refuses_unknown_active_time(
    tmp_path: Path,
    mode: str,
    expected: tuple[int, int, int],
    matrix_run: bool,
):
    path = tmp_path / "run"
    accepted, started, exit_code = expected
    child = subprocess.run(  # noqa: S603 - fixed fault-injection process
        [
            sys.executable,
            "-c",
            matrix_script() if matrix_run else SCRIPT,
            str(path),
            mode,
        ],
        capture_output=True,
        text=True,
        timeout=20,
        check=False,
    )
    assert child.returncode == exit_code, child.stderr
    # A normal SQLite writer open recovers an interrupted native transaction.
    with own_run(path) as connection:
        connection.execute("SELECT count(*) FROM commits").fetchone()
    before = da.inspect(path, verify=True)
    assert before.accepted == accepted
    assert before.counts["started"] == started
    assert before.counts["in_progress"] == (1 if mode == "before" else 0)
    assert before.resumable is False
    saved = (path / "run.sqlite3").read_bytes()
    with pytest.raises(ValueError, match="unknown active time"):
        da.run(resume=path)
    assert (path / "run.sqlite3").read_bytes() == saved
    assert da.inspect(path, verify=True).to_dict() == before.to_dict()


def matrix_script() -> str:
    """Exercise the same real process boundaries with two independently active cells."""
    return SCRIPT.replace(
        "try:\n    da.run",
        """request=planning.MatrixSpec(
    base=request.with_changes(target=planning.Target(count=1)),
    axes={'x': {'a': planning.Variant(), 'b': planning.Variant()}},
    allocation=planning.Allocation(per_cell=2), max_cells=2)
try:
    da.run""",
    )
