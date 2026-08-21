"""Import-boundary tests for optional playback visualization dependencies."""

from __future__ import annotations

import subprocess
import sys


def test_playback_contract_import_does_not_require_visualization_dependencies() -> None:
    """Core playback contracts remain usable without Matplotlib or NetworkX."""
    script = """
import builtins

original_import = builtins.__import__

def guarded_import(name, *args, **kwargs):
    if name.partition('.')[0] in {'matplotlib', 'networkx'}:
        raise AssertionError(f'unexpected optional import: {name}')
    return original_import(name, *args, **kwargs)

builtins.__import__ = guarded_import
import dense_arrays.playback
from dense_arrays.playback.html import PlaybackDocument
from dense_arrays.playback.positions import radial_path_positions

assert PlaybackDocument is not None
assert radial_path_positions(1) == ((0.5, 0.5),)
"""
    subprocess.run([sys.executable, "-c", script], check=True)
