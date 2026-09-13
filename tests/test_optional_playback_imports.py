"""Import-boundary tests for optional playback visualization dependencies."""

from __future__ import annotations

import subprocess
import sys


def test_playback_contract_import_does_not_require_visualization_dependencies() -> None:
    """Contracts and semantic projection need no solver or raster dependencies."""
    script = """
import builtins

original_import = builtins.__import__

def guarded_import(name, *args, **kwargs):
    if name.partition('.')[0] in {'matplotlib', 'networkx', 'ortools'}:
        raise AssertionError(f'unexpected optional import: {name}')
    return original_import(name, *args, **kwargs)

builtins.__import__ = guarded_import
import dense_arrays.playback
from dense_arrays.playback import PlaybackDocument
from dense_arrays.playback.presentation import PlaybackDocument as NeutralDocument
from dense_arrays.playback.graph.projection import project_explanation_graph
from dense_arrays.playback.reconstruction import reconstruct_playback
from dense_arrays.playback.serialization import loads_playback_plan
from dense_arrays.playback.positions import radial_path_positions
from dense_arrays.realized import RealizedArray

assert PlaybackDocument is not None
assert PlaybackDocument is NeutralDocument
assert callable(project_explanation_graph)
assert callable(reconstruct_playback)
assert callable(loads_playback_plan)
assert RealizedArray is not None
assert radial_path_positions(1) == ((0.5, 0.5),)
"""
    subprocess.run([sys.executable, "-c", script], check=True)
