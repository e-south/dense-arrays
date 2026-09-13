"""Topology selection and renderer dependency boundaries.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import subprocess
import sys
from typing import TYPE_CHECKING

from dense_arrays.playback import reconstruct_playback
from dense_arrays.playback.graph.isotropic_layout import NetworkXIsotropicLayout
from dense_arrays.playback.graph.layout import build_graph_scene
from dense_arrays.realized import Placement, PlacementKind, RealizedArray

if TYPE_CHECKING:
    import pytest

    from dense_arrays.playback.graph.model import (
        ExplanationGraph,
        GraphLayoutSpec,
        GraphPosition,
    )


def test_default_and_injected_engines_receive_same_pruned_topology(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    plan = reconstruct_playback(
        RealizedArray(
            source_id="test:overlap",
            sequence="AAAAAAAA",
            placements=tuple(
                Placement(str(i), str(i), PlacementKind.TFBS, "AAAA", i)
                for i in range(5)
            ),
        )
    )
    engine = NetworkXIsotropicLayout()
    original = type(engine).layout
    calls = []

    def record(
        self: NetworkXIsotropicLayout,
        graph: ExplanationGraph,
        spec: GraphLayoutSpec,
        *,
        seed: int,
    ) -> tuple[GraphPosition, ...]:
        calls.append(graph)
        return original(self, graph, spec, seed=seed)

    monkeypatch.setattr(type(engine), "layout", record)
    default = build_graph_scene(plan, max_context_edges_per_source=0, seed=17013)
    injected = build_graph_scene(
        plan, max_context_edges_per_source=0, engine=engine, seed=17013
    )
    assert calls[0].context_edges == calls[1].context_edges == ()
    assert default.positions == injected.positions


def test_semantic_graph_projection_imports_without_matplotlib() -> None:
    script = """
import sys
class RejectMatplotlib:
    def find_spec(self, fullname, *args):
        if fullname == "matplotlib" or fullname.startswith("matplotlib."):
            raise AssertionError("semantic projection imported matplotlib")
sys.meta_path.insert(0, RejectMatplotlib())
from dense_arrays.playback.graph.projection import project_explanation_graph
"""
    result = subprocess.run(  # noqa: S603 - isolated import check with fixed code
        [sys.executable, "-c", script], capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, result.stderr
