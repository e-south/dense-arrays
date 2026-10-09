"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/preparation/sets.py

Named independent recipes for one prepared collection.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from types import MappingProxyType

from dense_arrays._record_validation import integer, required_text
from dense_arrays.constraints import Length
from dense_arrays.parts.models import PartTable
from dense_arrays.parts.motifs.artifacts import PWMArtifact
from dense_arrays.parts.motifs.windows import MotifWindow
from dense_arrays.parts.preparation.requests import PreparationSpec


@dataclass(frozen=True)
class PreparationSet:
    """Prepare named sampled recipes with independent budgets and selection.

    Mapping order declares publication order. Each recipe retains its own seed;
    names do not implicitly seed or relabel a biological group.
    ``core_collisions='error'`` rejects equal observed, motif-oriented cores
    across recipes. Parts without a core are excluded from that comparison.
    """

    recipes: Mapping[str, PreparationSpec]
    sequence_collisions: str = "error"
    core_collisions: str = field(default="preserve", kw_only=True)

    def __post_init__(self) -> None:
        """Freeze a nonempty collection of complete, non-nested sampled requests."""
        for name in ("sequence_collisions", "core_collisions"):
            if getattr(self, name) not in ("error", "preserve"):
                msg = f"{name} requires error or preserve"
                raise ValueError(msg)
        if not isinstance(self.recipes, Mapping) or not self.recipes:
            msg = "preparation set requires a nonempty recipe mapping"
            raise ValueError(msg)
        for name, request in self.recipes.items():
            required_text(name, field_name="recipe ID")
            if not isinstance(request, PreparationSpec) or isinstance(
                request.source, PartTable
            ):
                msg = "preparation set recipes must be sampled PreparationSpec values"
                raise TypeError(msg)
        object.__setattr__(self, "recipes", MappingProxyType(dict(self.recipes)))

    def with_changes(self, **changes: object) -> "PreparationSet":
        """Validate an immutable recipe-set edit without running preparation."""
        return replace(self, **changes)

    @classmethod
    def from_windows(  # noqa: PLR0913 - explicit geometry, bounds and collision choices
        cls,
        base: PreparationSpec,
        *,
        windows: Mapping[str, MotifWindow],
        candidate_length: str,
        max_recipes: int,
        sequence_collisions: str = "error",
        core_collisions: str = "preserve",
    ) -> "PreparationSet":
        """Expand named source windows into independently scored recipes.

        ``candidate_length='base'`` preserves the declared sampling length;
        ``'window'`` makes each candidate exactly its selected window's width.
        Budgets, retained targets and seeds are copied per recipe. No files are
        opened, and source coordinates always refer to the unwindowed model.
        """
        integer(max_recipes, field_name="window max_recipes", minimum=1)
        if not isinstance(windows, Mapping) or not windows:
            msg = "windows requires a nonempty ordered mapping"
            raise ValueError(msg)
        if len(windows) > max_recipes:
            msg = "window count exceeds max_recipes"
            raise ValueError(msg)
        if not isinstance(base, PreparationSpec) or not isinstance(
            base.source, PWMArtifact
        ):
            msg = "window expansion requires a PWM PreparationSpec"
            raise TypeError(msg)
        if base.source.window is not None:
            msg = "window expansion requires an unwindowed source motif"
            raise ValueError(msg)
        if candidate_length not in ("base", "window"):
            msg = "candidate_length requires base or window"
            raise ValueError(msg)
        recipes = {}
        for name, window in windows.items():
            if not isinstance(window, MotifWindow):
                msg = "named windows must be MotifWindow values"
                raise TypeError(msg)
            if (
                candidate_length == "base"
                and base.sampling.minimum_length < window.length
            ):
                msg = "base candidate length must accommodate every declared window"
                raise ValueError(msg)
            sampling = (
                base.sampling
                if candidate_length == "base"
                else replace(base.sampling, length=Length(exact=window.length))
            )
            recipes[name] = replace(
                base, source=replace(base.source, window=window), sampling=sampling
            )
        return cls(recipes, sequence_collisions, core_collisions=core_collisions)
