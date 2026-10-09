"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/preparation/requests.py

Preparation request encoding shared by Python and CLI input files.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import os
from dataclasses import asdict, replace
from pathlib import Path

from dense_arrays._record_validation import integer, object_fields, required_text
from dense_arrays.constraints import Length
from dense_arrays.parts import (
    MMR,
    Background,
    CandidateBudget,
    ConditionalLimits,
    Eligibility,
    FimoScoring,
    LengthRange,
    MiningTarget,
    MotifWindow,
    PartFilter,
    PartTable,
    PreparationSet,
    PreparationSpec,
    PWMArtifact,
    PWMExclusion,
    Retention,
    Sampling,
    ScoreBands,
    ScoringLimits,
    Uniqueness,
)
from dense_arrays.parts.serialization import table_from_dict, table_to_dict

PREPARE_SCHEMA = "dense_arrays.prepare.v1"
PREPARE_SET_SCHEMA = "dense_arrays.preparation_set.v1"
PREPARE_WINDOWS_SCHEMA = "dense_arrays.preparation_windows.v1"


def preparation_to_dict(
    request: PreparationSpec | PreparationSet, *, base: Path | None = None
) -> dict[str, object]:
    """Encode a preparation request using the same fields as the CLI."""
    if isinstance(request, PreparationSet):
        return {
            "schema": PREPARE_SET_SCHEMA,
            "sequence_collisions": request.sequence_collisions,
            **(
                {"core_collisions": request.core_collisions}
                if request.core_collisions != "preserve"
                else {}
            ),
            "recipes": [
                {"id": name, "request": preparation_to_dict(item, base=base)}
                for name, item in request.recipes.items()
            ],
        }
    if not isinstance(request.source, PartTable):
        return sampled_to_dict(request, base=base)
    selection = (
        None if request.retain.select is None else request.retain.select.to_dict()
    )
    if selection is not None:
        selection.pop("schema")
    source = table_to_dict(request.source)
    if base is not None:
        source["table"] = os.path.relpath(request.source.table.absolute(), base)
    return {
        "schema": PREPARE_SCHEMA,
        "source": {"kind": "table", **source},
        "retain": {"select": selection},
    }


def preparation_from_dict(
    value: object, *, base: Path | None = None
) -> PreparationSpec | PreparationSet:
    """Parse a declared preparation request without importing optional scorers."""
    if isinstance(value, dict) and value.get("schema") == PREPARE_WINDOWS_SCHEMA:
        return _windows_from_dict(value, base=base)
    if isinstance(value, dict) and value.get("schema") == PREPARE_SET_SCHEMA:
        return _set_from_dict(value, base=base)
    if (
        isinstance(value, dict)
        and isinstance(value.get("source"), dict)
        and value["source"].get("kind") != "table"
    ):
        return sampled_from_dict(value, base=base)
    data = object_fields(value, {"schema", "source", "retain"}, "preparation")
    if data.pop("schema", None) != PREPARE_SCHEMA:
        msg = f"unsupported preparation schema; supported: {PREPARE_SCHEMA}"
        raise ValueError(msg)
    if not isinstance(data.get("source"), dict):
        msg = "preparation.source must be a table source object"
        raise TypeError(msg)
    source = dict(data.pop("source"))
    if source.pop("kind", None) != "table":
        msg = "unsupported preparation source; supported: table"
        raise ValueError(msg)
    source = table_from_dict(source)
    if base is not None and not source.table.is_absolute():
        source = replace(source, table=base / source.table)
    retention = object_fields(data.pop("retain", {}), {"select"}, "retain")
    selected = retention.get("select")
    return PreparationSpec(
        source,
        Retention(
            None if selected is None else PartFilter.from_dict(selected, declared=False)
        ),
    )


def _set_from_dict(value: object, *, base: Path | None) -> PreparationSet:
    """Decode a set of explicit complete preparation recipes."""
    data = object_fields(
        value,
        {"schema", "recipes", "sequence_collisions", "core_collisions"},
        "preparation set",
    )
    if not isinstance(data.get("recipes"), list):
        msg = "preparation set recipes must be an ordered array"
        raise TypeError(msg)
    recipes = {}
    for entry in data["recipes"]:
        item = object_fields(entry, {"id", "request"}, "preparation recipe")
        if set(item) != {"id", "request"} or item["id"] in recipes:
            msg = "preparation recipes require unique IDs and complete requests"
            raise ValueError(msg)
        recipes[item["id"]] = preparation_from_dict(item["request"], base=base)
    return PreparationSet(
        recipes,
        data.get("sequence_collisions", "error"),
        core_collisions=data.get("core_collisions", "preserve"),
    )


def _windows_from_dict(value: object, *, base: Path | None) -> PreparationSet:
    """Expand a bounded window declaration into the ordinary recipe-set contract."""
    required = {"schema", "base", "windows", "candidate_length", "max_recipes"}
    data = object_fields(
        value,
        required | {"sequence_collisions", "core_collisions"},
        "window preparation",
    )
    if not required <= data.keys():
        msg = "incomplete window preparation request"
        raise ValueError(msg)
    integer(data["max_recipes"], field_name="window max_recipes", minimum=1)
    declared = data["windows"]
    if not isinstance(declared, list) or not declared:
        msg = "windows must be a nonempty ordered array"
        raise TypeError(msg)
    if len(declared) > data["max_recipes"]:
        msg = "window count exceeds max_recipes"
        raise ValueError(msg)
    windows = {}
    for raw in declared:
        item = object_fields(raw, {"id", "window"}, "window choice")
        if set(item) != {"id", "window"}:
            msg = "window choices require an ID and window"
            raise ValueError(msg)
        name = required_text(item["id"], field_name="window choice ID")
        if name in windows:
            msg = "window preparation requires unique window IDs"
            raise ValueError(msg)
        windows[name] = MotifWindow.from_dict(item["window"])
    if (
        not isinstance(data["base"], dict)
        or data["base"].get("schema") != PREPARE_SCHEMA
    ):
        msg = "window base requires one complete preparation request"
        raise ValueError(msg)
    return PreparationSet.from_windows(
        preparation_from_dict(data["base"], base=base),
        windows=windows,
        candidate_length=data["candidate_length"],
        max_recipes=data["max_recipes"],
        sequence_collisions=data.get("sequence_collisions", "error"),
        core_collisions=data.get("core_collisions", "preserve"),
    )


def _locator(value: Path | str | None, base: Path | None) -> str | None:
    if value is None:
        return None
    return str(value) if base is None else os.path.relpath(Path(value).absolute(), base)


def _absolute(value: str | None, base: Path | None) -> Path | None:
    if value is None:
        return None
    path = Path(value)
    return base / path if base is not None and not path.is_absolute() else path


def sampled_to_dict(
    request: PreparationSpec, *, base: Path | None = None
) -> dict[str, object]:
    """Write explicit sampled policies and location-relative source references."""
    source = request.source
    encoded = (
        {
            "kind": "background",
            "base_probabilities": list(source.base_probabilities),
            "group": source.group,
        }
        if isinstance(source, Background)
        else {
            "kind": "pwm_artifact",
            "path": _locator(source.path, base),
            "motif_ids": list(source.motif_ids),
            **({"format": source.format} if source.format != "artifact" else {}),
            **({"window": source.window.to_dict()} if source.window else {}),
        }
    )
    scoring = scoring_to_dict(request.scoring, base=base)
    return {
        "schema": PREPARE_SCHEMA,
        "source": encoded,
        "sampling": {
            "strategy": request.sampling.strategy,
            **(
                {"limits": asdict(request.sampling.limits)}
                if request.sampling.limits is not None
                else {}
            ),
            **(
                {"base_probabilities": list(request.sampling.base_probabilities)}
                if request.sampling.base_probabilities is not None
                else {}
            ),
            "length": (
                asdict(request.sampling.length)
                if isinstance(request.sampling.length, LengthRange)
                else {"exact": request.sampling.length.exact}
            ),
        },
        "budget": request.budget.to_dict(),
        **(
            {"mining_target": asdict(request.mining_target)}
            if request.mining_target is not None
            else {}
        ),
        "scoring": scoring,
        **(
            {"score_bands": request.score_bands.to_dict()}
            if request.score_bands is not None
            else {}
        ),
        "eligibility": asdict(request.eligibility),
        "uniqueness": asdict(request.uniqueness),
        "retain": {
            "count": request.retain.count,
            "policy": request.retain.policy,
            "rank_by": request.retain.rank_by,
            **(
                {"mmr": request.retain.mmr.to_dict()}
                if request.retain.mmr is not None
                else {}
            ),
        },
        "screening": [screen_to_dict(rule, base=base) for rule in request.screening],
        "seed": request.seed,
    }


def sampled_from_dict(value: object, *, base: Path | None = None) -> PreparationSpec:
    """Resolve declared sampled fields; unknown options never degrade silently."""
    data = object_fields(
        value,
        {
            "schema",
            "source",
            "sampling",
            "budget",
            "mining_target",
            "scoring",
            "score_bands",
            "eligibility",
            "uniqueness",
            "retain",
            "screening",
            "seed",
        },
        "preparation",
    )
    if data.pop("schema", None) != PREPARE_SCHEMA:
        msg = "unsupported preparation schema"
        raise ValueError(msg)
    source = data.pop("source")
    if source.get("kind") == "background":
        fields = object_fields(
            source, {"kind", "base_probabilities", "group"}, "background source"
        )
        fields.pop("kind")
        resolved = Background(**fields)
    elif source.get("kind") == "pwm_artifact":
        fields = object_fields(
            source, {"kind", "path", "motif_ids", "window", "format"}, "PWM source"
        )
        if fields.get("path") is None:
            msg = "PWM source requires a path"
            raise ValueError(msg)
        resolved = PWMArtifact(
            _absolute(fields["path"], base),
            fields.get("motif_ids", ()),
            MotifWindow.from_dict(fields["window"]) if "window" in fields else None,
            fields.get("format", "artifact"),
        )
    else:
        msg = (
            "unsupported preparation source; supported: table, pwm_artifact, background"
        )
        raise ValueError(msg)
    sampling = object_fields(
        data.pop("sampling", {}),
        {"strategy", "length", "base_probabilities", "limits"},
        "sampling",
    )
    length = object_fields(
        sampling.pop("length", {}), {"exact", "minimum", "maximum"}, "sampling.length"
    )
    length_type = LengthRange if "minimum" in length else Length
    if "limits" in sampling:
        sampling["limits"] = ConditionalLimits(
            **object_fields(
                sampling["limits"],
                {"states", "automaton_states", "mass_bits", "seconds"},
                "sampling.limits",
            )
        )
    sampling = Sampling(length=length_type(**length), **sampling)
    budget = CandidateBudget(
        **object_fields(
            data.pop("budget", {}),
            {"candidates", "seconds", "batch_size", "batch_bases", "total_bases"},
            "budget",
        )
    )
    retained = object_fields(
        data.pop("retain", {}), {"count", "policy", "rank_by", "mmr"}, "retain"
    )
    if "mmr" in retained:
        retained["mmr"] = MMR.from_dict(retained["mmr"])
    retain = Retention(**retained)
    eligibility = Eligibility(
        **object_fields(
            data.pop("eligibility", {}), {"best_hit_score_min_exclusive"}, "eligibility"
        )
    )
    uniqueness = Uniqueness(
        **object_fields(data.pop("uniqueness", {}), {"key"}, "uniqueness")
    )
    scoring = scoring_from_dict(data.pop("scoring", None), base=base)
    if "score_bands" in data:
        data["score_bands"] = ScoreBands.from_dict(data["score_bands"])
    if "mining_target" in data:
        data["mining_target"] = MiningTarget(
            **object_fields(
                data["mining_target"],
                {"eligible_unique", "minimum_candidates", "max_retained_fraction"},
                "mining_target",
            )
        )
    screens = data.pop("screening", [])
    if not isinstance(screens, list):
        msg = "screening must be an ordered array"
        raise TypeError(msg)
    return PreparationSpec(
        resolved,
        retain,
        sampling,
        budget,
        scoring,
        eligibility,
        uniqueness,
        tuple(screen_from_dict(rule, base=base) for rule in screens),
        **data,
    )


def scoring_to_dict(
    scoring: FimoScoring | None, *, base: Path | None = None
) -> dict[str, object] | None:
    """Encode optional scorer settings once for sources and exclusion screens."""
    if scoring is None:
        return None
    result = asdict(scoring)
    result["background"] = _locator(scoring.background, base)
    result["executable"] = _locator(scoring.executable, base)
    result["backend"] = "fimo"
    return result


def scoring_from_dict(
    scoring: object, *, base: Path | None = None
) -> FimoScoring | None:
    """Parse complete named scorer settings with relative locators."""
    if scoring is not None:
        scoring = object_fields(
            scoring,
            {
                "backend",
                "hit_pvalue_max",
                "background",
                "strands",
                "pseudocount",
                "executable",
                "limits",
            },
            "scoring",
        )
        if scoring.pop("backend", None) != "fimo":
            msg = "supported scorer: fimo"
            raise ValueError(msg)
        for name in ("background", "executable"):
            scoring[name] = _absolute(scoring.get(name), base)
        scoring["limits"] = ScoringLimits(
            **object_fields(
                scoring.get("limits", {}),
                {"seconds", "windows", "output_bytes"},
                "scoring.limits",
            )
        )
        scoring = FimoScoring(**scoring)
    return scoring


def screen_to_dict(rule: object, *, base: Path | None = None) -> dict[str, object]:
    """Keep preparation-only PWM rules separate from packing requirements."""
    from dense_arrays.planning.serialization import requirement_to_dict  # noqa: PLC0415

    if not isinstance(rule, PWMExclusion):
        return requirement_to_dict(rule)
    return {
        "id": rule.id,
        "kind": "pwm_exclusion",
        "motifs": [
            {
                "path": _locator(m.path, base),
                "motif_ids": list(m.motif_ids),
                **({"format": m.format} if m.format != "artifact" else {}),
                **({"window": m.window.to_dict()} if m.window else {}),
            }
            for m in rule.motifs
        ],
        "scoring": scoring_to_dict(rule.scoring, base=base),
        "reject": rule.reject,
        "score_field": rule.score_field,
        "threshold": rule.threshold,
    }


def screen_from_dict(value: object, *, base: Path | None = None) -> object:
    """Read sequence requirements or an explicit preparation-only PWM screen."""
    from dense_arrays.planning.serialization import (  # noqa: PLC0415
        requirement_from_dict,
    )

    if not isinstance(value, dict) or value.get("kind") != "pwm_exclusion":
        return requirement_from_dict(value)
    data = object_fields(
        value,
        {"id", "kind", "motifs", "scoring", "reject", "score_field", "threshold"},
        "PWM exclusion",
    )
    data.pop("kind")
    motifs = data.pop("motifs", None)
    if not isinstance(motifs, list):
        msg = "PWM exclusion motifs must be an ordered array"
        raise TypeError(msg)
    sources = []
    for motif in motifs:
        item = object_fields(
            motif, {"path", "motif_ids", "window", "format"}, "screen motif"
        )
        if item.get("path") is None:
            msg = "screen motif requires path"
            raise ValueError(msg)
        sources.append(
            PWMArtifact(
                _absolute(item["path"], base),
                item.get("motif_ids", ()),
                MotifWindow.from_dict(item["window"]) if "window" in item else None,
                item.get("format", "artifact"),
            )
        )
    return PWMExclusion(
        motifs=tuple(sources),
        scoring=scoring_from_dict(data.pop("scoring", None), base=base),
        **data,
    )
