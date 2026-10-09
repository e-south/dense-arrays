"""Named windows expand to ordinary independently scored preparation recipes.

Author: Eric J. South.
"""

import json
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.cli import app
from dense_arrays.planning.preparation.requests import preparation_to_dict
from dense_arrays.workflow.inputs import read_source

from .test_fimo import HEADER, executable
from .test_windows import wide_source


def base_recipe() -> parts.PreparationSpec:
    return parts.PreparationSpec(
        parts.PWMArtifact("absent.json"),
        parts.Retention(count=2, policy="top_score", rank_by="best_hit_score"),
        sampling=parts.Sampling(parts.LengthRange(6, 8)),
        budget=parts.CandidateBudget(50),
        scoring=parts.FimoScoring(),
        seed=7,
    )


def test_named_windows_expand_without_io_and_preserve_unselected_policies():
    base = base_recipe()
    windows = {"short": parts.MotifWindow(3), "long": parts.MotifWindow(5, "uniform")}
    result = parts.PreparationSet.from_windows(
        base,
        windows=windows,
        candidate_length="base",
        max_recipes=2,
        sequence_collisions="preserve",
    )
    expected = parts.PreparationSet(
        {
            "short": replace(
                base, source=replace(base.source, window=windows["short"])
            ),
            "long": replace(base, source=replace(base.source, window=windows["long"])),
        },
        sequence_collisions="preserve",
    )
    assert result == expected
    assert tuple(result.recipes) == ("short", "long")
    assert base.source.window is None
    assert all(
        r.sampling.length == parts.LengthRange(6, 8) for r in result.recipes.values()
    )
    assert all(
        r.seed == 7 and r.budget.candidates == 50 and r.retain.count == 2
        for r in result.recipes.values()
    )


def test_matching_candidate_lengths_keeps_each_core_model_uniform():
    result = parts.PreparationSet.from_windows(
        base_recipe(),
        windows={"short": parts.MotifWindow(3), "long": parts.MotifWindow(5)},
        candidate_length="window",
        max_recipes=2,
    )
    assert [r.sampling.length for r in result.recipes.values()] == [
        planning.Length(exact=3),
        planning.Length(exact=5),
    ]


@pytest.mark.parametrize("maximum", [0, True, 1])
def test_window_expansion_is_bounded_before_source_resolution(maximum: object):
    with pytest.raises((ValueError, TypeError), match="max_recipes"):
        parts.PreparationSet.from_windows(
            base_recipe(),
            windows={"a": parts.MotifWindow(3), "b": parts.MotifWindow(5)},
            candidate_length="window",
            max_recipes=maximum,
        )


@pytest.mark.parametrize("mode", ["automatic", None, True])
def test_candidate_length_mode_requires_explicit_meaning(mode: object):
    with pytest.raises((ValueError, TypeError), match="candidate_length"):
        parts.PreparationSet.from_windows(
            base_recipe(),
            windows={"a": parts.MotifWindow(3)},
            candidate_length=mode,
            max_recipes=1,
        )


def test_window_expansion_rejects_nested_source_transformations():
    base = base_recipe().with_changes(
        source=parts.PWMArtifact("absent.json", window=parts.MotifWindow(5))
    )
    with pytest.raises(ValueError, match="unwindowed"):
        parts.PreparationSet.from_windows(
            base,
            windows={"a": parts.MotifWindow(3)},
            candidate_length="window",
            max_recipes=1,
        )


def test_window_request_cli_and_manual_set_share_plans_and_model_bound_records(
    tmp_path: Path,
):

    tool = executable(
        tmp_path,
        'if [ "$1" = "--version" ]; then printf "5.5.9\\n"; exit; fi\n'
        "printf '" + HEADER + "'\n"
        "for last; do :; done\n"
        '''awk '
        /^>/{id=substr($0,2);next}
        {printf "motif\\t\\t%s\\t1\\t%d\\t+\\t%d\\t0.01\\t\\t%s\\n",
        id,length($0),length($0),$0}
        ' "$last"''',
    )
    base = parts.PreparationSpec(
        wide_source(tmp_path),
        parts.Retention(
            count=1,
            policy="mmr",
            rank_by="best_hit_score",
            mmr=parts.MMR(3, 0.5, "score_percentile"),
        ),
        sampling=parts.Sampling(planning.Length(exact=6), strategy="consensus"),
        budget=parts.CandidateBudget(3),
        scoring=parts.FimoScoring(executable=tool),
        score_bands=parts.ScoreBands((0.5,)),
    )
    windows = {"short": parts.MotifWindow(2), "long": parts.MotifWindow(3)}
    request = parts.PreparationSet.from_windows(
        base, windows=windows, candidate_length="window", max_recipes=2
    )
    compact = {
        "schema": "dense_arrays.preparation_windows.v1",
        "base": preparation_to_dict(base, base=tmp_path),
        "windows": [
            {"id": name, "window": window.to_dict()} for name, window in windows.items()
        ],
        "candidate_length": "window",
        "max_recipes": 2,
    }
    raw = tmp_path / "windows.json"
    raw.write_text(json.dumps(compact))
    decoded = read_source(raw)
    assert decoded == request
    expected = parts.PreparationSet(
        {
            name: replace(
                base,
                source=replace(base.source, window=window),
                sampling=replace(
                    base.sampling, length=planning.Length(exact=window.length)
                ),
            )
            for name, window in windows.items()
        }
    )
    plan = da.plan(request)
    assert da.plan(expected).plan_id == plan.plan_id
    resolved = tuple(plan.resolved.recipes.values())
    assert [r.source.motif.width for r in resolved] == [2, 3]
    assert [(r.source.selection.start, r.source.selection.end) for r in resolved] == [
        (1, 3),
        (0, 3),
    ]
    assert len({r.source.scoring.binding_id for r in resolved}) == 2
    assert len({r.source.input.motif.model_id for r in resolved}) == 1
    assert plan.preview["candidate_budget"] == 6
    assert plan.preview["requested_retention"] == 2
    saved = tmp_path / "plan.json"
    plan.write(saved)
    cli_plan = tmp_path / "cli.plan.json"
    cli = CliRunner().invoke(app, ["plan", str(raw), "--out", str(cli_plan), "--json"])
    assert cli.exit_code == 0, cli.output
    assert read_source(cli_plan).plan_id == plan.plan_id
    human = CliRunner().invoke(app, ["plan", str(raw)])
    assert human.exit_code == 0, human.output
    assert "Motif window: [1, 3)" in human.stdout
    assert "Motif window: [0, 3)" in human.stdout
    pool = da.prepare(plan, out=tmp_path / "python")
    cli = CliRunner().invoke(
        app, ["prepare", str(raw), "--out", str(tmp_path / "cli"), "--json"]
    )
    assert cli.exit_code == 0, cli.output
    assert json.loads(cli.stdout)["pool_id"] == pool.pool_id
    with da.inspect(pool, view="candidates", all=True).records() as rows:
        candidates = [r.candidate for r in rows]
    assert [len(c.part.sequence) for c in candidates] == [2, 2, 2, 3, 3, 3]
    assert [len(c.score.core) for c in candidates] == [2, 2, 2, 3, 3, 3]
    exported = tmp_path / "expanded.json"
    da.export(decoded, view="request", out=exported)
    assert (
        json.loads(exported.read_text())["schema"] == "dense_arrays.preparation_set.v1"
    )
    assert da.plan(read_source(exported)).plan_id == plan.plan_id
    base.source.path.unlink()
    tool.unlink()
    with patch(
        "subprocess.Popen", side_effect=AssertionError("verification invoked a tool")
    ):
        assert da.inspect(pool, verify=True).retained_parts == 2
        assert read_source(saved).plan_id == plan.plan_id


@pytest.mark.parametrize("change", ["duplicate", "unknown", "excess", "nested"])
def test_compact_window_requests_reject_ambiguous_or_unbounded_expansion(
    tmp_path: Path, change: str
):

    value = {
        "schema": "dense_arrays.preparation_windows.v1",
        "base": preparation_to_dict(base_recipe()),
        "windows": [{"id": "a", "window": {"length": 3}}],
        "max_recipes": 1,
        "candidate_length": "window",
    }
    if change == "duplicate":
        value["max_recipes"] = 2
        value["windows"] *= 2
    elif change == "unknown":
        value["automatic"] = True
    elif change == "excess":
        value["windows"] *= 2
        value["base"] = {}
    elif change == "nested":
        value["base"] = {"schema": "dense_arrays.preparation_set.v1", "recipes": []}
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(value))
    expected = {
        "duplicate": "unique window IDs",
        "unknown": "window preparation: unknown fields",
        "excess": "max_recipes",
        "nested": "window base requires",
    }[change]
    with pytest.raises((ValueError, TypeError), match=expected):
        read_source(path)


def test_preserved_candidate_length_must_fit_every_declared_window():
    with pytest.raises(ValueError, match="candidate length"):
        parts.PreparationSet.from_windows(
            base_recipe(),
            windows={"too_wide": parts.MotifWindow(7)},
            candidate_length="base",
            max_recipes=1,
        )
