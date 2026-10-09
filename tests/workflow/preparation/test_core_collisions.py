"""Recipe-wide collision choices preserve observed core and sequence semantics.

Author: Eric J. South.
"""

import json
import sqlite3
from dataclasses import replace
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.cli import app
from dense_arrays.parts.candidates import Candidate
from dense_arrays.parts.retention.collisions import validate_collisions
from dense_arrays.planning.preparation.requests import (
    preparation_from_dict,
    preparation_to_dict,
)
from dense_arrays.workflow import preparation
from dense_arrays.workflow.inputs import read_source

from .test_fimo import scorer
from .test_motifs import source


def background():
    return parts.PreparationSpec(
        parts.Background(),
        parts.Retention(count=1, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=4)),
        budget=parts.CandidateBudget(2),
    )


def test_core_collision_policy_is_bound_without_changing_default_plans(tmp_path: Path):
    original = parts.PreparationSet({"one": background()})
    explicit = original.with_changes(core_collisions="preserve")
    previous = da.plan(original)
    frozen = json.loads(
        (
            Path(__file__).parents[2] / "fixtures/workflow/preparation-set-v1.json"
        ).read_text()
    )
    assert previous.to_dict() == frozen
    assert planning.PreparationPlan.from_dict(frozen).to_dict() == frozen
    assert da.plan(explicit).to_dict() == previous.to_dict()
    assert "core_collisions" not in preparation_to_dict(original)
    assert "core_collisions" not in previous.preview
    strict = original.with_changes(core_collisions="error")
    plan = da.plan(strict)
    assert plan.plan_id != previous.plan_id
    assert plan.preview["core_collisions"] == "error"
    assert preparation_from_dict(preparation_to_dict(strict)) == strict
    saved = tmp_path / "plan.json"
    plan.write(saved)
    assert read_source(saved).to_dict() == plan.to_dict()
    cli = CliRunner().invoke(app, ["plan", str(saved), "--json"])
    assert cli.exit_code == 0, cli.output
    assert json.loads(cli.stdout)["plan_id"] == plan.plan_id
    human = CliRunner().invoke(app, ["plan", str(saved)])
    assert human.exit_code == 0, human.output
    assert "observed core collisions error" in human.stdout
    assert "parts without a core are excluded" in human.stdout


@pytest.mark.parametrize("policy", [None, True, "allow", "deduplicate"])
def test_unsupported_core_collision_policies_fail_before_planning(policy: object):
    with pytest.raises((ValueError, TypeError), match="core_collisions"):
        parts.PreparationSet({"one": background()}, core_collisions=policy)
    encoded = preparation_to_dict(parts.PreparationSet({"one": background()}))
    encoded["core_collisions"] = policy
    with pytest.raises((ValueError, TypeError), match="core_collisions"):
        preparation_from_dict(encoded)


def test_distinct_sequences_with_shared_cores_require_explicit_preservation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    motif = source(tmp_path)
    tool = scorer(tmp_path, "motif\t\tcandidate_0\t2\t3\t+\t3\t0.0625\t\tAC\n")
    recipe = parts.PreparationSpec(
        motif,
        parts.Retention(count=1, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=4)),
        budget=parts.CandidateBudget(1),
        scoring=parts.FimoScoring(executable=tool),
        seed=1,
    )
    monkeypatch.setattr(
        preparation,
        "sample_sequence",
        lambda **kwargs: "TACT" if kwargs["seed"] == 1 else "GACA",
    )
    request = parts.PreparationSet(
        {"first": recipe, "second": recipe.with_changes(seed=2)},
        core_collisions="error",
    )
    with pytest.raises(ValueError, match="retained core") as conflict:
        da.prepare(request, out=tmp_path / "conflict")
    assert "first" in str(conflict.value)
    assert "second" in str(conflict.value)
    with pytest.raises(ValueError, match=r"manifest|commit"):
        da.inspect(tmp_path / "conflict")

    kept = request.with_changes(core_collisions="preserve")
    pool = da.prepare(kept, out=tmp_path / "preserved")
    assert da.inspect(pool, verify=True).retained_parts == 2
    with da.inspect(pool, view="parts", all=True).records() as rows:
        retained = [r.part for r in rows]
    assert [p.sequence for p in retained] == ["TACT", "GACA"]
    assert [p.core_sequence for p in retained] == ["AC", "AC"]

    saved = tmp_path / "strict.json"
    da.export(request, view="request", out=saved)
    cli = CliRunner().invoke(
        app, ["prepare", str(saved), "--out", str(tmp_path / "cli-conflict")]
    )
    assert cli.exit_code == 2, cli.output
    assert "retained core" in cli.stderr

    # A checksummed plan with a stricter policy must also reject saved collisions.
    wire = da.plan(request).to_dict(base=pool.path)
    with sqlite3.connect(pool.path / "pool.sqlite3") as connection:
        connection.execute(
            "UPDATE preparation SET payload=?,digest=? WHERE id=1",
            (canonical_json(wire), semantic_digest(wire)),
        )
    with pytest.raises(ValueError, match="retained core"):
        da.inspect(pool, verify=True)


def candidate(
    index: int,
    sequence: str,
    *,
    orientation: str | None = "forward",
    width: int = 2,
    recipe: str,
):
    return Candidate(
        index,
        parts.Part(
            f"{recipe}/candidate_{index}",
            sequence,
            group="shared_label",
            core_start=1 if orientation else None,
            core_end=1 + width if orientation else None,
            core_orientation=orientation,
        ),
        representative=index,
        rank=1,
        retained=True,
        recipe_id=recipe,
        recipe_index=1,
    )


@pytest.mark.parametrize(
    "left,right,strand,width,collides",
    [
        ("TACG", "CGTA", "reverse", 2, True),  # AC == reverse-complement(GT)
        ("TACG", "CGTA", "forward", 2, False),  # AC != GT
        ("CATG", "GATT", "reverse", 2, True),  # palindrome AT == AT
        ("TACG", "TACGT", "forward", 3, False),  # AC != ACG
    ],
)
def test_core_equivalence_uses_exact_observed_orientation(
    left: str, right: str, strand: str, width: int, collides: bool
):
    records = (
        candidate(1, left, recipe="one"),
        candidate(2, right, recipe="two", orientation=strand, width=width),
    )
    if collides:
        with pytest.raises(ValueError, match="retained core"):
            validate_collisions(records, sequences="preserve", cores="error")
    else:
        validate_collisions(records, sequences="preserve", cores="error")
    validate_collisions(records, sequences="preserve", cores="preserve")


@pytest.mark.parametrize("options", [{"recipe": "one"}, {"orientation": None}])
def test_core_comparison_excludes_other_populations(options: dict):
    records = (
        candidate(1, "TACG", recipe="one"),
        candidate(2, "CACA", **({"recipe": "two"} | options)),
    )
    validate_collisions(records, sequences="error", cores="error")


def test_unretained_candidates_do_not_create_collisions():
    records = (
        candidate(1, "TACG", recipe="one"),
        replace(candidate(2, "TACG", recipe="two"), retained=False),
    )
    validate_collisions(records, sequences="error", cores="error")


def test_sequence_collisions_remain_independent_of_core_policy():
    records = (
        candidate(1, "AAAA", recipe="one", orientation=None),
        candidate(2, "AAAA", recipe="two", orientation=None),
    )
    validate_collisions(records, sequences="preserve", cores="error")
    with pytest.raises(ValueError, match="retained sequence"):
        validate_collisions(records, sequences="error", cores="preserve")


def test_named_windows_preserve_explicit_core_policy_without_reading_inputs():
    base = background().with_changes(
        source=parts.PWMArtifact("motif.json"), scoring=parts.FimoScoring()
    )
    window = parts.MotifWindow(2)
    expected = parts.PreparationSet.from_windows(
        base,
        windows={"two_bases": window},
        candidate_length="window",
        max_recipes=1,
        core_collisions="error",
    )
    value = {
        "schema": "dense_arrays.preparation_windows.v1",
        "base": preparation_to_dict(base),
        "windows": [{"id": "two_bases", "window": window.to_dict()}],
        "candidate_length": "window",
        "max_recipes": 1,
        "core_collisions": "error",
    }
    assert preparation_from_dict(value) == expected
    value["core_collisions"] = "deduplicate"
    with pytest.raises(ValueError, match="core_collisions"):
        preparation_from_dict(value)
