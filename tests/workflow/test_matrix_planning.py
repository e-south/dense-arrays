"""Bounded matrix previews make every allocation and inactive cell explicit.

Author: Eric J. South.
"""

import json
import shutil
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays.cli import app
from dense_arrays.planning.matrices import Allocation, MatrixSpec, Variant, resolution


def test_variant_adds_requirements_only_to_its_selected_cell(tmp_path: Path):
    base = planning.DesignSpec(
        [parts.Part("up", "AAA"), parts.Part("down", "CCC")],
        planning.Length(maximum=6),
        strands="single",
    )
    added = (
        planning.Fixed("up-fixed", "up", "forward"),
        planning.Fixed("down-fixed", "down", "forward"),
        planning.Spacing("gap", "up", "down", 0, 0),
    )
    spec = MatrixSpec(
        base,
        {"mode": {"free": Variant(), "fixed": Variant(add_requirements=added)}},
        Allocation(per_cell=1),
        max_cells=2,
    )
    plan = da.plan(spec)
    assert plan.cells[0].plan.request.requirements == ()
    assert plan.cells[1].plan.request.requirements == added
    assert MatrixSpec.from_dict(spec.to_dict()) == spec
    assert "add_requirements" not in Variant().to_dict()
    da.export(spec, out=tmp_path / "request.json")
    result = CliRunner().invoke(app, ["plan", str(tmp_path / "request.json"), "--json"])
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout)["plan_id"] == plan.plan_id
    run = da.run(plan, out=tmp_path / "run")
    assert da.inspect(run, verify=True).accepted == 2
    with da.inspect(run, view="designs", all=True).records() as records:
        fixed = next(d for d in records if d.cell_id == "mode=fixed")
    assert fixed.realized.sequence == "AAACCC"
    comparison = da.inspect(
        plan,
        view="plan",
        compare=da.plan(
            spec.with_changes(
                axes={"mode": {"free": Variant(), "fixed": Variant()}},
            )
        ),
    )
    assert any(
        change.path[:3] == ("cells", "mode=fixed", "requirements")
        for change in comparison.changes
    )
    adjusted = spec.with_changes(
        axes={
            "mode": {
                "free": Variant(),
                "fixed": Variant(
                    add_requirements=(
                        *added[:2],
                        planning.Spacing("gap", "up", "down", -1, 0),
                    )
                ),
            }
        }
    )
    changes = {
        change.path: change
        for change in da.inspect(plan, view="plan", compare=da.plan(adjusted)).changes
    }
    assert (
        changes[
            ("matrix", "axes", "mode", "fixed", "add_requirements", "gap", "min")
        ].after
        == -1
    )


@pytest.mark.parametrize("conflict", ["base", "axes", "replace"])
def test_variant_additions_reject_identity_conflicts(conflict: str):
    rule = planning.Fixed("anchor", "up", "forward")
    first = Variant(add_requirements=(rule,))
    base = planning.DesignSpec(
        [parts.Part("up", "AAA")],
        planning.Length(maximum=3),
        requirements=(rule,) if conflict == "base" else (),
    )
    axes = {"first": {"one": first}}
    if conflict != "base":
        axes["second"] = {
            "one": first if conflict == "axes" else Variant(requirements=(rule,))
        }
    with pytest.raises(ValueError, match="requirement"):
        da.plan(MatrixSpec(base, axes, Allocation(per_cell=1), max_cells=1))


def test_variant_replacement_still_rejects_unknown_ids():
    base = planning.DesignSpec([parts.Part("up", "AAA")], planning.Length(maximum=3))
    with pytest.raises(ValueError, match="replace known"):
        da.plan(
            MatrixSpec(
                base,
                {
                    "mode": {
                        "fixed": Variant(
                            requirements=(planning.Fixed("new", "up", "forward"),)
                        )
                    }
                },
                Allocation(per_cell=1),
                max_cells=1,
            )
        )


def test_variant_addition_references_are_validated_after_combining_axes():
    base = planning.DesignSpec(
        [parts.Part("up", "AAA"), parts.Part("down", "CCC")],
        planning.Length(maximum=6),
    )
    axes = {
        "anchors": {
            "on": Variant(
                add_requirements=(
                    planning.Fixed("up-fixed", "up", "forward"),
                    planning.Fixed("down-fixed", "down", "forward"),
                )
            )
        },
        "spacing": {
            "zero": Variant(
                add_requirements=(planning.Spacing("gap", "up", "down", 0, 0),)
            )
        },
    }
    assert (
        len(
            da.plan(MatrixSpec(base, axes, Allocation(per_cell=1), max_cells=1))
            .cells[0]
            .plan.request.requirements
        )
        == 3
    )
    with pytest.raises(ValueError, match="fixed"):
        da.plan(
            MatrixSpec(
                base, {"spacing": axes["spacing"]}, Allocation(per_cell=1), max_cells=1
            )
        )


@pytest.mark.parametrize(
    "additions", [({"kind": "gc"},), (planning.Fixed("a", "up", "forward"),) * 2]
)
def test_variant_addition_syntax_is_strict(additions: tuple):
    with pytest.raises((ValueError, TypeError), match=r"requirement|identities"):
        Variant(add_requirements=additions)


def request(
    allocation: Allocation,
    *,
    pairing: str = "cross_product",
    pairs: tuple[dict[str, str], ...] = (),
):
    return MatrixSpec(
        base=planning.DesignSpec(
            parts=[parts.Part("up", "AAA"), parts.Part("down", "CCC")],
            length=planning.Length(maximum=6),
            strands="single",
            requirements=[
                planning.Fixed("up-fixed", "up", "forward"),
                planning.Fixed("down-fixed", "down", "forward"),
            ],
        ),
        axes={
            "up": {
                "one": Variant(parts=[parts.Part("up", "AAA")]),
                "two": Variant(parts=[parts.Part("up", "GGG")]),
            },
            "down": {
                "one": Variant(parts=[parts.Part("down", "CCC")]),
                "two": Variant(parts=[parts.Part("down", "TTT")]),
            },
        },
        allocation=allocation,
        max_cells=4,
        pairing=pairing,
        pairs=pairs,
    )


def test_cross_product_has_declared_order_and_balanced_remainders(
    monkeypatch: pytest.MonkeyPatch,
):
    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("matrix planning built a solver")

    monkeypatch.setattr(da.Optimizer, "build_model", forbidden)
    result = da.plan(request(Allocation(total=10, policy="balanced")))
    assert [dict(cell.choices) for cell in result.cells] == [
        {"up": "one", "down": "one"},
        {"up": "one", "down": "two"},
        {"up": "two", "down": "one"},
        {"up": "two", "down": "two"},
    ]
    assert [cell.target for cell in result.cells] == [3, 3, 2, 2]
    assert result.total == 10
    assert result.cells[-1].plan.request.parts[0].sequence == "GGG"
    assert result.cells[-1].plan.request.parts[1].sequence == "TTT"
    assert result.preview["allocation_policy"] == "ordered_balanced.v1"
    assert result.preview["execution_supported"] is True


def test_small_total_requires_named_inactive_cells_and_does_not_drop_them():
    with pytest.raises(ValueError, match="zero-target"):
        da.plan(request(Allocation(total=2, policy="balanced")))
    result = da.plan(
        request(
            Allocation(
                total=2,
                policy="balanced",
                zero_cells=(
                    "down=two,up=one",
                    "down=two,up=two",
                ),
            )
        )
    )
    assert [cell.target for cell in result.cells] == [1, 0, 1, 0]
    assert len(result.cells) == 4
    assert [cell.active for cell in result.cells] == [True, False, True, False]
    assert result.cells[1].plan.request.target.count == 0


def test_zip_pairs_by_choice_identity_and_explicit_pairs_preserve_order():
    zipped = da.plan(request(Allocation(per_cell=2), pairing="zip"))
    assert [cell.cell_id for cell in zipped.cells] == [
        "down=one,up=one",
        "down=two,up=two",
    ]
    explicit = da.plan(
        request(
            Allocation(per_cell=1),
            pairing="explicit",
            pairs=(
                {"up": "two", "down": "one"},
                {"up": "one", "down": "two"},
            ),
        )
    )
    assert [cell.cell_id for cell in explicit.cells] == [
        "down=one,up=two",
        "down=two,up=one",
    ]


def test_matrix_limit_precedes_input_read_and_bad_allocations_fail():
    spec = request(Allocation(per_cell=1))
    with pytest.raises(ValueError, match="max_cells"):
        da.plan(
            spec.with_changes(
                base=spec.base.with_changes(
                    parts=parts.PartTable(Path("absent.csv"), "csv")
                ),
                max_cells=3,
            )
        )
    with pytest.raises(ValueError, match="exactly one"):
        Allocation(total=4, per_cell=1, policy="balanced")
    with pytest.raises((TypeError, ValueError), match="integer"):
        Allocation(per_cell=True)
    with pytest.raises(ValueError, match="cell"):
        da.plan(request(Allocation(counts={"unknown": 4})))


def test_matrix_plan_python_cli_round_trip_and_execution(
    tmp_path: Path,
):
    spec = request(Allocation(per_cell=1))
    path = tmp_path / "matrix.json"
    path.write_text(json.dumps(spec.to_dict()))
    resolved = da.plan(spec)
    saved = tmp_path / "plan.json"
    response = CliRunner().invoke(
        app, ["plan", str(path), "--out", str(saved), "--json"]
    )
    assert response.exit_code == 0, response.output
    assert json.loads(response.stdout) == resolved.to_dict()
    assert (
        type(resolved).from_dict(json.loads(saved.read_text())).to_dict()
        == resolved.to_dict()
    )
    run = da.run(resolved, out=tmp_path / "run")
    assert da.inspect(run, verify=True).accepted == 4


def test_explicit_zero_target_publishes_an_empty_verified_run_without_solving(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("zero-target run built a solver")

    monkeypatch.setattr(da.Optimizer, "build_model", forbidden)
    spec = request(Allocation(per_cell=1)).base.with_changes(
        target=planning.Target(count=0)
    )
    run = da.run(spec, out=tmp_path / "empty")
    summary = da.inspect(run, verify=True)
    assert summary.state == "completed"
    assert summary.accepted == summary.target == summary.counts["started"] == 0
    quality = da.inspect(run, view="quality").to_dict()
    assert quality["selection"]["designs"] == 0
    selected = da.inspect(
        run,
        view="selection",
        select=reporting.LibrarySelection(take=reporting.Take(count=0)),
    )
    assert selected.selected == selected.shortfall == 0
    bundle = tmp_path / "empty-bundle"
    receipt = da.export(run, all=True, format="bundle", out=bundle)
    assert receipt.records == 0
    assert da.inspect(bundle, verify=True).designs == 0
    fasta = tmp_path / "empty.fasta"
    da.export(bundle, all=True, view="sequences", format="fasta", out=fasta)
    assert fasta.read_bytes() == b""


def test_axis_order_survives_canonical_json_and_preserves_cell_streams():
    spec = request(Allocation(per_cell=1))
    ordered = da.plan(spec)
    reversed_axes = {
        axis: dict(reversed(tuple(options.items())))
        for axis, options in reversed(tuple(spec.axes.items()))
    }
    reordered = da.plan(spec.with_changes(axes=reversed_axes))
    assert reordered.plan_id != ordered.plan_id
    assert {c.cell_id: c.plan.plan_id for c in reordered.cells} == {
        c.cell_id: c.plan.plan_id for c in ordered.cells
    }
    assert len({c.plan.request.seed for c in reordered.cells}) == 4
    wire = json.loads(json.dumps(reordered.to_dict(), sort_keys=True))
    restored = planning.MatrixPlan.from_dict(wire)
    assert [c.cell_id for c in restored.cells] == [
        "down=two,up=two",
        "down=two,up=one",
        "down=one,up=two",
        "down=one,up=one",
    ]
    assert restored.plan_id == reordered.plan_id


def test_matrix_saved_plan_moves_without_reopening_sources(tmp_path: Path):
    folder = tmp_path / "original"
    folder.mkdir()
    table = folder / "parts.csv"
    table.write_text("part_id,sequence\nup,AAA\ndown,CCC\n")
    spec = request(Allocation(per_cell=1))
    resolved = da.plan(
        spec.with_changes(
            base=spec.base.with_changes(parts=parts.PartTable(table, "csv"))
        )
    )
    resolved.write(folder / "plan.json")
    moved = tmp_path / "moved"
    shutil.move(folder, moved)
    restored = planning.MatrixPlan.from_dict(
        json.loads((moved / "plan.json").read_text()), base=moved
    )
    assert restored.plan_id == resolved.plan_id
    restored.base.verify_inputs()
    (moved / "parts.csv").unlink()
    assert (
        planning.MatrixPlan.from_dict(
            json.loads((moved / "plan.json").read_text()), base=moved
        ).plan_id
        == restored.plan_id
    )
    with pytest.raises(FileNotFoundError):
        restored.cells[0].plan.verify_inputs()


@pytest.mark.parametrize("change", ["target", "policy", "seed", "cell"])
def test_saved_matrix_rejects_changed_evidence(change: str):
    wire = da.plan(request(Allocation(per_cell=1))).to_dict()
    if change == "target":
        wire["cells"][0]["target"] = 2
    elif change == "policy":
        wire["preview"]["allocation_policy"] = "unknown.v2"
    elif change == "seed":
        wire["cells"][0]["plan"]["request"]["seed"] = 123
    else:
        wire["cells"][0]["cell_id"] = "up=missing"
    with pytest.raises(ValueError, match="identities, targets or policies"):
        planning.MatrixPlan.from_dict(wire)


@pytest.mark.parametrize("missing", ["base", "axes", "allocation", "max_cells"])
def test_incomplete_matrix_has_actionable_python_and_cli_error(
    tmp_path: Path, missing: str
):
    wire = request(Allocation(per_cell=1)).to_dict()
    wire.pop(missing)
    with pytest.raises(ValueError, match=missing):
        MatrixSpec.from_dict(wire)
    source = tmp_path / "invalid.json"
    source.write_text(json.dumps(wire))
    result = CliRunner().invoke(app, ["plan", str(source), "--json"])
    assert result.exit_code == 2, result.output
    assert json.loads(result.stdout)["code"] == "invalid_input"


def test_pairing_and_replacement_errors_do_not_silently_change_a_matrix():
    spec = request(Allocation(per_cell=1))
    with pytest.raises(ValueError, match="identical named choices"):
        da.plan(
            spec.with_changes(
                pairing="zip",
                axes={"up": spec.axes["up"], "down": {"one": spec.axes["down"]["one"]}},
            )
        )
    with pytest.raises(ValueError, match="repeat a cell"):
        da.plan(
            spec.with_changes(
                pairing="explicit",
                pairs=({"up": "one", "down": "one"}, {"up": "one", "down": "one"}),
            )
        )
    for identity, message in (("up", "both replace"), ("absent", "known part")):
        axes = dict(spec.axes)
        axes["down"] = {"one": Variant(parts=[parts.Part(identity, "CCC")])}
        with pytest.raises(ValueError, match=message):
            da.plan(spec.with_changes(axes=axes))


@pytest.mark.parametrize("field", ["choices", "variant"])
def test_incomplete_ordered_axis_records_fail_with_field_context(field: str):
    wire = request(Allocation(per_cell=1)).to_dict()
    record = wire["axes"][0] if field == "choices" else wire["axes"][0]["choices"][0]
    record.pop(field)
    with pytest.raises(ValueError, match=field):
        MatrixSpec.from_dict(wire)


def test_matrix_plan_rejects_unresolved_constructor_inputs():
    spec = request(Allocation(per_cell=1))
    with pytest.raises(TypeError, match="GenerationPlan"):
        planning.MatrixPlan(spec, spec.base)
    with pytest.raises(TypeError, match="MatrixSpec"):
        planning.MatrixPlan(None, da.plan(spec.base))


@pytest.mark.parametrize("pairing", ["cross_product", "zip", "explicit"])
def test_saved_matrix_rejects_understated_cell_count_before_expansion(
    pairing: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    base = planning.DesignSpec(
        [parts.Part("site", "ACGTTGCAAGTCCTGA")], planning.Length(maximum=16)
    )
    request = MatrixSpec(
        base,
        {axis: {str(i): Variant() for i in range(3)} for axis in ("x", "y")},
        Allocation(per_cell=1),
        max_cells=9,
        pairing=pairing,
        pairs=({"x": "0", "y": "0"},) if pairing == "explicit" else (),
    )
    plan = da.plan(request)
    value = plan.to_dict(base=tmp_path)
    assert planning.MatrixPlan.from_dict(value, base=tmp_path).plan_id == plan.plan_id
    value["cells"] = []
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(value))

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("malformed saved matrix reached cell expansion")

    monkeypatch.setattr(resolution, "expand", forbidden)
    monkeypatch.setattr(resolution, "_cell", forbidden)
    with pytest.raises(ValueError, match="cell count"):
        da.inspect(path, view="plan", read_limits=reporting.ReadLimits(identities=1))


def test_saved_matrix_rejects_understated_child_parts_before_compiling(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    base = planning.DesignSpec(
        [parts.Part("site", "ACGTTGCAAGTCCTGA")], planning.Length(maximum=16)
    )
    plan = da.plan(
        MatrixSpec(
            base,
            {"x": {str(i): Variant() for i in range(3)}},
            Allocation(per_cell=1),
            max_cells=3,
        )
    )
    value = plan.to_dict(base=tmp_path)
    for cell in value["cells"]:
        cell["plan"]["request"]["parts"] = []
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(value))

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("understated matrix dimensions reached child plan construction")

    monkeypatch.setattr(resolution, "_cell", forbidden)
    with pytest.raises(ValueError, match="inconsistent parts count"):
        da.inspect(path, view="plan", read_limits=reporting.ReadLimits(identities=4))
