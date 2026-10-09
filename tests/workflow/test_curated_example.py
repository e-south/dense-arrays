"""Exercise the maintained curated-library consumer through native operations.

Author: Eric J. South.
"""

import csv
import hashlib
import io
import json
import shutil
from collections import Counter
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import planning
from dense_arrays.cli import app

ROOT = Path(__file__).parents[2]
EXAMPLE = ROOT / "docs/examples/curated-library"
FIXTURE = ROOT / "tests/fixtures/workflow/densegen-curated-demo-v1.json"


def stage_example(destination: Path) -> None:
    """Copy only maintained inputs; every test owns its runtime destination."""
    shutil.copytree(EXAMPLE, destination, dirs_exist_ok=True)


def test_curated_example_preserves_source_rows_and_declares_two_targets(
    tmp_path: Path,
):
    stage_example(tmp_path)
    fixture = json.loads(FIXTURE.read_text())
    with (tmp_path / "parts.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 752
    assert Counter(row["group"] for row in rows) == {
        "TF_A": 250,
        "TF_B": 250,
        "TF_C": 250,
        "anchors": 2,
    }
    original = io.StringIO(newline="")
    writer = csv.writer(original, lineterminator="\n")
    writer.writerow(["tf", "tfbs"])
    writer.writerows((r["group"], r["sequence"]) for r in rows[:750])
    assert (
        hashlib.sha256(original.getvalue().encode()).hexdigest()
        == fixture["source_input_sha256"]
    )
    preparation = da.inspect(tmp_path / "prepare.yaml", view="request").request
    pool = da.prepare(preparation, out=tmp_path / "pool")
    assert da.inspect(pool, verify=True).retained_parts == 752
    plan = da.plan(da.inspect(tmp_path / "design.yaml", view="request").request)
    assert [cell.cell_id for cell in plan.cells] == [
        "architecture=baseline",
        "architecture=fixed_pair",
    ]
    assert [cell.target for cell in plan.cells] == [50, 50]
    assert [len(cell.plan.request.parts) for cell in plan.cells] == [750, 752]
    assert plan.request.base.length.exact == 100
    cli = CliRunner().invoke(app, ["plan", str(tmp_path / "design.yaml"), "--json"])
    assert cli.exit_code == 0, cli.output
    assert json.loads(cli.stdout)["plan_id"] == plan.plan_id


def test_curated_example_generation_preserves_geometry_and_portable_evidence(
    tmp_path: Path,
):
    stage_example(tmp_path)
    preparation = da.inspect(tmp_path / "prepare.yaml", view="request").request
    da.prepare(preparation, out=tmp_path / "pool")
    request = da.inspect(tmp_path / "design.yaml", view="request").request
    bounded = request.with_changes(
        allocation=planning.Allocation(per_cell=2),
        base=request.base.with_changes(limits=planning.Limits(40, 60, 10)),
    )
    run = da.run(bounded, out=tmp_path / "run")
    summary = da.inspect(run, verify=True)
    assert summary.state == "completed"
    assert summary.accepted == 4
    with da.inspect(run, view="designs", all=True).records() as records:
        for design in records:
            sequence = design.realized.sequence
            assert len(sequence) == 100
            assert all(rule["passed"] for rule in design.requirements)
            placements = {p.feature_id: p for p in design.realized.placements}
            pad_length = design.realized.provenance["assembly"]["padding_length"]
            if pad_length:
                pad = sequence[:pad_length]
                assert 0.4 <= (pad.count("G") + pad.count("C")) / len(pad) <= 0.6
            if design.cell_id == "architecture=fixed_pair":
                up, down = placements["anchor_up"], placements["anchor_down"]
                assert up.sequence == "TTGACA"
                assert down.sequence == "TATAAT"
                assert 16 <= down.start - up.end <= 18
            else:
                assert not {"anchor_up", "anchor_down"} & placements.keys()
    da.export(run, all=True, format="bundle", out=tmp_path / "library")
    shutil.rmtree(tmp_path / "run")
    shutil.rmtree(tmp_path / "pool")
    (tmp_path / "parts.csv").unlink()
    assert da.inspect(tmp_path / "library", verify=True).designs == 4


def test_curated_example_rejects_changed_input_before_preparation(tmp_path: Path):
    stage_example(tmp_path)
    plan = da.plan(da.inspect(tmp_path / "prepare.yaml", view="request").request)
    with (tmp_path / "parts.csv").open("a") as stream:
        stream.write("unexpected,AAA,TF_A\n")
    with pytest.raises(ValueError, match="changed"):
        da.prepare(plan, out=tmp_path / "pool")
    assert not (tmp_path / "pool").exists()
