"""Invalid tables report bounded row evidence before any output is published.

Module Author(s): Eric J. South
"""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.cli import app
from dense_arrays.parts.ingestion import read_parts


@pytest.mark.parametrize("format_name,separator", [("csv", ","), ("tsv", "\t")])
def test_invalid_row_totals_and_sample_are_bounded(
    tmp_path: Path, format_name: str, separator: str
):
    table = tmp_path / f"parts.{format_name}"
    rows = ["site,bases", "valid,ACGTTGCAAGTCCTGA"]
    rows.extend(f"bad-{index},ACGNTGCAAGTCCTGA" for index in range(25))
    table.write_text("\n".join(rows).replace(",", separator) + "\n")
    source = parts.PartTable(
        table, format_name, columns={"part_id": "site", "sequence": "bases"}
    )
    with pytest.raises(parts.TableImportError) as caught:
        read_parts(source)
    error = caught.value
    assert error.rows == 26
    assert error.invalid_rows == 25
    assert len(error.diagnostics) == error.sample_limit == 20
    assert [item.row for item in error.diagnostics] == list(range(2, 22))
    assert all(item.columns["sequence"] == "bases" for item in error.diagnostics)
    assert error.to_dict()["omitted_rows"] == 5
    assert "25 invalid rows" in str(error)
    assert "showing 20" in str(error)


def test_duplicate_ids_preserve_original_row_locations(tmp_path: Path):
    table = tmp_path / "parts.csv"
    table.write_text(
        "part_id,sequence\nsame,ACGTTGCAAGTCCTGA\nsame,TTGACCGATAGCTACG\nsame,ATGCTTAGGACGTTCA\n"
    )
    with pytest.raises(parts.TableImportError) as caught:
        read_parts(parts.PartTable(table, "csv"))
    error = caught.value
    assert error.rows == 3
    assert error.invalid_rows == 2
    assert [(item.row, item.related_row) for item in error.diagnostics] == [
        (2, 1),
        (3, 1),
    ]
    assert all(item.columns == {"part_id": "part_id"} for item in error.diagnostics)


def test_import_failure_is_shared_by_python_and_cli_without_partial_output(
    tmp_path: Path,
):
    table = tmp_path / "parts.csv"
    table.write_text(
        "site,bases\na,ACGNTGCAAGTCCTGA\nb,ACGTTGCAAGTCCTGA\nc,acgttgcaagtcctga\n"
    )
    source = parts.PartTable(
        table, "csv", columns={"part_id": "site", "sequence": "bases"}
    )
    request = planning.DesignSpec(parts=source, length=planning.Length(maximum=40))
    with pytest.raises(parts.TableImportError) as caught:
        da.run(request, out=tmp_path / "python-run")
    expected = caught.value.to_dict()
    assert not (tmp_path / "python-run").exists()
    recipe = tmp_path / "request.json"
    da.export(request, view="request", out=recipe)
    prepare_recipe = tmp_path / "prepare.json"
    da.export(parts.PreparationSpec(source), view="request", out=prepare_recipe)
    runner = CliRunner()
    for operation, recipe_path in (
        ("plan", recipe),
        ("run", recipe),
        ("prepare", prepare_recipe),
    ):
        out = tmp_path / f"cli-{operation}"
        result = runner.invoke(
            app, [operation, str(recipe_path), "--out", str(out), "--json"]
        )
        assert result.exit_code == 2, result.output
        error = json.loads(result.stdout)
        assert error["code"] == "invalid_input"
        assert error["import_report"] == expected
        assert not out.exists()
    with pytest.raises(parts.TableImportError) as prepared:
        da.prepare(parts.PreparationSpec(source=source), out=tmp_path / "pool")
    assert prepared.value.to_dict() == expected
    assert not (tmp_path / "pool").exists()


def test_invalid_identity_still_anchors_later_duplicate(tmp_path: Path):
    table = tmp_path / "parts.csv"
    table.write_text(
        "part_id,sequence,core_start,core_end,core_orientation\n"
        "same,ACGNTGCAAGTCCTGA,,,\n"
        "same,ACGTTGCAAGTCCTGA,,,\n"
        "other,ACGTTGCAAGTCCTGA,1,99,forward\n"
    )
    with pytest.raises(parts.TableImportError) as caught:
        read_parts(parts.PartTable(table, "csv"))
    error = caught.value
    assert error.rows == error.invalid_rows == 3
    assert error.diagnostics[1].related_row == 1
    assert error.diagnostics[2].columns["core_end"] == "core_end"
    assert "core interval" in error.diagnostics[2].message


def test_structural_failure_never_claims_complete_row_counts(tmp_path: Path):
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence\na,ACGN\nb,ACGTTGCAAGTCCTGA,extra\n")
    with pytest.raises(ValueError, match="one value per header column") as caught:
        read_parts(parts.PartTable(table, "csv"))
    assert not isinstance(caught.value, parts.TableImportError)


@pytest.mark.parametrize("format_name", ["parquet", "xlsx"])
def test_optional_formats_share_counted_row_diagnostics(
    tmp_path: Path, format_name: str
):
    table = tmp_path / f"parts.{format_name}"
    rows = [
        {"site": "a", "bases": "ACGNTGCAAGTCCTGA"},
        {"site": "b", "bases": "ACGTTGCAAGTCCTGA"},
        {"site": "c", "bases": "acgttgcaagtcctga"},
    ]
    if format_name == "parquet":
        arrow = pytest.importorskip("pyarrow")
        parquet = pytest.importorskip("pyarrow.parquet")
        parquet.write_table(arrow.Table.from_pylist(rows), table)
    else:
        excel = pytest.importorskip("openpyxl")
        book = excel.Workbook()
        book.active.append(["site", "bases"])
        for row in rows:
            book.active.append(list(row.values()))
        book.save(table)
        book.close()
    source = parts.PartTable(
        table, format_name, columns={"part_id": "site", "sequence": "bases"}
    )
    with pytest.raises(parts.TableImportError) as caught:
        da.plan(parts.PreparationSpec(source))
    assert caught.value.rows == 3
    assert caught.value.invalid_rows == 2
    assert [item.row for item in caught.value.diagnostics] == [1, 3]
    assert all(item.columns["sequence"] == "bases" for item in caught.value.diagnostics)
