"""Curated inputs preserve identities and require explicit normalization.

Author: Eric J. South.
"""

from pathlib import Path

import pytest

from dense_arrays.parts import Normalization, Part, PartSelector, PartTable
from dense_arrays.parts.ingestion import read_parts


def test_mapped_table_preserves_duplicate_sequences_and_metadata(tmp_path: Path):
    table = tmp_path / "parts.csv"
    table.write_text("site,bases,family,note,ignored\na,ACGT,A,one,x\nb,ACGT,B,two,y\n")
    result = read_parts(
        PartTable(
            table=table,
            format="csv",
            columns={"part_id": "site", "sequence": "bases", "group": "family"},
            metadata_columns={"caller.note": "note"},
        )
    )
    assert [p.part_id for p in result.parts] == ["a", "b"]
    assert [p.sequence for p in result.parts] == ["ACGT", "ACGT"]
    assert result.parts[1].metadata == {"caller.note": "two"}
    assert result.ignored_columns == ("ignored",)
    assert result.transformations == ()


def test_import_errors_identify_original_row_and_column(tmp_path: Path):
    table = tmp_path / "parts.csv"
    table.write_text("site,bases\na,ACGT\nb,acgt\n")
    with pytest.raises(ValueError, match=r"row 2.*bases.*sequence"):
        read_parts(
            PartTable(table, "csv", columns={"part_id": "site", "sequence": "bases"})
        )


def test_explicit_normalization_is_recorded(tmp_path: Path):
    table = tmp_path / "parts.tsv"
    table.write_text("sequence\n acgt \nACGT\n")
    result = read_parts(
        PartTable(
            table,
            "tsv",
            id_policy="row",
            normalization=Normalization(uppercase=True, trim_outer_whitespace=True),
        )
    )
    assert [p.part_id for p in result.parts] == ["row:1", "row:2"]
    assert result.transformations == ((1, "sequence", " acgt ", "ACGT"),)


@pytest.mark.parametrize(
    "contents",
    [
        "part_id,sequence\na,AAA\na,CCC\n",
        "part_id,sequence,sequence\na,AAA,CCC\n",
        "part_id,sequence\na,AAA,unexpected\n",
        "part_id,sequence\n,AAA\n",
        "part_id,sequence\na,A CG\n",
        "part_id,sequence\n",
    ],
)
def test_malformed_tables_do_not_return_partial_inputs(tmp_path: Path, contents: str):
    table = tmp_path / "parts.csv"
    table.write_text(contents)
    with pytest.raises(ValueError, match=r"part_id|header|column|bases|at least"):
        read_parts(PartTable(table, "csv"))


def test_typed_parts_freeze_metadata_and_validate_core():
    metadata = {"caller": {"values": [1]}}
    part = Part(
        "a",
        "ACGT",
        core_start=1,
        core_end=3,
        core_orientation="reverse",
        metadata=metadata,
    )
    metadata["caller"]["values"].append(2)
    assert part.metadata["caller"]["values"] == (1,)
    with pytest.raises(ValueError, match="core"):
        Part("b", "ACGT", core_start=1)
    with pytest.raises(ValueError, match="core"):
        Part("b", "ACGT", core_start=1, core_end=5, core_orientation="forward")


def test_selectors_do_not_infer_groups_or_accept_unknown_ids():
    pool = (Part("a", "AAA", group="A"), Part("b", "AAA", group="B"))
    assert PartSelector(groups=("A",)).indices(pool) == (0,)
    with pytest.raises(ValueError, match="unknown"):
        PartSelector(part_ids=("c",)).indices(pool)
    with pytest.raises(ValueError, match="either"):
        PartSelector(part_ids=("a",), groups=("A",))
