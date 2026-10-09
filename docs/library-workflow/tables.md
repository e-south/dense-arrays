---
title: Import part tables
description: Prepare parts from CSV, TSV, Parquet or Excel with explicit mappings and typed values.
author: Eric J. South
---

# Import part tables

Use `parts.PartTable` wherever a preparation or design request accepts parts.
CSV and TSV work with the base installation. Add the `tables` extra for Parquet
and Excel. Follow [optional-feature installation](../installation.md#optional-features)
to add the table readers to your environment.

Every format uses the same column mapping, DNA validation, occurrence identities
and normalization evidence. The format is explicit; filename extensions do not
select a parser. Parquet input is one file, not a partitioned dataset. Excel
input is one worksheet in an `.xlsx` workbook.

## Prepare a table through either interface

This example creates two small input files with identical typed contents.
Run it in a new directory with the `tables` extra installed:

```python
from pathlib import Path

import openpyxl
import pyarrow as pa
import pyarrow.parquet as pq

import dense_arrays as da

# Use the typed requests and operations needed by this example.
from dense_arrays import parts, planning

rows = [
    {"site": "upstream", "bases": "ACGTTGCAAGTCCTGA", "confidence": 0.8},
    {"site": "downstream", "bases": "GATCAGTACCTAGGTC", "confidence": 0.9},
]
pq.write_table(pa.Table.from_pylist(rows), "sites.parquet")
book = openpyxl.Workbook()
sheet = book.active
sheet.title = "binding_sites"
sheet.append(list(rows[0]))
for row in rows:
    sheet.append(list(row.values()))
book.save("sites.xlsx")
book.close()

for input_format in ("parquet", "xlsx"):
    source = parts.PartTable(
        f"sites.{input_format}",
        input_format,
        columns={"part_id": "site", "sequence": "bases"},
        metadata_columns={"confidence": "confidence"},
        **({"sheet": "binding_sites"} if input_format == "xlsx" else {}),
    )
    request = parts.PreparationSpec(source)
    da.export(request, view="request", out=f"{input_format}.json")
    prepared = da.prepare(request, out=f"pools/{input_format}-python")
    assert da.inspect(prepared, verify=True).retained_parts == 2
```

The same requests work through the CLI:

```bash
# Validate inputs and inspect or save the resolved plan.
dense-arrays plan parquet.json --out parquet.plan.json --json
# Prepare the declared pool or offered batch.
dense-arrays prepare parquet.plan.json --out pools/parquet-cli --json
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect pools/parquet-cli --view parts --all --json
# Validate inputs and inspect or save the resolved plan.
dense-arrays plan xlsx.json --out xlsx.plan.json --json
# Prepare the declared pool or offered batch.
dense-arrays prepare xlsx.plan.json --out pools/xlsx-cli --json
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect pools/xlsx-cli --view parts --all --json
```

Reuse a pool in a generation request. The optional readers are only needed when
reading the original table; persisted pools and resolved records contain native
parts:

```python
# Declare the part collection, sequence bounds and generation policy.
request = planning.DesignSpec(
    parts.PoolSource("pools/parquet-python"),
    planning.Length(maximum=32),
    strands="single",
)
# Generate under the declared bounds into a new output directory.
library = da.run(request, out="runs/from-table")
assert da.inspect(library, verify=True).accepted == 1
Path("sites.parquet").unlink()
assert da.inspect("pools/parquet-python", verify=True).retained_parts == 2
```

## Values and worksheets

| Input value | Contract |
| --- | --- |
| IDs and DNA | Stored text; numeric IDs are not converted into strings. DNA is uppercase A/C/G/T unless normalization is explicitly enabled. |
| Core coordinates | Integers, zero-based and half-open; fractional numbers and booleans fail. CSV/TSV integer text is parsed as before. |
| Optional fields | Empty text or null means absent. A core must still have a complete, valid coordinate and orientation declaration. |
| Mapped metadata | JSON-compatible values retain their types. Dates, binary values and nonfinite numbers require explicit conversion before import. CSV/TSV metadata remains text. |
| Unmapped columns | Named in the import report and excluded from parts. |

An Excel sheet must have unique, nonempty text headers in its first row. With
one worksheet, omitting `sheet` selects it. With several worksheets, supply its
exact name; hidden worksheets are included in this ambiguity check. Headers and
sheet names are case-sensitive. Completely empty worksheet rows are skipped;
reported row numbers and `id_policy="row"` use logical data-row order.
Values beyond the header fail. A formula or spreadsheet error in a mapped cell
fails with its worksheet and cell address. Save explicit values before importing;
DenseArrays does not evaluate formulas or infer whether cached results are current.

The importer captures one byte snapshot and fingerprints it, then parses mapped
columns through the shared validator. The original byte snapshot and resulting
parts are held in memory during import. Parquet decoding uses
[Arrow record batches](https://arrow.apache.org/docs/python/generated/pyarrow.parquet.ParquetFile.html#pyarrow.parquet.ParquetFile.iter_batches);
Excel uses an explicitly closed
[read-only workbook](https://openpyxl.readthedocs.io/en/stable/optimized.html#read-only-mode).
These reduce parser overhead; they do not make total import memory constant.

Missing optional readers fail before pool publication, with the required install
extra in the error. A saved executable plan still checks its input fingerprint
before running. Inspection of persisted parts uses their saved evidence.
