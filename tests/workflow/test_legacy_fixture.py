"""Compare a pinned DenseGen input/packing fixture without importing DNADesign.

Author: Eric J. South.
"""

import json
from pathlib import Path

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.parts.ingestion import read_parts

FIXTURE = Path(__file__).parents[1] / "fixtures/workflow/densegen-curated-v1.json"


def test_curated_import_and_objective_match_pinned_producer(tmp_path: Path):
    legacy = json.loads(FIXTURE.read_text())
    source = tmp_path / "parts.csv"
    source.write_text(legacy["input_csv"])
    imported = read_parts(parts.PartTable(source, "csv"))
    assert [(p.part_id, p.sequence, p.group, p.source) for p in imported.parts] == [
        (row["site_id"], row["tfbs"], row["tf"], row["source"])
        for row in legacy["normalized_rows"]
    ]
    assert imported.parts[0].part_id != imported.parts[3].part_id
    assert imported.parts[0].core_start is None  # Unknown, not a synthetic PWM hit.
    result = da.run(
        planning.DesignSpec(
            parts=imported.parts,
            length=planning.Length(maximum=legacy["packing"]["maximum"]),
            strands=legacy["packing"]["strands"],
            requirements=[planning.GroupCoverage(id="both", groups=("A", "B"), min=2)],
        ),
        out=tmp_path / "native",
    )
    assert da.inspect(result, verify=True).state == "completed"
    with da.inspect(result, view="designs").records() as records:
        design = next(records)
    assert len(design.realized.placements) == legacy["packing"]["optimal_occurrences"]
    assert {placement.label for placement in design.realized.placements} == {"A", "B"}


def test_matrix_pairing_and_allocations_match_pinned_functions():
    fixture = json.loads(FIXTURE.with_name("densegen-matrix-v1.json").read_text())
    axes = {
        axis: {
            name: planning.Variant(parts=[parts.Part(axis, sequence)])
            for name, sequence in options.items()
        }
        for axis, options in fixture["motif_sets"].items()
    }
    base = planning.DesignSpec(
        [parts.Part("up", "TTGACA"), parts.Part("down", "TATAAT")],
        planning.Length(maximum=12),
        strands="single",
    )
    for case in fixture["cases"]:
        pairing = "explicit" if case["pairing"] == "explicit_pairs" else case["pairing"]
        pairs = tuple({"up": c[0], "down": c[2]} for c in case["combinations"])
        zeros = tuple(
            f"down={c[2]},up={c[0]}"
            for c, count in zip(
                case["combinations"], case["quotas_before_zero_filter"], strict=True
            )
            if count == 0
        )
        result = da.plan(
            planning.MatrixSpec(
                base=base,
                axes=axes,
                max_cells=4,
                pairing=pairing,
                pairs=pairs if pairing == "explicit" else (),
                allocation=planning.Allocation(
                    total=case["total"],
                    policy="balanced",
                    zero_cells=zeros,
                ),
            )
        )
        assert [c.target for c in result.cells] == case["quotas_before_zero_filter"]
        assert [
            [
                c.choices["up"],
                c.plan.request.parts[0].sequence,
                c.choices["down"],
                c.plan.request.parts[1].sequence,
            ]
            for c in result.cells
        ] == case["combinations"]
