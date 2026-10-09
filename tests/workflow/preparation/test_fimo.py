"""Optional scoring binds its inputs before any candidate is evaluated.

Author: Eric J. South.
"""

import json
from pathlib import Path

import pytest

from dense_arrays import parts
from dense_arrays.parts.motifs import read_artifact

from .test_motifs import source


def executable(tmp_path: Path, body: str) -> Path:
    path = tmp_path / "fimo fixture"
    path.write_text("#!/bin/sh\n" + body + "\n")
    path.chmod(0o700)
    return path


def test_preflight_queries_only_version_and_binds_input_bytes(tmp_path: Path):
    from dense_arrays.parts.scoring import bind_fimo  # noqa: PLC0415

    log = tmp_path / "calls"
    tool = executable(
        tmp_path,
        f'printf "%s\\n" "$*" >> "{log}"\n'
        'if [ "$1" = "--version" ]; then printf "5.5.9\\n"; else exit 92; fi',
    )
    motif = read_artifact(source(tmp_path)).motif
    config = parts.FimoScoring(executable=tool, hit_pvalue_max=0.1)
    binding = bind_fimo(motif, config)
    assert binding.version == "5.5.9"
    assert binding.background == (0.25,) * 4
    assert binding.effective_background == (0.25,) * 4
    assert binding.motif == motif
    assert log.read_text() == "--version\n"
    assert json.loads(json.dumps(binding.to_dict())) == binding.to_dict()
    binding.verify()
    tool.write_text("#!/bin/sh\nexit 0\n")
    with pytest.raises(ValueError, match="changed"):
        binding.verify()


def test_missing_fimo_fails_during_preflight(tmp_path: Path):
    from dense_arrays.parts.scoring import ScoringError, bind_fimo  # noqa: PLC0415

    with pytest.raises(ScoringError) as caught:
        bind_fimo(
            read_artifact(source(tmp_path)).motif,
            parts.FimoScoring(executable=tmp_path / "absent"),
        )
    assert caught.value.reason == "unavailable"


@pytest.mark.parametrize(
    "fields",
    [
        {"hit_pvalue_max": 0},
        {"hit_pvalue_max": 1.1},
        {"hit_pvalue_max": True},
        {"hit_pvalue_max": float("nan")},
        {"pseudocount": -1},
        {"pseudocount": float("inf")},
        {"strands": "both"},
    ],
)
def test_scoring_options_reject_ambiguous_values(fields: dict):
    with pytest.raises((TypeError, ValueError)):
        parts.FimoScoring(**fields)


def test_background_is_explicit_bound_and_strand_policy_is_visible(tmp_path: Path):
    from dense_arrays.parts.scoring import bind_fimo  # noqa: PLC0415

    tool = executable(tmp_path, 'printf "5.5.9\\n"')
    background = tmp_path / "background.txt"
    background.write_text("# zero order\nA 0.1 C 0.2\nG 0.3 T 0.4\n")
    binding = bind_fimo(
        read_artifact(source(tmp_path)).motif,
        parts.FimoScoring(executable=tool, background=background, strands="double"),
    )
    assert binding.background == (0.1, 0.2, 0.3, 0.4)
    assert binding.effective_background == (0.25,) * 4
    assert binding.to_dict()["background_policy"] == "complement_average.v1"
    background.write_text("A 0.25 C 0.25 G 0.25 T 0.25\n")
    with pytest.raises(ValueError, match=r"background.*changed"):
        binding.verify()


@pytest.mark.parametrize(
    "text", ["A 0.5 C 0.5", "A 0.2 A 0.3 C 0.2 G 0.1 T 0.2", "AA 1", "A 0 C 0 G 0 T 1"]
)
def test_invalid_background_is_rejected_before_querying_tool(tmp_path: Path, text: str):
    from dense_arrays.parts.scoring import bind_fimo  # noqa: PLC0415

    background = tmp_path / "background.txt"
    background.write_text(text)
    with pytest.raises(ValueError, match="background"):
        bind_fimo(
            read_artifact(source(tmp_path)).motif,
            parts.FimoScoring(executable=tmp_path / "missing", background=background),
        )


@pytest.mark.parametrize("kind", ["timeout", "output_limit", "backend", "malformed"])
def test_preflight_execution_failures_remain_distinct(tmp_path: Path, kind: str):
    from dense_arrays.parts.scoring import (  # noqa: PLC0415
        ScoringError,
        ScoringLimits,
        bind_fimo,
    )

    body = {
        "timeout": "exec sleep 10",
        "output_limit": 'while true; do printf "excess output\\n"; done',
        "backend": 'printf "unavailable database\\n" >&2; exit 17',
        "malformed": 'printf "not a version\\n"',
    }[kind]
    tool = executable(tmp_path, body)
    settings = parts.FimoScoring(
        executable=tool,
        limits=ScoringLimits(
            seconds=0.1 if kind == "timeout" else 5, windows=100, output_bytes=1024
        ),
    )
    with pytest.raises(ScoringError) as caught:
        bind_fimo(read_artifact(source(tmp_path)).motif, settings)
    assert caught.value.reason == kind
    if kind == "backend":
        assert "unavailable database" in str(caught.value)


HEADER = (
    "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\t"
    "p-value\tq-value\tmatched_sequence\n"
)
REFERENCE = "motif\t\tcandidate_0\t1\t2\t+\t3\t0.0625\t\tAC\n"


def scorer(tmp_path: Path, rows: str) -> Path:
    return executable(
        tmp_path,
        'if [ "$1" = "--version" ]; then printf "5.5.9\\n"; exit; fi\n'
        'case "$*" in\n'
        "*reference.fa*) cat <<'EOF'\n" + HEADER + REFERENCE + "EOF\n;;\n"
        "*) cat <<'EOF'\n" + HEADER + rows + "EOF\n;;\nesac",
    )


def test_scoring_uses_validated_geometry_and_order_independent_ties(tmp_path: Path):
    from dense_arrays.parts.scoring import bind_fimo, scan_fimo  # noqa: PLC0415

    rows = (
        "motif\t\tcandidate_0\t5\t6\t+\t3\t0.0625\t\tAC\n"
        "motif\t\tcandidate_0\t3\t4\t-\t3\t0.0625\t\tAC\n"
    )
    binding = bind_fimo(
        read_artifact(source(tmp_path)).motif,
        parts.FimoScoring(executable=scorer(tmp_path, rows), hit_pvalue_max=0.1),
    )
    result = scan_fimo(binding, ("TTGTAC", "TTTTTT"))
    hit, absent = result.hits
    assert (hit.start, hit.end, hit.strand, hit.core) == (2, 4, "reverse", "AC")
    assert hit.raw == 3
    assert hit.per_base == 1.5
    assert hit.pvalue == 0.0625
    assert hit.theoretical_max == 3
    assert hit.fraction_of_max == 1
    assert hit.units == "fimo_log2_odds"
    assert absent is None
    assert result.processed == 2
    assert result.candidate_windows == 20
    assert result.calibration_windows == 2
    assert result.reported_hits == 2
    assert result.binding_id == binding.binding_id
    assert type(hit).from_dict(hit.to_dict()) == hit


@pytest.mark.parametrize(
    "replacement",
    [
        {"motif_id": "other"},
        {"sequence_name": "candidate_9"},
        {"start": "0"},
        {"stop": "20"},
        {"strand": "?"},
        {"score": "nan"},
        {"score": "4"},
        {"p-value": "-1"},
        {"matched_sequence": "GT"},
        {"q-value": "0.02"},
    ],
)
def test_invalid_backend_evidence_is_never_candidate_rejection(
    tmp_path: Path, replacement: dict
):
    from dense_arrays.parts.scoring import (  # noqa: PLC0415
        ScoringError,
        bind_fimo,
        scan_fimo,
    )

    columns = HEADER.strip().split("\t")
    row = dict(zip(columns, REFERENCE.rstrip("\n").split("\t"), strict=True))
    row.update(replacement)
    rows = "\t".join(row[k] for k in columns) + "\n"
    binding = bind_fimo(
        read_artifact(source(tmp_path)).motif,
        parts.FimoScoring(executable=scorer(tmp_path, rows), hit_pvalue_max=0.1),
    )
    with pytest.raises(ScoringError) as caught:
        scan_fimo(binding, ("AC",))
    assert caught.value.reason == "malformed"


def test_candidate_work_is_admitted_before_scoring(tmp_path: Path):
    from dense_arrays.parts.scoring import bind_fimo, scan_fimo  # noqa: PLC0415

    marker = tmp_path / "scored"
    tool = executable(
        tmp_path,
        f'if [ "$1" = "--version" ]; then printf "5.5.9\\n"; else touch "{marker}"; fi',
    )
    binding = bind_fimo(
        read_artifact(source(tmp_path)).motif,
        parts.FimoScoring(executable=tool, limits=parts.ScoringLimits(windows=4)),
    )
    for sequences in (("AAAA",), ("NAAA",), ("A",), ("aa",)):
        with pytest.raises(ValueError, match=r"window limit|uppercase|length"):
            scan_fimo(binding, sequences)
    assert not marker.exists()


def test_forward_wins_exact_strand_tie_independent_of_output_order(tmp_path: Path):
    from dense_arrays.parts.scoring import bind_fimo, scan_fimo  # noqa: PLC0415

    rows = (
        "motif\t\tcandidate_0\t1\t2\t-\t0\t0.0625\t\tAT\n"
        "motif\t\tcandidate_0\t1\t2\t+\t0\t0.0625\t\tAT\n"
    )
    binding = bind_fimo(
        read_artifact(source(tmp_path)).motif,
        parts.FimoScoring(executable=scorer(tmp_path, rows), hit_pvalue_max=0.1),
    )
    assert scan_fimo(binding, ("AT",)).hits[0].strand == "forward"


def test_binding_and_result_round_trip_without_invoking_tools(tmp_path: Path):
    from dense_arrays.parts.scoring import (  # noqa: PLC0415
        FimoBinding,
        FimoResult,
        bind_fimo,
        scan_fimo,
    )

    binding = bind_fimo(
        read_artifact(source(tmp_path)).motif,
        parts.FimoScoring(executable=scorer(tmp_path, REFERENCE), hit_pvalue_max=0.1),
    )
    result = scan_fimo(binding, ("AC",))
    binding.executable.unlink()
    assert FimoBinding.from_dict(binding.to_dict()) == binding
    assert FimoResult.from_dict(result.to_dict()) == result
    bad = result.to_dict()
    bad["processed"] = 99
    with pytest.raises(ValueError, match="disagree"):
        FimoResult.from_dict(bad)
    bad_binding = binding.to_dict()
    bad_binding["effective_background"] = [0.1, 0.2, 0.3, 0.4]
    with pytest.raises(ValueError, match="disagree"):
        FimoBinding.from_dict(bad_binding)


def test_repeated_backend_hit_is_malformed_not_extra_evidence(tmp_path: Path):
    from dense_arrays.parts.scoring import (  # noqa: PLC0415
        ScoringError,
        bind_fimo,
        scan_fimo,
    )

    binding = bind_fimo(
        read_artifact(source(tmp_path)).motif,
        parts.FimoScoring(
            executable=scorer(tmp_path, REFERENCE * 2), hit_pvalue_max=0.1
        ),
    )
    with pytest.raises(ScoringError, match="duplicate"):
        scan_fimo(binding, ("ACAC",))


def test_combined_output_limit_includes_calibration_stderr(tmp_path: Path):
    from dense_arrays.parts.scoring import (  # noqa: PLC0415
        ScoringError,
        bind_fimo,
        scan_fimo,
    )

    tool = scorer(tmp_path, REFERENCE)
    text = tool.read_text().replace(
        "*reference.fa*) cat", "*reference.fa*) printf '" + "e" * 500 + "' >&2; cat"
    )
    tool.write_text(text)
    settings = parts.FimoScoring(
        executable=tool,
        limits=parts.ScoringLimits(output_bytes=750),
        hit_pvalue_max=0.1,
    )
    binding = bind_fimo(read_artifact(source(tmp_path)).motif, settings)
    with pytest.raises(ScoringError) as caught:
        scan_fimo(binding, ("AC",))
    assert caught.value.reason == "output_limit"


def test_pinned_fimo_fixture_preserves_scores_and_earliest_position(tmp_path: Path):
    from dense_arrays.parts.scoring import bind_fimo, scan_fimo  # noqa: PLC0415

    fixture_path = Path(__file__).parents[2] / "fixtures/workflow/densegen-fimo-v1.json"
    fixture = json.loads(fixture_path.read_text())
    artifact_path = tmp_path / "motif.json"
    artifact_path.write_text(json.dumps(fixture["artifact"]))
    tool = scorer(tmp_path, fixture["tsv"].split("\n", 1)[1])
    binding = bind_fimo(
        read_artifact(parts.PWMArtifact(artifact_path)).motif,
        parts.FimoScoring(executable=tool, hit_pvalue_max=fixture["hit_pvalue_max"]),
    )
    result = scan_fimo(binding, fixture["sequences"])
    for hit, expected in zip(result.hits, fixture["expected"], strict=True):
        if expected is None:
            assert hit is None
        else:
            assert hit.start == expected["start"] - 1
            assert hit.end == expected["stop"]
            assert hit.strand == {"+": "forward", "-": "reverse"}[expected["strand"]]
            assert hit.core == expected["matched_sequence"]
            assert hit.raw == expected["score"]


def test_maximum_reference_uses_ratios_without_overflow(tmp_path: Path):
    from dataclasses import replace  # noqa: PLC0415

    from dense_arrays.parts.scoring import bind_fimo, scan_fimo  # noqa: PLC0415

    tool = scorer(tmp_path, "")
    text = tool.read_text().replace(REFERENCE, REFERENCE.replace("\tAC\n", "\tCA\n"))
    tool.write_text(text)
    motif = replace(
        read_artifact(source(tmp_path)).motif,
        background=(1e-309, 1e-320, 0.5, 0.5),
        probabilities=((0.5, 0.5, 0, 0), (1, 0, 0, 0)),
    )
    binding = bind_fimo(motif, parts.FimoScoring(executable=tool, strands="single"))
    assert scan_fimo(binding, ("CA",)).theoretical_max == 3
