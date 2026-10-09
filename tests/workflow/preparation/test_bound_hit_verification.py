"""Persisted FIMO observations must fit their declared motif and strand scope.

Author: Eric J. South.
"""

from pathlib import Path

import pytest

from dense_arrays import parts, planning
from dense_arrays.artifacts.preparation.records import recount
from dense_arrays.artifacts.preparation.verification import verify_decisions
from dense_arrays.parts.candidates import Candidate
from dense_arrays.parts.motifs import Motif, MotifInput
from dense_arrays.parts.retention.selection import select_candidates
from dense_arrays.parts.scoring import FimoBinding, FimoHit
from dense_arrays.parts.screening import ScreenObservation
from dense_arrays.planning.preparation.motifs import MotifSource
from dense_arrays.planning.preparation.sampled import SampledPreparation
from dense_arrays.planning.preparation.screening import BoundExclusion


def _verify_saved_hit(role: str, *, strands: str, strand: str, width: int = 1) -> None:
    """Verify complete synthetic decisions without motif files or a scorer."""
    motif = Motif("m", ((0.7, 0.1, 0.1, 0.1),), (0.25,) * 4, None, "test")
    scoring = parts.FimoScoring(executable=Path("fimo"), strands=strands)
    artifact = parts.PWMArtifact("motif.json")
    binding = FimoBinding(
        motif, scoring, Path("fimo"), "0" * 64, "5.5.9", (0.25,) * 4, None, None
    )
    source = MotifSource(MotifInput(motif, artifact.path, "0" * 64), binding)
    # Reported p-values are rounded; the backend decides whether a hit qualifies.
    hit = FimoHit(
        0,
        width,
        strand,
        ("A" if strand == "reverse" else "T") * width,
        1,
        scoring.hit_pvalue_max + 1e-10,
        2,
    )
    exclusion = parts.PWMExclusion("exclude", (artifact,), scoring)
    request = parts.PreparationSpec(
        source=artifact if role == "primary" else parts.Background(),
        sampling=parts.Sampling(planning.Length(exact=2)),
        budget=parts.CandidateBudget(1),
        scoring=scoring if role == "primary" else None,
        retain=parts.Retention(count=1, policy="first_eligible"),
        screening=() if role == "primary" else (exclusion,),
    )
    resolved = SampledPreparation(
        request,
        source if role == "primary" else request.source,
        () if role == "primary" else (BoundExclusion("exclude", (source,)),),
    )
    part = parts.Part(
        "candidate_1",
        "TT",
        group="m" if role == "primary" else "background",
        source="pwm_artifact" if role == "primary" else "background",
        **(
            {
                "core_start": hit.start,
                "core_end": hit.end,
                "core_orientation": hit.strand,
                "metadata": {"score": hit.to_dict()},
            }
            if role == "primary"
            else {}
        ),
    )
    candidate = Candidate(
        1,
        part,
        reasons=() if role == "primary" else ("exclude",),
        screening=()
        if role == "primary"
        else (ScreenObservation("exclude", binding.binding_id, hit),),
    )
    selected = select_candidates(
        (candidate,), request.uniqueness, request.retain, motif=motif
    )
    recorded = recount(selected, target=1, budget=1, stop_reason="candidate_budget")
    verify_decisions(selected, resolved, recorded)


@pytest.mark.parametrize("role", ["primary", "exclusion"])
def test_single_strand_binding_rejects_saved_reverse_hits(role: str):
    with pytest.raises(ValueError, match="strand"):
        _verify_saved_hit(role, strands="single", strand="reverse")


@pytest.mark.parametrize("role", ["primary", "exclusion"])
@pytest.mark.parametrize(
    "strands,strand",
    [("single", "forward"), ("double", "forward"), ("double", "reverse")],
)
def test_saved_hits_preserve_strand_scope_and_backend_threshold_authority(
    role: str, strands: str, strand: str
):
    _verify_saved_hit(role, strands=strands, strand=strand)


@pytest.mark.parametrize("role", ["primary", "exclusion"])
def test_saved_hit_width_must_fit_its_bound_motif(role: str):
    with pytest.raises(ValueError, match="motif"):
        _verify_saved_hit(role, strands="double", strand="reverse", width=2)
