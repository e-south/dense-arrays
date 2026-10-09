"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/preparation/records.py

Stage accounting for completed or incomplete sampled preparation.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    immutable_json_mapping,
    integer,
    mutable_json,
    object_fields,
)
from dense_arrays.artifacts.preparation.bands import (
    summarize_bands,
    validate_band_report,
)
from dense_arrays.artifacts.preparation.retention import (
    retention_summary,
    validate_retention,
)
from dense_arrays.parts.background import ConstructionReport

if TYPE_CHECKING:
    from collections.abc import Mapping

    from dense_arrays.parts.candidates import Candidate
    from dense_arrays.parts.retention.bands import ScoreBands
    from dense_arrays.parts.retention.pool import PoolSize

COUNTS = {
    "processed",
    "eligibility_rejected",
    "eligible",
    "execution_error",
    "duplicate_discarded",
    "eligible_unique",
    "retained",
    "not_selected",
}


@dataclass(frozen=True)
class PoolAccounting:
    """Keep effort, eligibility, uniqueness and retention counts independently."""

    counts: Mapping[str, int]
    requested_retention: int
    candidate_budget: int
    stop_reason: str
    rejections: Mapping[str, int]
    retention: Mapping[str, object] | None = None
    mining_target: Mapping[str, object] | None = None
    score_bands: Mapping[str, object] | None = None
    construction: ConstructionReport | None = None

    def __post_init__(self) -> None:
        """Reconcile every stage without equating effort exhaustion with success."""
        counts = object_fields(self.counts, COUNTS, "pool counts")
        if set(counts) != COUNTS:
            msg = "incomplete pool stage counts"
            raise ValueError(msg)
        for key, value in counts.items():
            integer(value, field_name=f"counts.{key}", minimum=0)
        integer(self.requested_retention, field_name="requested_retention", minimum=0)
        integer(self.candidate_budget, field_name="candidate_budget", minimum=1)
        if (
            counts["processed"]
            != counts["eligibility_rejected"]
            + counts["eligible"]
            + counts["execution_error"]
            or counts["eligible"]
            != counts["duplicate_discarded"] + counts["eligible_unique"]
            or counts["eligible_unique"] != counts["retained"] + counts["not_selected"]
            or counts["retained"] > self.requested_retention
            or counts["processed"] > self.candidate_budget
        ):
            msg = "pool preparation counts do not reconcile"
            raise ValueError(msg)
        if self.stop_reason not in {
            "candidate_budget",
            "time_budget",
            "execution_error",
            "mining_target",
            "construction_limited",
            "construction_infeasible",
        }:
            msg = "unsupported preparation stop reason"
            raise ValueError(msg)
        if (self.stop_reason == "execution_error") != (counts["execution_error"] > 0):
            msg = "execution error stopping reason and counts must agree"
            raise ValueError(msg)
        if (
            self.stop_reason == "candidate_budget"
            and counts["processed"] != self.candidate_budget
        ):
            msg = "candidate budget stop requires all candidates accounted"
            raise ValueError(msg)
        for value in self.rejections.values():
            integer(value, field_name="rejection count", minimum=1)
            if value > counts["eligibility_rejected"]:
                msg = "rejection count exceeds rejected candidates"
                raise ValueError(msg)
        object.__setattr__(self, "counts", immutable_json_mapping(counts))
        object.__setattr__(self, "rejections", immutable_json_mapping(self.rejections))
        self._validate_retention()
        self._validate_mining_target()
        self._validate_construction()
        if self.score_bands is not None:
            validate_band_report(
                self.score_bands,
                total=counts["eligible_unique"],
                retained=counts["retained"],
            )
            object.__setattr__(
                self, "score_bands", immutable_json_mapping(self.score_bands)
            )

    def _validate_construction(self) -> None:
        if self.construction is not None and not isinstance(
            self.construction, ConstructionReport
        ):
            msg = "conditional construction requires ConstructionReport"
            raise TypeError(msg)
        stopped = self.stop_reason.startswith("construction_")
        failed = (
            self.construction is not None and self.construction.status != "feasible"
        )
        if stopped != failed or (
            failed
            and (
                self.stop_reason != f"construction_{self.construction.status}"
                or self.counts["processed"] != 0
            )
        ):
            msg = "construction outcome disagrees with candidate counts or stop reason"
            raise ValueError(msg)

    def _validate_mining_target(self) -> None:
        if self.mining_target is None:
            if self.stop_reason == "mining_target":
                msg = "mining target stop requires its declared target"
                raise ValueError(msg)
            return
        target = object_fields(
            self.mining_target,
            {"eligible_unique", "minimum_candidates", "met"},
            "mining target",
        )
        if set(target) != {"eligible_unique", "minimum_candidates", "met"}:
            msg = "incomplete mining target outcome"
            raise ValueError(msg)
        for name in ("eligible_unique", "minimum_candidates"):
            integer(target[name], field_name=f"mining_target.{name}", minimum=0)
        met = (
            self.counts["eligible_unique"] >= target["eligible_unique"]
            and self.counts["processed"] >= target["minimum_candidates"]
        )
        if not isinstance(target["met"], bool) or target["met"] != met:
            msg = "mining target outcome disagrees with observed counts"
            raise ValueError(msg)
        if self.stop_reason == "mining_target" and not met:
            msg = "mining target stop requires the target to be met"
            raise ValueError(msg)
        object.__setattr__(self, "mining_target", immutable_json_mapping(target))

    def _validate_retention(self) -> None:
        if self.retention is not None:
            selection = validate_retention(
                self.retention, counts=self.counts, target=self.requested_retention
            )
            object.__setattr__(self, "retention", immutable_json_mapping(selection))

    @property
    def state(self) -> str:
        """Completion requires retention, declared supply and no execution errors."""
        return (
            "completed"
            if self.counts["retained"] == self.requested_retention
            and self.counts["execution_error"] == 0
            and (self.mining_target is None or self.mining_target["met"])
            and (self.construction is None or self.construction.status == "feasible")
            else "incomplete"
        )

    def to_dict(self) -> dict[str, object]:
        """Expose reconciled population counts and the exact stopping condition."""
        return {
            "schema": "dense_arrays.pool_accounting.v1",
            "counts": mutable_json(self.counts),
            "requested_retention": self.requested_retention,
            "candidate_budget": self.candidate_budget,
            "stop_reason": self.stop_reason,
            "rejections": mutable_json(self.rejections),
            **(
                {"construction": self.construction.to_dict()}
                if self.construction is not None
                else {}
            ),
            **(
                {"score_bands": mutable_json(self.score_bands)}
                if self.score_bands is not None
                else {}
            ),
            **(
                {"mining_target": mutable_json(self.mining_target)}
                if self.mining_target is not None
                else {}
            ),
            **(
                {"retention": mutable_json(self.retention)}
                if self.retention is not None
                else {}
            ),
        }

    @classmethod
    def from_dict(cls, value: object) -> PoolAccounting:
        """Read complete counts without inventing omitted stage evidence."""
        keys = {
            "schema",
            "counts",
            "requested_retention",
            "candidate_budget",
            "stop_reason",
            "rejections",
        }
        data = object_fields(
            value,
            keys | {"retention", "mining_target", "score_bands", "construction"},
            "pool accounting",
        )
        retention = data.pop("retention", None)
        mining_target = data.pop("mining_target", None)
        score_bands = data.pop("score_bands", None)
        construction = data.pop("construction", None)
        if set(data) != keys or data.pop("schema") != "dense_arrays.pool_accounting.v1":
            msg = "unsupported or incomplete pool accounting"
            raise ValueError(msg)
        return cls(
            **data,
            retention=retention,
            mining_target=mining_target,
            score_bands=score_bands,
            construction=None
            if construction is None
            else ConstructionReport.from_dict(construction),
        )


def recount(  # noqa: PLR0913 - independent effort, retention and supply contracts
    candidates: tuple[Candidate, ...],
    *,
    target: int,
    budget: int,
    stop_reason: str,
    mmr: bool = False,
    pool_sizing: PoolSize | None = None,
    mining_target: Mapping[str, int] | None = None,
    score_bands: ScoreBands | None = None,
    scoring_id: str | None = None,
    construction: ConstructionReport | None = None,
) -> PoolAccounting:
    """Recount persisted decisions without invoking sampling or scoring."""
    if pool_sizing is not None and not mmr:
        msg = "pool sizing requires MMR retention"
        raise ValueError(msg)
    rejected = sum(bool(c.reasons) and c.error is None for c in candidates)
    errors = sum(c.error is not None for c in candidates)
    eligible = len(candidates) - rejected - errors
    unique = sum(c.representative == c.index for c in candidates)
    retained = sum(c.retained for c in candidates)
    return PoolAccounting(
        {
            "processed": len(candidates),
            "eligibility_rejected": rejected,
            "eligible": eligible,
            "execution_error": errors,
            "duplicate_discarded": eligible - unique,
            "eligible_unique": unique,
            "retained": retained,
            "not_selected": unique - retained,
        },
        target,
        budget,
        stop_reason,
        dict(
            Counter(
                reason for c in candidates if c.error is None for reason in c.reasons
            )
        ),
        retention_summary(candidates, target=target, sizing=pool_sizing)
        if mmr
        else None,
        None
        if mining_target is None
        else {
            **mining_target,
            "met": unique >= mining_target["eligible_unique"]
            and len(candidates) >= mining_target["minimum_candidates"],
        },
        None
        if score_bands is None
        else summarize_bands(candidates, score_bands, scoring_id=scoring_id),
        construction,
    )
