"""Parse the pgsc_calc match summary (`{sampleset}_summary.csv`).

Each row is one (score, match_status, flag combination) with a variant count; the counts for a
score sum to its variant total.
"""

import csv
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from pgs_mcp.qc.errors import QCInputError
from pgs_mcp.qc.models import VariantMetrics

REQUIRED_COLUMNS = {
    "dataset",
    "accession",
    "score_pass",
    "match_status",
    "ambiguous",
    "is_multiallelic",
    "duplicate_best_match",
    "duplicate_ID",
    "match_flipped",
    "count",
}
EXCLUSION_FLAGS = ("ambiguous", "is_multiallelic", "duplicate_best_match", "duplicate_ID")

_PGS_ID = re.compile(r"^(PGS\d{6})")
_BUILD = re.compile(r"(GRCh3[78])")


def parse_accession(accession: str) -> tuple[str, str | None]:
    """`PGS000004_hmPOS_GRCh37` -> (`PGS000004`, `GRCh37`). Custom scores keep their name."""
    pgs_id = _PGS_ID.match(accession)
    build = _BUILD.search(accession)
    return (pgs_id.group(1) if pgs_id else accession, build.group(1) if build else None)


@dataclass
class SummaryScore:
    dataset: str
    accession: str
    pipeline_pass: bool
    n_matched: int = 0
    n_matched_ambiguous: int = 0
    n_ambiguous: int = 0
    n_flipped: int = 0
    n_excluded: int = 0
    n_unmatched: int = 0
    excluded_by_reason: Counter[str] = field(default_factory=Counter)
    matched_multiallelic: bool = False
    excluded_multiallelic: bool = False
    excluded_ambiguous: bool = False

    @property
    def n_variants(self) -> int:
        return self.n_matched + self.n_excluded + self.n_unmatched

    def metrics(self) -> VariantMetrics:
        def pct(n: int) -> float:
            return 100 * n / self.n_variants if self.n_variants else 0.0

        return VariantMetrics(
            matched_pct=pct(self.n_matched),
            matched_unambiguous_pct=pct(self.n_matched - self.n_matched_ambiguous),
            ambiguous_pct=pct(self.n_ambiguous),
            flipped_pct=pct(self.n_flipped),
            excluded_pct=pct(self.n_excluded),
            unmatched_pct=pct(self.n_unmatched),
            excluded_by_reason=dict(self.excluded_by_reason),
        )


def read_summary(path: Path) -> list[SummaryScore]:
    """Read a match summary into one SummaryScore per accession, in file order."""
    try:
        with path.open(newline="") as f:
            reader = csv.DictReader(f)
            missing = REQUIRED_COLUMNS - set(reader.fieldnames or ())
            if missing:
                raise QCInputError(
                    f"{path} is not a pgsc_calc match summary; missing columns: "
                    f"{', '.join(sorted(missing))}"
                )
            scores: dict[str, SummaryScore] = {}
            for row in reader:
                _add_row(scores, row)
    except OSError as e:
        raise QCInputError(f"Cannot read match summary {path}: {e}") from e

    if not scores:
        raise QCInputError(f"Match summary {path} contains no scores.")
    datasets = {s.dataset for s in scores.values()}
    if len(datasets) > 1:
        raise QCInputError(
            f"Match summary {path} mixes datasets ({', '.join(sorted(datasets))}); "
            "pass one sampleset's summary at a time."
        )
    return list(scores.values())


def _add_row(scores: dict[str, SummaryScore], row: dict[str, str]) -> None:
    accession = row["accession"]
    score = scores.get(accession)
    if score is None:
        score = scores[accession] = SummaryScore(
            dataset=row["dataset"],
            accession=accession,
            pipeline_pass=row["score_pass"].lower() == "true",
        )

    try:
        count = int(row["count"])
    except ValueError as e:
        raise QCInputError(f"Non-integer count {row['count']!r} for {accession}") from e

    flags = {k: row[k].lower() == "true" for k in EXCLUSION_FLAGS + ("match_flipped",)}
    status = row["match_status"]
    if status == "matched":
        score.n_matched += count
        score.n_flipped += count * flags["match_flipped"]
        if flags["ambiguous"]:
            score.n_matched_ambiguous += count
            score.n_ambiguous += count
        score.matched_multiallelic |= flags["is_multiallelic"]
    elif status == "excluded":
        score.n_excluded += count
        score.n_ambiguous += count * flags["ambiguous"]
        for flag in EXCLUSION_FLAGS:
            if flags[flag]:
                score.excluded_by_reason[flag] += count
        score.excluded_ambiguous |= flags["ambiguous"]
        score.excluded_multiallelic |= flags["is_multiallelic"]
    elif status == "unmatched":
        score.n_unmatched += count
    else:
        raise QCInputError(f"Unexpected match_status {status!r} for {accession}")
