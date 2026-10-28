"""Output schema for match QC.

Every metric is reported raw. Warnings flag structural problems (mostly categorical, no
cutoff); observations are mechanical, descriptive statements. Neither is a verdict.
Percentages are 0-100 and always relative to the score's own variant count.
"""

from typing import Any

from pydantic import BaseModel, Field


class QCWarning(BaseModel):
    code: str = Field(description="Stable machine-readable code, e.g. CHROM_ABSENT.")
    message: str
    evidence: dict[str, Any] = Field(default_factory=dict)


class ChromCoverage(BaseModel):
    chrom: str
    n_total: int
    n_matched: int
    n_excluded: int = Field(description="Found in the target but excluded from scoring.")


class VariantMetrics(BaseModel):
    """From the match summary: share of the score's variants in each state."""

    matched_pct: float
    matched_unambiguous_pct: float
    ambiguous_pct: float = Field(description="Strand-ambiguous variants, matched or excluded.")
    flipped_pct: float = Field(description="Matched only after flipping strand.")
    excluded_pct: float
    unmatched_pct: float
    excluded_by_reason: dict[str, int] = Field(
        default_factory=dict,
        description="Excluded variant counts per flag; a variant can carry several flags.",
    )


class WeightMetrics(BaseModel):
    weight_coverage_pct: float | None = Field(
        description="Sum of |effect_weight| over matched variants / sum over all variants."
    )
    max_abs_weight: float | None
    largest_unmatched_abs_weight: float | None = Field(
        description="Largest |effect_weight| among unmatched variants; null if none unmatched."
    )


class LogMetrics(BaseModel):
    """From the variant-level match log; absent when no log is supplied."""

    weights: WeightMetrics
    match_types: dict[str, int] = Field(description="Matched variant counts per match_type.")
    n_not_best_candidates: int = Field(
        description="Extra candidate matches discarded in favour of a better one."
    )
    snv_total: int
    snv_matched: int
    indel_total: int = Field(description="Indels and other non-SNV variants.")
    indel_matched: int
    n_missing_position: int
    chromosomes: list[ChromCoverage]


class ScoreQC(BaseModel):
    pgs_id: str
    accession: str
    genome_build: str | None
    n_variants: int
    pipeline_pass: bool = Field(description="pgsc_calc's own score_pass verdict.")
    variants: VariantMetrics
    log: LogMetrics | None
    warnings: list[QCWarning]
    observations: list[str]


class ScoreIndexRow(BaseModel):
    """Compact one-line view of a score, always returned for every score."""

    pgs_id: str
    n_variants: int
    pipeline_pass: bool
    matched_pct: float
    weight_coverage_pct: float | None
    warning_codes: list[str]


class InferredSettings(BaseModel):
    min_overlap_lower: float | None = Field(
        description="Highest match rate (fraction) among failing scores: min_overlap is above it."
    )
    min_overlap_upper: float | None = Field(
        description="Lowest match rate (fraction) among passing scores: min_overlap is at most it."
    )
    keep_ambiguous: bool | None
    keep_multiallelic: bool | None


class RunQC(BaseModel):
    dataset: str
    genome_builds: list[str]
    n_scores: int
    n_pipeline_pass: int
    n_pipeline_fail: int
    inferred_settings: InferredSettings
    log_used: bool
    chromosomes: list[ChromCoverage] | None = Field(
        description="Run-wide coverage over score-variant pairs; a shared variant counts per score."
    )
    warnings: list[QCWarning]
    score_index: list[ScoreIndexRow]
    scores: list[ScoreQC] = Field(description="Full detail for selected scores.")
    detail_note: str | None = Field(
        default=None, description="Set when detail was limited; explains how to get more."
    )
