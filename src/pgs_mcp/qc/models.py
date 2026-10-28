"""Output schema for match QC.

Every metric is reported raw. Warnings flag structural problems (mostly categorical, no
cutoff); observations are mechanical, descriptive statements. Neither is a verdict.
Percentages are 0-100 and always relative to the score's own variant count.

Values are computed at full precision; rounding happens only when serialized.
"""

from typing import Annotated, Any, NamedTuple

from pydantic import BaseModel, Field, PlainSerializer

Pct = Annotated[float, PlainSerializer(lambda v: round(v, 2), return_type=float)]
Fraction = Annotated[float, PlainSerializer(lambda v: round(v, 4), return_type=float)]
# Effect weights span orders of magnitude, so keep significant figures, not decimals.
Weight = Annotated[float, PlainSerializer(lambda v: float(f"{v:.4g}"), return_type=float)]

N_TOP_UNMATCHED = 3


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

    matched_pct: Pct
    matched_unambiguous_pct: Pct
    ambiguous_pct: Pct = Field(description="Strand-ambiguous variants, matched or excluded.")
    flipped_pct: Pct = Field(description="Matched only after flipping strand.")
    excluded_pct: Pct
    unmatched_pct: Pct
    excluded_by_reason: dict[str, int] = Field(
        default_factory=dict,
        description="Excluded variant counts per flag; a variant can carry several flags.",
    )


class UnmatchedVariant(BaseModel):
    chrom: str | None = Field(description="Null when the scoring file gives no position.")
    position: int | None
    effect_allele: str | None
    other_allele: str | None
    abs_weight: Weight

    @property
    def label(self) -> str:
        """`6:32658079 A/G` (effect/other), `no position, A/G`, or `no position, effect HLA-...`."""
        if self.other_allele:
            alleles = f"{self.effect_allele}/{self.other_allele}"
        else:  # effect-only entries such as HLA haplotypes may contain slashes themselves
            alleles = f"effect {self.effect_allele}"
        if self.chrom is None or self.position is None:
            return f"no position, {alleles}"
        return f"{self.chrom}:{self.position} {alleles}"


class WeightMetrics(BaseModel):
    weight_coverage_pct: Pct | None = Field(
        description="Sum of |effect_weight| over matched variants / sum over all variants."
    )
    total_abs_weight: Weight = Field(description="Sum of |effect_weight| over all variants.")
    max_abs_weight: Weight | None
    largest_unmatched: list[UnmatchedVariant] = Field(
        description=f"Up to {N_TOP_UNMATCHED} unmatched variants with the largest |effect_weight|, "
        "largest first."
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
    indel_total: int = Field(
        description="Indels and other non-SNV entries (e.g. HLA alleles, presence/absence codes)."
    )
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


class ScoreIndexRow(NamedTuple):
    """One score in the index. Serialized as an array, in this field order."""

    pgs_id: str
    n_variants: int
    pipeline_pass: bool
    matched_pct: Pct
    weight_coverage_pct: Pct | None
    warning_codes: list[str]


class ScoreIndex(BaseModel):
    """Every reported score, one row each, sent as a table to keep large runs small."""

    columns: list[str] = Field(default=list(ScoreIndexRow._fields))
    rows: list[ScoreIndexRow]


class InferredSettings(BaseModel):
    min_overlap_lower: Fraction | None = Field(
        description="Highest match rate (fraction) among failing scores: min_overlap is above it."
    )
    min_overlap_upper: Fraction | None = Field(
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
        description="Run-wide coverage over score-variant pairs; a shared variant counts per "
        "score. Null without a log, or when specific scores were requested."
    )
    warnings: list[QCWarning]
    score_index: ScoreIndex
    scores: list[ScoreQC] = Field(description="Full detail for selected scores.")
    detail_note: str | None = Field(
        default=None, description="Set when detail was limited; explains how to get more."
    )
