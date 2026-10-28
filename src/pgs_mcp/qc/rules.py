"""Warnings (structural problems) and observations (descriptive statements) from raw metrics.

No rule here grades a score. Warnings fire on categorical conditions except the two named in
thresholds.py, whose cutoff is always included in the evidence.
"""

import statistics

from pgs_mcp.qc import thresholds
from pgs_mcp.qc.log import LogSummary
from pgs_mcp.qc.models import ChromCoverage, LogMetrics, QCWarning, VariantMetrics


def absent_chromosomes(chromosomes: list[ChromCoverage]) -> list[ChromCoverage]:
    """Chromosomes where no score variant was found in the target (matched or excluded)."""
    return [c for c in chromosomes if c.n_total > 0 and c.n_matched + c.n_excluded == 0]


def _fmt_chroms(chroms: list[ChromCoverage]) -> str:
    return ", ".join(f"{c.chrom} ({c.n_total:,})" for c in chroms)


def score_warnings(
    n_variants: int,
    variants: VariantMetrics,
    log: LogMetrics | None,
    keep_ambiguous: bool | None,
    run_absent: set[str],
) -> list[QCWarning]:
    warnings: list[QCWarning] = []

    if n_variants < thresholds.SMALL_SCORE_MAX_VARIANTS:
        warnings.append(
            QCWarning(
                code="SMALL_SCORE",
                message=f"Score has only {n_variants} variants; each one carries substantial "
                "weight, so single missing variants shift scores noticeably.",
                evidence={
                    "n_variants": n_variants,
                    "threshold": thresholds.SMALL_SCORE_MAX_VARIANTS,
                },
            )
        )

    duplicates = {
        k: v
        for k, v in variants.excluded_by_reason.items()
        if k in ("duplicate_ID", "duplicate_best_match")
    }
    if duplicates:
        warnings.append(
            QCWarning(
                code="DUPLICATES_EXCLUDED",
                message="Variants were excluded because they matched duplicated target IDs or "
                "had duplicate best matches.",
                evidence=duplicates,
            )
        )

    if keep_ambiguous and variants.ambiguous_pct > 0:
        warnings.append(
            QCWarning(
                code="AMBIGUOUS_RETAINED",
                message=f"{variants.ambiguous_pct:.1f}% of variants are strand-ambiguous "
                "(A/T or C/G) and were kept; their effect allele cannot be confirmed by strand.",
                evidence={"ambiguous_pct": round(variants.ambiguous_pct, 2)},
            )
        )

    if log is None:
        return warnings

    absent = absent_chromosomes(log.chromosomes)
    if absent:
        warnings.append(
            QCWarning(
                code="CHROM_ABSENT",
                message=f"No variants found in the target on chromosome(s) "
                f"{_fmt_chroms(absent)} (score variants in parentheses).",
                evidence={
                    "chromosomes": {c.chrom: c.n_total for c in absent},
                    "absent_run_wide": sorted(c.chrom for c in absent if c.chrom in run_absent),
                },
            )
        )

    if log.indel_total > 0 and log.indel_matched == 0:
        warnings.append(
            QCWarning(
                code="INDELS_UNMATCHED",
                message=f"None of the score's {log.indel_total:,} indel/non-SNV variants matched.",
                evidence={"indel_total": log.indel_total},
            )
        )

    if log.n_missing_position:
        warnings.append(
            QCWarning(
                code="MISSING_POSITION",
                message="Variants with no chromosome or position in the scoring file, which "
                f"cannot be matched by position: {log.n_missing_position:,}.",
                evidence={"n_missing_position": log.n_missing_position},
            )
        )

    return warnings


def score_observations(variants: VariantMetrics, log: LogMetrics | None) -> list[str]:
    if log is None or variants.unmatched_pct + variants.excluded_pct == 0:
        return []
    observations: list[str] = []
    w = log.weights

    if w.weight_coverage_pct is not None:
        gap = w.weight_coverage_pct - variants.matched_pct
        share = "more" if gap < 0 else "less"
        observations.append(
            f"Weight coverage is {w.weight_coverage_pct:.1f}% against a variant match rate of "
            f"{variants.matched_pct:.1f}% ({gap:+.1f} points): the variants that did not match "
            f"carry {share} than their share of total |effect weight|."
        )

    if w.largest_unmatched_abs_weight is not None and w.largest_unmatched_abs_weight == (
        w.max_abs_weight
    ):
        observations.append(
            f"The largest-effect variant in the score (|weight| {w.max_abs_weight:.4g}) "
            "is unmatched."
        )

    n_unmatched_indels = log.indel_total - log.indel_matched
    n_unmatched = (log.snv_total + log.indel_total) - (log.snv_matched + log.indel_matched)
    if n_unmatched and n_unmatched_indels:
        observations.append(
            f"{n_unmatched_indels:,} of {n_unmatched:,} unmatched or excluded variants are "
            "indels/non-SNVs."
        )
    return observations


def run_warnings(log: LogSummary | None, matched_pcts: list[float]) -> list[QCWarning]:
    warnings: list[QCWarning] = []

    median = statistics.median(matched_pcts) if matched_pcts else None
    if median is not None and median < thresholds.BUILD_MISMATCH_MEDIAN_MATCHED_PCT:
        warnings.append(
            QCWarning(
                code="BUILD_MISMATCH_SUSPECTED",
                message=f"Median variant match rate across scores is {median:.1f}%. Uniformly "
                "low matching usually means the scoring files and target data use different "
                "genome builds or variant ID conventions.",
                evidence={
                    "median_matched_pct": round(median, 2),
                    "threshold": thresholds.BUILD_MISMATCH_MEDIAN_MATCHED_PCT,
                },
            )
        )

    if log is None:
        return warnings

    absent = absent_chromosomes(log.chromosomes)
    if absent:
        warnings.append(
            QCWarning(
                code="CHROM_ABSENT",
                message=f"No score variant found in the target on chromosome(s) "
                f"{_fmt_chroms(absent)} "
                "(score-variant pairs in parentheses). The target data likely lacks these "
                "chromosomes or labels them differently.",
                evidence={"chromosomes": {c.chrom: c.n_total for c in absent}},
            )
        )

    if log.indel_total > 0 and log.indel_matched == 0:
        warnings.append(
            QCWarning(
                code="TARGET_LACKS_INDELS",
                message=f"0 of {log.indel_total:,} indel/non-SNV score variants matched while "
                f"{log.snv_matched:,} SNVs did. The target data likely contains no indels "
                "(array data, or indels filtered or not normalised).",
                evidence={
                    "indel_total": log.indel_total,
                    "snv_total": log.snv_total,
                    "snv_matched": log.snv_matched,
                },
            )
        )
    return warnings
