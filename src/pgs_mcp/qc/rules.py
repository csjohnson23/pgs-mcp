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


def _n(count: int, noun: str) -> str:
    """`_n(1, "variant")` -> "1 variant"; `_n(2, "variant")` -> "2 variants"."""
    return f"{count:,} {noun}" if count == 1 else f"{count:,} {noun}s"


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
                message=f"Score has only {_n(n_variants, 'variant')}; each one carries substantial "
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
                message=f"None of the score's {_n(log.indel_total, 'indel/non-SNV variant')} "
                "matched.",
                evidence={"indel_total": log.indel_total},
            )
        )

    if log.n_missing_position:
        warnings.append(
            QCWarning(
                code="MISSING_POSITION",
                message="Variants missing a chromosome or position in the scoring file, which "
                f"cannot be matched by position: {log.n_missing_position:,}.",
                evidence={"n_missing_position": log.n_missing_position},
            )
        )

    return warnings


def _autosome_range(chromosomes: list[ChromCoverage]) -> str | None:
    """Spread of per-autosome match rates: narrow means even loss, wide means regional loss."""
    rates = [
        (100 * c.n_matched / c.n_total, c.chrom)
        for c in chromosomes
        if c.n_total and c.chrom.isdigit() and 1 <= int(c.chrom) <= 22
    ]
    if len(rates) < 2:
        return None
    (low, low_chrom), (high, high_chrom) = min(rates), max(rates)
    if f"{low:.1f}" == f"{high:.1f}":
        return f"Variant match rate is {low:.1f}% on each of the {len(rates)} autosomes."
    return (
        f"Variant match rate per autosome ranges from {low:.1f}% (chr{low_chrom}) to "
        f"{high:.1f}% (chr{high_chrom}) across {len(rates)} autosomes."
    )


def score_observations(
    n_variants: int, variants: VariantMetrics, log: LogMetrics | None
) -> list[str]:
    if log is None or variants.unmatched_pct + variants.excluded_pct == 0:
        return []
    if variants.matched_pct == 0:
        return ["No variant matched, so this score cannot be computed in the target data."]
    observations: list[str] = []
    w = log.weights

    if w.weight_coverage_pct is not None:
        gap = w.weight_coverage_pct - variants.matched_pct
        if abs(gap) < 0.05:
            share = "about their share"
        else:
            share = ("more" if gap < 0 else "less") + " than their share"
        observations.append(
            f"Weight coverage is {w.weight_coverage_pct:.1f}% against a variant match rate of "
            f"{variants.matched_pct:.1f}% ({gap:+.1f} points): the variants that did not match "
            f"carry {share} of total |effect weight|."
        )

    top = w.largest_unmatched
    if top and top[0].abs_weight == w.max_abs_weight:
        share = 100 * top[0].abs_weight / w.total_abs_weight if w.total_abs_weight else 0
        observations.append(
            f"The largest-effect variant in the score ({top[0].label}, |weight| "
            f"{top[0].abs_weight:.4g}, {share:.2g}% of total |weight|) is unmatched."
        )

    autosome_range = _autosome_range(log.chromosomes)
    if autosome_range and n_variants >= thresholds.SMALL_SCORE_MAX_VARIANTS:
        observations.append(autosome_range)

    n_lost_indels = log.indel_total - log.indel_matched
    n_lost = (log.snv_total + log.indel_total) - (log.snv_matched + log.indel_matched)
    if n_lost_indels == n_lost == 1:
        observations.append("The one unmatched or excluded variant is an indel/non-SNV.")
    elif n_lost_indels == n_lost:
        observations.append(f"All {n_lost:,} unmatched or excluded variants are indels/non-SNVs.")
    elif n_lost_indels:
        observations.append(
            f"{n_lost_indels:,} of {n_lost:,} unmatched or excluded variants are indels/non-SNVs."
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
        warnings.append(
            QCWarning(
                code="LOG_NOT_SUPPLIED",
                message="No variant-level match log was supplied, so weight coverage, the "
                "largest unmatched effects, per-chromosome coverage, indel matching and missing "
                "positions are unknown. A score can pass on variant match rate while losing most "
                "of its effect weight; supply the match log to check.",
            )
        )
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
