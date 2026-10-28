"""Assemble a RunQC from a match summary and, optionally, the variant-level match log."""

from pathlib import Path

from pgs_mcp.qc import thresholds
from pgs_mcp.qc.errors import QCInputError
from pgs_mcp.qc.log import LogSummary, ProgressFn, aggregate_log, summarize_log
from pgs_mcp.qc.models import (
    InferredSettings,
    QCWarning,
    RunQC,
    ScoreIndex,
    ScoreIndexRow,
    ScoreQC,
)
from pgs_mcp.qc.rules import absent_chromosomes, run_warnings, score_observations, score_warnings
from pgs_mcp.qc.summary import SummaryScore, parse_accession, read_summary


def assess_match_qc(
    summary_path: Path,
    log_path: Path | None = None,
    pgs_ids: list[str] | None = None,
    progress: ProgressFn | None = None,
) -> RunQC:
    """Run-wide context (inferred settings, run warnings) always uses every score in the file;
    `pgs_ids` limits which scores are reported, in the order given."""
    summary = read_summary(_existing(summary_path))
    selected = _select(summary, pgs_ids)
    inferred = _infer_settings(summary)

    log: LogSummary | None = None
    integrity: list[QCWarning] = []
    if log_path is not None:
        log = summarize_log(aggregate_log(_existing(log_path), progress))
        integrity = _check_log_matches_summary(summary, log)

    run_absent = {c.chrom for c in absent_chromosomes(log.chromosomes)} if log else set()
    scores = [_score_qc(s, log, inferred.keep_ambiguous, run_absent) for s in selected]

    if pgs_ids or len(scores) <= thresholds.DETAIL_SCORE_LIMIT:
        detailed, note = scores, None
    else:
        detailed = []
        note = (
            f"{len(scores)} scores: full detail omitted to keep the result small. Call again "
            "with pgs_ids for the scores you want in detail."
        )

    return RunQC(
        dataset=summary[0].dataset,
        genome_builds=sorted({s.genome_build for s in scores if s.genome_build}),
        n_scores=len(scores),
        n_pipeline_pass=sum(s.pipeline_pass for s in scores),
        n_pipeline_fail=sum(not s.pipeline_pass for s in scores),
        inferred_settings=inferred,
        log_used=log is not None,
        chromosomes=log.chromosomes if log and not pgs_ids else None,
        warnings=integrity + run_warnings(log, [s.metrics().matched_pct for s in summary]),
        score_index=ScoreIndex(rows=[_index_row(s) for s in scores]),
        scores=detailed,
        detail_note=note,
    )


def _existing(path: Path) -> Path:
    path = path.expanduser()
    if not path.is_file():
        raise QCInputError(f"File not found: {path}")
    return path


def _select(summary: list[SummaryScore], pgs_ids: list[str] | None) -> list[SummaryScore]:
    """Match PGS IDs or accessions case-insensitively, keeping the requested order."""
    if not pgs_ids:
        return summary
    by_key: dict[str, SummaryScore] = {}
    for s in summary:
        by_key[s.accession.lower()] = s
        by_key.setdefault(parse_accession(s.accession)[0].lower(), s)
    selected: list[SummaryScore] = []
    unknown: list[str] = []
    for requested in pgs_ids:
        match = by_key.get(requested.strip().lower())
        if match is None:
            unknown.append(requested)
        elif match not in selected:
            selected.append(match)
    if unknown:
        raise QCInputError(
            f"Not in this match summary: {', '.join(unknown)}. "
            f"It contains {len(summary)} scores, e.g. {summary[0].accession}."
        )
    return selected


def _infer_settings(summary: list[SummaryScore]) -> InferredSettings:
    rates = [(s.n_matched / s.n_variants, s.pipeline_pass) for s in summary if s.n_variants]
    failing = [r for r, passed in rates if not passed]
    passing = [r for r, passed in rates if passed]

    def flag(kept: bool, dropped: bool) -> bool | None:
        return True if kept else False if dropped else None

    return InferredSettings(
        min_overlap_lower=max(failing) if failing else None,
        min_overlap_upper=min(passing) if passing else None,
        keep_ambiguous=flag(
            any(s.n_matched_ambiguous for s in summary), any(s.excluded_ambiguous for s in summary)
        ),
        keep_multiallelic=flag(
            any(s.matched_multiallelic for s in summary),
            any(s.excluded_multiallelic for s in summary),
        ),
    )


def _check_log_matches_summary(summary: list[SummaryScore], log: LogSummary) -> list[QCWarning]:
    """Catch a log and summary from different runs before reporting numbers from both."""
    missing = [s.accession for s in summary if s.accession not in log.n_variants]
    differing = {
        s.accession: {"summary": s.n_variants, "log": log.n_variants[s.accession]}
        for s in summary
        if s.accession in log.n_variants and log.n_variants[s.accession] != s.n_variants
    }
    if not missing and not differing:
        return []
    return [
        QCWarning(
            code="LOG_SUMMARY_MISMATCH",
            message=f"The match log does not agree with the summary: {len(missing)} scores "
            f"missing from the log, {len(differing)} with different variant counts. They may be "
            "from different runs; log-derived metrics for these scores are unreliable.",
            evidence={
                "missing_from_log": missing[:20],
                "count_differences": dict(list(differing.items())[:20]),
            },
        )
    ]


def _score_qc(
    s: SummaryScore, log: LogSummary | None, keep_ambiguous: bool | None, run_absent: set[str]
) -> ScoreQC:
    pgs_id, build = parse_accession(s.accession)
    variants = s.metrics()
    log_metrics = log.scores.get(s.accession) if log else None
    return ScoreQC(
        pgs_id=pgs_id,
        accession=s.accession,
        genome_build=build,
        n_variants=s.n_variants,
        pipeline_pass=s.pipeline_pass,
        variants=variants,
        log=log_metrics,
        warnings=score_warnings(s.n_variants, variants, log_metrics, keep_ambiguous, run_absent),
        observations=score_observations(variants, log_metrics),
    )


def _index_row(s: ScoreQC) -> ScoreIndexRow:
    return ScoreIndexRow(
        pgs_id=s.pgs_id,
        n_variants=s.n_variants,
        pipeline_pass=s.pipeline_pass,
        matched_pct=s.variants.matched_pct,
        weight_coverage_pct=s.log.weights.weight_coverage_pct if s.log else None,
        warning_codes=[w.code for w in s.warnings],
    )
