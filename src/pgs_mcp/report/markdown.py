"""Render a RunQC as a Markdown report.

One score gets a single-score report; several get a run report with an overview, a table of
every score, and full sections for the scores returned in detail. An interpretation written by
the client model (guided by the interpretation guide) is placed near the top when supplied.
"""

import re
import statistics
from collections import Counter

from pgs_mcp.qc import thresholds
from pgs_mcp.qc.models import ChromCoverage, QCWarning, RunQC, ScoreQC


def render_report(run: RunQC, interpretation: str | None = None) -> str:
    if run.n_scores == 1 and run.scores:
        blocks = _single_score(run, run.scores[0], interpretation)
    else:
        blocks = _multi_score(run, interpretation)
    return "\n\n".join(b for b in blocks if b) + "\n"


def _single_score(run: RunQC, s: ScoreQC, interpretation: str | None) -> list[str]:
    return [
        f"# Match QC: {s.pgs_id}",
        _facts(
            ("Dataset", run.dataset),
            ("Scoring file", f"`{s.accession}`"),
            ("Pipeline verdict", "pass" if s.pipeline_pass else "**fail**"),
            ("Variant log", "used" if run.log_used else "not supplied (summary metrics only)"),
        ),
        _interpretation(interpretation),
        *_score_sections(s, level=2),
        "## Run context",
        _settings(run),
        _warnings(run.warnings, empty="No run-wide warnings."),
    ]


def _multi_score(run: RunQC, interpretation: str | None) -> list[str]:
    blocks = [
        f"# Match QC: {run.dataset}",
        _facts(
            ("Scores", f"{run.n_scores:,}"),
            ("Pipeline verdict", f"{run.n_pipeline_pass:,} pass, {run.n_pipeline_fail:,} fail"),
            ("Genome build", ", ".join(run.genome_builds) or "unknown"),
            ("Variant log", "used" if run.log_used else "not supplied (summary metrics only)"),
        ),
        _interpretation(interpretation),
        "## Run-wide warnings",
        _warnings(run.warnings, empty="No run-wide warnings."),
        "## Inferred pipeline settings",
        _settings(run),
        "## Distribution across scores",
        _distribution(run),
    ]
    if run.chromosomes:
        blocks += ["## Chromosome coverage (score-variant pairs)", _chrom_table(run.chromosomes)]
    blocks += ["## Warnings by score", _warning_counts(run), "## All scores", _score_table(run)]
    if run.scores:
        blocks.append("## Score details")
        for s in run.scores:
            blocks += [f"### {s.pgs_id}", *_score_sections(s, level=4)]
    elif run.detail_note:
        blocks.append(
            "_Per-score detail is omitted for runs of more than "
            f"{thresholds.DETAIL_SCORE_LIMIT} scores; generate the report for specific scores "
            "to include it._"
        )
    return blocks


def _score_sections(s: ScoreQC, level: int) -> list[str]:
    h = "#" * level
    v = s.variants
    rows = [
        ("Variants in score", f"{s.n_variants:,}"),
        ("Matched", _pct(v.matched_pct)),
        ("Matched, unambiguous", _pct(v.matched_unambiguous_pct)),
        ("Strand-ambiguous", _pct(v.ambiguous_pct)),
        ("Matched after strand flip", _pct(v.flipped_pct)),
        ("Excluded", _pct(v.excluded_pct)),
        ("Unmatched", _pct(v.unmatched_pct)),
    ]
    if s.log:
        w = s.log.weights
        rows += [
            ("Weight coverage (share of total abs. weight matched)", _pct(w.weight_coverage_pct)),
            ("Largest abs. weight in score", _num(w.max_abs_weight)),
            (
                "Largest unmatched (effect/other, abs. weight)",
                "; ".join(f"{v.label}, {_num(v.abs_weight)}" for v in w.largest_unmatched)
                or "none",
            ),
            ("SNVs matched", f"{s.log.snv_matched:,} / {s.log.snv_total:,}"),
            ("Indels/non-SNVs matched", f"{s.log.indel_matched:,} / {s.log.indel_total:,}"),
            ("Variants missing a position", f"{s.log.n_missing_position:,}"),
        ]
    if v.excluded_by_reason:
        reasons = ", ".join(f"{k}: {n:,}" for k, n in v.excluded_by_reason.items())
        rows.append(("Exclusion flags", reasons))

    sections = [f"{h} Metrics", _table(("Metric", "Value"), rows)]
    sections += [f"{h} Warnings", _warnings(s.warnings, empty="No score-level warnings.")]
    if s.observations:
        sections += [f"{h} Observations", "\n".join(f"- {o}" for o in s.observations)]
    if s.log:
        types = [(t, f"{n:,}") for t, n in s.log.match_types.items()]
        types.append(("discarded candidate matches (not_best)", f"{s.log.n_not_best_candidates:,}"))
        sections += [f"{h} Match types", _table(("Match type", "Variants"), types)]
        sections += [f"{h} Chromosome coverage", _chrom_table(s.log.chromosomes)]
    return sections


def _settings(run: RunQC) -> str:
    i = run.inferred_settings
    lo, hi = i.min_overlap_lower, i.min_overlap_upper
    if lo is not None and hi is not None:
        overlap = f"between {lo:.3f} (best failing score) and {hi:.3f} (worst passing score)"
    elif hi is not None:
        overlap = f"at most {hi:.3f} (no score failed)"
    elif lo is not None:
        overlap = f"above {lo:.3f} (no score passed)"
    else:
        overlap = "unknown"

    def flag(value: bool | None) -> str:
        return "unknown" if value is None else "yes" if value else "no"

    return "\n".join(
        [
            "Inferred from the match summary, not read from the run configuration:",
            "",
            f"- `min_overlap`: {overlap}",
            f"- `keep_ambiguous`: {flag(i.keep_ambiguous)}",
            f"- `keep_multiallelic`: {flag(i.keep_multiallelic)}",
        ]
    )


def _distribution(run: RunQC) -> str:
    metrics = {"Variants matched (%)": [r.matched_pct for r in run.score_index.rows]}
    weights = [
        r.weight_coverage_pct for r in run.score_index.rows if r.weight_coverage_pct is not None
    ]
    if weights:
        metrics["Weight coverage (%)"] = weights
        gaps = [
            r.weight_coverage_pct - r.matched_pct
            for r in run.score_index.rows
            if r.weight_coverage_pct is not None
        ]
        metrics["Weight minus variant coverage (points)"] = gaps
    metrics["Variants per score"] = [r.n_variants for r in run.score_index.rows]
    rows = [
        (name, _num(min(vals)), _num(statistics.median(vals)), _num(max(vals)))
        for name, vals in metrics.items()
    ]
    return _table(("Metric", "Min", "Median", "Max"), rows)


def _warning_counts(run: RunQC) -> str:
    counts = Counter(code for r in run.score_index.rows for code in r.warning_codes)
    if not counts:
        return "No score-level warnings."
    return _table(("Warning", "Scores"), [(f"`{c}`", f"{n:,}") for c, n in counts.most_common()])


def _score_table(run: RunQC) -> str:
    if run.log_used:
        # Weight coverage surfaces scores that passed on variant count but lost their signal.
        order = "weight coverage"
        rows = sorted(
            run.score_index.rows,
            key=lambda r: (r.weight_coverage_pct is None, r.weight_coverage_pct or 0),
        )
    else:
        order = "variant match rate"
        rows = sorted(run.score_index.rows, key=lambda r: r.matched_pct)
    return f"Sorted by {order}, lowest first.\n\n" + _table(
        ("Score", "Variants", "Pipeline", "Matched", "Weight coverage", "Warnings"),
        [
            (
                r.pgs_id,
                f"{r.n_variants:,}",
                "pass" if r.pipeline_pass else "**fail**",
                _pct(r.matched_pct),
                _pct(r.weight_coverage_pct),
                ", ".join(r.warning_codes) or "",
            )
            for r in rows
        ],
    )


def _chrom_table(chroms: list[ChromCoverage]) -> str:
    return _table(
        ("Chromosome", "Variants", "Matched", "Matched %"),
        [
            (c.chrom, f"{c.n_total:,}", f"{c.n_matched:,}", _pct(100 * c.n_matched / c.n_total))
            for c in chroms
        ],
    )


def _warnings(warnings: list[QCWarning], empty: str) -> str:
    if not warnings:
        return empty
    return "\n".join(f"- **`{w.code}`**: {w.message}" for w in warnings)


_OWN_HEADING = re.compile(r"^#+\s*interpretation\b[^\n]*\n+", re.IGNORECASE)


def _interpretation(text: str | None) -> str:
    """Place the model's interpretation under one heading, dropping a duplicate of its own."""
    body = _OWN_HEADING.sub("", (text or "").strip(), count=1).strip()
    return f"## Interpretation\n\n{body}" if body else ""


def _facts(*pairs: tuple[str, str]) -> str:
    return "\n".join(f"- **{k}:** {v}" for k, v in pairs)


def _table(header: tuple[str, ...], rows: list[tuple[str, ...]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(str(c).replace("|", "\\|") for c in row) + " |" for row in rows]
    return "\n".join(lines)


def _pct(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.1f}%"


def _num(value: float | None) -> str:
    if value is None:
        return "n/a"
    if float(value).is_integer() and abs(value) >= 1:
        return f"{int(value):,}"
    return f"{value:.4g}"
