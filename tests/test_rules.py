"""Edge cases for observations and wording, built directly from metrics."""

from pgs_mcp.qc.models import LogMetrics, UnmatchedVariant, VariantMetrics, WeightMetrics
from pgs_mcp.qc.rules import score_observations, score_warnings


def _variants(matched_pct: float) -> VariantMetrics:
    return VariantMetrics(
        matched_pct=matched_pct,
        matched_unambiguous_pct=matched_pct,
        ambiguous_pct=0,
        flipped_pct=0,
        excluded_pct=0,
        unmatched_pct=100 - matched_pct,
    )


def _log(weight_coverage_pct: float, indel_total: int = 0) -> LogMetrics:
    return LogMetrics(
        weights=WeightMetrics(
            weight_coverage_pct=weight_coverage_pct,
            max_abs_weight=1.0,
            largest_unmatched=[
                UnmatchedVariant(
                    chrom=None, position=None, effect_allele="A", other_allele="G", abs_weight=0.5
                )
            ],
        ),
        match_types={},
        n_not_best_candidates=0,
        snv_total=1 - indel_total,
        snv_matched=0,
        indel_total=indel_total,
        indel_matched=0,
        n_missing_position=1,
        chromosomes=[],
    )


def test_nothing_matched_is_stated_plainly():
    assert score_observations(_variants(0), _log(0)) == [
        "No variant matched, so this score cannot be computed in the target data."
    ]


def test_equal_coverage_is_not_called_more_or_less():
    [obs] = score_observations(_variants(50), _log(50.01))[:1]
    assert "about their share" in obs


def test_singular_wording():
    warnings = score_warnings(1, _variants(0), _log(0, indel_total=1), None, set())
    messages = " ".join(w.message for w in warnings)
    assert "only 1 variant;" in messages
    assert "1 indel/non-SNV variant matched" in messages


def test_label_for_effect_only_entry():
    v = UnmatchedVariant(
        chrom=None, position=None, effect_allele="HLA-DQ2.5", other_allele=None, abs_weight=1
    )
    assert v.label == "no position, effect HLA-DQ2.5"
