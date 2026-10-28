"""Edge cases for observations and wording, built directly from metrics."""

from pgs_mcp.qc.models import (
    ChromCoverage,
    LogMetrics,
    UnmatchedVariant,
    VariantMetrics,
    WeightMetrics,
)
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
            total_abs_weight=2.0,
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
    assert score_observations(1, _variants(0), _log(0)) == [
        "No variant matched, so this score cannot be computed in the target data."
    ]


def test_equal_coverage_is_not_called_more_or_less():
    [obs] = score_observations(1, _variants(50), _log(50.01))[:1]
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


def test_autosome_range_only_for_larger_scores():
    log = _log(40)
    log.chromosomes = [
        ChromCoverage(chrom="1", n_total=100, n_matched=33, n_excluded=0),
        ChromCoverage(chrom="2", n_total=100, n_matched=37, n_excluded=0),
        ChromCoverage(chrom="X", n_total=10, n_matched=0, n_excluded=0),  # not an autosome
    ]
    [obs] = [o for o in score_observations(210, _variants(35), log) if "per autosome" in o]
    assert obs == (
        "Variant match rate per autosome ranges from 33.0% (chr1) to 37.0% (chr2) "
        "across 2 autosomes."
    )
    assert not [o for o in score_observations(5, _variants(35), log) if "per autosome" in o]


def test_largest_effect_share_of_total_weight():
    log = _log(10)
    log.weights.max_abs_weight = 0.5
    [obs] = [o for o in score_observations(1, _variants(50), log) if "largest-effect" in o]
    assert "25% of total |weight|" in obs


def test_autosome_range_when_all_equal():
    log = _log(100)
    log.chromosomes = [
        ChromCoverage(chrom=c, n_total=10, n_matched=10, n_excluded=0) for c in ("1", "6")
    ]
    obs = [o for o in score_observations(53, _variants(90), log) if "autosome" in o]
    assert obs == ["Variant match rate is 100.0% on each of the 2 autosomes."]
