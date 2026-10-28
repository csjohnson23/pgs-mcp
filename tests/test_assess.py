import pytest

from pgs_mcp.qc import QCInputError, assess_match_qc, thresholds


def _codes(warnings):
    return {w.code for w in warnings}


def test_summary_only(summary):
    run = assess_match_qc(summary)
    assert not run.log_used
    assert run.chromosomes is None
    assert _codes(run.warnings) == {"LOG_NOT_SUPPLIED"}
    assert (run.n_scores, run.n_pipeline_pass, run.n_pipeline_fail) == (3, 2, 1)
    assert all(s.log is None and not s.observations for s in run.scores)


def test_inferred_settings(summary):
    i = assess_match_qc(summary).inferred_settings
    assert (i.min_overlap_lower, i.min_overlap_upper) == (0.5, 0.75)
    assert i.keep_ambiguous is True
    assert i.keep_multiallelic is None


def test_run_warnings(summary, log):
    run = assess_match_qc(summary, log)
    assert _codes(run.warnings) == {"CHROM_ABSENT", "TARGET_LACKS_INDELS"}
    absent = next(w for w in run.warnings if w.code == "CHROM_ABSENT")
    assert absent.evidence["chromosomes"] == {"7": 1, "X": 2}  # chr9 was found, then excluded


def test_score_warnings(summary, log):
    scores = {s.pgs_id: s for s in assess_match_qc(summary, log).scores}
    assert _codes(scores["PGS000001"].warnings) == {
        "SMALL_SCORE",
        "AMBIGUOUS_RETAINED",
        "CHROM_ABSENT",
        "INDELS_UNMATCHED",
    }
    assert _codes(scores["PGS000002"].warnings) == {
        "SMALL_SCORE",
        "DUPLICATES_EXCLUDED",
        "CHROM_ABSENT",
    }
    assert _codes(scores["PGS000003"].warnings) == {"SMALL_SCORE", "MISSING_POSITION"}
    small = next(w for w in scores["PGS000003"].warnings if w.code == "SMALL_SCORE")
    assert small.evidence["threshold"] == thresholds.SMALL_SCORE_MAX_VARIANTS


def test_observations_describe_weight_loss(summary, log):
    one = assess_match_qc(summary, log, ["PGS000001"]).scores[0]
    assert any("-44.9 points" in o for o in one.observations)
    assert any("largest-effect variant" in o for o in one.observations)


def test_select_by_pgs_id_or_accession(summary):
    assert [s.pgs_id for s in assess_match_qc(summary, pgs_ids=["PGS000002"]).scores] == [
        "PGS000002"
    ]
    run = assess_match_qc(summary, pgs_ids=["PGS000003_hmPOS_GRCh38"])
    assert run.n_scores == 1
    # Run-wide context still comes from every score in the file.
    assert run.inferred_settings.min_overlap_lower == 0.5


def test_selection_is_case_insensitive_and_ordered(summary, log):
    run = assess_match_qc(summary, log, ["pgs000003", "PGS000001", "pgs000003"])
    assert [s.pgs_id for s in run.scores] == ["PGS000003", "PGS000001"]
    assert run.chromosomes is None  # run-wide table only for whole-run calls
    assert {w.code for w in run.warnings} == {"CHROM_ABSENT", "TARGET_LACKS_INDELS"}


def test_unknown_pgs_id(summary):
    with pytest.raises(QCInputError, match="PGS999999"):
        assess_match_qc(summary, pgs_ids=["PGS999999"])


def test_missing_file(tmp_path):
    with pytest.raises(QCInputError, match="File not found"):
        assess_match_qc(tmp_path / "nope.csv")


def test_detail_limited_for_large_runs(summary, monkeypatch):
    monkeypatch.setattr(thresholds, "DETAIL_SCORE_LIMIT", 2)
    run = assess_match_qc(summary)
    assert run.scores == []
    assert len(run.score_index.rows) == 3
    assert "pgs_ids" in run.detail_note


def test_log_from_another_run_is_flagged(summary, log, tmp_path):
    other = tmp_path / "other_log.csv"
    other.write_text(
        "".join(line for line in log.open() if "PGS000003" not in line and "X," not in line)
    )
    run = assess_match_qc(summary, other)
    mismatch = next(w for w in run.warnings if w.code == "LOG_SUMMARY_MISMATCH")
    assert mismatch.evidence["missing_from_log"] == ["PGS000003_hmPOS_GRCh38"]
    assert "PGS000002_hmPOS_GRCh38" in mismatch.evidence["count_differences"]
