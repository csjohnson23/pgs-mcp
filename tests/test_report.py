from pgs_mcp.qc import assess_match_qc, thresholds
from pgs_mcp.report import render_report


def test_single_score_report(single_summary, log):
    report = render_report(assess_match_qc(single_summary, log))
    assert report.startswith("# Match QC: PGS000001\n")
    assert "## Metrics" in report
    assert "## Run context" in report
    assert "## All scores" not in report


def test_single_score_selected_from_run(summary, log):
    report = render_report(assess_match_qc(summary, log, ["PGS000002"]))
    assert report.startswith("# Match QC: PGS000002\n")
    assert "**fail**" in report


def test_multi_score_report(summary, log):
    report = render_report(assess_match_qc(summary, log))
    assert report.startswith("# Match QC: TOY\n")
    for section in ("Run-wide warnings", "Distribution across scores", "All scores"):
        assert f"## {section}" in report
    table = report.split("## All scores")[1]
    assert table.index("PGS000002") < table.index("PGS000001")  # lowest weight coverage first
    assert "### PGS000003" in report


def test_model_heading_not_duplicated(summary, log):
    report = render_report(assess_match_qc(summary, log), "## Interpretation\n\nAll fine.")
    assert report.count("Interpretation") == 1


def test_table_sorted_by_weight_coverage_with_log(summary, log):
    table = render_report(assess_match_qc(summary, log)).split("## All scores")[1]
    assert "Sorted by weight coverage" in table
    # PGS000001 matches the most variants but keeps the least weight among passing scores.
    assert table.index("PGS000002") < table.index("PGS000001") < table.index("PGS000003")


def test_interpretation_placed_before_metrics(summary, log):
    report = render_report(assess_match_qc(summary, log), "Most scores lost indel weight.")
    assert report.index("## Interpretation") < report.index("## Run-wide warnings")
    assert "Most scores lost indel weight." in report


def test_summary_only_report_omits_log_sections(summary):
    report = render_report(assess_match_qc(summary))
    assert "not supplied" in report
    assert "Chromosome coverage" not in report


def test_large_run_report_points_to_detail(summary, monkeypatch):
    monkeypatch.setattr(thresholds, "DETAIL_SCORE_LIMIT", 2)
    report = render_report(assess_match_qc(summary))
    assert "## Score details" not in report
    assert "specific scores" in report
    assert "pgs_ids" not in report  # tool vocabulary stays out of the human report
