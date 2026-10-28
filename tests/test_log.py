import gzip

import pytest
from pydantic import TypeAdapter

from pgs_mcp.qc import QCInputError
from pgs_mcp.qc import log as log_module
from pgs_mcp.qc.log import LogSummary, aggregate_log, summarize_log


def _summary(path):
    return summarize_log(aggregate_log(path))


def test_weights_and_variant_classes(log):
    one = _summary(log).scores["PGS000001_hmPOS_GRCh38"]
    assert one.weights.weight_coverage_pct == pytest.approx(100 * 0.41 / 0.91)
    assert one.weights.max_abs_weight == pytest.approx(0.5)
    assert one.weights.largest_unmatched_abs_weight == pytest.approx(0.5)
    assert (one.snv_total, one.snv_matched, one.indel_total, one.indel_matched) == (9, 9, 1, 0)
    assert one.n_not_best_candidates == 1
    assert one.match_types["no_oa_alt"] == 1


def test_not_best_rows_do_not_count_as_variants(log):
    assert _summary(log).n_variants["PGS000001_hmPOS_GRCh38"] == 10


def test_chromosome_coverage(log):
    s = _summary(log)
    by_chrom = {c.chrom: c for c in s.chromosomes}
    assert (by_chrom["X"].n_total, by_chrom["X"].n_matched) == (2, 0)
    assert by_chrom["9"].n_excluded == 1
    assert [c.chrom for c in s.chromosomes][-1] == "X"  # natural order, not lexical
    assert s.scores["PGS000003_hmPOS_GRCh38"].n_missing_position == 1


def _rounded(value):
    """Chunked float sums can differ in the last bit; compare at a sane precision."""
    if isinstance(value, float):
        return round(value, 9)
    if isinstance(value, dict):
        return {k: _rounded(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_rounded(v) for v in value]
    return value


def test_gzip_and_small_chunks_match_plain(log, tmp_path, monkeypatch):
    expected = _summary(log)
    gz = tmp_path / "toy_log.csv.gz"
    gz.write_bytes(gzip.compress(log.read_bytes()))
    monkeypatch.setattr(log_module, "CHUNK_BYTES", 64)  # forces lines to split across chunks
    actual = _summary(gz)
    dump = TypeAdapter(LogSummary).dump_python
    assert _rounded(dump(actual)) == _rounded(dump(expected))


def test_cache_reused_until_file_changes(log, tmp_path, isolated_cache):
    copy = tmp_path / "log.csv"
    copy.write_bytes(log.read_bytes())
    aggregate_log(copy)
    assert len(list(isolated_cache.glob("*.parquet"))) == 1
    aggregate_log(copy)
    assert len(list(isolated_cache.glob("*.parquet"))) == 1
    copy.write_bytes(log.read_bytes() + log.read_bytes().splitlines(keepends=True)[1])
    aggregate_log(copy)
    assert len(list(isolated_cache.glob("*.parquet"))) == 2


def test_rejects_non_log(summary):
    with pytest.raises(QCInputError, match="not a pgsc_calc match log"):
        aggregate_log(summary)
