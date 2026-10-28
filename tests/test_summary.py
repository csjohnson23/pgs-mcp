import pytest

from pgs_mcp.qc import QCInputError
from pgs_mcp.qc.summary import parse_accession, read_summary


def test_variant_metrics(summary):
    scores = {s.accession: s for s in read_summary(summary)}
    one = scores["PGS000001_hmPOS_GRCh38"]
    assert one.n_variants == 10
    assert one.pipeline_pass
    m = one.metrics()
    assert m.matched_pct == pytest.approx(90)
    assert m.matched_unambiguous_pct == pytest.approx(70)
    assert m.ambiguous_pct == pytest.approx(20)
    assert m.flipped_pct == pytest.approx(10)
    assert m.unmatched_pct == pytest.approx(10)

    two = scores["PGS000002_hmPOS_GRCh38"].metrics()
    assert two.excluded_pct == pytest.approx(100 / 6)
    assert two.excluded_by_reason == {"duplicate_ID": 1}


@pytest.mark.parametrize(
    ("accession", "expected"),
    [
        ("PGS000004_hmPOS_GRCh37", ("PGS000004", "GRCh37")),
        ("PGS000004", ("PGS000004", None)),
        ("my_custom_score", ("my_custom_score", None)),
    ],
)
def test_parse_accession(accession, expected):
    assert parse_accession(accession) == expected


def test_rejects_non_summary(log):
    with pytest.raises(QCInputError, match="missing columns"):
        read_summary(log)


def test_rejects_mixed_datasets(tmp_path, summary):
    lines = summary.read_text().splitlines()
    lines.append(lines[1].replace("TOY,", "OTHER,", 1).replace("PGS000001", "PGS000009"))
    mixed = tmp_path / "mixed.csv"
    mixed.write_text("\n".join(lines) + "\n")
    with pytest.raises(QCInputError, match="mixes datasets"):
        read_summary(mixed)
