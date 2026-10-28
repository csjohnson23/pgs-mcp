"""Aggregate the pgsc_calc variant-level match log (`{sampleset}_log.csv.gz`).

The log has one row per candidate match and can run to hundreds of millions of rows, so it is
streamed in chunks and reduced to a small aggregate table (one row per score x chromosome x
variant class x match outcome). Everything downstream is derived from that table, which is
cached on disk keyed by the log's path, size and modification time.
"""

import gzip
import hashlib
import os
from collections import Counter, defaultdict
from collections.abc import Callable, Iterator
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

from pgs_mcp.qc.errors import QCInputError
from pgs_mcp.qc.models import ChromCoverage, LogMetrics, WeightMetrics

ProgressFn = Callable[[int, int], None]

COLUMNS = (
    "accession",
    "chr_name",
    "chr_position",
    "effect_allele",
    "other_allele",
    "effect_weight",
    "match_type",
    "match_status",
)
GROUP_KEYS = ["accession", "chr_name", "is_snv", "missing_position", "match_status", "match_type"]
CHUNK_BYTES = 256 * 1024 * 1024
CACHE_VERSION = 1  # bump when the aggregate table's shape or meaning changes


def aggregate_log(path: Path, progress: ProgressFn | None = None) -> pl.DataFrame:
    """Return the aggregate table for a match log, from cache when the file is unchanged."""
    cache = _cache_path(path)
    if cache.exists():
        return pl.read_parquet(cache)

    partials = [_aggregate_chunk(chunk) for chunk in _iter_csv_chunks(path, progress)]
    if not partials:
        raise QCInputError(f"Match log {path} contains no rows.")
    agg = (
        pl.concat(partials)
        .group_by(GROUP_KEYS)
        .agg(pl.col("n").sum(), pl.col("abs_w_sum").sum(), pl.col("abs_w_max").max())
    )

    cache.parent.mkdir(parents=True, exist_ok=True)
    tmp = cache.with_suffix(".tmp")
    agg.write_parquet(tmp)
    tmp.replace(cache)
    return agg


def _cache_path(path: Path) -> Path:
    stat = path.stat()
    key = f"{path.resolve()}|{stat.st_size}|{stat.st_mtime_ns}|v{CACHE_VERSION}"
    root = Path(
        os.environ.get("PGS_MCP_CACHE_DIR")
        or Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "pgs-mcp"
    )
    return root / f"log-{hashlib.sha256(key.encode()).hexdigest()[:20]}.parquet"


def _iter_csv_chunks(path: Path, progress: ProgressFn | None) -> Iterator[bytes]:
    """Yield the header plus a block of whole lines, for each ~CHUNK_BYTES of decompressed data."""
    total = path.stat().st_size
    try:
        with path.open("rb") as raw:
            is_gzip = raw.read(2) == b"\x1f\x8b"
            raw.seek(0)
            stream = gzip.GzipFile(fileobj=raw) if is_gzip else raw
            header = stream.readline()
            _check_header(path, header)
            carry = b""
            while block := stream.read(CHUNK_BYTES):
                block = carry + block
                cut = block.rfind(b"\n") + 1
                carry = block[cut:]
                if cut:
                    yield header + block[:cut]
                if progress:
                    progress(raw.tell(), total)
            if carry.strip():
                yield header + carry
    except (OSError, EOFError) as e:
        raise QCInputError(f"Cannot read match log {path}: {e}") from e


def _check_header(path: Path, header: bytes) -> None:
    names = set(header.decode(errors="replace").strip().split(","))
    missing = set(COLUMNS) - names
    if missing:
        raise QCInputError(
            f"{path} is not a pgsc_calc match log; missing columns: {', '.join(sorted(missing))}"
        )


def _aggregate_chunk(chunk: bytes) -> pl.DataFrame:
    abs_w = pl.col("effect_weight").cast(pl.Float64, strict=False).abs()
    return (
        pl.read_csv(chunk, columns=list(COLUMNS), infer_schema=False)
        .lazy()
        .with_columns(
            is_snv=(pl.col("effect_allele").str.len_chars().fill_null(0) <= 1)
            & (pl.col("other_allele").str.len_chars().fill_null(0) <= 1),
            missing_position=pl.col("chr_name").is_null() | pl.col("chr_position").is_null(),
        )
        .group_by(GROUP_KEYS)
        .agg(n=pl.len(), abs_w_sum=abs_w.sum(), abs_w_max=abs_w.max())
        .collect()
    )


@dataclass
class _ScoreAccumulator:
    n_variants: int = 0
    w_total: float = 0.0
    w_matched: float = 0.0
    w_max: float | None = None
    w_max_unmatched: float | None = None
    match_types: Counter[str] = field(default_factory=Counter)
    n_not_best: int = 0
    snv_total: int = 0
    snv_matched: int = 0
    indel_total: int = 0
    indel_matched: int = 0
    n_missing_position: int = 0
    chrom_total: Counter[str] = field(default_factory=Counter)
    chrom_matched: Counter[str] = field(default_factory=Counter)
    chrom_excluded: Counter[str] = field(default_factory=Counter)

    def metrics(self) -> LogMetrics:
        return LogMetrics(
            weights=WeightMetrics(
                weight_coverage_pct=100 * self.w_matched / self.w_total if self.w_total else None,
                max_abs_weight=self.w_max,
                largest_unmatched_abs_weight=self.w_max_unmatched,
            ),
            match_types=dict(self.match_types.most_common()),
            n_not_best_candidates=self.n_not_best,
            snv_total=self.snv_total,
            snv_matched=self.snv_matched,
            indel_total=self.indel_total,
            indel_matched=self.indel_matched,
            n_missing_position=self.n_missing_position,
            chromosomes=_chrom_coverage(self.chrom_total, self.chrom_matched, self.chrom_excluded),
        )


def _max(a: float | None, b: float | None) -> float | None:
    return b if a is None else a if b is None else max(a, b)


@dataclass
class LogSummary:
    scores: dict[str, LogMetrics]
    n_variants: dict[str, int]
    chromosomes: list[ChromCoverage]
    snv_total: int
    snv_matched: int
    indel_total: int
    indel_matched: int


def summarize_log(agg: pl.DataFrame) -> LogSummary:
    """Derive per-score and run-wide metrics from the aggregate table."""
    acc: dict[str, _ScoreAccumulator] = defaultdict(_ScoreAccumulator)
    for row in agg.iter_rows(named=True):
        s = acc[row["accession"]]
        n, status = row["n"], row["match_status"]
        if status == "not_best":
            s.n_not_best += n
            continue
        matched = status == "matched"
        s.n_variants += n
        s.w_total += row["abs_w_sum"] or 0.0
        s.w_max = _max(s.w_max, row["abs_w_max"])
        if matched:
            s.w_matched += row["abs_w_sum"] or 0.0
            s.match_types[row["match_type"] or "unknown"] += n
        elif status == "unmatched":
            s.w_max_unmatched = _max(s.w_max_unmatched, row["abs_w_max"])
        if row["is_snv"]:
            s.snv_total += n
            s.snv_matched += n * matched
        else:
            s.indel_total += n
            s.indel_matched += n * matched
        if row["missing_position"]:
            s.n_missing_position += n
        else:
            s.chrom_total[row["chr_name"]] += n
            s.chrom_matched[row["chr_name"]] += n * matched
            s.chrom_excluded[row["chr_name"]] += n * (status == "excluded")

    run_total: Counter[str] = Counter()
    run_matched: Counter[str] = Counter()
    run_excluded: Counter[str] = Counter()
    for s in acc.values():
        run_total.update(s.chrom_total)
        run_matched.update(s.chrom_matched)
        run_excluded.update(s.chrom_excluded)

    return LogSummary(
        scores={a: s.metrics() for a, s in acc.items()},
        n_variants={a: s.n_variants for a, s in acc.items()},
        chromosomes=_chrom_coverage(run_total, run_matched, run_excluded),
        snv_total=sum(s.snv_total for s in acc.values()),
        snv_matched=sum(s.snv_matched for s in acc.values()),
        indel_total=sum(s.indel_total for s in acc.values()),
        indel_matched=sum(s.indel_matched for s in acc.values()),
    )


_CHROM_ORDER = {str(i): i for i in range(1, 23)} | {"X": 23, "Y": 24, "XY": 25, "MT": 26}


def _chrom_coverage(
    total: Counter[str], matched: Counter[str], excluded: Counter[str]
) -> list[ChromCoverage]:
    order = sorted(total, key=lambda c: (_CHROM_ORDER.get(c, 99), c))
    return [
        ChromCoverage(chrom=c, n_total=total[c], n_matched=matched[c], n_excluded=excluded[c])
        for c in order
    ]
