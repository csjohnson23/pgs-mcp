"""Aggregate the pgsc_calc variant-level match log (`{sampleset}_log.csv.gz`).

The log has one row per candidate match and can run to hundreds of millions of rows, so it is
streamed in chunks and reduced to two small tables: counts and weight sums per score x
chromosome x variant class x match outcome, and each score's largest-effect unmatched variants.
Everything downstream is derived from those, cached on disk keyed by the log's path, size and
modification time.
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
from pgs_mcp.qc.models import (
    N_TOP_UNMATCHED,
    ChromCoverage,
    LogMetrics,
    UnmatchedVariant,
    WeightMetrics,
)

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
VARIANT_COLUMNS = ["accession", "chr_name", "chr_position", "effect_allele", "other_allele"]
CACHE_VERSION = 2  # bump when a cached table's shape or meaning changes


@dataclass
class LogAggregate:
    groups: pl.DataFrame  # GROUP_KEYS + n, abs_w_sum, abs_w_max
    top_unmatched: pl.DataFrame  # VARIANT_COLUMNS + abs_w, up to N_TOP_UNMATCHED per score


def aggregate_log(path: Path, progress: ProgressFn | None = None) -> LogAggregate:
    """Return the aggregate tables for a match log, from cache when the file is unchanged."""
    groups_path, top_path = _cache_paths(path)
    if groups_path.exists() and top_path.exists():
        return LogAggregate(pl.read_parquet(groups_path), pl.read_parquet(top_path))

    partials = [_aggregate_chunk(chunk) for chunk in _iter_csv_chunks(path, progress)]
    if not partials:
        raise QCInputError(f"Match log {path} contains no rows.")
    result = LogAggregate(
        groups=pl.concat(p.groups for p in partials)
        .group_by(GROUP_KEYS)
        .agg(pl.col("n").sum(), pl.col("abs_w_sum").sum(), pl.col("abs_w_max").max()),
        top_unmatched=_top_per_score(pl.concat(p.top_unmatched for p in partials)),
    )

    groups_path.parent.mkdir(parents=True, exist_ok=True)
    for frame, target in ((result.groups, groups_path), (result.top_unmatched, top_path)):
        tmp = target.with_suffix(".tmp")
        frame.write_parquet(tmp)
        tmp.replace(target)
    return result


def _cache_paths(path: Path) -> tuple[Path, Path]:
    stat = path.stat()
    key = f"{path.resolve()}|{stat.st_size}|{stat.st_mtime_ns}|v{CACHE_VERSION}"
    root = Path(
        os.environ.get("PGS_MCP_CACHE_DIR")
        or Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "pgs-mcp"
    )
    stem = f"log-{hashlib.sha256(key.encode()).hexdigest()[:20]}"
    return root / f"{stem}-groups.parquet", root / f"{stem}-top-unmatched.parquet"


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


def _aggregate_chunk(chunk: bytes) -> LogAggregate:
    rows = (
        pl.read_csv(chunk, columns=list(COLUMNS), infer_schema=False)
        .lazy()
        .with_columns(abs_w=pl.col("effect_weight").cast(pl.Float64, strict=False).abs())
    )
    groups = (
        rows.with_columns(
            is_snv=(pl.col("effect_allele").str.len_chars().fill_null(0) <= 1)
            & (pl.col("other_allele").str.len_chars().fill_null(0) <= 1),
            missing_position=pl.col("chr_name").is_null() | pl.col("chr_position").is_null(),
        )
        .group_by(GROUP_KEYS)
        .agg(n=pl.len(), abs_w_sum=pl.col("abs_w").sum(), abs_w_max=pl.col("abs_w").max())
    )
    unmatched = rows.filter(pl.col("match_status") == "unmatched").select(*VARIANT_COLUMNS, "abs_w")
    groups_df, unmatched_df = pl.collect_all([groups, unmatched])
    return LogAggregate(groups_df, _top_per_score(unmatched_df))


def _top_per_score(variants: pl.DataFrame) -> pl.DataFrame:
    return (
        variants.sort("abs_w", descending=True, nulls_last=True)
        .group_by("accession", maintain_order=True)
        .head(N_TOP_UNMATCHED)
    )


@dataclass
class _ScoreAccumulator:
    n_variants: int = 0
    w_total: float = 0.0
    w_matched: float = 0.0
    w_max: float | None = None
    top_unmatched: list[UnmatchedVariant] = field(default_factory=list)
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
                largest_unmatched=self.top_unmatched,
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


def summarize_log(agg: LogAggregate) -> LogSummary:
    """Derive per-score and run-wide metrics from the aggregate tables."""
    acc: dict[str, _ScoreAccumulator] = defaultdict(_ScoreAccumulator)
    for row in agg.top_unmatched.iter_rows(named=True):
        acc[row["accession"]].top_unmatched.append(
            UnmatchedVariant(
                chrom=row["chr_name"],
                position=int(pos) if (pos := row["chr_position"]) and pos.isdigit() else None,
                effect_allele=row["effect_allele"],
                other_allele=row["other_allele"],
                abs_weight=row["abs_w"] or 0.0,
            )
        )
    for row in agg.groups.iter_rows(named=True):
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
