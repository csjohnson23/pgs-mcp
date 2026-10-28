"""Fixtures are synthetic pgsc_calc outputs for dataset TOY (GRCh38), built so each warning fires:

PGS000001  10 variants, pass. 2 ambiguous kept, 1 strand flip, 1 no_oa match, 1 not_best
           candidate. Its one unmatched variant is an indel on chr7 carrying the largest weight.
PGS000002  6 variants, fail (50%). Both chrX variants unmatched, 1 duplicate_ID exclusion.
PGS000003  4 variants, pass (75%). 1 variant with no chromosome/position.
"""

from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture(autouse=True)
def isolated_cache(tmp_path, monkeypatch):
    cache = tmp_path / "cache"
    monkeypatch.setenv("PGS_MCP_CACHE_DIR", str(cache))
    return cache


@pytest.fixture
def summary() -> Path:
    return FIXTURES / "toy_summary.csv"


@pytest.fixture
def single_summary() -> Path:
    return FIXTURES / "toy_single_summary.csv"


@pytest.fixture
def log() -> Path:
    return FIXTURES / "toy_log.csv"


@pytest.fixture
def anyio_backend():
    return "asyncio"
