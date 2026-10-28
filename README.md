# pgs-mcp

An [MCP](https://modelcontextprotocol.io) server that tells you how much to trust a polygenic
score before you use it.

The first release covers **variant matching QC for [pgsc_calc](https://github.com/PGScatalog/pgsc_calc)**:
how faithfully each score was reproduced in your genotype data, what was lost, and why. Score
context from the PGS Catalog and source publications is next.

## Why

pgsc_calc gives each score a pass/fail based on one number, the share of variants matched. That
number hides the things that matter:

- **Which variants were lost.** A score can match 72% of its variants and pass while keeping only
  16% of its total effect weight, because the variants that failed to match were the
  large-effect ones.
- **Why they were lost.** Missing chromosomes, no indels in the target, unharmonised positions
  and duplicate IDs each leave a distinct signature, and most are fixable upstream.
- **What was assumed.** Strand-ambiguous variants retained, strand flips, single-allele matches.

pgs-mcp reports all of this as raw numbers plus structural warnings. It never grades a score.
The interpretation is written by your MCP client's model, guided by a reviewable
[interpretation guide](src/pgs_mcp/guides/match_qc.md) shipped with the server.

## What it exposes

| Kind | Name | Purpose |
|---|---|---|
| Tool | `assess_match_qc` | Raw match metrics, warnings and observations for one score or a whole run |
| Tool | `render_qc_report` | Markdown report for one score or a whole run, with optional interpretation |
| Resource | `pgs://guides/match-qc` | What each metric and warning means, likely causes, how to interpret |
| Prompt | `interpret_match_qc` | Runs the assessment and writes a guided interpretation |

### Inputs

Both tools take explicit paths to pgsc_calc outputs:

- `summary_path` (required): the match summary, `{sampleset}_summary.csv`. Fast.
- `log_path` (optional): the variant-level match log, `{sampleset}_log.csv.gz`. Adds effect
  weight coverage, per-chromosome coverage, indel/SNV breakdown and match types. Logs can be
  several GB; the first read streams the file (about a minute for a 3.5 GB log) and caches a
  small aggregate, so later calls are instant.
- `pgs_ids` (optional): limit the report to specific scores. Large runs return a compact index
  of every score, and full detail for the scores you name.

### What you get back

Per score: variant match rate, unambiguous match rate, strand-ambiguous and flipped shares,
exclusions by reason, effect weight coverage, the largest-effect unmatched variants (position
and alleles), SNV and indel match
counts, match types and per-chromosome coverage. Per run: the inferred `min_overlap`,
`keep_ambiguous` and `keep_multiallelic` settings and run-wide coverage.

Warnings flag structural problems and almost all are categorical, with no tunable cutoff:

| Warning | Fires when |
|---|---|
| `CHROM_ABSENT` | No variant on a chromosome was found in the target (run-wide and per score) |
| `TARGET_LACKS_INDELS` | Not a single indel matched anywhere in the run |
| `INDELS_UNMATCHED` | All of a score's indels failed to match |
| `MISSING_POSITION` | Score variants have no chromosome or position |
| `DUPLICATES_EXCLUDED` | Variants were dropped for duplicate IDs or tied best matches |
| `AMBIGUOUS_RETAINED` | Strand-ambiguous variants were kept |
| `LOG_SUMMARY_MISMATCH` | The log and summary appear to come from different runs |
| `LOG_NOT_SUPPLIED` | Only the summary was given, so weight-based metrics are unavailable |
| `SMALL_SCORE` | Fewer than 20 variants (cutoff reported in evidence) |
| `BUILD_MISMATCH_SUSPECTED` | Median match rate below 20% (cutoff reported in evidence) |

## Install

Requires [uv](https://docs.astral.sh/uv/).

### Claude Code

```bash
claude mcp add pgs-mcp -- uvx --from git+https://github.com/csjohnson23/pgs-mcp pgs-mcp
```

### Claude Desktop and other MCP clients

```json
{
  "mcpServers": {
    "pgs-mcp": {
      "command": "uvx",
      "args": ["--from", "git+https://github.com/csjohnson23/pgs-mcp", "pgs-mcp"]
    }
  }
}
```

Everything runs locally over stdio. Your files never leave your machine and the server makes no
network calls.

## Use

Ask your assistant, for example:

> Assess the match QC for `/data/run1/PROFILE_summary.csv` with the log
> `/data/run1/PROFILE_log.csv.gz`, and explain which scores I should be careful with.

Or run the `interpret_match_qc` prompt. MCP prompts are started by the user, not the model: in
Claude Code type `/mcp__pgs-mcp__interpret_match_qc`; in Claude Desktop it is listed under the
connected server in the message box. Plain requests work just as well, since the tools point
the model to the interpretation guide.

The same QC runs from the shell, which is also a convenient way to warm the cache for a very
large log:

```bash
uvx --from git+https://github.com/csjohnson23/pgs-mcp pgs-mcp qc summary.csv --log log.csv.gz --report
```

## Development

```bash
uv sync
uv run pytest
uv run ruff check . && uv run ruff format --check .
```

```
src/pgs_mcp/
├── server.py          MCP adapter: tools, resource, prompt. No logic.
├── qc/
│   ├── summary.py     match summary -> per-score variant metrics
│   ├── log.py         streaming log aggregation + on-disk cache
│   ├── rules.py       warnings and observations
│   ├── assess.py      assembles the result
│   ├── models.py      output schema (published to clients as JSON Schema)
│   └── thresholds.py  the only numeric cutoffs
├── report/markdown.py single-score and run reports
└── guides/match_qc.md interpretation guide (served as a resource)
```

Tests use small synthetic fixtures in `tests/fixtures/`, each built so specific warnings fire.
No real genotype-derived data is committed.

The cache lives in `~/.cache/pgs-mcp` (override with `PGS_MCP_CACHE_DIR`).

## Roadmap

- PGS Catalog metadata with explicit gap reporting (ancestry, validation, benchmarks)
- Publication context: score construction method, study purpose, GWAS / training / testing
  cohorts, reported and missing benchmarks
- Per-sample coverage from pgsc_calc scoring outputs
- Discovery of inputs from a pgsc_calc output directory

## License

MIT
