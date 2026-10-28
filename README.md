# pgs-mcp

An MCP server that evaluates polygenic score matching from [pgsc_calc](https://github.com/PGScatalog/pgsc_calc) in a target population, flags matching warnings and silent errors, and contextualizes findings into a user report.  

The pgsc_calc package passes or fails each polygenic score based on a single number: variant match rate. However, the percentage of variants matched is an incomplete picture of matching accuracy. 

For example: A score can match a threshold-passing percentage of its variants, say 72%, while keeping
only 16% of the total effect weight, because the unmatched variants were the large-effect ones. 

This server reads pgsc_calc's matching outputs and gives the full picture: weight
coverage, which variants were lost, where, and the likely upstream cause, as well as analysis of strand flips, ambiguous alleles, and missing data. It gives recommendations for scores to drop from analysis even if they passed the threshold, as well as upstream data checks to reduce unmatched variants. The tools are adaptable to one score or a run of hundreds. 

## Requirements

- [uv](https://docs.astral.sh/uv/) (`brew install uv`, or see the uv docs)
- An MCP client, such as Claude Code
- Match log and match summary files from a pgsc_calc run (specifically, the match_variants step ahead of pgs calculation)

## Install

**Claude Code**

```bash
claude mcp add pgs-mcp -- uvx --from git+https://github.com/csjohnson23/pgs-mcp pgs-mcp
```

**Other MCP clients:** add this to your client's MCP configuration (for Claude
Desktop, `claude_desktop_config.json`) and restart the client.

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

The server runs locally.

## Find your inputs

pgsc_calc writes both files to the `match/` folder for each sampleset:

```
results/<sampleset>/match/<sampleset>_summary.csv   match summary (required)
results/<sampleset>/match/<sampleset>_log.csv.gz    variant-level match log (recommended)
```

The log is required for effect weight coverage, the
largest lost variants, per-chromosome coverage and indel matching, which is where most problems
show up. The first log read is cached, but expect initial runtime of up to a minute for larger logs (>3GB). 

## Use

Ask your assistant for pgs matching analysis in plain language, with the paths to your files. For example:

> Check the match QC for `results/cohort/match/cohort_summary.csv` with the log
> `results/cohort/match/cohort_log.csv.gz`. Which scores should I be careful with, and why?

> Look at PGS000316 in detail. What was lost, and does it matter?

> Write a QC report for the whole run to `qc_report.md`, with your interpretation at the top.

For a run with many scores you get a run-level picture first (problems shared by every score,
such as missing chromosomes or no indels in your data), then scores grouped by pattern. Ask about
specific scores for full detail.

You can also start the guided interpretation directly: in Claude Code type
`/mcp__pgs-mcp__interpret_match_qc`; in Claude Desktop it is listed under the connected server in
the message box.

## What you get back

The server reports measurements, never grades. Your assistant interprets them using an
[interpretation guide](src/pgs_mcp/guides/match_qc.md) that ships with the server, and is told
to separate what the data shows from what it infers.

**For each score**
- variant match rate, and the rate excluding strand-ambiguous variants
- effect weight coverage: the share of the score's total |effect weight| that matched
- the largest-effect unmatched variants, with position and alleles
- strand-ambiguous, strand-flipped and excluded shares, with exclusion reasons
- SNV and indel matching, match types, and match rate by chromosome
- pgsc_calc's own pass/fail verdict, alongside the above

**For the run:** the pgsc_calc settings implied by the results (`min_overlap`, whether ambiguous
variants were kept), coverage by chromosome across all scores, and how scores are distributed.

**Warnings** flag structural problems:

| Warning | Meaning |
|---|---|
| `CHROM_ABSENT` | No variant on a chromosome was found in your data (run-wide or for one score) |
| `TARGET_LACKS_INDELS` | Not a single indel matched anywhere in the run |
| `INDELS_UNMATCHED` | All of a score's indels failed to match |
| `MISSING_POSITION` | The score contains entries with no genomic position, such as HLA alleles |
| `DUPLICATES_EXCLUDED` | Variants were dropped for duplicate IDs or tied matches |
| `AMBIGUOUS_RETAINED` | Strand-ambiguous (A/T, C/G) variants were kept and scored as same-strand |
| `SMALL_SCORE` | Fewer than 20 variants, so each one matters |
| `BUILD_MISMATCH_SUSPECTED` | Matching is uniformly low (median under 20%) |
| `LOG_SUMMARY_MISMATCH` | The log and summary appear to come from different runs |
| `LOG_NOT_SUPPLIED` | Only the summary was given; weight-based checks were not possible |

**Reports** are Markdown: a focused page for one score, or for a run an overview, distributions,
chromosome coverage and a table of every score sorted by weight coverage.

## Command line

The same checks run without an assistant, which is also a quick way to pre-read a large log:

```bash
uvx --from git+https://github.com/csjohnson23/pgs-mcp pgs-mcp qc cohort_summary.csv --log cohort_log.csv.gz --report > qc_report.md
```

Add `--pgs-id PGS000316` (repeatable) for specific scores, or drop `--report` for JSON.

## Good to know

- **Match QC measures coverage, not validity.** A perfectly matched score can still be poorly
  suited to your population; these checks cannot tell you that.
- **Weight coverage is a proxy.** It uses effect sizes without allele frequencies, so it
  approximates, rather than measures, how much of a score's signal was retained.
- **Your cohort is unknown to the server.** If ancestry, genotyping platform or imputation details
  matter to an interpretation, tell your assistant, or it will ask.
- **Cache:** aggregated log data is stored in `~/.cache/pgs-mcp` (set `PGS_MCP_CACHE_DIR` to
  move it). It is small and safe to delete.

## Contributing

Issues and pull requests are welcome, especially improvements to the
[interpretation guide](src/pgs_mcp/guides/match_qc.md) from PGS practitioners.

```bash
uv sync
uv run pytest
uv run ruff check . && uv run ruff format --check .
```

The QC logic in `src/pgs_mcp/qc/` is plain Python, independent of MCP; `server.py` only exposes
it. Tests use small synthetic pgsc_calc outputs in `tests/fixtures/`; no real genotype-derived data
is committed.

## License

MIT
