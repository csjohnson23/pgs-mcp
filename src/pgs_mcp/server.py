"""MCP server: registers tools, the interpretation guide resource and the interpretation prompt.

All logic lives in `pgs_mcp.qc` and `pgs_mcp.report`; this module only adapts it to MCP.
"""

from importlib.resources import files
from pathlib import Path
from typing import Annotated

import anyio.from_thread
import anyio.to_thread
from mcp.server.mcpserver import Context, MCPServer
from mcp.server.mcpserver.exceptions import ToolError
from mcp.types import CallToolResult, TextContent, ToolAnnotations
from pydantic import Field

from pgs_mcp.qc import QCInputError
from pgs_mcp.qc import assess_match_qc as run_assessment
from pgs_mcp.qc.models import RunQC
from pgs_mcp.report import render_report

GUIDE_URI = "pgs://guides/match-qc"

mcp = MCPServer(
    name="pgs-mcp",
    instructions=(
        "Tools for polygenic score quality control. assess_match_qc reads pgsc_calc variant "
        "matching outputs and returns raw metrics plus structural warnings, never a pass/fail "
        f"grade. Before interpreting its output, read the guide at {GUIDE_URI}, or use the "
        "interpret_match_qc prompt. render_qc_report produces a Markdown report."
    ),
)

SummaryPath = Annotated[
    str,
    Field(description="Path to the pgsc_calc match summary, e.g. `{sampleset}_summary.csv`."),
]
LogPath = Annotated[
    str | None,
    Field(
        description="Optional path to the variant-level match log, e.g. "
        "`{sampleset}_log.csv.gz`. Enables weight, chromosome, indel and match-type metrics. "
        "The first read of a large log takes about a minute per few GB; it is then cached."
    ),
]
PgsIds = Annotated[
    list[str] | None,
    Field(
        description="Optional PGS IDs (e.g. PGS000004) or full accessions to report on. "
        "Omit for every score in the file."
    ),
]

READ_ONLY = ToolAnnotations(read_only_hint=True, open_world_hint=False, idempotent_hint=True)


async def _assess(
    summary_path: str, log_path: str | None, pgs_ids: list[str] | None, ctx: Context
) -> RunQC:
    def progress(done: int, total: int) -> None:
        message = f"Reading match log: {_size(done)} of {_size(total)}"
        anyio.from_thread.run(ctx.report_progress, done, total, message)

    try:
        return await anyio.to_thread.run_sync(
            lambda: run_assessment(
                Path(summary_path), Path(log_path) if log_path else None, pgs_ids, progress
            )
        )
    except QCInputError as e:
        raise ToolError(str(e)) from e  # input problems go back to the model; crashes do not


def _size(n_bytes: int) -> str:
    return f"{n_bytes / 1e9:.1f} GB" if n_bytes >= 1e9 else f"{n_bytes / 1e6:.0f} MB"


@mcp.tool(annotations=READ_ONLY)
async def assess_match_qc(
    summary_path: SummaryPath,
    ctx: Context,
    log_path: LogPath = None,
    pgs_ids: PgsIds = None,
) -> Annotated[CallToolResult, RunQC]:
    """Assess how well polygenic scores matched a target genotype dataset in pgsc_calc.

    Works for a single score or a whole run. Returns raw metrics (variant match rate, weight
    coverage, strand handling, exclusions, per-chromosome coverage), structural warnings and
    descriptive observations. Large runs return a compact index of every score; request
    specific pgs_ids for full detail. Interpret with the guide at pgs://guides/match-qc.
    """
    run = await _assess(summary_path, log_path, pgs_ids, ctx)
    # The text block is what most clients show the model: send it compact, not indented.
    return CallToolResult(
        content=[TextContent(type="text", text=run.model_dump_json())],
        structured_content=run.model_dump(mode="json"),
    )


@mcp.tool(annotations=ToolAnnotations(read_only_hint=False, open_world_hint=False))
async def render_qc_report(
    summary_path: SummaryPath,
    ctx: Context,
    log_path: LogPath = None,
    pgs_ids: PgsIds = None,
    interpretation: Annotated[
        str | None,
        Field(
            description="Optional Markdown interpretation to place at the top of the report, "
            "written after reading the assess_match_qc output and the interpretation guide. "
            "The report adds its own Interpretation heading."
        ),
    ] = None,
    output_path: Annotated[
        str | None,
        Field(description="Optional .md path to write the report to instead of returning it."),
    ] = None,
) -> str:
    """Build a Markdown match QC report for one score or a whole run.

    A single score gets a focused report; several scores get a run overview, distributions,
    chromosome coverage, a table of every score and detail sections. For large runs, pass
    output_path to write the report to disk rather than returning it.
    """
    report = render_report(await _assess(summary_path, log_path, pgs_ids, ctx), interpretation)
    if output_path is None:
        return report
    out = Path(output_path).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(report)
    return f"Report written to {out} ({len(report):,} characters)."


@mcp.resource(
    GUIDE_URI,
    name="match-qc-guide",
    title="Interpreting pgsc_calc match QC",
    description="What each match QC metric and warning means, likely causes, and how to write "
    "an interpretation for one score or a full run.",
    mime_type="text/markdown",
)
def match_qc_guide() -> str:
    return files("pgs_mcp.guides").joinpath("match_qc.md").read_text()


@mcp.prompt(title="Interpret pgsc_calc match QC")
def interpret_match_qc(
    summary_path: Annotated[
        str, Field(description="Path to the pgsc_calc match summary (`*_summary.csv`).")
    ],
    log_path: Annotated[
        str, Field(description="Optional path to the variant-level match log (`*_log.csv.gz`).")
    ] = "",
    pgs_ids: Annotated[
        str, Field(description="Optional comma-separated PGS IDs; blank for every score.")
    ] = "",
) -> str:
    """Assess pgsc_calc match QC and write a guided interpretation of the results."""
    args = [f"summary_path={summary_path!r}"]
    if log_path:
        args.append(f"log_path={log_path!r}")
    ids = [i.strip() for i in pgs_ids.split(",") if i.strip()]
    if ids:
        args.append(f"pgs_ids={ids!r}")
    return f"""Interpret the pgsc_calc variant matching QC for these results.

1. Call `assess_match_qc` with {", ".join(args)}.
2. Read the resource `{GUIDE_URI}` and follow it.
3. Write the interpretation the guide describes: a single-score interpretation if one score is
   reported, otherwise a run-level interpretation that groups scores by pattern. Quote the
   numbers you rely on, keep observed facts separate from inferred causes, and do not invent
   pass/fail cutoffs.
4. Offer to produce a report with `render_qc_report`, passing your interpretation as
   `interpretation`."""
