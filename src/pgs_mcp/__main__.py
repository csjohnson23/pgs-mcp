"""Entry point. `pgs-mcp` serves MCP over stdio; `pgs-mcp qc ...` runs match QC from the shell."""

import argparse
import sys
from pathlib import Path

from pgs_mcp.qc import QCInputError, assess_match_qc
from pgs_mcp.report import render_report


def main(argv: list[str] | None = None) -> None:
    argv = sys.argv[1:] if argv is None else argv
    if not argv or argv[0] != "qc":
        from pgs_mcp.server import mcp

        mcp.run()
        return

    parser = argparse.ArgumentParser(prog="pgs-mcp qc", description="Run pgsc_calc match QC.")
    parser.add_argument("summary", type=Path, help="pgsc_calc match summary CSV")
    parser.add_argument("--log", type=Path, help="variant-level match log (.csv or .csv.gz)")
    parser.add_argument("--pgs-id", action="append", dest="pgs_ids", help="repeatable")
    parser.add_argument("--report", action="store_true", help="print Markdown instead of JSON")
    args = parser.parse_args(argv[1:])

    def progress(done: int, total: int) -> None:
        print(f"\rReading match log: {100 * done / total:5.1f}%", end="", file=sys.stderr)

    try:
        run = assess_match_qc(args.summary, args.log, args.pgs_ids, progress)
    except QCInputError as e:
        parser.exit(2, f"error: {e}\n")
    if args.log:
        print(file=sys.stderr)
    print(render_report(run) if args.report else run.model_dump_json(indent=2))


if __name__ == "__main__":
    main()
