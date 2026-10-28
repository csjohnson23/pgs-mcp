"""End-to-end through the MCP protocol, using the SDK's in-process client."""

import json

import pytest
from mcp.client import Client

from pgs_mcp.server import GUIDE_URI, mcp

pytestmark = pytest.mark.anyio


async def test_lists_tools_resource_and_prompt():
    async with Client(mcp) as client:
        tools = {t.name: t for t in (await client.list_tools()).tools}
        assert set(tools) == {"assess_match_qc", "render_qc_report"}
        assert tools["assess_match_qc"].annotations.read_only_hint
        assert tools["assess_match_qc"].output_schema is not None
        resources = (await client.list_resources()).resources
        assert [str(r.uri) for r in resources] == [GUIDE_URI]
        prompts = (await client.list_prompts()).prompts
        assert [p.name for p in prompts] == ["interpret_match_qc"]
        assert all(arg.description for arg in prompts[0].arguments)


async def test_assess_returns_structured_result(summary, log):
    async with Client(mcp) as client:
        result = await client.call_tool(
            "assess_match_qc", {"summary_path": str(summary), "log_path": str(log)}
        )
    assert not result.is_error
    assert result.structured_content["n_scores"] == 3
    assert result.structured_content["log_used"] is True
    text = result.content[0].text
    assert "\n" not in text  # compact: this is the block clients hand the model
    assert json.loads(text) == result.structured_content
    # Rounded on the wire: 100 * 0.41 / 0.91 = 45.0549...
    detail = {sc["pgs_id"]: sc for sc in result.structured_content["scores"]}
    assert detail["PGS000001"]["log"]["weights"]["weight_coverage_pct"] == 45.05
    assert result.structured_content["inferred_settings"]["min_overlap_upper"] == 0.75


async def test_bad_path_is_a_tool_error(tmp_path):
    async with Client(mcp, raise_exceptions=False) as client:
        result = await client.call_tool(
            "assess_match_qc", {"summary_path": str(tmp_path / "missing.csv")}
        )
    assert result.is_error
    assert "File not found" in result.content[0].text


async def test_report_written_to_file(summary, tmp_path):
    out = tmp_path / "reports" / "qc.md"
    async with Client(mcp) as client:
        result = await client.call_tool(
            "render_qc_report",
            {
                "summary_path": str(summary),
                "interpretation": "Looks fine.",
                "output_path": str(out),
            },
        )
    assert "Report written" in result.content[0].text
    assert "Looks fine." in out.read_text()


async def test_guide_and_prompt():
    async with Client(mcp) as client:
        guide = await client.read_resource(GUIDE_URI)
        assert "Weight coverage" in guide.contents[0].text
        prompt = await client.get_prompt(
            "interpret_match_qc", {"summary_path": "s.csv", "pgs_ids": "PGS000001, PGS000002"}
        )
        text = prompt.messages[0].content.text
        assert "pgs_ids=['PGS000001', 'PGS000002']" in text
        assert GUIDE_URI in text
