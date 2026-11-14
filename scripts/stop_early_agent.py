"""Stop-early baseline as a Harbor custom agent (no model involved).

It connects to the task's first HTTP MCP server from inside the agent
container, performs a fixed, short list of tool calls and quits. Default: one
`get_brief` call and nothing else, i.e. "looked at the task and stopped".

    PYTHONPATH=scripts harbor run -p tasks/_template \
        -a stop_early_agent:StopEarlyAgent -o jobs

Other call lists through an agent kwarg (JSON list of [tool, arguments]):

    --ak 'calls=[["get_brief",{}],["submit_regimen",{"regimen":{"phases":[{"days":84,"doses":{"MEKi":1}}]}}]]'     # blind early submit

The calls run with the Python `mcp` client that the task image ships (the same
one solution/solve.py uses), as the task's agent user. Output and tool results
land in <trial>/agent/stop-early.txt.
"""

from __future__ import annotations

import json
import shlex

from harbor.agents.base import BaseAgent
from harbor.agents.installed.base import NonZeroAgentExitCodeError
from harbor.environments.base import BaseEnvironment
from harbor.models.agent.context import AgentContext

_CLIENT = r"""
import asyncio, json, sys
from mcp.client.session import ClientSession
try:
    from mcp.client.streamable_http import streamablehttp_client
except ImportError:
    from mcp.client.streamable_http import streamable_http_client as streamablehttp_client

async def main(url, calls):
    async with streamablehttp_client(url) as streams:
        async with ClientSession(streams[0], streams[1]) as session:
            await session.initialize()
            print("tools:", [t.name for t in (await session.list_tools()).tools], flush=True)
            for tool, args in calls:
                r = await session.call_tool(tool, args)
                text = "".join(c.text for c in r.content if getattr(c, "type", "") == "text")
                print(json.dumps({"tool": tool, "args": args, "is_error": r.isError, "result": text[:2000]}), flush=True)

asyncio.run(main(sys.argv[1], json.loads(sys.argv[2])))
"""


class StopEarlyAgent(BaseAgent):
    def __init__(self, *args, calls: list | str | None = None, **kwargs):
        super().__init__(*args, **kwargs)
        if isinstance(calls, str):
            calls = json.loads(calls)
        self._calls = calls if calls is not None else [["get_brief", {}]]

    @staticmethod
    def name() -> str:
        return "stop-early"

    def version(self) -> str:
        return "1.0.0"

    async def setup(self, environment: BaseEnvironment) -> None:
        return

    async def run(self, instruction: str, environment: BaseEnvironment, context: AgentContext) -> None:
        urls = [s.url for s in self.mcp_servers if s.url]
        if not urls:
            raise RuntimeError("task declares no HTTP MCP server in [[environment.mcp_servers]]")
        command = (
            f"python3 -c {shlex.quote(_CLIENT)} {shlex.quote(urls[0])} "
            f"{shlex.quote(json.dumps(self._calls))} 2>&1 | tee /logs/agent/stop-early.txt"
        )
        result = await environment.exec(command=command, timeout_sec=120)
        if result.return_code != 0:
            # NonZeroAgentExitCodeError keeps the verifier running; any other
            # exception type makes Harbor skip verification (reward = None).
            raise NonZeroAgentExitCodeError(f"stop-early client failed: {result.stdout or ''}{result.stderr or ''}")
