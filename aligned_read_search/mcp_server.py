"""Minimal MCP server stub exposing search_phenotype as one tool.

This is intentionally thin (v1): it wraps the same ``search_phenotype`` used by
the Python API and CLI. Requires the optional ``mcp`` extra:

    pip install "aligned-read-search[mcp]"
    python -m aligned_read_search.mcp_server
"""

from __future__ import annotations

from .search import search_phenotype


def build_server():
    try:
        from mcp.server.fastmcp import FastMCP
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise SystemExit(
            "The 'mcp' package is required. Install with: "
            "pip install 'aligned-read-search[mcp]'"
        ) from exc

    server = FastMCP("aligned-read-search")

    @server.tool()
    def search_datasets(
        phenotype: str,
        archives: str = "ena,sra",
        limit: int = 50,
        exact: bool = False,
        include_unaligned: bool = False,
        library_strategy: str = "",
    ) -> list[dict]:
        """Find SRA/ENA sequencing datasets matching a phenotype or MONDO id.

        Returns one record per sequencing run, including whether a BAM/CRAM
        alignment file is available and its download URLs. By default only runs
        with an alignment file are returned; set ``include_unaligned`` to include
        all matches. ``library_strategy`` optionally filters by strategy
        (comma-separated, e.g. "WGS,WXS").
        """
        df = search_phenotype(
            phenotype,
            archives=[a.strip() for a in archives.split(",") if a.strip()],
            limit=limit,
            exact=exact,
            include_unaligned=include_unaligned,
            library_strategy=library_strategy or None,
        )
        return df.to_dict("records")

    return server


def main():  # pragma: no cover - entry point
    build_server().run()


if __name__ == "__main__":  # pragma: no cover
    main()
