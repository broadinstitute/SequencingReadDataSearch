"""Command-line interface for aligned-read-search."""

from __future__ import annotations

import json
import logging

import click

from .search import search_phenotype


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("phenotype")
@click.option(
    "-a",
    "--archive",
    "archives",
    default="ena,sra",
    help="Comma-separated archives to query (ena,sra).",
)
@click.option(
    "-l", "--limit", default=50, show_default=True, help="Max results per archive."
)
@click.option(
    "-e",
    "--exact",
    is_flag=True,
    help="Match the literal term (disable MONDO ontology expansion).",
)
@click.option(
    "-u",
    "--include-unaligned",
    is_flag=True,
    help="Also include runs without a BAM/CRAM alignment file (default: aligned only).",
)
@click.option(
    "-s",
    "--library-strategy",
    default=None,
    help="Filter by library strategy, comma-separated (e.g. WGS,WXS).",
)
@click.option(
    "-j", "--json", "as_json", is_flag=True, help="Emit JSON instead of a table."
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Print progress (what's being expanded/searched/enriched) to stderr.",
)
@click.option(
    "--no-enrich",
    is_flag=True,
    help="Skip the ENA alignment-status lookup for SRA rows.",
)
def main(
    phenotype,
    archives,
    limit,
    exact,
    include_unaligned,
    library_strategy,
    as_json,
    verbose,
    no_enrich,
):
    """Find sequencing datasets matching PHENOTYPE on SRA and ENA.

    By default only runs with a BAM/CRAM alignment file are returned; pass
    --include-unaligned to see all matches. PHENOTYPE may be a disease name
    (e.g. "Ataxia") or a MONDO id (e.g. MONDO:0100254).
    """
    logging.basicConfig(
        level=logging.INFO if verbose else logging.WARNING,
        format="%(levelname)s: %(message)s",
    )
    archive_list = [a.strip() for a in archives.split(",") if a.strip()]
    df = search_phenotype(
        phenotype,
        archives=archive_list,
        exact=exact,
        limit=limit,
        enrich_alignment=not no_enrich,
        include_unaligned=include_unaligned,
        library_strategy=library_strategy,
    )

    if as_json:
        click.echo(json.dumps(df.to_dict("records"), indent=2, default=str))
        return

    if df.empty:
        click.echo("No datasets found.")
        return

    n_aligned = int(df["has_alignment"].sum())
    display_cols = [
        "run_accession",
        "source_archive",
        "organism",
        "library_strategy",
        "has_alignment",
        "phenotype_match",
        "title",
    ]
    table = df[display_cols].copy()
    table["title"] = table["title"].str.slice(0, 60)
    with pandas_display_opts():
        click.echo(table.to_string(index=False))
    click.echo(
        f"\n{len(df)} datasets ({n_aligned} with alignment files) "
        f"from {', '.join(archive_list)}."
    )


def pandas_display_opts():
    import pandas as pd

    return pd.option_context(
        "display.max_rows", None, "display.max_colwidth", 60, "display.width", 200
    )


if __name__ == "__main__":
    main()
