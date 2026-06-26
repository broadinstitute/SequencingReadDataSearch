"""Top-level orchestrator: phenotype -> list of datasets across SRA + ENA."""

from __future__ import annotations

import logging
from typing import List, Optional, Sequence

import pandas as pd

from .archives.base import ArchiveClient
from .archives.ena import EnaClient
from .archives.sra import SraClient
from .models import Dataset
from .query import IdentityExpander, OntologyExpander, QueryExpander

logger = logging.getLogger(__name__)

_CLIENTS = {"ena": EnaClient, "sra": SraClient}

# Column order for the returned DataFrame.
_COLUMNS = [
    "run_accession",
    "source_archive",
    "organism",
    "library_strategy",
    "platform",
    "instrument_model",
    "has_alignment",
    "phenotype_match",
    "title",
    "study_accession",
    "sample_accession",
    "experiment_accession",
    "read_count",
    "base_count",
    "alignment_urls",
    "read_urls",
]


def _dedup(datasets: List[Dataset]) -> List[Dataset]:
    """Dedup by run accession; prefer the record that has alignment info."""
    by_run: dict[str, Dataset] = {}
    for ds in datasets:
        key = ds.run_accession
        existing = by_run.get(key)
        if existing is None:
            by_run[key] = ds
            continue
        # Prefer a record that already knows about an alignment / has more URLs.
        better = ds.has_alignment and not existing.has_alignment
        more_urls = len(ds.alignment_urls) > len(existing.alignment_urls)
        if better or more_urls:
            by_run[key] = ds
    return list(by_run.values())


def search_phenotype(
    term: str,
    archives: Sequence[str] = ("ena", "sra"),
    expander: Optional[QueryExpander] = None,
    exact: bool = False,
    limit: int = 50,
    enrich_alignment: bool = True,
    include_unaligned: bool = False,
    library_strategy: Optional[str] = None,
) -> pd.DataFrame:
    """Search SRA and ENA for datasets matching ``term``.

    Parameters
    ----------
    term : phenotype string (e.g. "Ataxia") or a MONDO id (e.g. "MONDO:0100254").
    archives : which archives to query.
    expander : a QueryExpander; defaults to OntologyExpander (or IdentityExpander
        when ``exact=True``).
    exact : opt out of ontology expansion and match the literal term.
    limit : max results requested per archive.
    enrich_alignment : look up alignment files for runs whose status is unknown.
    include_unaligned : by default only runs with a BAM/CRAM alignment file are
        returned; set True to also include unaligned runs.
    library_strategy : optional comma-separated filter (case-insensitive), e.g.
        "WGS" or "WGS,WXS", keeping only runs with a matching library_strategy.

    Notes
    -----
    A single archive that fails (network/import error) is logged and skipped; if
    *all* requested archives fail a RuntimeError is raised. Filters
    (``include_unaligned``, ``library_strategy``) are applied *after* the
    per-archive ``limit``, so the final count can be smaller than ``limit``.
    """
    if expander is None:
        expander = IdentityExpander() if exact else OntologyExpander()

    tokens = expander.expand(term)

    results: List[Dataset] = []
    failures: List[tuple] = []
    for name in archives:
        client_cls = _CLIENTS.get(name)
        if client_cls is None:
            raise ValueError(f"Unknown archive: {name!r} (choose from {list(_CLIENTS)})")
        client: ArchiveClient = client_cls()
        logger.info("Searching %s (limit %d)", name, limit)
        try:
            found = client.search(tokens, limit=limit)
        except Exception as exc:  # noqa: BLE001 - isolate per-archive backend failures
            logger.warning("Archive %r search failed: %s", name, exc)
            failures.append((name, exc))
            continue
        logger.info("%s returned %d runs", name, len(found))
        results.extend(found)
    if failures and len(failures) == len(archives):
        detail = "; ".join(f"{n}: {e}" for n, e in failures)
        raise RuntimeError(f"All archives failed: {detail}")

    deduped = _dedup(results)
    logger.info("Deduplicated %d -> %d runs", len(results), len(deduped))
    results = deduped

    # Aligned-only filtering depends on alignment status, which for SRA rows is
    # only known after enrichment. Disabling enrichment while filtering to
    # aligned-only would silently drop every SRA row, so force it back on.
    if not include_unaligned and not enrich_alignment:
        logger.warning(
            "aligned-only results require alignment enrichment; "
            "ignoring enrich_alignment=False"
        )
        enrich_alignment = True

    if enrich_alignment:
        _enrich(results)

    if library_strategy:
        wanted = {s.strip().upper() for s in library_strategy.split(",") if s.strip()}
        kept = [d for d in results if d.library_strategy.upper() in wanted]
        logger.info(
            "library_strategy filter %s: %d -> %d runs",
            sorted(wanted),
            len(results),
            len(kept),
        )
        results = kept

    if not include_unaligned:
        kept = [d for d in results if d.has_alignment]
        logger.info("aligned-only filter: %d -> %d runs", len(results), len(kept))
        results = kept

    return to_dataframe(results)


def _enrich(datasets: List[Dataset]) -> None:
    """Fill in alignment info for runs that don't yet have it (e.g. SRA rows)."""
    from .alignment import fetch_alignment_info

    unknown = [d for d in datasets if not d.has_alignment]
    if not unknown:
        return
    logger.info("Enriching alignment status for %d runs via ENA", len(unknown))
    info = fetch_alignment_info(d.run_accession for d in unknown)
    for d in unknown:
        hit = info.get(d.run_accession)
        if not hit:
            continue
        if hit["has_alignment"]:
            d.has_alignment = True
            d.alignment_urls = hit["alignment_urls"]
        if not d.read_urls and hit.get("read_urls"):
            d.read_urls = hit["read_urls"]


def to_dataframe(datasets: List[Dataset]) -> pd.DataFrame:
    if not datasets:
        return pd.DataFrame(columns=_COLUMNS)
    df = pd.DataFrame([d.to_dict() for d in datasets])
    return df.reindex(columns=_COLUMNS)
