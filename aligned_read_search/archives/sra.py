"""SRA client - wraps pysradb's SraSearch.

pysradb passes the query straight into NCBI's ``term`` (All Fields search) and
honours boolean ``OR``, so a single request covers every expansion token with
good recall. Alignment-file presence is not exposed by SRA's metadata here; it
is filled in later by ``alignment.fetch_alignment_info`` via ENA (which mirrors
INSDC runs).
"""

from __future__ import annotations

import logging
from typing import List, Optional

from ..models import Dataset
from .base import ArchiveClient, all_token_matches

logger = logging.getLogger(__name__)


def _pick(row, *names) -> str:
    for n in names:
        if n in row and row[n] not in (None, ""):
            return str(row[n])
    return ""


def _pick_int(row, *names) -> Optional[int]:
    raw = _pick(row, *names)
    try:
        return int(float(raw)) if raw else None
    except (TypeError, ValueError):
        return None


class SraClient(ArchiveClient):
    name = "sra"

    def __init__(self, verbosity: int = 2):
        self.verbosity = verbosity

    def _run_search(self, query: str, limit: int):
        from pysradb.search import SraSearch

        instance = SraSearch(verbosity=self.verbosity, return_max=limit, query=query)
        instance.search()
        return instance.get_df()

    def search(self, tokens: List[str], limit: int = 50) -> List[Dataset]:
        tokens = [t for t in tokens if t]
        if not tokens:
            return []
        # NCBI All-Fields search; quote multi-char tokens, OR them together.
        query = " OR ".join(tokens)
        logger.info("Querying NCBI SRA via pysradb (limit %d): %s", limit, query)
        df = self._run_search(query, limit)
        if df is None or df.empty:
            return []

        datasets = []
        for row in df.to_dict("records"):
            run_acc = _pick(row, "run_accession", "run_1_accession")
            if not run_acc:
                continue
            title = _pick(row, "experiment_title", "study_study_title")
            organism = _pick(row, "sample_scientific_name", "organism_name", "common_name")
            match_text = " ".join(
                _pick(row, c)
                for c in (
                    "experiment_title",
                    "study_study_title",
                    "study_study_abstract",
                )
            )
            datasets.append(
                Dataset(
                    run_accession=run_acc,
                    source_archive=self.name,
                    study_accession=_pick(row, "study_accession"),
                    sample_accession=_pick(row, "sample_accession"),
                    experiment_accession=_pick(row, "experiment_accession"),
                    organism=organism,
                    platform=_pick(row, "experiment_platform"),
                    instrument_model=_pick(row, "experiment_instrument_model"),
                    library_strategy=_pick(row, "experiment_library_strategy"),
                    title=title,
                    phenotype_match=all_token_matches(match_text, tokens),
                    has_alignment=False,  # enriched via ENA later
                    alignment_urls=[],
                    read_urls=[],
                    read_count=_pick_int(row, "run_total_spots", "total_spots"),
                    base_count=_pick_int(row, "run_total_bases", "total_bases"),
                )
            )
        return datasets
