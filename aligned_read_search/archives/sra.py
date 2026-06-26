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

# Descriptive columns where a disease/phenotype legitimately appears. We avoid
# scanning every column (urls/filenames/md5s) so short acronym tokens don't
# match spuriously. ``sample_attributes_*_value`` columns are added dynamically.
_MATCH_FIELDS = (
    "experiment_title",
    "study_study_title",
    "study_study_abstract",
    "sample_title",
    "experiment_design_description",
    "experiment_library_name",
    "sample_alias",
)


def _match_text(row) -> str:
    parts = [_pick(row, c) for c in _MATCH_FIELDS]
    parts += [
        str(v)
        for k, v in row.items()
        if k.startswith("sample_attributes_") and k.endswith("_value") and v
    ]
    return " ".join(p for p in parts if p)


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

    def __init__(self, verbosity: int = 3):
        # verbosity=3 makes SraSearch return the descriptive columns
        # (study_study_title, study_study_abstract, sample_title,
        # sample_attributes_*) where the phenotype actually lives; lower
        # verbosity omits them entirely.
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
        skipped = 0
        for row in df.to_dict("records"):
            run_acc = _pick(row, "run_accession", "run_1_accession")
            if not run_acc:
                continue
            # Require a visible justification: NCBI All-Fields can match on
            # fields we don't surface, so drop runs with no token in our
            # descriptive text rather than return an empty phenotype_match.
            matched = all_token_matches(_match_text(row), tokens)
            if not matched:
                skipped += 1
                continue
            title = _pick(row, "experiment_title", "study_study_title")
            organism = _pick(row, "sample_scientific_name", "organism_name", "common_name")
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
                    phenotype_match=matched,
                    has_alignment=False,  # enriched via ENA later
                    alignment_urls=[],
                    read_urls=[],
                    read_count=_pick_int(row, "run_1_total_spots", "run_total_spots", "total_spots"),
                    base_count=_pick_int(row, "run_1_total_bases", "run_total_bases", "total_bases"),
                )
            )
        if skipped:
            logger.info("Dropped %d SRA run(s) with no phenotype match in descriptive fields", skipped)
        return datasets
