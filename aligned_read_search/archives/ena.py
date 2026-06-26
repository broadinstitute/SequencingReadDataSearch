"""ENA client - queries the ENA Portal API directly.

We do NOT use pysradb's EnaSearch here: it only matches the query against
``experiment_title``, but the disease name overwhelmingly lives in
``study_title`` (3220 hits vs 112 for ``sample_description`` in live tests).
Querying the Portal API directly lets us OR the tokens across several text
fields AND request the rich field set (bam links, submitted format, sizes,
disease, sample attributes) in a single call.
"""

from __future__ import annotations

import logging
from typing import List

import requests

from ..models import Dataset
from .base import ArchiveClient, all_token_matches

logger = logging.getLogger(__name__)

PORTAL = "https://www.ebi.ac.uk/ena/portal/api/search"

# Free-text fields the phenotype tokens are matched against (OR'd).
TEXT_FIELDS = ["study_title", "sample_title", "sample_description", "experiment_title"]

# ENA's Portal API silently returns no rows once a query has too many OR
# clauses (empirically it breaks between ~144 and ~152). We chunk well under
# that so a large expanded token set still searches every field.
MAX_CLAUSES_PER_REQUEST = 100

# Fields requested back from the Portal API (all valid read_run returnFields).
RETURN_FIELDS = [
    "run_accession",
    "study_accession",
    "sample_accession",
    "experiment_accession",
    "scientific_name",
    "instrument_platform",
    "instrument_model",
    "library_strategy",
    "library_source",
    "library_selection",
    "study_title",
    "experiment_title",
    "sample_title",
    "sample_description",
    "read_count",
    "base_count",
    "fastq_ftp",
    "submitted_ftp",
    "submitted_format",
    "bam_ftp",
]

ALIGNMENT_FORMATS = ("BAM", "CRAM")


def _split_urls(value: str) -> List[str]:
    if not value:
        return []
    out = []
    for part in value.split(";"):
        part = part.strip()
        if not part:
            continue
        out.append(part if "://" in part else "https://" + part)
    return out


def _to_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def build_query(tokens: List[str]) -> str:
    """OR every token across every text field, e.g.
    ``study_title="*ataxia*" OR ... OR experiment_title="*friedreich*"``.
    """
    clauses = []
    for tok in tokens:
        safe = tok.replace('"', "").replace("*", "")
        if not safe:
            continue
        for field in TEXT_FIELDS:
            clauses.append(f'{field}="*{safe}*"')
    return " OR ".join(clauses)


class EnaClient(ArchiveClient):
    name = "ena"

    def __init__(self, timeout: float = 60.0, session: requests.Session | None = None):
        self.timeout = timeout
        self.session = session or requests.Session()

    def _request(self, query: str, limit: int) -> List[dict]:
        payload = {
            "result": "read_run",
            "query": query,
            "fields": ",".join(RETURN_FIELDS),
            "format": "json",
            "limit": limit,
        }
        resp = self.session.post(PORTAL, data=payload, timeout=self.timeout)
        resp.raise_for_status()
        if not resp.text.strip():
            return []
        return resp.json()

    def _row_to_dataset(self, row: dict, tokens: List[str]) -> Dataset:
        submitted_fmt = (row.get("submitted_format") or "").upper()
        bam_urls = _split_urls(row.get("bam_ftp", ""))
        submitted_urls = _split_urls(row.get("submitted_ftp", ""))
        aln_urls = list(bam_urls)
        if any(fmt in submitted_fmt for fmt in ALIGNMENT_FORMATS):
            aln_urls.extend(submitted_urls)
        has_alignment = bool(aln_urls)

        match_text = " ".join(row.get(f, "") for f in TEXT_FIELDS)
        matched = all_token_matches(match_text, tokens)

        return Dataset(
            run_accession=row.get("run_accession", ""),
            source_archive=self.name,
            study_accession=row.get("study_accession", ""),
            sample_accession=row.get("sample_accession", ""),
            experiment_accession=row.get("experiment_accession", ""),
            organism=row.get("scientific_name", ""),
            platform=row.get("instrument_platform", ""),
            instrument_model=row.get("instrument_model", ""),
            library_strategy=row.get("library_strategy", ""),
            title=row.get("study_title") or row.get("experiment_title", ""),
            phenotype_match=matched,
            has_alignment=has_alignment,
            alignment_urls=aln_urls,
            read_urls=_split_urls(row.get("fastq_ftp", "")),
            read_count=_to_int(row.get("read_count")),
            base_count=_to_int(row.get("base_count")),
        )

    def search(self, tokens: List[str], limit: int = 50) -> List[Dataset]:
        tokens = [t for t in tokens if t]
        if not tokens:
            return []

        # Chunk tokens so each request stays under ENA's OR-clause limit, then
        # union the runs (deduped) up to `limit`.
        per_chunk = max(1, MAX_CLAUSES_PER_REQUEST // len(TEXT_FIELDS))
        n_chunks = (len(tokens) + per_chunk - 1) // per_chunk
        logger.info(
            "Querying ENA Portal API: %d token(s) across %d field(s) in %d chunk(s)",
            len(tokens),
            len(TEXT_FIELDS),
            n_chunks,
        )
        seen: set[str] = set()
        datasets: List[Dataset] = []
        for i, start in enumerate(range(0, len(tokens), per_chunk), 1):
            if len(datasets) >= limit:
                break
            chunk = tokens[start : start + per_chunk]
            rows = self._request(build_query(chunk), limit)
            logger.info("ENA chunk %d/%d -> %d rows", i, n_chunks, len(rows))
            for row in rows:
                acc = row.get("run_accession")
                if not acc or acc in seen:
                    continue
                seen.add(acc)
                # Match against the full token set, not just this chunk.
                datasets.append(self._row_to_dataset(row, tokens))
                if len(datasets) >= limit:
                    break
        return datasets
