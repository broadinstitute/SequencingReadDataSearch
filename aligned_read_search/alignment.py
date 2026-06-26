"""Alignment-file detection via ENA (which mirrors INSDC runs from all archives).

Used to enrich SRA-sourced rows (whose alignment status SRA metadata does not
expose) and as a standalone ``has_alignment_for`` helper.
"""

from __future__ import annotations

import logging
from typing import Dict, Iterable, List

import requests

from .archives.ena import ALIGNMENT_FORMATS, PORTAL, _split_urls

logger = logging.getLogger(__name__)

_FIELDS = "run_accession,bam_ftp,submitted_ftp,submitted_format,fastq_ftp"


def _alignment_from_row(row: dict) -> dict:
    submitted_fmt = (row.get("submitted_format") or "").upper()
    urls = _split_urls(row.get("bam_ftp", ""))
    if any(fmt in submitted_fmt for fmt in ALIGNMENT_FORMATS):
        urls.extend(_split_urls(row.get("submitted_ftp", "")))
    return {
        "has_alignment": bool(urls),
        "alignment_urls": urls,
        "read_urls": _split_urls(row.get("fastq_ftp", "")),
    }


def fetch_alignment_info(
    run_accessions: Iterable[str],
    timeout: float = 60.0,
    batch_size: int = 100,
    session: requests.Session | None = None,
) -> Dict[str, dict]:
    """Map each run accession -> {has_alignment, alignment_urls} via ENA."""
    session = session or requests.Session()
    accs = [a for a in dict.fromkeys(run_accessions) if a]  # dedup, keep order
    info: Dict[str, dict] = {}
    if accs:
        n_batches = (len(accs) + batch_size - 1) // batch_size
        logger.info(
            "Looking up alignment files for %d run(s) via ENA (%d batch(es))",
            len(accs),
            n_batches,
        )
    for start in range(0, len(accs), batch_size):
        batch = accs[start : start + batch_size]
        query = " OR ".join(f'run_accession="{a}"' for a in batch)
        payload = {
            "result": "read_run",
            "query": query,
            "fields": _FIELDS,
            "format": "json",
            "limit": len(batch),
        }
        resp = session.post(PORTAL, data=payload, timeout=timeout)
        resp.raise_for_status()
        rows = resp.json() if resp.text.strip() else []
        for row in rows:
            acc = row.get("run_accession")
            if acc:
                info[acc] = _alignment_from_row(row)
    return info


def has_alignment_for(run_accession: str, **kwargs) -> dict:
    """Convenience wrapper for a single run accession."""
    return fetch_alignment_info([run_accession], **kwargs).get(
        run_accession, {"has_alignment": False, "alignment_urls": [], "read_urls": []}
    )
