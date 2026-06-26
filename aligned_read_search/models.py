"""Normalized dataset record shared across archives."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import List, Optional


@dataclass
class Dataset:
    """A single sequencing run, normalized across SRA and ENA.

    The v1 fields below are populated for every result. The v2 enrichment
    fields (see plan) are optional and default to empty/None so the record
    can grow without breaking callers.
    """

    # --- identity ---
    run_accession: str
    source_archive: str  # "ena" | "sra"
    study_accession: str = ""
    sample_accession: str = ""
    experiment_accession: str = ""

    # --- core descriptive ---
    organism: str = ""
    platform: str = ""
    instrument_model: str = ""
    library_strategy: str = ""  # WGS, WXS (exome), RNA-Seq, ...
    title: str = ""  # study or experiment title
    phenotype_match: str = ""  # which expansion token(s) matched, + field

    # --- aligned-reads focus ---
    has_alignment: bool = False
    alignment_urls: List[str] = field(default_factory=list)  # bam / cram
    read_urls: List[str] = field(default_factory=list)  # fastq / sra

    # --- light volume metrics (free from ENA) ---
    read_count: Optional[int] = None
    base_count: Optional[int] = None

    def to_dict(self) -> dict:
        return asdict(self)
