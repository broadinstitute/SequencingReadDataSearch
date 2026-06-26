"""Offline unit tests (no network)."""

import logging

import pandas as pd
import pytest

from aligned_read_search.alignment import _alignment_from_row
from aligned_read_search.archives.base import all_token_matches
from aligned_read_search.archives.ena import (
    EnaClient,
    _split_urls,
    build_query,
)
from aligned_read_search.archives.sra import SraClient
from aligned_read_search.models import Dataset
from aligned_read_search.query import IdentityExpander, tokenize
from aligned_read_search import search as search_mod
from aligned_read_search.search import _dedup, _enrich, search_phenotype, to_dataframe


def test_identity_expander():
    assert IdentityExpander().expand("  Ataxia ") == ["Ataxia"]
    assert IdentityExpander().expand("   ") == []


def test_tokenize_drops_stopwords_and_numbers():
    toks = tokenize("autosomal recessive spinocerebellar ataxia 12")
    assert "spinocerebellar" in toks
    assert "ataxia" in toks
    assert "autosomal" not in toks
    assert "recessive" not in toks
    assert "12" not in toks


def test_tokenize_keeps_known_short():
    assert "sca" in tokenize("SCA")
    assert "frda" in tokenize("FRDA")


def test_all_token_matches():
    assert (
        all_token_matches("Spinocerebellar ataxia cohort", ["foo", "ataxia", "cerebellar"])
        == "ataxia, cerebellar"
    )
    assert all_token_matches("nothing here", ["ataxia"]) == ""
    # De-duplicates case-insensitively, preserves token order.
    assert all_token_matches("ATAXIA ataxia", ["ataxia", "Ataxia"]) == "ataxia"


def test_split_urls():
    assert _split_urls("") == []
    out = _split_urls("ftp.x/a.bam;ftp.x/b.bam")
    assert out == ["https://ftp.x/a.bam", "https://ftp.x/b.bam"]


def test_build_query_multifield_singletoken():
    q = build_query(["ataxia", "friedreich"])
    assert 'study_title="*ataxia*"' in q
    assert 'experiment_title="*friedreich*"' in q
    assert " OR " in q


def test_ena_row_to_dataset_detects_bam():
    row = {
        "run_accession": "ERR1",
        "study_accession": "ERP1",
        "scientific_name": "Homo sapiens",
        "library_strategy": "WGS",
        "study_title": "Friedreich ataxia genomes",
        "submitted_format": "BAM",
        "submitted_ftp": "ftp.sra.ebi.ac.uk/x.bam",
        "bam_ftp": "",
        "fastq_ftp": "ftp.sra.ebi.ac.uk/x.fastq.gz",
        "read_count": "1000",
        "base_count": "150000",
    }
    ds = EnaClient()._row_to_dataset(row, ["ataxia", "friedreich"])
    assert ds.has_alignment is True
    assert ds.alignment_urls == ["https://ftp.sra.ebi.ac.uk/x.bam"]
    assert ds.read_urls == ["https://ftp.sra.ebi.ac.uk/x.fastq.gz"]
    # Title "Friedreich ataxia genomes" matches both tokens, comma-joined.
    assert ds.phenotype_match == "ataxia, friedreich"
    assert ds.read_count == 1000


def test_ena_row_no_alignment():
    row = {"run_accession": "ERR2", "submitted_format": "FASTQ", "study_title": "ataxia"}
    ds = EnaClient()._row_to_dataset(row, ["ataxia"])
    assert ds.has_alignment is False
    assert ds.alignment_urls == []


def test_sra_mapper_via_fake_df():
    import pandas as pd

    class FakeSra(SraClient):
        def _run_search(self, query, limit):
            return pd.DataFrame(
                [
                    {
                        "run_accession": "SRR1",
                        "study_accession": "SRP1",
                        "sample_scientific_name": "Homo sapiens",
                        "experiment_library_strategy": "WXS",
                        "experiment_title": "Exome of ataxia patient",
                        "run_total_spots": "500",
                        "run_total_bases": "75000",
                    }
                ]
            )

    out = FakeSra().search(["ataxia"], limit=5)
    assert len(out) == 1
    assert out[0].run_accession == "SRR1"
    assert out[0].library_strategy == "WXS"
    assert out[0].read_count == 500
    assert out[0].has_alignment is False
    assert out[0].phenotype_match == "ataxia"


def test_dedup_prefers_alignment():
    a = Dataset(run_accession="R", source_archive="sra", has_alignment=False)
    b = Dataset(
        run_accession="R",
        source_archive="ena",
        has_alignment=True,
        alignment_urls=["x.bam"],
    )
    out = _dedup([a, b])
    assert len(out) == 1
    assert out[0].source_archive == "ena"
    assert out[0].has_alignment is True


def test_to_dataframe_empty_has_columns():
    df = to_dataframe([])
    assert "run_accession" in df.columns
    assert "has_alignment" in df.columns
    assert df.empty


# --- alignment read_urls back-fill -------------------------------------------


def test_alignment_from_row_read_urls():
    info = _alignment_from_row({"fastq_ftp": "ftp.x/a.fastq.gz", "submitted_format": "FASTQ"})
    assert info["read_urls"] == ["https://ftp.x/a.fastq.gz"]
    assert info["has_alignment"] is False


def test_enrich_backfills_read_urls(monkeypatch):
    import aligned_read_search.alignment as aln

    def fake_info(accs, **kw):
        return {
            "SRR1": {
                "has_alignment": True,
                "alignment_urls": ["https://x/a.bam"],
                "read_urls": ["https://x/a.fastq.gz"],
            }
        }

    monkeypatch.setattr(aln, "fetch_alignment_info", fake_info)
    d = Dataset(run_accession="SRR1", source_archive="sra", has_alignment=False)
    _enrich([d])
    assert d.has_alignment is True
    assert d.alignment_urls == ["https://x/a.bam"]
    assert d.read_urls == ["https://x/a.fastq.gz"]


# --- search_phenotype filtering & failure isolation --------------------------


def _make_client(datasets=(), fail=False):
    """Build a fake archive client class for search._CLIENTS monkeypatching."""
    rows = list(datasets)

    class _Fake:
        def search(self, tokens, limit=50):
            if fail:
                raise RuntimeError("backend down")
            return [
                Dataset(**{k: v for k, v in d.to_dict().items()}) for d in rows
            ]

    return _Fake


def _ds(run, aligned=False, strategy="WGS"):
    return Dataset(
        run_accession=run,
        source_archive="ena",
        has_alignment=aligned,
        library_strategy=strategy,
    )


def test_aligned_only_is_default(monkeypatch):
    clients = {
        "ena": _make_client([_ds("ERR1", aligned=True), _ds("ERR2", aligned=False)])
    }
    monkeypatch.setattr(search_mod, "_CLIENTS", clients)
    monkeypatch.setattr(search_mod, "_enrich", lambda results: None)

    df = search_phenotype("ataxia", archives=("ena",), exact=True)
    assert list(df["run_accession"]) == ["ERR1"]

    df_all = search_phenotype("ataxia", archives=("ena",), exact=True, include_unaligned=True)
    assert set(df_all["run_accession"]) == {"ERR1", "ERR2"}


def test_library_strategy_filter(monkeypatch):
    clients = {
        "ena": _make_client(
            [
                _ds("ERR1", aligned=True, strategy="WGS"),
                _ds("ERR2", aligned=True, strategy="WXS"),
                _ds("ERR3", aligned=True, strategy="RNA-Seq"),
            ]
        )
    }
    monkeypatch.setattr(search_mod, "_CLIENTS", clients)
    monkeypatch.setattr(search_mod, "_enrich", lambda results: None)

    df = search_phenotype("ataxia", archives=("ena",), exact=True, library_strategy="wgs")
    assert list(df["run_accession"]) == ["ERR1"]

    df2 = search_phenotype(
        "ataxia", archives=("ena",), exact=True, library_strategy="WGS,WXS"
    )
    assert set(df2["run_accession"]) == {"ERR1", "ERR2"}


def test_no_enrich_forced_on_for_aligned_only(monkeypatch, caplog):
    calls = []
    monkeypatch.setattr(search_mod, "_CLIENTS", {"ena": _make_client([_ds("ERR1", aligned=True)])})
    monkeypatch.setattr(search_mod, "_enrich", lambda results: calls.append(True))

    with caplog.at_level(logging.WARNING):
        search_phenotype("ataxia", archives=("ena",), exact=True, enrich_alignment=False)
    assert calls == [True]  # enrichment was forced back on
    assert any("enrichment" in r.message for r in caplog.records)


def test_archive_failure_is_isolated(monkeypatch, caplog):
    clients = {
        "ena": _make_client([_ds("ERR1", aligned=True)]),
        "sra": _make_client(fail=True),
    }
    monkeypatch.setattr(search_mod, "_CLIENTS", clients)
    monkeypatch.setattr(search_mod, "_enrich", lambda results: None)

    with caplog.at_level(logging.WARNING):
        df = search_phenotype("ataxia", archives=("ena", "sra"), exact=True)
    assert list(df["run_accession"]) == ["ERR1"]
    assert any("sra" in r.message and "failed" in r.message for r in caplog.records)


def test_all_archives_failing_raises(monkeypatch):
    clients = {"ena": _make_client(fail=True), "sra": _make_client(fail=True)}
    monkeypatch.setattr(search_mod, "_CLIENTS", clients)
    with pytest.raises(RuntimeError, match="All archives failed"):
        search_phenotype("ataxia", archives=("ena", "sra"), exact=True)
