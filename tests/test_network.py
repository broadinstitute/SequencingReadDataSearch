"""Live integration tests (hit SRA/ENA/OLS). Run with: pytest -m network."""

import pytest

from aligned_read_search.alignment import has_alignment_for
from aligned_read_search.query import OntologyExpander
from aligned_read_search.search import search_phenotype

pytestmark = pytest.mark.network


def test_ontology_expander_ataxia():
    toks = OntologyExpander().expand("Ataxia")
    assert "ataxia" in toks
    # Descendant subtypes should pull in distinctive tokens.
    assert any(t in toks for t in ("spinocerebellar", "friedreich", "cerebellar"))


def test_ontology_expander_accepts_mondo_id():
    # MONDO:0100254 = hereditary ataxia
    toks = OntologyExpander().expand("MONDO:0100254")
    assert "ataxia" in toks


def test_search_ena_returns_rows():
    df = search_phenotype("ataxia", archives=("ena",), exact=True, limit=10,
                          enrich_alignment=False, include_unaligned=True)
    assert not df.empty
    assert {"run_accession", "has_alignment", "library_strategy"} <= set(df.columns)


def test_expansion_increases_hits():
    exact = search_phenotype("ataxia", archives=("ena",), exact=True, limit=200,
                             enrich_alignment=False, include_unaligned=True)
    expanded = search_phenotype("ataxia", archives=("ena",), limit=200,
                                enrich_alignment=False, include_unaligned=True)
    # Expanded should find at least as many distinct runs as the literal term.
    assert len(expanded) >= len(exact)


def test_alignment_lookup_smoke():
    # Just verify the helper returns the expected shape against a real run.
    info = has_alignment_for("ERR1160846")
    assert set(info) == {"has_alignment", "alignment_urls", "read_urls"}
