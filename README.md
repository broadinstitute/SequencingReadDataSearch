# aligned-read-search

Find sequencing datasets on **NCBI SRA** and **EMBL-EBI ENA** by **phenotype**
(e.g. *Ataxia*), with a focus on **aligned reads**. Usable as a Python API, a
CLI, or an MCP tool.

## Why this exists

Neither archive indexes "phenotype" directly — disease lives in free-text
metadata. This tool:

1. **Expands** the phenotype to synonyms + disease subtypes via the **MONDO**
   ontology (EBI OLS), because a literal match misses most data
   (`study_title` "*ataxia*" ≈ 3220 ENA runs; adding "*Friedreich*" alone adds
   ~1500 more). Use `--exact` to opt out.
2. **Searches** both archives — SRA via `pysradb` (NCBI all-fields), ENA via the
   Portal API across `study_title`/`sample_title`/`sample_description`/
   `experiment_title` (highest-recall fields).
3. **Tags** each run with whether a **BAM/CRAM alignment file** exists and its
   download URLs. Because the focus is aligned reads, results are **restricted to
   runs with an alignment file by default**; pass `--include-unaligned` to see all
   matches flagged instead.

## Install

```bash
pip install -e .
# optional MCP server:
pip install -e ".[mcp]"
```

> **Env note:** `pysradb`'s compiled deps (pandas/pyarrow) in some conda
> environments are built against NumPy 1.x. If you see
> `ImportError: numpy.core.multiarray failed to import`, pin `numpy<2`.

## CLI

```bash
aligned-read-search ataxia                      # both archives, aligned reads only
aligned-read-search ataxia --archive ena --limit 20
aligned-read-search ataxia --exact              # literal match, no expansion
aligned-read-search ataxia --include-unaligned  # also show runs without a BAM/CRAM file
aligned-read-search ataxia --library-strategy WGS,WXS   # whole-genome / exome only
aligned-read-search ataxia --json               # machine-readable
aligned-read-search MONDO:0100254               # a MONDO id also works
```

## Python API

```python
from aligned_read_search import search_phenotype

df = search_phenotype("ataxia", archives=("ena", "sra"), limit=50)
df[df.has_alignment][["run_accession", "library_strategy", "alignment_urls"]]
```

`search_phenotype` returns a pandas DataFrame. Swap the matching strategy with
`expander=`:

```python
from aligned_read_search import OntologyExpander, IdentityExpander
search_phenotype("ataxia", expander=IdentityExpander())          # exact
search_phenotype("ataxia", expander=OntologyExpander(max_terms=30))
```

## MCP

```bash
python -m aligned_read_search.mcp_server
```
Exposes one tool,
`search_datasets(phenotype, archives, limit, exact, include_unaligned, library_strategy)`.

## Architecture

| Module | Role |
|---|---|
| `query.py` | `QueryExpander` interface; `OntologyExpander` (default, MONDO/OLS) + `IdentityExpander` |
| `archives/ena.py` | ENA Portal API client (multi-field OR, chunked under the ~145-clause limit) |
| `archives/sra.py` | SRA client wrapping `pysradb.search.SraSearch` |
| `alignment.py` | BAM/CRAM detection via ENA (mirrors INSDC runs from all archives) |
| `search.py` | Orchestrator: expand → query → dedup by run → enrich alignment |
| `cli.py` / `mcp_server.py` | CLI and MCP surfaces over `search_phenotype` |

See `aligned_read_search` plan notes for the v2 roadmap (access tier, file
format, WGS/exome, cloud/https/ftp access methods, publications, sample
attributes, reference genome, coverage).

## Tests

```bash
pytest -m "not network"     # offline unit tests
pytest -m network           # live SRA/ENA/OLS integration
```
