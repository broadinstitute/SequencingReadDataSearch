---
name: search-datasets
description: Find sequencing datasets (NCBI SRA / EMBL-EBI ENA runs) by phenotype or disease, focused on aligned BAM/CRAM reads. Use when the user wants to locate sequencing data, aligned reads, BAM/CRAM files, FASTQ runs, or datasets for a disease, phenotype, or MONDO id.
when_to_use: Trigger on requests like "find aligned reads for ataxia", "what SRA/ENA datasets exist for <disease>", "find BAM/CRAM files for <phenotype>", "search sequencing runs for <MONDO id>", or "WGS/exome datasets for <disease>".
allowed-tools: Bash(aligned-read-search *) Bash(python -m aligned_read_search *) Bash(pip install *) Bash(pip show *) Bash(command -v *)
---

# Find sequencing datasets by phenotype

This skill drives the `aligned-read-search` CLI, which searches **NCBI SRA** and
**EMBL-EBI ENA** for sequencing runs matching a phenotype/disease and flags which
runs have an aligned **BAM/CRAM** file. It expands the phenotype to synonyms and
disease subtypes via the MONDO ontology, so a generic term like "ataxia" also
finds its subtypes.

## Step 1 — Ensure the tool is installed

Check whether the CLI is available:

```bash
command -v aligned-read-search
```

If that prints nothing, install the package. Prefer the copy bundled with this
plugin when available (`$CLAUDE_PLUGIN_ROOT`), otherwise install from GitHub:

```bash
pip install -e "${CLAUDE_PLUGIN_ROOT}" 2>/dev/null \
  || pip install "git+https://github.com/broadinstitute/SequencingReadDataSearch.git"
```

If a later command fails with `numpy.core.multiarray failed to import` or
`_ARRAY_API not found`, the environment has NumPy 2.x but pandas/pyarrow were
built against 1.x. Fix it and retry:

```bash
pip install 'numpy<2'
```

## Step 2 — Run the search

Always pass `--json` so the output is machine-readable:

```bash
aligned-read-search "<phenotype>" --json
```

Useful flags (short / long):

- `-a, --archive ena,sra` — which archives to query (default: both).
- `-l, --limit N` — max results per archive (default: 50).
- `-e, --exact` — literal match; skip MONDO ontology expansion.
- `-u, --include-unaligned` — also include runs **without** a BAM/CRAM file.
- `-s, --library-strategy WGS,WXS` — filter by strategy (e.g. WGS, WXS, RNA-Seq).
- `-v, --verbose` — print progress to stderr (handy when a query is slow).

Important defaults and tips:

- By **default only runs with an aligned BAM/CRAM file are returned**. Pass
  `-u/--include-unaligned` to see everything.
- Aligned files are sparse for many phenotypes. If the result is empty, retry
  with a larger limit (e.g. `-l 500`) before concluding there is nothing.
- A MONDO id works in place of a name, e.g. `aligned-read-search MONDO:0100254`.

## Step 3 — Report results

Parse the JSON and summarize for the user:

- How many runs matched, and how many have alignment files.
- For aligned runs, list `run_accession`, `library_strategy`, `organism`, the
  matched `phenotype_match` tokens, and the `alignment_urls` (BAM/CRAM download
  links).
- Surface `read_urls` (FASTQ) when the user wants the raw reads.

Example:

```bash
aligned-read-search ataxia -a ena -l 200 --json
```
