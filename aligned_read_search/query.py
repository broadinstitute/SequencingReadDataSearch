"""Phenotype term -> archive search terms.

Two expanders implement the same ``QueryExpander`` interface so call sites are
identical:

* ``IdentityExpander``  - literal exact match (opt-out via ``--exact``).
* ``OntologyExpander``  - default; resolves the phenotype against MONDO (via the
  EBI OLS4 API), pulls synonyms + descendant subtypes, and reduces them to
  distinctive single-word tokens. Single tokens are required because ENA
  wildcards do not match multi-word phrases (see plan grounding facts).
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import tempfile
from pathlib import Path
from typing import List, Protocol, runtime_checkable

import requests

logger = logging.getLogger(__name__)

OLS_BASE = "https://www.ebi.ac.uk/ols4/api"
MONDO_CURIE_RE = re.compile(r"^MONDO[:_]\d+$", re.IGNORECASE)

# Generic words that carry no disease signal once a phrase is tokenized.
_STOPWORDS = {
    "disease",
    "disorder",
    "syndrome",
    "type",
    "form",
    "autosomal",
    "recessive",
    "dominant",
    "hereditary",
    "familial",
    "congenital",
    "idiopathic",
    "acquired",
    "chronic",
    "acute",
    "early",
    "late",
    "onset",
    "with",
    "without",
    "and",
    "the",
    "due",
    "to",
    "of",
    "non",
    "like",
    "associated",
    "related",
    "linked",
    "deficiency",
    "complex",
    "primary",
    "secondary",
}

# Generic anatomy / body-region / finding / age-descriptor words that show up in
# disease *subtype* labels (e.g. "infantile liver failure", "intellectual
# disability") but carry no disease-specific signal. These are filtered out of
# tokens derived from descendant subtypes ONLY -- never from the user's own term
# or the primary disease label -- so a direct search for e.g. "anemia" still
# keeps that word. Heuristic and tunable.
_GENERIC_SUBTYPE_STOPWORDS = {
    # anatomy / body regions
    "liver", "kidney", "renal", "hepatic", "cardiac", "cardiovascular", "brain",
    "cerebral", "spinal", "muscular", "bone", "skeletal", "optic", "corneal",
    "retinal", "choroidal", "macular", "ocular", "eye", "blood", "nerve",
    "vision", "hearing", "pulmonary", "respiratory", "gastrointestinal",
    "peripheral", "central", "proximal", "distal", "motor", "sensory",
    # generic findings / symptoms
    "failure", "anemia", "disability", "weakness", "atrophy", "dystrophy",
    "hypotonia", "hypogonadism", "seizure", "seizures", "intellectual",
    "developmental", "cognitive", "growth", "retardation", "intrusion",
    "saccadic", "spastic", "spasticity", "aniridia", "progressive", "generalized",
    # age / onset descriptors
    "infantile", "juvenile", "neonatal", "adult", "childhood", "prenatal",
    "perinatal",
}

# Short tokens worth keeping despite being < 3 chars or numeric-ish.
_KEEP_SHORT = {"sca", "ea1", "ea2", "frda", "ftd", "als", "ad", "pd"}

# Acronym-like tokens in a synonym/label, e.g. SCA, FRDA, EA1, AOA2, SCAR12.
_ACRONYM_RE = re.compile(r"\b[A-Z][A-Z0-9]{1,6}\b")


def _acronyms(text: str) -> List[str]:
    return [m.group(0).lower() for m in _ACRONYM_RE.finditer(text or "")]


def _keep_acronym(acronym: str) -> bool:
    """Keep only acronym tokens specific enough to avoid false positives.

    Short pure-alpha codes (``scan``, ``capa``, ``asat``) collide with common
    words even under whole-word matching, so we keep an acronym only if it
    carries a digit (``ea1``, ``scar12``, ``kcna1``), is reasonably long
    (>= 5 chars, e.g. ``scasi``), or is in the curated short allowlist.
    """
    return (
        any(c.isdigit() for c in acronym)
        or len(acronym) >= 5
        or acronym in _KEEP_SHORT
    )


@runtime_checkable
class QueryExpander(Protocol):
    """Turns a phenotype term into a list of single-token search terms."""

    def expand(self, term: str) -> List[str]: ...


class IdentityExpander:
    """Literal, exact-match expander. Returns the term unchanged."""

    def expand(self, term: str) -> List[str]:
        term = term.strip()
        logger.info("Exact match: searching for the literal term %r", term)
        return [term] if term else []


def tokenize(phrase: str, extra_stop: set | None = None) -> List[str]:
    """Split a label/synonym into distinctive single-word search tokens.

    ``extra_stop`` adds further stopwords for this call (used to strip generic
    anatomy/finding words from descendant *subtype* labels without affecting the
    primary term).
    """
    stop = _STOPWORDS if not extra_stop else _STOPWORDS | extra_stop
    tokens = []
    for raw in re.split(r"[^A-Za-z0-9]+", phrase.lower()):
        if not raw:
            continue
        if raw in _KEEP_SHORT:
            tokens.append(raw)
            continue
        if raw in stop:
            continue
        if raw.isdigit():
            continue
        if len(raw) < 3:
            continue
        tokens.append(raw)
    return tokens


class _DiskCache:
    """Tiny JSON disk cache keyed by URL (OLS is the only network dependency)."""

    def __init__(self, namespace: str = "aligned_read_search"):
        base = os.environ.get("XDG_CACHE_HOME") or os.path.join(
            tempfile.gettempdir(), ".cache"
        )
        self.dir = Path(base) / namespace
        self.dir.mkdir(parents=True, exist_ok=True)

    def _path(self, key: str) -> Path:
        digest = hashlib.sha256(key.encode()).hexdigest()[:16]
        return self.dir / f"{digest}.json"

    def get(self, key: str):
        path = self._path(key)
        if path.exists():
            try:
                return json.loads(path.read_text())
            except (ValueError, OSError):
                return None
        return None

    def set(self, key: str, value) -> None:
        try:
            self._path(key).write_text(json.dumps(value))
        except OSError:
            pass


class OntologyExpander:
    """Default expander: MONDO-backed synonym + subtype expansion via EBI OLS4."""

    def __init__(
        self,
        ontology: str = "mondo",
        include_descendants: bool = True,
        max_terms: int = 50,
        search_rows: int = 20,
        timeout: float = 30.0,
        session: requests.Session | None = None,
        cache: _DiskCache | None = None,
    ):
        self.ontology = ontology
        self.include_descendants = include_descendants
        self.max_terms = max_terms
        self.search_rows = search_rows
        self.timeout = timeout
        self.session = session or requests.Session()
        self.cache = cache if cache is not None else _DiskCache()

    # -- HTTP helper with caching ------------------------------------------
    def _get_json(self, url: str, params: dict | None = None):
        key = url + "?" + json.dumps(params or {}, sort_keys=True)
        cached = self.cache.get(key)
        if cached is not None:
            return cached
        resp = self.session.get(url, params=params, timeout=self.timeout)
        resp.raise_for_status()
        data = resp.json()
        self.cache.set(key, data)
        return data

    # -- resolve a term/CURIE to OLS term docs -----------------------------
    def _search_docs(self, term: str) -> List[dict]:
        """Return the top matching OLS term docs (single doc for a CURIE)."""
        field_list = "iri,label,obo_id,synonym"
        if MONDO_CURIE_RE.match(term):
            obo_id = term.upper().replace("_", ":")
            data = self._get_json(
                f"{OLS_BASE}/search",
                {
                    "q": obo_id,
                    "ontology": self.ontology,
                    "queryFields": "obo_id",
                    "fieldList": field_list,
                },
            )
            return ((data.get("response") or {}).get("docs") or [])[:1]

        data = self._get_json(
            f"{OLS_BASE}/search",
            {
                "q": term,
                "ontology": self.ontology,
                "rows": self.search_rows,
                "fieldList": field_list,
            },
        )
        return (data.get("response") or {}).get("docs") or []

    @staticmethod
    def _label_of(doc: dict) -> str:
        return doc.get("label") or ""

    @staticmethod
    def _synonyms_of(doc: dict) -> List[str]:
        out: List[str] = []
        for key in ("synonym", "obo_synonym", "synonyms"):
            val = doc.get(key)
            if isinstance(val, list):
                out.extend(s for s in val if isinstance(s, str))
        return out

    def _collect(self, doc: dict, labels: List[str], acronyms: List[str]) -> None:
        """Add a doc's label (full tokens) and its acronym-like synonym tokens.

        Labels are clean disease names; synonyms are noisy (eponyms, "caused by
        mutation in ...") so from them we keep only acronyms like SCA/FRDA/EA1,
        which are exactly the abbreviations datasets use.
        """
        label = self._label_of(doc)
        if label:
            labels.append(label)
            acronyms.extend(_acronyms(label))
        for syn in self._synonyms_of(doc):
            acronyms.extend(_acronyms(syn))

    def _descendant_docs(self, iri: str) -> List[dict]:
        # OLS requires the IRI to be double URL-encoded in the path.
        from urllib.parse import quote

        enc = quote(quote(iri, safe=""), safe="")
        url = f"{OLS_BASE}/ontologies/{self.ontology}/terms/{enc}/hierarchicalDescendants"
        try:
            data = self._get_json(url, {"size": 200})
        except requests.RequestException:
            return []
        return (data.get("_embedded") or {}).get("terms") or []

    def expand(self, term: str) -> List[str]:
        term = term.strip()
        if not term:
            return []
        logger.info("Resolving %r against the %s ontology (EBI OLS)", term, self.ontology)
        docs = self._search_docs(term)
        if not docs:
            # No ontology hit -> fall back to the literal term.
            logger.info("No %s match for %r; falling back to the literal term", self.ontology, term)
            return [term]
        logger.info("Matched %d %s term(s) for %r", len(docs), self.ontology, term)

        # The user's term and the single primary disease label define the
        # disease's identity and are tokenized loosely. Every other ontology
        # label -- sibling top-matches and descendant subtypes, whose multi-word
        # names carry generic anatomy/finding words like "liver failure" or
        # "intellectual disability" -- is tokenized with the generic stoplist so
        # that noise doesn't become a search token.
        primary_labels: List[str] = []
        other_labels: List[str] = []
        acronyms: List[str] = []
        primary: dict | None = None

        if MONDO_CURIE_RE.match(term):
            primary = docs[0]
            self._collect(primary, primary_labels, acronyms)
        else:
            # Aggregate across top docs that are actually about this disease
            # (label shares a token with the query, e.g. "ataxia"), so a generic
            # term like "Ataxia" picks up its many subtypes.
            orig = set(tokenize(term))
            for doc in docs:
                label = (doc.get("label") or "").lower()
                if orig and any(tok in label for tok in orig):
                    if primary is None:
                        primary = doc
                        self._collect(doc, primary_labels, acronyms)
                    else:
                        self._collect(doc, other_labels, acronyms)
            if primary is None:  # no label overlap; fall back to top hit
                primary = docs[0]
                self._collect(primary, primary_labels, acronyms)

        if self.include_descendants and primary and primary.get("iri"):
            logger.info("Fetching descendant subtypes of %s", primary.get("label") or primary["iri"])
            descendants = self._descendant_docs(primary["iri"])
            logger.info("Found %d descendant subtype(s)", len(descendants))
            for doc in descendants:
                self._collect(doc, other_labels, acronyms)

        # Order: primary label + user term (loose), then sibling/subtype labels
        # (generic-filtered), then acronyms. De-dup preserving order.
        seen: set[str] = set()
        tokens: List[str] = []

        def _add(source: str, extra_stop: set | None = None) -> None:
            for tok in tokenize(source, extra_stop=extra_stop):
                if tok not in seen:
                    seen.add(tok)
                    tokens.append(tok)

        for source in (*primary_labels, term):
            _add(source)
        for source in other_labels:
            _add(source, extra_stop=_GENERIC_SUBTYPE_STOPWORDS)
        for tok in acronyms:
            if tok not in seen and _keep_acronym(tok):
                seen.add(tok)
                tokens.append(tok)
        final = tokens[: self.max_terms]
        logger.info("Expanded %r to %d search token(s): %s", term, len(final), ", ".join(final))
        return final
