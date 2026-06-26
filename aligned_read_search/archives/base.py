"""Archive client interface."""

from __future__ import annotations

import re
from abc import ABC, abstractmethod
from typing import List

from ..models import Dataset


class ArchiveClient(ABC):
    """A source of datasets for a set of (already-expanded) search tokens."""

    name: str = "archive"

    @abstractmethod
    def search(self, tokens: List[str], limit: int = 50) -> List[Dataset]:
        """Return datasets whose metadata matches any of ``tokens``."""
        raise NotImplementedError


def all_token_matches(text: str, tokens: List[str]) -> str:
    """Return all tokens found in ``text`` as a comma-joined string, preserving
    token order and de-duplicating (case-insensitive).

    Matching is **whole-word** on alphanumeric boundaries, not substring: a token
    like ``ea1`` matches a standalone ``EA1`` but not ``area123``, and
    ``cerebellar`` does not match inside ``spinocerebellar``. This avoids
    spurious hits from short acronyms / semi-generic words substring-matching
    unrelated records.
    """
    text = text or ""
    seen: set[str] = set()
    out: List[str] = []
    for tok in tokens:
        key = tok.lower()
        if key in seen or not key:
            continue
        pattern = r"(?<![A-Za-z0-9])" + re.escape(key) + r"(?![A-Za-z0-9])"
        if re.search(pattern, text, re.IGNORECASE):
            seen.add(key)
            out.append(tok)
    return ", ".join(out)
