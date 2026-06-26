"""Archive client interface."""

from __future__ import annotations

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
    """Return all tokens found (case-insensitively) in ``text`` as a
    comma-joined string, preserving token order and de-duplicating.
    """
    low = (text or "").lower()
    seen: set[str] = set()
    out: List[str] = []
    for tok in tokens:
        key = tok.lower()
        if key in low and key not in seen:
            seen.add(key)
            out.append(tok)
    return ", ".join(out)
