"""Utility functions for GOProbioticProfiler."""

import re
from pathlib import Path
from typing import Optional, Union


def ensure_dir(path: Union[str, Path]) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def safe_id(text: str) -> str:
    """Convert any string into a safe snake_case identifier."""
    text = re.sub(r'[^\w\s-]', '', text).strip().lower()
    return re.sub(r'[-\s]+', '_', text)


def fdr_correction(pvalues: list[float]) -> list[float]:
    """Benjamini-Hochberg FDR correction."""
    from statsmodels.stats.multitest import multipletests
    if not pvalues:
        return []
    _, fdr, _, _ = multipletests(pvalues, alpha=0.05, method='fdr_bh')
    return list(fdr)
