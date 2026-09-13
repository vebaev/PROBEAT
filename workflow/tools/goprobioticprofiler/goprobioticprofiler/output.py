"""Output writers for TSV and JSON files."""

import json
from pathlib import Path

import pandas as pd

from .utils import ensure_dir


def write_tsv(df: pd.DataFrame, path: Path | str, float_format: str = "%.4g") -> None:
    path = Path(path)
    ensure_dir(path.parent)
    df.to_csv(path, sep="\t", index=False, float_format=float_format)


def write_json(data: dict, path: Path | str) -> None:
    path = Path(path)
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, ensure_ascii=False)
