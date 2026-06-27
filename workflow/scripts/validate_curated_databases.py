#!/usr/bin/env python3
"""Validate the small curated interpretation tables distributed by PROBEAT."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys


def validate_pgmdb(path: Path) -> list[str]:
    errors = []
    seen = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(row) != 3:
                errors.append(f"{path}:{line_number}: expected 3 columns")
                continue
            key = tuple(value.strip() for value in row)
            if not all(key):
                errors.append(f"{path}:{line_number}: empty required value")
            if key in seen:
                errors.append(f"{path}:{line_number}: duplicate row")
            seen.add(key)
    return errors


def validate_header_table(path: Path, required: set[str], unique: str | None = None) -> list[str]:
    errors = []
    seen = set()
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            return [f"{path}: missing columns: {', '.join(sorted(missing))}"]
        for line_number, row in enumerate(reader, 2):
            if unique:
                value = (row.get(unique) or "").strip()
                if not value:
                    errors.append(f"{path}:{line_number}: empty {unique}")
                elif value in seen:
                    errors.append(f"{path}:{line_number}: duplicate {unique}: {value}")
                seen.add(value)
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database-root", type=Path, default=Path("data/databases"))
    args = parser.parse_args()
    root = args.database_root
    errors = []
    errors.extend(validate_pgmdb(root / "PGMDB.tab"))
    errors.extend(
        validate_header_table(
            root / "gut_brain_axis_marker_database.tsv",
            {"marker_id", "metabolite_or_trait", "primary_markers", "marker_logic", "interpretation"},
            "marker_id",
        )
    )
    errors.extend(
        validate_header_table(
            root / "probiotic_cazy.txt",
            {"Prebiotic substrate / group", "Main dbCAN / CAZy families", "Typical enzymes / functions"},
        )
    )
    for error in errors:
        print(f"[FAIL] {error}")
    if errors:
        return 1
    print("[OK] Curated database structure and duplicate checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
