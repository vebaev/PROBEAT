#!/usr/bin/env python3
"""Compare existing BV7/DPL5 core metrics with the publication baseline."""

from __future__ import annotations

import re
from pathlib import Path
import sys

import yaml


ROOT = Path(__file__).resolve().parents[2]


def table(path: Path) -> dict[str, str]:
    rows = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                rows[parts[0]] = parts[1]
    return rows


def metrics(sample: str) -> dict[str, float]:
    result = ROOT / "data" / "results" / sample
    quast = table(result / "quast" / "report.tsv")
    bakta = (result / "bakta" / f"{sample}.txt").read_text(encoding="utf-8")
    busco = (result / "busco" / "busco_summary.txt").read_text(encoding="utf-8")
    pgmdb = table(result / "pgmdb" / "summary.tsv")
    busco_match = re.search(r"C:([\d.]+)%\[S:[\d.]+%,D:([\d.]+)%\]", busco)
    cds_match = re.search(r"CDSs:\s+(\d+)", bakta)
    if not busco_match or not cds_match:
        raise ValueError(f"Could not parse baseline metrics for {sample}")
    return {
        "contigs_quast": int(quast["# contigs"]),
        "genome_size_bp_quast": int(quast["Total length"]),
        "gc_percent": float(quast["GC (%)"]),
        "n50_bp": int(quast["N50"]),
        "busco_complete_percent": float(busco_match.group(1)),
        "busco_duplicated_percent": float(busco_match.group(2)),
        "bakta_cds": int(cds_match.group(1)),
        "pgmdb_total_hits": int(pgmdb["total_hits"]),
        "pgmdb_unique_markers": int(
            next(line.split("\t")[2] for line in (result / "pgmdb" / "summary.tsv").read_text().splitlines() if line.startswith("total_hits\t"))
        ),
    }


def main() -> int:
    baseline = yaml.safe_load((ROOT / "manifests" / "release-baseline.yaml").read_text())["baseline"]
    failed = 0
    for sample, expected in baseline["samples"].items():
        observed = metrics(sample)
        for key, expected_value in expected.items():
            if key == "layout":
                continue
            tolerance = 0.01 if key == "gc_percent" else 0.1 if key.startswith("busco_") else 0
            ok = abs(observed[key] - expected_value) <= tolerance
            print(f"{'[OK]' if ok else '[FAIL]'} {sample} {key}: {observed[key]} expected {expected_value}")
            failed += int(not ok)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
