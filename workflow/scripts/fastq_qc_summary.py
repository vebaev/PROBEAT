#!/usr/bin/env python3
"""Write a fastp-compatible summary JSON for raw and cleaned ONT FASTQ reads."""

import gzip
import json
import sys
from pathlib import Path


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def summarize(path):
    reads = bases = q30_bases = gc_bases = 0
    with open_text(path) as handle:
        while True:
            name = handle.readline()
            if not name:
                break
            sequence = handle.readline().strip()
            handle.readline()
            quality = handle.readline().strip()
            reads += 1
            bases += len(sequence)
            gc_bases += sum(base in "GgCc" for base in sequence)
            q30_bases += sum((ord(char) - 33) >= 30 for char in quality)
    return {
        "total_reads": reads,
        "total_bases": bases,
        "q30_rate": q30_bases / bases if bases else 0,
        "gc_content": gc_bases / bases if bases else 0,
    }


def main():
    raw_path, cleaned_path, output_path = map(Path, sys.argv[1:4])
    before = summarize(raw_path)
    after = summarize(cleaned_path)
    payload = {
        "summary": {
            "before_filtering": before,
            "after_filtering": after,
        },
        "filtering_result": {
            "passed_filter_reads": after["total_reads"],
            "low_quality_reads": 0,
            "too_many_N_reads": 0,
        },
        "adapter_cutting": {
            "adapter_trimmed_reads": max(0, before["total_reads"] - after["total_reads"]),
            "adapter_trimmed_bases": max(0, before["total_bases"] - after["total_bases"]),
        },
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
