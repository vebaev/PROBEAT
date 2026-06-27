#!/usr/bin/env python3
import csv
import sys
from pathlib import Path


def score(row):
    category = row.get("refseq_category", "").lower()
    level = row.get("assembly_level", "").lower()
    type_material = row.get("relation_to_type_material", "").lower()
    version = row.get("version_status", "").lower()
    return (
        {"reference genome": 3, "representative genome": 2}.get(category, 1),
        {"complete genome": 4, "chromosome": 3, "scaffold": 2, "contig": 1}.get(level, 0),
        int(bool(type_material and type_material != "na")),
        int(version == "latest"),
        row.get("seq_rel_date", ""),
    )


def main():
    source, output = map(Path, sys.argv[1:3])
    header = None
    representatives = {}
    with source.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("#assembly_accession") or line.startswith("# assembly_accession"):
                header = line.lstrip("# ").rstrip("\n").split("\t")
                break
        if not header:
            raise SystemExit("Could not find NCBI assembly summary header")
        reader = csv.DictReader(handle, fieldnames=header, delimiter="\t")
        for row in reader:
            species_taxid = row.get("species_taxid", "").strip()
            if not species_taxid or species_taxid == "na":
                continue
            current = representatives.get(species_taxid)
            if current is None or score(row) > score(current):
                representatives[species_taxid] = row

    wanted = [
        "species_taxid", "organism_name", "assembly_accession", "refseq_category",
        "assembly_level", "relation_to_type_material", "genome_size", "gc_percent",
        "contig_count", "scaffold_count", "ftp_path", "seq_rel_date",
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=wanted, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in sorted(representatives.values(), key=lambda item: item.get("organism_name", "")):
            writer.writerow({key: row.get(key, "") for key in wanted})


if __name__ == "__main__":
    main()
