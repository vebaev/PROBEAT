#!/usr/bin/env python3
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


def aliases(symbol):
    return {
        value.strip().lower()
        for value in re.split(r"[/;,]", symbol)
        if value.strip()
    }


def categories(text):
    normalized = []
    replacements = {
        "acid stress": "Acid stress",
        "bile resistance": "Bile resistance",
        "temperature": "Temperature stress",
        "temperature stress": "Temperature stress",
        "oxidative stress": "Oxidative stress",
        "antioxidant": "Antioxidant",
    }
    for value in re.split(r"[/;]", text):
        value = value.strip().replace("Asid stress", "Acid stress")
        if value:
            normalized.append(replacements.get(value.lower(), value))
    return normalized


def load_bakta(path):
    rows = []
    by_locus = {}
    with path.open(newline="", errors="replace") as handle:
        for row in csv.reader((line for line in handle if not line.startswith("#")), delimiter="\t"):
            if len(row) < 9 or row[1].lower() != "cds":
                continue
            record = {
                "contig": row[0],
                "start": row[2],
                "stop": row[3],
                "strand": row[4],
                "locus_tag": row[5],
                "gene": row[6].strip(),
                "product": row[7],
            }
            rows.append(record)
            if record["locus_tag"]:
                by_locus[record["locus_tag"]] = record
    return rows, by_locus


def load_eggnog(path):
    if not path.exists():
        return []
    header = None
    rows = []
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("#query\t"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if not header or not line.strip() or line.startswith("#"):
                continue
            values = line.rstrip("\n").split("\t")
            values.extend([""] * (len(header) - len(values)))
            rows.append(dict(zip(header, values)))
    return rows


def main():
    database, bakta_tsv, eggnog_tsv, hits_out, summary_out = map(Path, sys.argv[1:6])
    markers = defaultdict(list)
    entries = []
    with database.open(newline="", errors="replace") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 3:
                continue
            symbol, description, function = row[:3]
            entry = {
                "marker": symbol,
                "description": description,
                "function": function.replace("Asid stress", "Acid stress"),
            }
            entry["id"] = (entry["marker"], entry["description"], entry["function"])
            entries.append(entry)
            for alias in aliases(symbol):
                markers[alias].append(entry)

    hits = []
    seen_hits = set()
    found_in_bakta = set()
    bakta_rows, bakta_by_locus = load_bakta(bakta_tsv)

    def add_hit(record, gene, marker, source, product=None):
        key = (
            marker["marker"].strip().lower(),
            marker["function"].strip().lower(),
            record.get("locus_tag", ""),
        )
        if key in seen_hits:
            return
        seen_hits.add(key)
        hits.append({
            "contig": record.get("contig", ""),
            "start": record.get("start", ""),
            "stop": record.get("stop", ""),
            "strand": record.get("strand", ""),
            "locus_tag": record.get("locus_tag", ""),
            "gene": gene,
            "product": product if product is not None else record.get("product", ""),
            "marker": marker["marker"],
            "description": marker["description"],
            "function": marker["function"],
            "source": source,
        })

    for record in bakta_rows:
        gene = record["gene"]
        if not gene:
            continue
        alias = gene.lower()
        for marker in markers.get(alias, []):
            add_hit(record, gene, marker, "Bakta")
            found_in_bakta.add((marker["id"], alias))

    remaining = {
        (entry["id"], alias)
        for entry in entries
        for alias in aliases(entry["marker"])
        if (entry["id"], alias) not in found_in_bakta
    }
    for row in load_eggnog(eggnog_tsv):
        preferred_name = row.get("Preferred_name", "").strip()
        if not preferred_name or preferred_name == "-":
            continue
        locus = row.get("query", "").strip()
        record = bakta_by_locus.get(locus, {"locus_tag": locus})
        for alias in aliases(preferred_name):
            for marker in markers.get(alias, []):
                if (marker["id"], alias) not in remaining:
                    continue
                add_hit(
                    record,
                    preferred_name,
                    marker,
                    "eggNOG fallback",
                    row.get("Description", "") or record.get("product", ""),
                )

    hits_out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["contig", "start", "stop", "strand", "locus_tag", "gene", "product", "marker", "description", "function", "source"]
    with hits_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(hits)

    unique_markers = {hit["marker"] for hit in hits}
    category_loci = defaultdict(set)
    category_markers = defaultdict(set)
    for hit in hits:
        for category in categories(hit["function"]):
            category_loci[category].add(hit["locus_tag"])
            category_markers[category].add(hit["marker"])

    with summary_out.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value", "unique_markers"])
        writer.writerow(["total_hits", len(hits), len(unique_markers)])
        for category in sorted(category_loci, key=lambda key: (-len(category_loci[key]), key)):
            writer.writerow([category, len(category_loci[category]), len(category_markers[category])])


if __name__ == "__main__":
    main()
