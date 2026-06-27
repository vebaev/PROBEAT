#!/usr/bin/env python3
"""Normalize the gut-brain marker database to one row per unique gene."""
import csv
import sys
from collections import OrderedDict
from pathlib import Path


FIELDS = [
    "marker_id", "category", "subcategory", "metabolite_or_trait",
    "primary_markers", "secondary_markers", "ko_or_ec_hints",
    "marker_logic", "confidence_if_present", "interpretation",
    "risk_or_benefit", "report_label", "notes",
]


def values(text, separator=";"):
    return [value.strip() for value in (text or "").split(separator) if value.strip()]


def add_unique(target, value):
    if value and value not in target:
        target.append(value)


def main():
    source = Path(sys.argv[1])
    output = Path(sys.argv[2])
    genes = OrderedDict()

    with source.open(newline="", errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            marker_roles = [
                ("primary", marker) for marker in values(row.get("primary_markers"))
            ] + [
                ("secondary", marker) for marker in values(row.get("secondary_markers"))
            ]
            for role, marker in marker_roles:
                key = marker.lower()
                record = genes.setdefault(key, {
                    "gene": marker,
                    "roles": [],
                    "marker_ids": [],
                    **{field: [] for field in FIELDS if field not in ("marker_id", "primary_markers", "secondary_markers")},
                })
                add_unique(record["roles"], role)
                add_unique(record["marker_ids"], row.get("marker_id", ""))
                for field in (
                    "category", "subcategory", "metabolite_or_trait",
                    "marker_logic", "confidence_if_present", "interpretation",
                    "risk_or_benefit", "report_label", "notes",
                ):
                    add_unique(record[field], row.get(field, ""))
                for hint in values(row.get("ko_or_ec_hints")):
                    add_unique(record["ko_or_ec_hints"], hint)

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        for index, record in enumerate(genes.values(), start=1):
            notes = list(record["notes"])
            add_unique(notes, f"Original roles: {', '.join(record['roles'])}")
            add_unique(notes, f"Original marker IDs: {', '.join(record['marker_ids'])}")
            writer.writerow({
                "marker_id": f"GBA_GENE_{index:03d}",
                "category": " | ".join(record["category"]),
                "subcategory": " | ".join(record["subcategory"]),
                "metabolite_or_trait": " | ".join(record["metabolite_or_trait"]),
                "primary_markers": record["gene"],
                "secondary_markers": "",
                "ko_or_ec_hints": "; ".join(record["ko_or_ec_hints"]),
                "marker_logic": " | ".join(record["marker_logic"]),
                "confidence_if_present": " | ".join(record["confidence_if_present"]),
                "interpretation": " | ".join(record["interpretation"]),
                "risk_or_benefit": " | ".join(record["risk_or_benefit"]),
                "report_label": " | ".join(record["report_label"]),
                "notes": " | ".join(notes),
            })

    print(f"Normalized {source} to {len(genes)} unique gene rows in {output}")


if __name__ == "__main__":
    main()
