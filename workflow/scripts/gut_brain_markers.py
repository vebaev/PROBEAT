#!/usr/bin/env python3
"""Gut-brain axis marker detection from Bakta + eggNOG annotations."""
import csv, sys, re
from pathlib import Path
from collections import defaultdict


def product_matches_marker(product, marker):
    """Match descriptive markers as complete words/phrases, never substrings."""
    phrase = re.sub(r"[_-]+", " ", marker.lower()).strip()
    if not phrase or phrase.endswith(" like") or len(phrase.split()) < 2:
        return False
    pattern = r"(?<![a-z0-9])" + r"[\s_-]+".join(
        re.escape(word) for word in phrase.split()
    ) + r"(?![a-z0-9])"
    return re.search(pattern, product.lower()) is not None


def main():
    db_path = Path(sys.argv[1])
    bakta_path = Path(sys.argv[2])
    eggnog_path = Path(sys.argv[3])
    out_path = Path(sys.argv[4])
    fieldnames = [
        "metabolite_or_trait",
        "marker",
        "bakta_locus",
        "ko_or_ec_hints",
        "marker_logic",
        "confidence_if_present",
        "interpretation",
    ]

    # Read Bakta annotations: gene name -> set of locus_tags
    bakta_genes = defaultdict(set)
    bakta_products = defaultdict(set)
    with bakta_path.open(errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 8:
                continue
            locus = cols[5]
            gene = cols[6]
            product = cols[7]
            if gene:
                bakta_genes[gene.lower()].add(locus)
            bakta_products[product.lower()].add(locus)

    # Read eggNOG annotations: preferred gene name -> query IDs
    eggnog_genes = defaultdict(set)
    if eggnog_path.exists():
        with eggnog_path.open(errors="replace") as f:
            header = None
            for line in f:
                if line.startswith("#query\t"):
                    header = line.lstrip("#").rstrip("\n").split("\t")
                    continue
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if not header:
                    continue
                cols.extend([""] * (len(header) - len(cols)))
                row = dict(zip(header, cols))
                query = row.get("query", "")
                preferred_name = row.get("Preferred_name", "")
                if query and preferred_name and preferred_name != "-":
                    for gene in re.split(r"[/;,]", preferred_name):
                        gene = gene.strip().lower()
                        if gene:
                            eggnog_genes[gene].add(query)

    # Search gut-brain axis markers
    results = []
    with db_path.open(errors="replace") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            roles = re.search(r"Original roles:\s*([^|]+)", row.get("notes", ""), re.I)
            if roles and "primary" not in roles.group(1).lower():
                continue
            trait = row.get("metabolite_or_trait", "")
            primary = [g.strip() for g in row.get("primary_markers", "").split(";") if g.strip()]
            secondary = [g.strip() for g in row.get("secondary_markers", "").split(";") if g.strip()]
            hints = row.get("ko_or_ec_hints", "")
            logic = row.get("marker_logic", "")
            interpretation = row.get("interpretation", "")

            all_markers = primary + secondary
            matched_genes = set()
            found_loci = set()

            # Search gene names in Bakta
            for marker in all_markers:
                ml = marker.lower()
                # Direct gene name match
                loci = bakta_genes.get(ml, set())
                if not loci:
                    # Try product description match for genes > 3 chars
                    if len(marker) > 3:
                        for prod, ploci in bakta_products.items():
                            if product_matches_marker(prod, marker):
                                loci |= ploci
                if not loci:
                    loci |= eggnog_genes.get(ml, set())
                if loci:
                    matched_genes.add(marker)
                    found_loci |= loci

            if found_loci:
                markers_str = "; ".join(sorted(matched_genes)) if matched_genes else "-"
                loci_str = ", ".join(sorted(found_loci))
                results.append({
                    "metabolite_or_trait": trait,
                    "marker": markers_str,
                    "bakta_locus": loci_str,
                    "ko_or_ec_hints": hints,
                    "marker_logic": logic,
                    "confidence_if_present": row.get("confidence_if_present", ""),
                    "interpretation": interpretation,
                })

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        if results:
            w.writerows(results)
    if results:
        print("Found %d gut-brain axis markers -> %s" % (len(results), out_path))
    else:
        print("Found 0 gut-brain axis markers -> %s" % out_path)


if __name__ == "__main__":
    main()
