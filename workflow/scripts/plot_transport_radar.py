#!/usr/bin/env python3
"""Export transport TC table from eggNOG KEGG_TC annotations as TSV."""
import sys, csv
from pathlib import Path
from collections import defaultdict

CATEGORIES = [
    ("Carbohydrate uptake", ["4.A.1", "3.A.1", "2.A.1"]),
    ("Oligosaccharide uptake", ["3.A.1.5"]),
    ("Peptide uptake", ["3.A.1.17"]),
    ("Amino acid & polyamine uptake", ["2.A.3"]),
    ("Stress resistance", ["2.A.12", "2.A.66", "2.A.6"]),
    ("Ion homeostasis", ["3.A.3"]),
]

CATEGORY_TABLE = {
    "Carbohydrate uptake": [
        ("4.A.1", "PTS systems", "Glucose, fructose, mannose, cellobiose uptake; prebiotic utilization"),
        ("3.A.1", "ABC sugar transporters", "Oligosaccharide and sugar uptake"),
        ("2.A.1", "MFS transporters", "Lactose, galactose and carbohydrate transport"),
    ],
    "Oligosaccharide uptake": [
        ("3.A.1.5", "Oligosaccharide ABC transporters", "FOS, GOS, XOS, raffinose-family oligosaccharide utilization"),
    ],
    "Peptide uptake": [
        ("3.A.1.17", "Opp/Dpp peptide transporters", "Peptide nutrition, gut adaptation and colonization"),
    ],
    "Amino acid & polyamine uptake": [
        ("2.A.3", "APC family transporters", "Amino acid uptake, polyamine metabolism"),
    ],
    "Stress resistance": [
        ("2.A.12", "Glycine betaine/choline transporters", "Osmotic stress tolerance"),
        ("2.A.66", "MATE transporters", "Bile tolerance and multidrug efflux"),
        ("2.A.6", "RND efflux systems", "Toxin and stress resistance"),
    ],
    "Ion homeostasis": [
        ("3.A.3", "ATP-driven ion transporters", "pH regulation and ion balance"),
    ],
}


def main():
    eggnog = Path(sys.argv[1])
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else eggnog.with_name(eggnog.stem + "_transport.tsv")

    seen = defaultdict(set)
    with open(eggnog) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) >= 18 and cols[17] not in ("", "-"):
                for tc in cols[17].split(","):
                    tc = tc.strip()
                    if not tc:
                        continue
                    for cat_name, prefixes in CATEGORIES:
                        for prefix in prefixes:
                            if tc == prefix or tc.startswith(prefix + "."):
                                seen[cat_name].add(tc)

    with open(out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Category", "TC class", "Transporter family", "Probiotic relevance", "Unique TCs"])
        for cat_name, rows in CATEGORY_TABLE.items():
            detected_set = seen.get(cat_name, set())
            for i, (tc, tc_name, desc) in enumerate(rows):
                count = sum(1 for t in detected_set if t == tc or t.startswith(tc + "."))
                w.writerow([cat_name if i == 0 else "", tc, tc_name, desc, str(count)])

    print("Saved: %s" % out)

if __name__ == "__main__":
    main()
