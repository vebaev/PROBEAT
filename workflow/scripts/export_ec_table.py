#!/usr/bin/env python3
"""Export detailed EC classification with per-subcategory counts from eggNOG annotations."""
import sys, csv
from pathlib import Path
from collections import defaultdict

EC_TABLE = [
    ("EC 1", "Dehydrogenases", ["1.1"], "Acting on CH-OH groups"),
    ("EC 1", "Aldehyde dehydrogenases", ["1.2"], "Acting on aldehydes"),
    ("EC 1", "Oxidases", ["1.3", "1.4"], "Oxidation reactions"),
    ("EC 1", "Peroxidases/Catalases", ["1.11"], "Oxidative stress defense"),
    ("EC 1", "Superoxide dismutases", ["1.15.1"], "Reactive oxygen detoxification"),
    ("EC 2", "Kinases", ["2.7"], "Phosphoryl transfer"),
    ("EC 2", "Glycosyltransferases", ["2.4"], "EPS and cell wall synthesis"),
    ("EC 2", "Acyltransferases", ["2.3"], "Lipid and SCFA metabolism"),
    ("EC 2", "Aminotransferases", ["2.6"], "Amino acid metabolism"),
    ("EC 2", "Methyltransferases", ["2.1"], "Methylation reactions"),
    ("EC 3", "Glycosidases", ["3.2"], "Carbohydrate degradation"),
    ("EC 3", "Peptidases", ["3.4"], "Protein degradation"),
    ("EC 3", "Amidases", ["3.5"], "Nitrogen metabolism"),
    ("EC 3", "Phosphatases", ["3.1.3", "3.6"], "Phosphate metabolism"),
    ("EC 3", "Esterases/Lipases", ["3.1.1"], "Lipid degradation"),
    ("EC 4", "Decarboxylases", ["4.1"], "GABA and polyamine pathways"),
    ("EC 4", "Aldolases", ["4.1.2"], "Central carbon metabolism"),
    ("EC 4", "Ammonia-lyases", ["4.3"], "Amino acid metabolism"),
    ("EC 5", "Isomerases", ["5.3"], "Sugar interconversion"),
    ("EC 5", "Epimerases/Racemases", ["5.1"], "Stereochemical conversions"),
    ("EC 5", "Mutases", ["5.4"], "Metabolic rearrangements"),
    ("EC 6", "C-O bond ligases", ["6.1"], "Aminoacyl-tRNA synthetases"),
    ("EC 6", "C-N bond ligases", ["6.3"], "Peptide/cofactor biosynthesis"),
    ("EC 6", "DNA ligases", ["6.5"], "Genome maintenance"),
    ("EC 7", "Hydrogen ion translocases", ["7.1"], "Proton motive force"),
    ("EC 7", "Inorganic cation translocases", ["7.2"], "Ion homeostasis"),
    ("EC 7", "ABC transport translocases", ["7.6"], "ATP-dependent transport"),
]


def ec_matches(ec, prefixes):
    for prefix in prefixes:
        if ec == prefix:
            return True
        if ec.startswith(prefix + "."):
            return True
        # Handle prefix like "1.1" matching "1.1.2.3"
        parts = ec.split(".")
        prefix_parts = prefix.split(".")
        if len(parts) >= len(prefix_parts):
            if all(parts[i] == prefix_parts[i] for i in range(len(prefix_parts))):
                return True
    return False


def main():
    eggnog = Path(sys.argv[1])
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else eggnog.with_name(eggnog.stem + "_ec.tsv")

    ec_by_subcat = defaultdict(set)
    total_ec = set()

    with open(eggnog) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) >= 11 and cols[10] not in ("", "-"):
                for ec in cols[10].split(","):
                    ec = ec.strip()
                    if not ec:
                        continue
                    total_ec.add(ec)
                    for cat, subcat, prefixes, desc in EC_TABLE:
                        if ec_matches(ec, prefixes):
                            ec_by_subcat[(cat, subcat)].add(ec)

    with open(out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["EC class", "Subcategory", "EC prefix", "Description", "Unique ECs"])
        seen_cats = set()
        for cat, subcat, prefixes, desc in EC_TABLE:
            count = len(ec_by_subcat.get((cat, subcat), set()))
            prefix_str = ",".join(p + ".*.*" if len(p.split(".")) == 1 else p + (".*" if len(p.split(".")) == 2 else "") for p in prefixes)
            # Format prefix nicely
            formatted_prefix = []
            for p in prefixes:
                parts = p.split(".")
                if len(parts) == 1:
                    formatted_prefix.append(p + ".*.*")
                elif len(parts) == 2:
                    formatted_prefix.append(p + ".*")
                else:
                    formatted_prefix.append(p + ".*")
            prefix_display = ", ".join(formatted_prefix)
            if cat not in seen_cats:
                w.writerow([cat, subcat, prefix_display, desc, str(count)])
                seen_cats.add(cat)
            else:
                w.writerow(["", subcat, prefix_display, desc, str(count)])

    print("Saved: %s (%d unique ECs)" % (out, len(total_ec)))

if __name__ == "__main__":
    main()
