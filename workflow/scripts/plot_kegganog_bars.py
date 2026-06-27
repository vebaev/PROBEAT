#!/usr/bin/env python3
"""Bar chart from KEGGaNOG filtered pathways, sorted by category then value."""
import csv, sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CATEGORIES = {
    "Carbon fixation": ["3-Hydroxypropionate Bicycle", "4-Hydroxybutyrate/3-hydroxypropionate", "CBB Cycle", "gluconeogenesis", "rTCA Cycle", "RuBisCo", "Wood-Ljungdahl"],
    "Carbohydrate metabolism": ["Entner-Doudoroff Pathway", "glycolysis", "Sulfolipid biosynthesis", "TCA Cycle", "Glyoxylate shunt", "Mixed acid: Acetate", "Mixed acid: Ethanol, Acetate to Acetylaldehyde", "Mixed acid: Ethanol, Acetyl-CoA to Acetylaldehyde (reversible)", "Mixed acid: Ethanol, Acetylaldehyde to Ethanol", "Mixed acid: Formate", "Mixed acid: Formate to CO2 & H2", "Mixed acid: Lactate", "Mixed acid: PEP to Succinate via OAA, malate & fumarate", "alpha-amylase", "polyhydroxybutyrate synthesis", "starch/glycogen degradation", "starch/glycogen synthesis", "Anaplerotic genes"],
    "Carbon degradation": ["beta-glucosidase", "beta-N-acetylhexosaminidase", "chitinase", "D-galacturonate epimerase", "D-galacturonate isomerase", "diacetylchitobiose deacetylase", "glucoamylase", "pullulanase", "DMS dehydrogenase", "DMSP demethylation", "Naphthalene degradation to salicylate", "alcohol oxidase", "basic endochitinase B", "bifunctional chitinase/lysozyme", "cellulase", "dimethylamine/trimethylamine dehydrogenase", "oligogalacturonide lyase", "pectinesterase", "soluble methane monooxygenase", "exopolygalacturonase", "PET degradation"],
    "Nitrogen metabolism": ["dissim nitrate reduction", "DNRA", "nitric oxide reduction", "nitrite oxidation", "nitrite reduction", "nitrogen fixation", "nitrous-oxide reduction", "ammonia oxidation (amo/pmmo)", "hydrazine dehydrogenase", "hydrazine synthase", "hydroxylamine oxidation"],
    "Sulfur metabolism": ["alt thiosulfate oxidation tsdA", "dissimilatory sulfate < > APS", "dissimilatory sulfite < > APS", "dissimilatory sulfite < > sulfide", "DMSO reductase", "sulfide oxidation", "sulfite dehydrogenase", "sulfite dehydrogenase (quinone)", "sulfur dioxygenase", "thiosulfate oxidation", "thiosulfate/polysulfide reductase", "alt thiosulfate oxidation doxAD", "sulfhydrogenase", "sulfur assimilation", "sulfur disproportionation", "sulfur reductase sreABC"],
    "Oxidative phosphorylation": ["Cytochrome bd complex", "Cytochrome c oxidase", "Cytochrome c oxidase, cbb3-type", "F-type ATPase", "Na-NADH-ubiquinone oxidoreductase", "NADH-quinone oxidoreductase", "Ubiquinol-cytochrome c reductase", "V-type ATPase", "Cytochrome aa3-600 menaquinol oxidase", "Cytochrome b6/f complex", "Cytochrome o ubiquinol oxidase", "NAD(P)H-quinone oxidoreductase"],
    "Hydrogen redox": ["NAD-reducing hydrogenase", "NiFe hydrogenase Hyd-1", "Coenzyme B/Coenzyme M regeneration", "Coenzyme M reduction to methane", "NADP-reducing hydrogenase", "NiFe hydrogenase", "ferredoxin hydrogenase", "hydrogen:quinone oxidoreductase", "membrane-bound hydrogenase"],
    "Amino acid metabolism": ["arginine", "asparagine", "glutamine", "histidine", "lysine", "serine", "threonine", "Serine pathway/formaldehyde assimilation", "alanine", "aspartate", "cysteine", "glutamate", "glycine", "isoleucine", "leucine", "methionine", "phenylalanine", "proline", "tryptophan", "tyrosine", "valine"],
    "Vitamin biosynthesis": ["cobalamin biosynthesis", "riboflavin biosynthesis", "thiamin biosynthesis", "MEP-DOXP pathway", "Retinal biosynthesis", "Retinal from apo-carotenals", "carotenoids backbone biosynthesis", "end-product astaxanthin", "end-product myxoxanthophylls", "end-product nostoxanthin", "end-product zeaxanthin diglucoside", "mevalonate pathway", "staphyloaxanthin biosynthesis"],
    "Cell mobility": ["Chemotaxis", "Flagellum", "Adhesion"],
    "Biofilm formation": ["Biofilm PGA Synthesis protein", "Biofilm regulator BssS", "Colanic acid and Biofilm protein A", "Colanic acid and Biofilm transcriptional regulator", "Curli fimbriae biosynthesis"],
    "Bacterial secretion systems": ["Sec-SRP", "Twin Arginine Targeting", "Type I Secretion", "Type II Secretion", "Type III Secretion", "Type IV Secretion", "Type Vabc Secretion", "Type VI Secretion"],
    "Genetic competence": ["competence factors", "competence-related core components", "competence-related related components"],
    "Metal transporters": ["Cobalt transporter CbiMQ", "Cobalt transporter CbtA", "Cobalt transporter CorA", "Copper transporter CopA", "Fe-Mn transporter MntH", "Ferric iron ABC-type substrate-binding AfuA", "Ferrous iron transporter FeoB", "Nickel ABC-type substrate-binding NikA"],
    "Methanogenesis": ["Methanogenesis via CO2", "Methanogenesis via acetate", "Methanogenesis via dimethylamine", "Methanogenesis via dimethylsulfide, methanethiol, methylpropanoate", "Methanogenesis via methanol", "Methanogenesis via methylamine", "Methanogenesis via trimethylamine", "methanol dehydrogenase"],
    "Photosynthesis": ["anoxygenic type-I reaction center", "anoxygenic type-II reaction center", "Photosystem I", "Photosystem II"],
    "Arsenic reduction": ["Arsenic reduction"],
    "Transporters": ["transporter: phosphate", "transporter: phosphonate", "transporter: thiamin", "transporter: urea", "transporter: vitamin B12"],
    "Miscellaneous": ["bidirectional polyphosphate", "C-P lyase cleavage PhnJ", "CP-lyase complex", "CP-lyase operon", "Anaplerotic genes"],
}

CAT_COLORS = {
    "Carbon fixation": "#059669",
    "Carbohydrate metabolism": "#0891b2",
    "Carbon degradation": "#7c3aed",
    "Nitrogen metabolism": "#2563eb",
    "Sulfur metabolism": "#d97706",
    "Oxidative phosphorylation": "#dc2626",
    "Hydrogen redox": "#db2777",
    "Amino acid metabolism": "#4f46e5",
    "Vitamin biosynthesis": "#ea580c",
    "Cell mobility": "#9333ea",
    "Biofilm formation": "#16a34a",
    "Bacterial secretion systems": "#c026d3",
    "Genetic competence": "#65a30d",
    "Metal transporters": "#ca8a04",
    "Methanogenesis": "#0d9488",
    "Photosynthesis": "#0284c7",
    "Arsenic reduction": "#be123c",
    "Transporters": "#475569",
    "Miscellaneous": "#94a3b8",
}


def main():
    tsv = Path(sys.argv[1])
    out = tsv.with_suffix(".png").as_posix()
    if len(sys.argv) > 2:
        out = sys.argv[2]

    with open(tsv) as f:
        rdr = list(csv.reader(f, delimiter="\t"))
    hdrs = rdr[0][1:]
    vals = [float(v) for v in rdr[1][1:]]
    sample_name = rdr[1][0]
    if not hdrs:
        fig, ax = plt.subplots(figsize=(7, 3))
        ax.axis("off")
        ax.text(0.5, 0.58, "KEGGaNOG unavailable", ha="center", va="center",
                fontsize=15, fontweight="bold", color="#1d2a57")
        ax.text(0.5, 0.38, "No usable KO annotations were available for this sample.",
                ha="center", va="center", fontsize=10, color="#5b6472")
        fig.savefig(out, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print("Saved: %s" % out)
        return

    # Assign categories
    name_to_pathway = {k.lower(): k for k in hdrs}
    pathway_to_category = {}
    for cat, pathways in CATEGORIES.items():
        for p in pathways:
            key = p.lower()
            if key in name_to_pathway:
                pathway_to_category[name_to_pathway[key]] = cat

    # Build rows with spacer between categories: (category, pathway, value)
    rows = []
    for name, val in zip(hdrs, vals):
        cat = pathway_to_category.get(name, "Miscellaneous")
        rows.append((cat, name, val))

    # Sort: by category order, then by value descending within category
    cat_order = {cat: i for i, cat in enumerate(CATEGORIES.keys())}
    rows.sort(key=lambda r: (-r[2], cat_order.get(r[0], 999)))
    rows.sort(key=lambda r: cat_order.get(r[0], 999))

    # Insert blank spacers between categories
    spaced_rows = []
    prev_cat = None
    for row in rows:
        if prev_cat is not None and row[0] != prev_cat:
            spaced_rows.append(("", "", -1))  # spacer
        spaced_rows.append(row)
        prev_cat = row[0]
    rows = spaced_rows

    n = len(rows)
    fig, ax = plt.subplots(figsize=(8, max(8, n * 0.28)))
    plt.subplots_adjust(left=0.38, right=0.95, top=0.97, bottom=0.03)

    # Draw bars: skip spacers (val < 0)
    y_positions = list(range(n))
    colors = [CAT_COLORS.get(r[0], "#94a3b8") if r[0] else "#ffffff00" for r in rows]
    vals_bars = [max(0, r[2]) for r in rows]
    bars = ax.barh(y_positions, vals_bars, color=colors, height=0.7)

    # Category labels on spacers
    tick_labels = []
    for i, (cat, name, val) in enumerate(rows):
        if val < 0:
            # Spacer row: show next category label here
            next_cat = None
            for j in range(i + 1, len(rows)):
                if rows[j][2] >= 0:
                    next_cat = rows[j][0]
                    break
            if next_cat:
                ax.annotate(next_cat, xy=(0, i), xycoords=("axes fraction", "data"),
                            xytext=(-28, 0), textcoords="offset points",
                            ha="right", va="center", fontsize=7, fontweight="bold",
                            color=CAT_COLORS.get(next_cat, "#1f2937"),
                            annotation_clip=False)
            tick_labels.append("")
        else:
            tick_labels.append(name)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(tick_labels, fontsize=8)
    ax.invert_yaxis()

    ax.set_xlabel("Pathway completeness", fontsize=11)
    real_n = sum(1 for r in rows if r[2] >= 0)
    ax.set_title(f"KEGG Pathway Completeness — {sample_name}\n({real_n} pathways with non-zero values)", fontsize=13, fontweight="bold")
    ax.set_xlim(0, 1.10)
    ax.grid(axis="x", alpha=0.3)

    for bar, val in zip(bars, [r[2] for r in rows]):
        if val >= 0:
            ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height() / 2,
                    f"{val:.2f}", va="center", fontsize=7, color="#1f2937")

    plt.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    print("Saved: %s" % out)

if __name__ == "__main__":
    main()
