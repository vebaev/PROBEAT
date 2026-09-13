#!/usr/bin/env /home/vebaev/codex/mamba_envs/goprobioticprofiler/bin/python
"""Generate PNG visualizations from GOProbioticProfiler JSON outputs."""

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag

# Probiotic GO Slim mapping (copied to avoid goatools import)
PROBIOTIC_GO_SLIM = {
    "Gut Survival": {
        "GO:0009408": "response_to_heat",
        "GO:0034605": "cellular_response_to_heat",
        "GO:0006970": "response_to_osmotic_stress",
        "GO:0034606": "cellular_response_to_osmotic_stress",
        "GO:0006979": "response_to_oxidative_stress",
        "GO:0034599": "cellular_response_to_oxidative_stress",
        "GO:0006972": "hyperosmotic_response",
        "GO:0010447": "response_to_acidic_pH",
        "GO:0071468": "cellular_response_to_acidic_pH",
        "GO:1990451": "cellular_stress_response_to_acidic_pH",
        "GO:0015721": "bile_acid_and_bile_salt_transport",
        "GO:0015125": "bile_acid_transmembrane_transporter_activity",
    },
    "Carbohydrate Utilization": {
        "GO:0005975": "carbohydrate_metabolic_process",
        "GO:0004553": "hydrolase_activity_hydrolyzing_O_glycosyl_compounds",
        "GO:0016798": "hydrolase_activity_acting_on_glycosyl_bonds",
        "GO:0000016": "lactose_metabolic_process",
    },
    "Host Interaction": {
        "GO:0030159": "cellular_adhesion",
        "GO:0044406": "adhesion_of_symbiont_to_host",
        "GO:0044036": "cell wall organization",
        "GO:0009275": "gram_positive_bacterium-type_cell_wall",
        "GO:0044277": "cellular_response_to_bacterium",
        "GO:0048870": "cell_motility",
        "GO:0042710": "biofilm_formation",
        "GO:0009372": "quorum_sensing",
    },
    "EPS Production": {
        "GO:0000271": "polysaccharide_biosynthetic_process",
        "GO:0046353": "exopolysaccharide_biosynthetic_process",
        "GO:0045227": "capsule_polysaccharide_biosynthetic_process",
    },
    "Vitamin Production": {
        "GO:0009110": "vitamin_biosynthetic_process",
        "GO:0051188": "cofactor_biosynthetic_process",
        "GO:0006775": "fat-soluble_vitamin_biosynthetic_process",
        "GO:0009235": "cobalamin_biosynthetic_process",
    },
    "Gut-Brain Axis": {
        "GO:0006563": "L-serine_metabolic_process",
        "GO:0006565": "L-tryptophan_metabolic_process",
        "GO:0006540": "GABA_shunt",
        "GO:0006596": "polyamine_catabolic_process",
        "GO:0006598": "polyamine_biosynthetic_process",
    },
    "SCFA Potential": {
        "GO:0015949": "nucleobase-containing_small_molecule_interconversion",
        "GO:0019654": "acetate_metabolic_process",
        "GO:0006085": "acetyl-CoA_metabolic_process",
        "GO:0006089": "lactate_metabolic_process",
        "GO:0006090": "pyruvate_metabolic_process",
        "GO:0019543": "propionate_metabolic_process",
        "GO:0046353": "exopolysaccharide_biosynthetic_process",
        "GO:0005975": "carbohydrate_metabolic_process",
        "GO:0006113": "fermentation",
    },
    "Antioxidant Capacity": {
        "GO:0016209": "antioxidant_activity",
        "GO:0004791": "thioredoxin_disulfide_reductase_activity",
        "GO:0015035": "protein_disulfide_oxidoreductase_activity",
        "GO:0006749": "glutathione_metabolic_process",
    },
    "Immunomodulation": {
        "GO:0002376": "immune_system_process",
        "GO:0044267": "cellular_protein_metabolic_process",
        "GO:0006477": "protein_amino_acid_lipidation",
        "GO:0000032": "cell_wall_mannoprotein_biosynthetic_process",
    },
    "Antimicrobial Potential": {
        "GO:0042742": "defense_response_to_bacterium",
        "GO:0017000": "antibiotic_biosynthetic_process",
        "GO:0009253": "peptidoglycan_catabolic_process",
    },
    "Metal and Environmental Stress": {
        "GO:0046688": "response_to_copper_ion",
        "GO:0010038": "response_to_metal_ion",
        "GO:0010044": "response_to_metal_ion_stress",
        "GO:0006826": "iron_ion_transport",
        "GO:0000041": "transition_metal_ion_transport",
        "GO:0042221": "response_to_chemical",
    },
}

NAMESPACE_LABELS = {
    "biological_process": "Biological Process",
    "molecular_function": "Molecular Function",
    "cellular_component": "Cellular Component",
}
NAMESPACE_COLORS = {
    "biological_process": "#2A9D8F",
    "molecular_function": "#E9C46A",
    "cellular_component": "#457B9D",
}


def plot_general_go_bubble(tsv_path: Path, out_png: Path) -> None:
    """Genome-wide GO landscape: coverage, gene count, and GO branch."""
    df = pd.read_csv(tsv_path, sep="\t")
    if df.empty:
        return
    df = df.sort_values(["Namespace", "Genome_Coverage_Percent"], ascending=[True, True])
    labels = [f"{name} ({go})" for name, go in zip(df["GO_Name"], df["GO_ID"])]
    y = np.arange(len(df))
    sizes = 35 + 14 * np.sqrt(df["Gene_Count"].astype(float))
    colors = [NAMESPACE_COLORS.get(ns, "#8D99AE") for ns in df["Namespace"]]
    fig, ax = plt.subplots(figsize=(12, max(8, len(df) * 0.28)))
    ax.scatter(df["Genome_Coverage_Percent"], y, s=sizes, c=colors, alpha=0.82,
               edgecolors="white", linewidths=0.8)
    ax.set_yticks(y, labels, fontsize=7)
    ax.set_xlabel("GO-annotated genome coverage (%)")
    ax.set_title("Genome GO Functional Landscape", fontsize=13, fontweight="bold")
    ax.grid(axis="x", color="#D9E2EC", linewidth=0.6)
    handles = [
        plt.Line2D([], [], marker="o", linestyle="", color=color, label=NAMESPACE_LABELS[ns])
        for ns, color in NAMESPACE_COLORS.items()
    ]
    ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_probiotic_go_bubble(tsv_path: Path, out_png: Path) -> None:
    """Category-level probiotic GO potential with explicit evidence tiers."""
    df = pd.read_csv(tsv_path, sep="\t")
    if df.empty:
        return
    tier_order = {"Significant": 0, "Supportive": 1, "Detected": 2}
    df["_tier_order"] = df["Evidence_Tier"].map(tier_order).fillna(3)
    df = df.sort_values(
        ["_tier_order", "Category_Coverage_Percent", "PGMDB_Gene_Count"],
        ascending=[False, True, True],
    )
    labels = list(df["Probiotic_Category"])
    y = np.arange(len(df))
    sizes = 100 + 55 * np.sqrt(df["PGMDB_Gene_Count"].clip(lower=1).astype(float))
    tier_colors = {
        "Significant": "#178A55",
        "Supportive": "#E69F00",
        "Detected": "#5B7FA3",
    }
    colors = [tier_colors.get(tier, "#94A3B8") for tier in df["Evidence_Tier"]]
    edges = ["#0F5132" if tier == "Significant" else "#334155" for tier in df["Evidence_Tier"]]
    widths = [1.8 if tier == "Significant" else 0.8 for tier in df["Evidence_Tier"]]
    fig, ax = plt.subplots(figsize=(11.5, 7.2))
    ax.scatter(
        df["Fold_Enrichment"], y, s=sizes, c=colors, alpha=0.90,
        edgecolors=edges, linewidths=widths, zorder=3,
    )
    ax.axvline(1.0, color="#64748B", linestyle="--", linewidth=0.8)
    ax.set_yticks(y, labels, fontsize=8.5)
    ax.set_xlabel("Fold enrichment in PGMDB-supported genes")
    ax.set_title("Probiotic GO Functional Potential", fontsize=13, fontweight="bold")
    ax.grid(axis="x", color="#D9E2EC", linewidth=0.6, zorder=0)
    x_max = max(float(df["Fold_Enrichment"].max()), 1.0)
    ax.set_xlim(left=0, right=x_max * 1.45 + 0.25)
    handles = [
        plt.Line2D(
            [], [], marker="o", linestyle="", markersize=8,
            markerfacecolor=color, markeredgecolor="#334155", label=tier,
        )
        for tier, color in tier_colors.items()
    ]
    ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=8)
    ax.text(
        0.0, -0.13,
        "Significant: FDR <= 0.05 | Supportive: FDR <= 0.20, fold >= 2, and >= 2 PGMDB genes",
        transform=ax.transAxes, fontsize=7, color="#475569",
    )
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_go_enrichment(json_path: Path, out_png: Path, tsv_path: Path, top_n: int = 20) -> None:
    with json_path.open() as fh:
        data = json.load(fh)
    go_ids = data["x"][:top_n]
    fold = data["y"][:top_n]

    # Read GO names from the TSV
    go_names = {}
    with tsv_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            go_names[row["GO_ID"]] = row["GO_Name"]

    labels = [go_names.get(go_id, go_id) for go_id in go_ids]

    plt.figure(figsize=(12, 8))
    colors = plt.cm.RdYlGn(np.linspace(0.2, 0.8, len(go_ids)))[::-1]
    bars = plt.barh(range(len(go_ids)), fold, color=colors)
    plt.yticks(range(len(go_ids)), labels, fontsize=8)
    plt.xlabel("Fold Enrichment", fontsize=10)
    plt.title("Top Enriched GO Terms", fontsize=12, fontweight="bold")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()


def plot_go_enrichment_by_category(tsv_path: Path, out_png: Path, obo_path: Path) -> None:
    """Aggregate enriched GO terms into probiotic GO Slim categories using keyword matching and GO ancestors."""
    # Load GO DAG
    godag = GODag(obo_path)
    # Build keyword -> category mapping
    keyword_map = {
        "Gut Survival": ["stress", "heat", "osmotic", "oxidative", "pH", "bile", "acid"],
        "Carbohydrate Utilization": ["carbohydrate", "glycosyl", "lactose", "sugar", "starch"],
        "Host Interaction": ["adhesion", "cell wall", "biofilm", "motility", "symbiont", "external encapsulating structure"],
        "EPS Production": ["polysaccharide", "exopolysaccharide", "capsule"],
        "Vitamin Production": ["vitamin", "cobalamin", "cofactor"],
        "Gut-Brain Axis": ["serine", "tryptophan", "GABA", "polyamine", "gamma-aminobutyric"],
        "SCFA Potential": ["acetate", "propionate", "butyrate", "SCFA", "fermentation", "monocarboxylic"],
        "Antioxidant Capacity": ["antioxidant", "thioredoxin", "glutathione", "disulfide"],
        "Immunomodulation": ["immune", "lipidation", "mannoprotein", "proteolysis", "protein metabolic"],
        "Metal and Environmental Stress": ["copper", "iron", "metal", "transition metal", "chemical"],
    }

    # Read enriched GO terms
    enriched = []
    with tsv_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            fdr = row.get("FDR", "")
            try:
                fdr_val = float(fdr) if fdr else 1.0
            except ValueError:
                fdr_val = 1.0
            if fdr_val < 0.05:
                fold = float(row.get("Fold_Enrichment", 0))
                enriched.append((row["GO_ID"], fold, row.get("GO_Name", ""), row.get("Namespace", "")))

    # Map each enriched GO term to probiotic categories via keywords and GO ancestors
    category_fold = defaultdict(float)
    category_count = defaultdict(int)
    for go_id, fold, name, ns in enriched:
        go = godag.get(go_id)
        ancestors = set()
        if go is not None:
            ancestors = go.get_all_parents() | {go_id}
        # Get names for ancestors
        ancestor_names = []
        for a in ancestors:
            a_go = godag.get(a)
            if a_go is not None:
                ancestor_names.append(a_go.name.lower())
        name_lower = name.lower()
        matched = False
        seen_cats = set()
        for category, keywords in keyword_map.items():
            for kw in keywords:
                if kw in name_lower or any(kw in an for an in ancestor_names):
                    if category not in seen_cats:
                        seen_cats.add(category)
                        category_fold[category] += fold
                        category_count[category] += 1
                    matched = True
                    break
        if not matched:
            category_fold["Other"] += fold
            category_count["Other"] += 1

    # Sort by fold enrichment descending
    cats = sorted(category_fold.keys(), key=lambda c: category_fold[c], reverse=True)
    folds = [category_fold[c] for c in cats]
    counts = [category_count[c] for c in cats]

    # Create labels with counts
    labels = [f"{cat}\n({counts[i]} GO terms)" for i, cat in enumerate(cats)]

    plt.figure(figsize=(10, 6))
    colors = plt.cm.Set2(np.linspace(0, 1, len(cats)))
    bars = plt.barh(range(len(cats)), folds, color=colors)
    plt.yticks(range(len(cats)), labels, fontsize=8)
    plt.xlabel("Total Fold Enrichment", fontsize=10)
    plt.title("Enriched GO Terms by Probiotic Category", fontsize=12, fontweight="bold")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()


def plot_go_slim(json_path: Path, out_png: Path) -> None:
    with json_path.open() as fh:
        data = json.load(fh)
    cats = data["x"]
    counts = data["y"]
    plt.figure(figsize=(10, 6))
    colors = plt.cm.Set3(np.linspace(0, 1, len(cats)))
    bars = plt.barh(range(len(cats)), counts, color=colors)
    plt.yticks(range(len(cats)), cats, fontsize=8)
    plt.xlabel("Gene Count", fontsize=10)
    plt.title("GO Slim Abundance", fontsize=12, fontweight="bold")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()


def plot_probiotic_heatmap(json_path: Path, out_png: Path) -> None:
    with json_path.open() as fh:
        data = json.load(fh)
    cats = data["categories"]
    values = data["values"]

    # Aggregate duplicates (same category can appear with different probiotic functions)
    agg = defaultdict(int)
    for cat, val in zip(cats, values):
        agg[cat] += val
    cats = list(agg.keys())
    values = list(agg.values())

    # Sort by value descending
    pairs = sorted(zip(cats, values), key=lambda x: x[1], reverse=True)
    cats, values = zip(*pairs) if pairs else ([], [])

    plt.figure(figsize=(8, 6))
    colors = plt.cm.YlOrRd(np.linspace(0.2, 0.8, len(values)))
    bars = plt.barh(range(len(cats)), values, color=colors)
    plt.yticks(range(len(cats)), cats, fontsize=8)
    plt.xlabel("Gene Count", fontsize=10)
    plt.title("Probiotic Category Gene Count", fontsize=12, fontweight="bold")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()


def plot_sankey(json_path: Path, out_png: Path) -> None:
    with json_path.open() as fh:
        data = json.load(fh)
    nodes = data["nodes"]
    sources = data["sources"]
    targets = data["targets"]
    values = data["values"]

    # Sankey via matplotlib is complex; use a simple layered graph instead
    layers = defaultdict(set)
    node_layer = {}
    for n in nodes:
        if n in ["Gut Survival", "Carbohydrate Utilization", "Host Interaction", "EPS Production",
                 "Vitamin Production", "Gut-Brain Axis", "SCFA Potential", "Antioxidant Capacity",
                 "Immunomodulation", "Metal and Environmental Stress"]:
            node_layer[n] = 2
            layers[2].add(n)
        else:
            # crude heuristic: if target in a category, it's a function (layer 1)
            pass
    for i, n in enumerate(nodes):
        if n not in node_layer:
            node_layer[n] = 0
            layers[0].add(n)
    for i, n in enumerate(nodes):
        if n not in node_layer:
            node_layer[n] = 1
            layers[1].add(n)

    fig, ax = plt.subplots(figsize=(14, 10))
    layer_x = {0: 0.1, 1: 0.5, 2: 0.9}
    node_y = {}
    for layer_id, x in layer_x.items():
        ns = sorted(layers[layer_id])
        n = len(ns)
        for i, node in enumerate(ns):
            y = 0.1 + (0.8 * i / max(n - 1, 1))
            node_y[node] = (x, y)
            ax.text(x, y, node, fontsize=6, ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", edgecolor="gray"))

    for s, t, v in zip(sources, targets, values):
        src = nodes[s]
        tgt = nodes[t]
        if src in node_y and tgt in node_y:
            x1, y1 = node_y[src]
            x2, y2 = node_y[tgt]
            ax.plot([x1, x2], [y1, y2], "k-", alpha=0.1, linewidth=0.5)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    plt.title("Gene → Function → GO Slim Category", fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()


def main():
    outdir = Path(sys.argv[1])
    obo_path = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("data/databases/go/go-basic.obo")

    plot_go_enrichment(outdir / "plot_go_enrichment.json", outdir / "go_enrichment.png", outdir / "go_enrichment.tsv")
    plot_go_enrichment_by_category(outdir / "go_enrichment.tsv", outdir / "go_enrichment_by_category.png", obo_path)
    plot_go_slim(outdir / "plot_go_slim.json", outdir / "go_slim.png")
    plot_probiotic_heatmap(outdir / "plot_probiotic_heatmap.json", outdir / "probiotic_heatmap.png")
    plot_sankey(outdir / "plot_sankey.json", outdir / "sankey.png")
    plot_general_go_bubble(outdir / "general_go_profile.tsv", outdir / "general_go_bubble.png")
    plot_probiotic_go_bubble(outdir / "probiotic_go_categories.tsv", outdir / "probiotic_go_bubble.png")


if __name__ == "__main__":
    main()
